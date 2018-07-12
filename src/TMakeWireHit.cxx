#include "TMakeWireHit.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TFADCHit.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <THandle.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>
#include <CaptGeomId.hxx>
#include <TChannelCalib.hxx>

#include <TSpectrum.h>
#include <TH1F.h>
#include <TH2F.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <sstream>
#include <iostream>

namespace {
    // A least squares fit to the shape of a set of samples to get the
    // curvature of the data.  This assumes that [begin, end] includes the
    // peak of the function.  This is just encapsulating the fit so that the
    // code is cleaner.  It's based of of the maxima output from
    // polynomial_peak_fit.mac.  The peak is fitted to a polynomial
    //
    // y ~ (a-b*x[i]-c*x[i]^2)
    //
    // This only returns the peak position and curvature (the height is
    // ignored).
    std::pair<double,double>
    FitPeakCurvature(double Xoff, double step,
                     std::vector<double>::const_iterator begin,
                     std::vector<double>::const_iterator peak,
                     std::vector<double>::const_iterator end) {
        double Xpeak = Xoff + (peak-begin)*step;

        // Only ever fit right around the peak.
        if (peak-begin >= 3) begin = peak - 3;
        else if (peak-begin >= 2) begin = peak - 2;
        else if (peak-begin >= 1) begin = peak - 1;

        if (end-peak >= 3) end = peak + 3;
        else if (end-peak >= 2) end = peak + 2;
        else if (end-peak >= 1) end = peak + 1;
        
        Xoff = Xpeak - step*(peak-begin);

        //   sYi = sum(Y[i],i,0,m)
        double sYi = 0;
        //   sXiYi = sum(X[i]*Y[i],i,0,m)
        double sXiYi = 0;
        //   sX2iYi = sum(X[i]^2*Y[i],i,0,m)
        double sX2iYi = 0;
        //   sXi = sum(X[i],i,0,m)
        double sXi = 0;
        //   sX2i = sum(X[i]^2,i,0,m)
        double sX2i = 0;
        //   sX3i = sum(X[i]^3,i,0,m)
        double sX3i = 0;
        //   sX4i = sum(X[i]^4,i,0,m)
        double sX4i = 0;
        double Xi = 0.0;
        double m = 0;
        do {
            double Yi = (*begin)/(*peak);
            m += 1.0;
            sYi += Yi;
            sXiYi += Xi*Yi;
            sX2iYi += Xi*Xi*Yi;
            sXi += Xi;
            sX2i += Xi*Xi;
            sX3i += Xi*Xi*Xi;
            sX4i += Xi*Xi*Xi*Xi;
            if (begin == end) break;
            Xi += 1.0;
            ++begin;
        } while (true);
        // Find the general solution...  a can be calculated, but not used.
#ifdef CALCULATE_CONSTANT_FACTOR
        double a = -std::pow(sX4i*std::pow(sXi,2)-2*sX2i*sX3i*sXi-m*sX2i*sX4i+m*std::pow(sX3i,2)+std::pow(sX2i,3),-1)*((sX2i*sX4i-std::pow(sX3i,2))*sYi+sXi*(sX2iYi*sX3i-sX4i*sXiYi)+sX2i*sX3i*sXiYi-std::pow(sX2i,2)*sX2iYi);
#endif
        double b = std::pow(sX4i*std::pow(sXi,2)-2*sX2i*sX3i*sXi-m*sX2i*sX4i+m*std::pow(sX3i,2)+std::pow(sX2i,3),-1)*((sX4i*sXi-sX2i*sX3i)*sYi+m*(sX2iYi*sX3i-sX4i*sXiYi)+std::pow(sX2i,2)*sXiYi-sX2i*sX2iYi*sXi);
        double c = -std::pow(sX4i*std::pow(sXi,2)-2*sX2i*sX3i*sXi-m*sX2i*sX4i+m*std::pow(sX3i,2)+std::pow(sX2i,3),-1)*((sX3i*sXi-std::pow(sX2i,2))*sYi+sX2i*sXi*sXiYi-m*sX3i*sXiYi-sX2iYi*std::pow(sXi,2)+m*sX2i*sX2iYi);
        // Peak position...
        double fitPeak = step*(-0.5*b/c);
        // curvature...
        double Xcurv = 2.0*c*step*step/(*peak);

        return std::make_pair(Xoff+fitPeak,Xcurv);
    }
};


CP::TMakeWireHit::TMakeWireHit(bool correctEfficiency) {
    fCorrectCollectionEfficiency = correctEfficiency;
}
CP::TMakeWireHit::~TMakeWireHit() { }

CP::THandle<CP::THit>
CP::TMakeWireHit::operator()(const CP::TCalibPulseDigit& digit, double step,
                             std::size_t beginIndex, std::size_t endIndex,
                             std::size_t extraSamples) {
    double charge = 0.0;
    int samples = 1;

    CP::TGeometryId geomId
        = CP::TChannelInfo::Get().GetGeometry(digit.GetChannelId());

    if (!geomId.IsValid()) {
        CaptWarn("Making hits for an invalid geometry id");
    }

    /// Protect against out of bounds problems.
    if (digit.GetSampleCount() <= endIndex) endIndex = digit.GetSampleCount()-1;
    if (endIndex <= beginIndex) return CP::THandle<THit>();

    TChannelCalib calib;
    
    // Copy samples that make up the pulse to be integrated into a vector.
    std::vector<double> sampleCharge(endIndex-beginIndex+1);
    for (std::size_t i = beginIndex; i<=endIndex; ++i) {
        sampleCharge[i-beginIndex] = digit.GetSample(i);
    }

    // Find the lower and upper bound time for the pulse.
    double startPulse = beginIndex*step + digit.GetFirstSample();
    double stopPulse = (sampleCharge.size()-1)*step + startPulse;

    // Find the peak charge.
    for (std::size_t j=beginIndex; j<=endIndex; ++j) {
        double v = digit.GetSample(j);
        if (v <= 0) continue;
        charge += v;
        ++samples;
    }

    // Protect against empty regions.
    if (charge<1) return CP::THandle<CP::THit>();

    // Find the baseline uncertainty.  This is based on the range outside of
    // the hit range.
    std::size_t baselineRange = 500;
    std::size_t baseLowIndex = beginIndex-baselineRange;
    if (beginIndex<=baselineRange) baseLowIndex = 1;
    std::size_t baseHiIndex = endIndex+baselineRange;
    if (baseHiIndex>digit.GetSampleCount()) {
        baseHiIndex = digit.GetSampleCount();
    }
    double baseline = 0.0;
    double baseline2 = 0.0;
    double baselineCount = 0.0;
    std::vector<double> baseVect;
    double baseDiff = 0.0;
    double baseDiff2 = 0.0;
    double baseDiffCount = 0.0;
    for (std::size_t j=baseLowIndex; j<=baseHiIndex; ++j) {
        if (beginIndex<=j && j<=endIndex) continue;
        double v = digit.GetSample(j);
        double dv = digit.GetSample(j)-digit.GetSample(j-1);
        baseDiff += dv;
        baseDiff2 += dv*dv;
        baseDiffCount += 1.0;
        if (std::abs(digit.GetSample(j-1))>std::abs(v)) continue;
        if (std::abs(digit.GetSample(j+1))>std::abs(v)) continue;
        baseVect.push_back(std::abs(v));
        baseline += v;
        baseline2 += v*v;
        baselineCount += 1.0;
    }
    double baselineUnc = 0.5*charge;
    if (baselineCount>5) {
        baseline /= baselineCount;
        baseline2 /= baselineCount;
        double baselineRMS=std::sqrt(baseline2-baseline*baseline);
        if (baseVect.size() > 10) {
            std::sort(baseVect.begin(),baseVect.end());
            // Exclude the biggest, and then take the 95% percentile.
            double vectRMS = baseVect[0.7*(baseVect.size()-3)];
            if (baselineRMS < vectRMS) baselineRMS = vectRMS;
        }
        baseDiff /= baseDiffCount;
        baseDiff2 /= baseDiffCount;
        double sampleUnc = std::sqrt(baseDiff2-baseDiff*baseDiff)/1.414;
        baselineUnc = sampleUnc*std::sqrt(samples);
        // Correct for correlations between samples.
        double b = std::abs(baseline) + baselineRMS;
        b = b*b - sampleUnc*sampleUnc;
        if (b > 0.0) b = std::sqrt(b);
        else b = 0.0;
        baselineUnc += b*samples;
    }
    
    double fwhm = step;
    // Find the peak, the bins that bracket the peak and the FWHM.
    std::vector<double>::iterator maxBin = sampleCharge.end();
    std::vector<double>::iterator lowBin = sampleCharge.end();
    std::vector<double>::iterator hiBin = sampleCharge.end();
    for (std::vector<double>::iterator s = sampleCharge.begin();
         s != sampleCharge.end(); ++s) {
        if (maxBin==sampleCharge.end() || *maxBin < *s) maxBin = s;
    }
    if (maxBin != sampleCharge.end()) {
        lowBin = maxBin;
        while (lowBin != sampleCharge.begin()) {
            --lowBin;
            if (*lowBin < 0.5*(*maxBin)) break;
        }
        hiBin = maxBin;
        while (hiBin != sampleCharge.end()) {
            if (*hiBin < 0.5*(*maxBin)) break;
            if (hiBin+1 == sampleCharge.end()) break;
            ++hiBin;
        }
        if (hiBin-lowBin<2) return CP::THandle<CP::THit>();
        double lowVal=(0.5*(*maxBin)-(*lowBin))/(*(lowBin+1)-(*lowBin));
        double hiVal = (0.5*(*maxBin)-(*hiBin))/(*(hiBin-1)-(*hiBin));
        int diff = hiBin-lowBin;
        // Convert FWHM into an RMS.  The 2.36 factor is the ratio between
        // the rms and the FWHM for a Gaussian peak.
        fwhm = 1.0*(diff - lowVal - hiVal)*step;
        if (fwhm < step) fwhm = step;
    }

    double timeRMS = 0.0;
    double timeAvg = 0.0;
    for (std::size_t j=0; j<sampleCharge.size(); ++j) {
        double v = sampleCharge[j];
        if (v <= 0) continue;
        ++samples;
        double t = step*j + startPulse;
        timeAvg += v*t;
        timeRMS += v*t*t;
    }
    timeAvg /= charge;
    timeRMS /= charge;
    timeRMS = timeRMS-timeAvg*timeAvg;
    if (timeRMS > 0.0) timeRMS = std::sqrt(timeRMS);
    else timeRMS = 0.0;

    // Base the uncertainty on the curvature at the peak.
    //
    // Nota Bene: If my reasoning, and algebra, is right then the time
    // uncertainty is just the one over curvature at the peak.  Since the peak
    // shape already includes the normalization of the peak, I think the
    // curvature correctly includes the usual (1/sqrt(n)) term.  This comment
    // is to note that I'm not terribly confident about this conclusion.
    double timeUnc = timeRMS;   // Have a default in case of numeric problems.
    std::pair<double,double> timePeak
        = FitPeakCurvature(step*(lowBin-sampleCharge.begin())+startPulse,step,
                           lowBin,maxBin,hiBin);
    // Protect against peaks with almost no curvature.  The uncertainty should
    // never be more than the peak RMS.
    if (std::isfinite(timePeak.second) && timePeak.second < 0.0) {
        timeUnc = -(*maxBin)/timePeak.second;
    }
    timeUnc = std::min(timeUnc, timeRMS);

    double timeHit = timePeak.first;

    // Check that the average and peak times are close together.  If they are
    // to far apart, move to the average and take the separation as the
    // uncertainty.
    if (std::abs(timeAvg-timeHit) > timeUnc) {
        double newUnc = std::abs(timeAvg-timeHit);
        double arg = newUnc/timeUnc;
        arg = std::exp(1.0-arg*arg);
        double newTime = ((1.0-arg)*timeAvg + arg*timeHit);
        timeUnc = newUnc;
        timeHit = newTime;
    }
   
    // Calculate the charge uncertainty (the calculate starts with chargeUnc
    // being the variance, and then takes the sqrt).  The fundamental
    // uncertainty is in the number of electrons (Poisson distributed).
    double chargeUnc = charge/unit::eplus;

    // Add the RMS for the baseline uncertainty.  This is analogous
    // to the deviation of the baseline from a "bestfit" straight line (the
    // s_plane parameter in the PDG notation).  The correlations are
    // introduced by the electronics shaping and the amount of correlation
    // depends on the number of samples being integrated over.
    double sigC = baselineUnc;
    chargeUnc += sigC*sigC;

    // Now take the sqrt of the variance to get the uncertainty.
    chargeUnc = std::sqrt(chargeUnc);

    CaptNamedInfo("TMakeWireHit", digit.GetChannelId() << " Q: " << charge
                  << " Samples: " << samples
                  << " Sigma: " << chargeUnc
                  << " (" << sigC << ")");
    
    // Correct for the wire collection efficiency.  Before this correction,
    // the wire is calibrated in terms of the charged collected (or induced)
    // on the wire.  After this correction, the wire is calibrated in terms of
    // "sensed" electrons which is the number of electrons passing in the
    // vicinity of the wire.  For the collection wires, the measured electrons
    // and sensed electrons are the same thing.
    if (fCorrectCollectionEfficiency) {
        charge /= calib.GetCollectionEfficiency(digit.GetChannelId());
        chargeUnc /= calib.GetCollectionEfficiency(digit.GetChannelId());
    }

    // Check the hit validity.
    bool isValid = true;
    if (!std::isfinite(charge)) {
        CaptError("Invalid hit charge for " << digit.GetChannelId());
        isValid = false;
    }

    if (!std::isfinite(timeRMS) || timeRMS <= 0.0) {
        CaptError("Time RMS for " << digit.GetChannelId()
                  << " is not positive and finite ("
                  << timeRMS << " out of " << sampleCharge.size() << ")");
        timeRMS = 1*unit::second;
        isValid = false;
    }

    if (!std::isfinite(timeUnc) || timeUnc <= 0.0) {
        CaptError("Time uncertainty for " << digit.GetChannelId()
                  << " is not positive and finite (" << timeUnc << ")");
        timeUnc = 1*unit::second;
        isValid = false;
    }
    
    if (!std::isfinite(chargeUnc) || chargeUnc <= 0.0) {
        CaptError("Charge uncertainty for " << digit.GetChannelId()
                  << " is not positive and finite " << chargeUnc);
        chargeUnc = 1*unit::coulomb;
        isValid = false;
    }

    if (!std::isfinite(startPulse)) {
        CaptError("Invalid hit start time for "
                  << digit.GetChannelId());
        isValid = false;
    }

    if (!std::isfinite(stopPulse)) {
        CaptError("Invalid hit stop time for "
                  << digit.GetChannelId());
        isValid = false;
    }

    if (sampleCharge.size()<1) {
        CaptError("Invalid hit charge samples "
                  << digit.GetChannelId()
                  << " (samples == "<< sampleCharge.size() << ")");
        isValid = false;
    }

    if (!isValid) return CP::THandle<CP::TFADCHit>();

    // Build the hit.
    CP::TWritableFADCHit hit;
    hit.SetGeomId(geomId);
    hit.SetDigit(digit.GetParent());
    hit.SetCharge(charge);
    hit.SetChargeUncertainty(chargeUnc);
    hit.SetTime(timeHit);
    hit.SetTimeUncertainty(timeUnc);
    hit.SetTimeRMS(timeRMS);
    hit.SetTimeLowerBound(startPulse);
    hit.SetTimeUpperBound(stopPulse);
    if (!calib.IsGoodChannel(digit.GetChannelId())) {
        hit.SetTimeValidity(false);
        hit.SetChargeValidity(false);
    }
    if (!calib.IsGoodWire(geomId)) {
        hit.SetTimeValidity(false);
        hit.SetChargeValidity(false);
    }

    // Find the full width of the peak.  A region equal to one peak full width
    // is added as a side band to the peak shape.  If the full width is less
    // than 4 (i.e. twice the typical pulse rise time), then force it to be
    // four.
    int fullWidth = endIndex - beginIndex;
    if (fullWidth<4) fullWidth = 4;

    // Find the start and stop of the samples to be added to the hit.  The
    // stopIndex is one past the end of the samples.
    int startIndex = beginIndex  - fullWidth - extraSamples;
    if (startIndex<0) startIndex = 0;
    std::size_t stopIndex = endIndex + fullWidth + extraSamples + 1;
    if (stopIndex > digit.GetSampleCount()) stopIndex = digit.GetSampleCount();
        
    // Find the lower and upper bound time for the pulse.
    double startTime = startIndex*step + digit.GetFirstSample();
    double stopTime = stopIndex*step + digit.GetFirstSample();
    hit.SetTimeStart(startTime);
    hit.SetTimeStop(stopTime);
    hit.SetTimeSamples(digit.GetSamples().begin()+startIndex,
                       digit.GetSamples().begin()+stopIndex);

    return CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit));
}
