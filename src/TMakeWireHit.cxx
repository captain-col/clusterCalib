#include "TMakeWireHit.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TFADCHit.hxx>
#include <HEPUnits.hxx>
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

CP::TMakeWireHit::TMakeWireHit(bool correctEfficiency) {
    fCorrectCollectionEfficiency = correctEfficiency;
}
CP::TMakeWireHit::~TMakeWireHit() { }

CP::THandle<CP::THit>
CP::TMakeWireHit::operator()(const CP::TCalibPulseDigit& digit, double step,
                             std::size_t beginIndex, std::size_t endIndex) {
    double charge = 0.0;
    int samples = 1;

    /// Protect against out of bounds problems.
    if (digit.GetSampleCount() <= endIndex) endIndex = digit.GetSampleCount()-1;
    if (endIndex < beginIndex) return CP::THandle<THit>();

    TChannelCalib calib;
    
    // Find the start and stop time.
    double startTime = beginIndex*step + digit.GetFirstSample();
    double stopTime = (endIndex+1)*step + digit.GetFirstSample();

    // Copy the samples into the hit.
    std::vector<double> sampleCharge(endIndex-beginIndex+1);
    for (std::size_t i = beginIndex; i<=endIndex; ++i) {
        sampleCharge[i-beginIndex] = digit.GetSample(i);
    }

    // Find the peak charge.
    for (std::size_t j=beginIndex; j<=endIndex; ++j) {
        double v = digit.GetSample(j);
        if (v <= 0) continue;
        charge += v;
        ++samples;
    }

    // Protect against empty regions.
    if (charge<1) return CP::THandle<CP::THit>();

    double fwhm = step;
    // Estimate the RMS using the FWHM.
    std::vector<double>::iterator maxBin = sampleCharge.end();
    for (std::vector<double>::iterator s = sampleCharge.begin();
         s != sampleCharge.end(); ++s) {
        if (maxBin==sampleCharge.end() || *maxBin < *s) maxBin = s;
    }
    if (maxBin != sampleCharge.end()) {
        std::vector<double>::iterator lowBin = maxBin;
        while (lowBin != sampleCharge.begin()) {
            --lowBin;
            if (*lowBin < 0.5*(*maxBin)) break;
        }
        std::vector<double>::iterator hiBin = maxBin;
        while (hiBin != sampleCharge.end()) {
            if (*hiBin < 0.5*(*maxBin)) break;
            if (hiBin+1 == sampleCharge.end()) break;
            ++hiBin;
        }
        double lowVal=(0.5*(*maxBin)-(*lowBin))/(*(lowBin+1)-(*lowBin));
        double hiVal = (0.5*(*maxBin)-(*hiBin))/(*(hiBin-1)-(*hiBin));
        int diff = hiBin-lowBin;
        // Convert FWHM into an RMS.  The 2.36 factor is the ratio between
        // the rms and the FWHM for a Gaussian peak.
        fwhm = 1.0*(diff - lowVal - hiVal)*step;
        if (fwhm > 1) fwhm *= step;
    }
    double rms = fwhm/2.36;
    
    double time = 0.0;
    // Find the time by looking at samples near the peak.
    maxBin = sampleCharge.end();
    for (std::vector<double>::iterator s = sampleCharge.begin();
         s != sampleCharge.end(); ++s) {
        if (maxBin==sampleCharge.end() || *maxBin < *s) maxBin = s;
    }
    if (maxBin != sampleCharge.end()) {
        std::vector<double>::iterator lowBin = maxBin;
        while (lowBin != sampleCharge.begin()) {
            --lowBin;
            if (2.0*rms < step*(maxBin-lowBin)) break;
        }
        std::vector<double>::iterator hiBin = maxBin;
        while (hiBin != sampleCharge.end()) {
            if (2.0*rms < step*(hiBin-maxBin)) break;
            ++hiBin;
        }
        double w = 0.0;
        int bin = 0;
        for (std::vector<double>::iterator b = lowBin; b != hiBin; ++b) {
            time += (*b)*bin;
            w += (*b);
            ++bin;
        }
        time *= step/w;
        time += step*(lowBin-sampleCharge.begin()) + startTime + step/2.0;
    }
    
    // Base the uncertainty in the time on the number of samples used to find
    // the RMS.
    double timeUnc = rms/std::sqrt(samples);

    // Calculate the charge uncertainty (the calculate starts with chargeUnc
    // being the variance, and then takes the sqrt).  The fundamental
    // uncertainty is in the number of electrons (Poisson distributed).
    double chargeUnc = charge/unit::eplus;

    // Add the RMS for the baseline uncertainty.  This is analogous
    // to the deviation of the baseline from a "bestfit" straight line (the
    // s_plane parameter in the PDG notation).  The correlations are
    // introduced by the electronics shaping and the amount of correlation
    // depends on the number of samples being integrated over.
    double sigC = 0.0;   // Not yet implemented!!!!
    chargeUnc += sigC*sigC;

    // Now take the sqrt of the variance to get the uncertainty.
    chargeUnc = std::sqrt(chargeUnc);

    CaptNamedInfo("TMakeWireHit", digit.GetChannelId() << " Q: " << charge
                  << " Samples: " << samples
                  << " Sigma: " << chargeUnc
                  << " (" << sigC << ")");
    
    CP::TGeometryId geomId
        = CP::TChannelInfo::Get().GetGeometry(digit.GetChannelId());
    
    if (!geomId.IsValid()) {
        // CaptError("Making hits for an invalid geometry id");
        // return CP::THandle<CP::THit>();
    }

    if (!std::isfinite(rms) || rms <= 0.0) {
        CaptError("Time RMS for " << digit.GetChannelId()
                  << " is not positive and finite ("
                  << rms << " out of " << sampleCharge.size() << ")");
        timeUnc = 1*unit::second;
    }

    if (!std::isfinite(timeUnc) || timeUnc <= 0.0) {
        CaptError("Time uncertainty for " << digit.GetChannelId()
                  << " is not positive and finite (" << timeUnc << ")");
        timeUnc = 1*unit::second;
    }

        if (!std::isfinite(chargeUnc) || chargeUnc <= 0.0) {
        CaptError("Charge uncertainty for " << digit.GetChannelId()
                  << " is not positive and finite " << samples);
        chargeUnc = 1*unit::coulomb;
    }

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
    if (!std::isfinite(chargeUnc) || chargeUnc <= 0.0) {
        CaptError("Invalid hit charge uncertainty for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(time)) {
        CaptError("Invalid hit time for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(timeUnc)) {
        CaptError("Invalid hit time uncertainty for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(rms)) {
        CaptError("Invalid hit time RMS for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(startTime)) {
        CaptError("Invalid hit start time for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(stopTime)) {
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
    hit.SetTime(time);
    hit.SetTimeRMS(rms);
    hit.SetTimeStart(startTime);
    hit.SetTimeStop(stopTime);
    hit.SetTimeSamples(sampleCharge.begin(), sampleCharge.end());
    hit.SetTimeUncertainty(timeUnc);

    return CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit));
}
