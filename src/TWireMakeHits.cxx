#include "TWireMakeHits.hxx"
#include "TChannelCalib.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TFADCHit.hxx>
#include <HEPUnits.hxx>
#include <THandle.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>

#include <TSpectrum.h>
#include <TH1F.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <sstream>
#include <iostream>

CP::TWireMakeHits::TWireMakeHits(bool correctLifetime,
                                 bool correctEfficiency) {
    fCorrectElectronLifetime = correctLifetime;
    fCorrectCollectionEfficiency = correctEfficiency;
    fNSource = 0;
    fSource = NULL;
    fDest = NULL;
    fWork = NULL;

    fPeakMaximumCol
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.charge.collection");

    fPeakMaximumInd
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.charge.induction");

    fPeakAreaCol
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.area.collection");

    fPeakAreaInd
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.area.induction");

    fPeakWidthCol
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.width.collection");

    fPeakWidthInd
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.width.induction");

    fNoiseThresholdCut
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.noise");

    fPeakRMSLimit
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.rmsLimit");

    fDigitEndSkip
        = CP::TRuntimeParameters::Get().GetParameterI(
            "clusterCalib.peakSearch.endSkip");

    fIntegrationChargeThreshold = 0.001;

    fIntegrationNoiseThreshold = 0.001;

}
CP::TWireMakeHits::~TWireMakeHits() {
    if (fSource) delete[] fSource;
    if (fDest) delete[] fDest;
    if (fWork) delete[] fWork;
}

CP::THandle<CP::THit>
CP::TWireMakeHits::MakeHit(const CP::TCalibPulseDigit& digit,
                           double step, double t0,
                           double baselineSigma, double sampleSigma,
                           int beginIndex, int endIndex, bool split) {
    double charge = 0.0;
    double sample = 0.0;
    double sampleSquared = 0.0;
    double samples = 1.0;

    // Find the start and stop time.
    double startTime = beginIndex*step + digit.GetFirstSample();
    double stopTime = endIndex*step + digit.GetFirstSample();

    // Copy the samples into the hit.
    std::vector<double> sampleCharge(endIndex-beginIndex);
    for (std::size_t i = beginIndex; i<endIndex; ++i) {
        sampleCharge[i-beginIndex] = digit.GetSample(i);
    }

    // Find the peak charge, average sample index and sample index squared for
    // the peak.  The average is charge weighted.
    double b1 = digit.GetSample(beginIndex+1); 
    double b2 = digit.GetSample(endIndex-1);
    for (int j=beginIndex+1; j<endIndex; ++j) {
        double v = digit.GetSample(j);
        double baseline = b1*(endIndex-j-1) + b2*(j-beginIndex-1);
        baseline /= (endIndex - beginIndex - 2);
        if (!split) v -= baseline;
        if (v <= 0) continue;
        charge += v;
        sample += v*j;
        sampleSquared += v*j*j;
        ++samples;
    }

    // Protect against empty regions.
    if (charge<1) return CP::THandle<CP::THit>();

    sample /= charge;
    sampleSquared /= charge;

    // Convert the sample index into a time.
    double time = (sample + 0.5)*step + digit.GetFirstSample();
    
    // Find the sample RMS, and then convert into a time RMS.
    double rms = step*std::sqrt(sampleSquared - sample*sample + 1.0);

    // The ifdef code controls how the RMS is handled when the pulse is being
    // split.  If the code is in, then the RMS is artificially set to be the
    // half with of the time slice.  Otherwise, the real RMS is used.  The
    // default is to use the real RMS.
#ifdef RESET_SPLIT_RMS
    // If this is from a pulse that is being split, then the RMS is defined by
    // the half width of the split part of the pulse.
    if (split) {
        rms = 0.5*step*(endIndex-beginIndex);
    }
#endif

    // Base the uncertainty in the time on the number of samples used to find
    // the RMS.
    double timeUnc = rms/std::sqrt(samples);

    // The charge uncertainty is calculated assuming the limit of Poisson
    // statistics, but assumes that the background uncertainties are
    // correlated.
    double sig = sampleSigma*std::sqrt(samples);
    double chargeUnc = std::sqrt(charge + sig*sig) + baselineSigma*samples;

    CP::TGeometryId geomId
        = CP::TChannelInfo::Get().GetGeometry(digit.GetChannelId());

    if (!geomId.IsValid()) return CP::THandle<CP::THit>();

    if (!std::isfinite(timeUnc) || timeUnc <= 0.0) {
        CaptError("Time uncertainty for " << digit.GetChannelId()
                  << " is not positive or and finite ");
        timeUnc = 1*unit::second;
    }
    if (!std::isfinite(chargeUnc) || chargeUnc <= 0.0) {
        CaptError("Charge uncertainty for " << digit.GetChannelId()
                  << " is not positive and finite " << samples);
        chargeUnc = 1*unit::coulomb;
    }

    TChannelCalib calib;

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

    // Correct for the electron drift lifetime.
    double deltaT = time - t0;
    if (fCorrectElectronLifetime && deltaT > 0.0) {
        charge *= std::exp(deltaT/calib.GetElectronLifetime());
        chargeUnc *= std::exp(deltaT/calib.GetElectronLifetime());
    }

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

std::pair<int, int> CP::TWireMakeHits::HitExtent(
    int peakIndex,
    const CP::TCalibPulseDigit& deconv,
    double baselineSigma, double sampleSigma) {
    double peak = deconv.GetSample(peakIndex);
    int beginIndex = peakIndex;
    while (beginIndex>0) {
        double v = deconv.GetSample(beginIndex);
        if (v < fIntegrationChargeThreshold*peak) break;
        --beginIndex;
    }
    int endIndex = peakIndex;
    while (endIndex<deconv.GetSampleCount()) {
        double v = deconv.GetSample(endIndex);
        if (v < fIntegrationChargeThreshold*peak) break;
        ++endIndex;
    }
    std::cout << beginIndex << " " << endIndex << std::endl;
    return std::make_pair(beginIndex,endIndex);
}

double CP::TWireMakeHits::operator() (CP::THitSelection& hits,
                                      const CP::TCalibPulseDigit& calib,
                                      const CP::TCalibPulseDigit& deconv,
                                      double t0,
                                      double baselineSigma,
                                      double sampleSigma) {
    double wireCharge = 0.0;

    // Find the time per sample in the digit.
    double digitStep = deconv.GetLastSample()-deconv.GetFirstSample();
    digitStep /= deconv.GetSampleCount();

    // Make sure we have enough memory allocated for the spectrum.
    if (fNSource < (int) deconv.GetSampleCount()) {
        if (fSource) delete[] fSource;
        if (fDest) delete[] fDest;
        if (fWork) delete[] fWork;
        fNSource = 2*deconv.GetSampleCount();
        fSource = new float[fNSource];
        fDest = new float[fNSource];
        fWork = new float[fNSource];
    }

    CP::TChannelCalib channelCalib;
    bool wireIsBipolar = channelCalib.IsBipolarSignal(deconv.GetChannelId());

    double peakMaximumCut = fPeakMaximumCol;
    double peakAreaCut = fPeakAreaCol;
    double peakWidth = fPeakWidthCol;
    if (wireIsBipolar) {
        peakMaximumCut = fPeakMaximumInd;
        peakAreaCut = fPeakAreaInd;
        peakWidth = fPeakWidthInd;
    }
    
    std::ostringstream histName;
    histName << "wire-"
             << std::setw(3)
             << std::setfill('0')
             << CP::TChannelInfo::Get().GetWireNumber(deconv.GetChannelId());

    std::ostringstream histTitle ;
    histTitle << "wire "
             << CP::TChannelInfo::Get().GetWireNumber(deconv.GetChannelId())
             << " (" << deconv.GetChannelId().AsString() << ")";

    // Fill the spectrum and reset the destination.
    double signalOffset = 100000.0;
    for (std::size_t i = 0; i<deconv.GetSampleCount(); ++i) {
        double p = deconv.GetSample(i);
        if (!std::isfinite(p)) {
            CaptError("Channel " << deconv.GetChannelId()
                      << " w/ invalid sample " << i << " " << p);
        }
        fSource[i] = deconv.GetSample(i);
        fDest[i] = 0.0;
        signalOffset = std::max(signalOffset, std::abs(fSource[i])+100.0);
    }

    // Add an artificial baseline to the source so that the values are always
    // positive definite (a TSpectrum requirement).
    for (std::size_t i = 0; i<deconv.GetSampleCount(); ++i) {
        fSource[i] += signalOffset;
    }

    std::auto_ptr<TSpectrum> spectrum(new TSpectrum(fNSource));

    // Set the parameters to use for the peak search.  These should be in the
    // parameters file, but the search isn't extremely sensitive to them.  The
    // most sensitive parameter is the expected peak width.  The width only
    // affects the search since the actual peak width is calculated after the
    // peaks are found.
    double sigma = peakWidth/digitStep;
    double threshold = 1;
    bool removeBkg = false;
    int iterations = 15;
    bool useMarkov = true;
    int window = 3;
    int found = spectrum->SearchHighRes(fSource,fDest,deconv.GetSampleCount(),
                                        sigma, threshold,
                                        removeBkg, iterations,
                                        useMarkov, window);

    // Remove the artificial baseline from the destination.
    for (std::size_t i = 0; i<deconv.GetSampleCount(); ++i) {
        fDest[i] -= signalOffset;
    }

    // Find the new baseline.
    for (std::size_t i = 0; i<deconv.GetSampleCount(); ++i) {
        fWork[i] = fDest[i];
    }
    std::sort(&fWork[0],&fWork[deconv.GetSampleCount()]);
    double baseline = fWork[deconv.GetSampleCount()/2];

    // Find the magnitude of noise for this channel.
    for (std::size_t i = 0; i<deconv.GetSampleCount(); ++i) {
        fWork[i] = std::abs(fDest[i] - baseline);
    }
    std::sort(&fWork[0],&fWork[deconv.GetSampleCount()]);
    int inoise = 0.7*deconv.GetSampleCount();
    // A threshold will be set in terms of standard deviations of the noise.
    // Peaks less than this are rejected as noise.
    double noise = fWork[inoise];

    // Protect against a "zero" channel.
    if (fWork[deconv.GetSampleCount()-1] < 100) {
        CaptLog("Wire with no signal: " << deconv.GetChannelId()
                << " noise: " << noise
                << " max: " << fWork[deconv.GetSampleCount()-1]);
        return wireCharge;
    }
    
    // Don't bother with channels that have crazy big noise.  This says if the
    // noise is bigger than the peak selection threshold, don't save the wire.
    if (noise > peakAreaCut) {
        CaptInfo("Wire with large peak noise: " << deconv.GetChannelId()
                << "  noise: " << noise
                << "  base: " << baseline );
        // return wireCharge;
    }

#define STANDARD_HISTOGRAM
#ifdef STANDARD_HISTOGRAM
#undef STANDARD_HISTOGRAM
    static TH1F* gNoiseHistogram = NULL;
    if (!gNoiseHistogram) {
        gNoiseHistogram
            = new TH1F("peakSearchNoise",
                       "Peak area sigma for all wires",
                       100, 0.0, 5000);
    }
    gNoiseHistogram->Fill(noise);
#endif

    // Get the peak positions, and determine which ones are above threshold.
    float* xx = spectrum->GetPositionX();
    std::vector<float> peaks;

#define STANDARD_HISTOGRAM
#ifdef STANDARD_HISTOGRAM
#undef STANDARD_HISTOGRAM
    static TH1F* gPeakHistogram = NULL;
    if (!gPeakHistogram) {
        gPeakHistogram
            = new TH1F("peakArea",
                       "Estimated area of found peaks on each wires",
                       100, 0.0, 50000);
    }
    for (int i=0; i<std::min(10,found); ++i) {
        int index = (int) (xx[i] + 0.5);
        double peak = fDest[index];
        gPeakHistogram->Fill(std::min(peak,49000.0));
    }
#endif

#define STANDARD_HISTOGRAM
#ifdef STANDARD_HISTOGRAM
#undef STANDARD_HISTOGRAM
    static TH1F* gHeightHistogram = NULL;
    if (!gHeightHistogram) {
        gHeightHistogram
            = new TH1F("peakHeight",
                       "Height of found peaks on each wires",
                       100, 0.0, 50000);
    }
    for (int i=0; i<std::min(10,found); ++i) {
        int index = (int) (xx[i] + 0.5);
        double peak = deconv.GetSample(index);
        gHeightHistogram->Fill(std::min(peak,49000.0));
    }
#endif

    // Now use fWork to mask out peaks that are part of another peak.
    std::fill(&fWork[0],&fWork[deconv.GetSampleCount()], 0.0);
    
    double noiseThreshold
        = std::sqrt(baselineSigma*baselineSigma + sampleSigma*sampleSigma);
    noiseThreshold *= fIntegrationNoiseThreshold;

    for (int i=0; i<found; ++i)  {
        // Check the peak size and deconvolution power.
        int index = (int) (xx[i] + 0.5);
        // No peaks at the ends of the digit.
        if (index < fDigitEndSkip) continue;
        if (index > (int) deconv.GetSampleCount()-fDigitEndSkip - 1) continue;
        // Apply a cut to the overall peak size. (The digit has had the
        // baseline remove and is in units of charge).
        if (deconv.GetSample(index) < peakMaximumCut) continue;
        // Apply a cut to the minimum value of the found peak.
        if (fDest[index] < peakAreaCut) continue;
        // Apply a cut based on the noise in the peak search.
        if (fDest[index] < noise*fNoiseThresholdCut + baseline) continue;
        // Skip peaks that close to other peaks.
        if (fWork[index] > 0.5) continue;
        std::pair<int, int> extent = HitExtent(index,
                                               deconv,
                                               baselineSigma,
                                               sampleSigma);
        for (int j=extent.first; j<extent.second; ++j) {
            fWork[j] = 1.0;
        }
        peaks.push_back(xx[i]);
    }

    // Sort the peaks by time.
    std::sort(peaks.begin(), peaks.end());

    // This ugly bit of code is taking a found peak and making sure that if
    // it's too wide (i.e. the RMS is too big), it's split into smaller hits.
    // The heuristic is that a peak contains all charge for bins that are more
    // than a charge and noise threshold.  If more than one peak is
    // found, then the digit is split at the halfway point between the peaks.
    // If the peak is too wide, then the peak is split.
    int hitCount = 0;
    for (std::vector<float>::iterator p = peaks.begin();
         p != peaks.end(); ++p) {
        int iPeak = (int) (*p + 0.5);
        int beginIndex = 0;
        int endIndex = deconv.GetSampleCount();
        double peak = deconv.GetSample(iPeak);

        // For the current peak, find the upper and lower bounds of the
        // integration region.  If the peak is close to another, the bound is
        // halfway between the two peaks.  After this is finished, the
        // "beginIndex" and "endIndex" will be the indices to be looking at.
        if (peaks.size() > 1) {
            for (std::vector<float>::iterator o = peaks.begin();
                 o != peaks.end(); ++o) {
                if (o == p) continue;
                if (*o < *p) {
                    int t = (int) ((*p+*o)/2);
                    if (beginIndex<t) beginIndex = t;
                }
                else {
                    int t = (int) ((*p+*o)/2);
                    if (endIndex>t) endIndex = t;
                }
            }
        }

        std::pair<int, int> extent = HitExtent(iPeak,
                                               deconv,
                                               baselineSigma,
                                               sampleSigma);

        if (beginIndex < extent.first) beginIndex = extent.first;
        
        if (extent.second < endIndex) endIndex = extent.second;
        
        double charge = 0.0;
        double time = 0.0;
        double timeSquared = 0.0;
        double samples = 1.0;

        // Find the raw RMS for the peak.
        for (int j=beginIndex; j<endIndex; ++j) {
            double v = deconv.GetSample(j);
            if (v<0.1) v = 0.1;
            charge += v;
            time += v*j;
            timeSquared += v*j*j;
            ++samples;
        }
        time /= charge;
        timeSquared /= charge;

        double rawRMS = (timeSquared - time*time);
        rawRMS = digitStep*std::sqrt(rawRMS+1.0);

        int split = (int) (rawRMS/fPeakRMSLimit + 1.0);
        double step = 1.0*(endIndex-beginIndex)/split;

        ++hitCount;
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        std::ostringstream hitName;
        hitName << histName.str()
                << "-hit"
                << std::setw(3)
                << std::setfill('0')
                << hitCount;
        TH1F* hitHist =
            new TH1F(hitName.str().c_str(),
                     ("Hit for " + hitName.str()).c_str(),
                     endIndex - beginIndex,
                     beginIndex*digitStep, endIndex*digitStep);
        for (int baseIndex = beginIndex;
             baseIndex<endIndex; ++baseIndex) {
            hitHist->SetBinContent(baseIndex-beginIndex+1,
                                   deconv.GetSample(baseIndex));
        }
#endif

        /// The loop looks a odd since it's being done with doubles,
        /// but this means that the "digitized" split size gets distributed
        /// over the entire range instead of all at one end.
        for (double baseIndex = beginIndex;
             baseIndex<endIndex;
             baseIndex += step) {
            int i = (int) baseIndex;
            int j = (int) (baseIndex + step);
            CP::THandle<CP::THit> newHit = MakeHit(deconv, digitStep, t0,
                                                   baselineSigma,
                                                   sampleSigma,
                                                   i, j,
                                                   (split > 1));
            if (!newHit) continue;
            wireCharge += newHit->GetCharge();
            hits.push_back(newHit);
        }
    }

    if (hitCount < 1) return wireCharge;
    
#define FILL_HISTOGRAM
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* calibHist
        = new TH1F((histName.str()+"-calib").c_str(),
                   ("Calibration for "
                    + histTitle.str()).c_str(),
                   calib.GetSampleCount(),
                   calib.GetFirstSample(), calib.GetLastSample());
    for (std::size_t i = 0; i<calib.GetSampleCount(); ++i) {
        calibHist->SetBinContent(i+1,calib.GetSample(i));
    }
#endif

#define FILL_HISTOGRAM
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* sourceHist
        = new TH1F((histName.str()+"-deconv").c_str(),
                   ("Deconvolution for "
                    + histTitle.str()).c_str(),
                   deconv.GetSampleCount(),
                   deconv.GetFirstSample(), deconv.GetLastSample());
    for (std::size_t i = 0; i<deconv.GetSampleCount(); ++i) {
        sourceHist->SetBinContent(i+1,deconv.GetSample(i));
    }
#endif

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* destHist
        = new TH1F((histName.str()+"-tspectrum").c_str(),
                   ("TSpectrum result for "
                    + histTitle.str()).c_str(),
                   deconv.GetSampleCount(),
                   deconv.GetFirstSample(), deconv.GetLastSample());
    for (std::size_t i = 0; i<deconv.GetSampleCount(); ++i) {
        destHist->SetBinContent(i+1,fDest[i]);
    }
#endif

    return wireCharge;
}
