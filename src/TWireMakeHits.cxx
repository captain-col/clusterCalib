#include "TWireMakeHits.hxx"
#include "TChannelCalib.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TDataHit.hxx>
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

CP::TWireMakeHits::TWireMakeHits(bool correctLifetime, 
                                 bool correctEfficiency) {
    fCorrectElectronLifetime = correctLifetime;
    fCorrectCollectionEfficiency = correctEfficiency;
    fNSource = 0;
    fSource = NULL;
    fDest = NULL;
    fWork = NULL;

    fPeakMaximumCut 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.charge");

    fPeakDeconvolutionCut 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.deconvolution");

    fNoiseThresholdCut 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.noise");

    fPeakRMSLimit 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.rmsLimit");

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
    
    // Find the peak charge, average sample index and sample index squared for
    // the peak.  The average is charge weighted.
    for (int j=beginIndex; j<endIndex; ++j) {
        double v = digit.GetSample(j);
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
                  << " is not positive and finite ");
        timeUnc = 1*unit::second;
    }
    if (!std::isfinite(chargeUnc) || chargeUnc <= 0.0) {
        CaptError("Charge uncertainty for " << digit.GetChannelId() 
                  << " is not positive and finite" << samples);
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
    }

    // Correct for the electron drift lifetime.  
    double deltaT = time - t0;
    if (fCorrectElectronLifetime && deltaT > 0.0) {
        charge *= std::exp(deltaT/calib.GetElectronLifetime());
    }

    // Build the hit.
    CP::TWritableDataHit hit;
    hit.SetGeomId(geomId);
    hit.SetDigit(digit.GetParent());
    hit.SetCharge(charge);
    hit.SetChargeUncertainty(chargeUnc);
    hit.SetTime(time);
    hit.SetTimeRMS(rms);
    hit.SetTimeUncertainty(timeUnc);

    return CP::THandle<CP::TDataHit>(new CP::TDataHit(hit));
}

void CP::TWireMakeHits::operator() (CP::THitSelection& hits, 
                                    const CP::TCalibPulseDigit& digit,
                                    double t0,
                                    double baselineSigma,
                                    double sampleSigma) {

    // Make sure we have enough memory allocated for the spectrum.
    if (fNSource < (int) digit.GetSampleCount()) {
        if (fSource) delete[] fSource;
        if (fDest) delete[] fDest;
        if (fWork) delete[] fWork;
        fNSource = 2*digit.GetSampleCount();
        fSource = new float[fNSource];
        fDest = new float[fNSource];
        fWork = new float[fNSource];
    }

    // Fill the spectrum and reset the destination.
    double signalOffset = 100000.0;
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        double p = digit.GetSample(i);
        if (!std::isfinite(p)) {
            CaptError("Channel " << digit.GetChannelId() 
                      << " w/ invalid sample " << i);
        }
        fSource[i] = digit.GetSample(i);
        fDest[i] = 0.0;
        signalOffset = std::max(signalOffset, std::abs(fSource[i])+100.0);
    }

    // Add an artificial baseline to the source so that the values are always
    // positive definite (a TSpectrum requirement).
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        fSource[i] += signalOffset;
    }
    
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* sourceHist 
        = new TH1F((digit.GetChannelId().AsString()+"-source").c_str(),
                   ("Source for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   digit.GetFirstSample(), digit.GetLastSample());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        sourceHist->SetBinContent(i+1,fSource[i]);
    }
#endif
        
    std::auto_ptr<TSpectrum> spectrum(new TSpectrum(fNSource));

    // Set the parameters to use for the peak search.  These should be in the
    // parameters file, but the search isn't extremely sensitive to them.  The
    // most sensitive parameter is the expected peak width.  The width only
    // affects the search since the actual peak width is calculated after the
    // peaks are found.
    double sigma = 2.0;
    double threshold = 1;
    bool removeBkg = false;
    int iterations = 15;
    bool useMarkov = true;
    int window = 2;
    int found = spectrum->SearchHighRes(fSource,fDest,digit.GetSampleCount(),
                                        sigma, threshold,
                                        removeBkg, iterations, 
                                        useMarkov, window);

    // Remove the artificial baseline from the destination.
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        fDest[i] -= signalOffset;
    }

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* destHist 
        = new TH1F((digit.GetChannelId().AsString()+"-peaks").c_str(),
                   ("Peaks found for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   digit.GetFirstSample(), digit.GetLastSample());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        destHist->SetBinContent(i+1,fDest[i]);
    }
#endif

    // Find the new baseline.
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        fWork[i] = fDest[i];
    }
    std::sort(&fWork[0],&fWork[digit.GetSampleCount()]);
    double baseline = fWork[digit.GetSampleCount()/2];

    // Find the magnitude of noise for this channel.
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        fWork[i] = std::abs(fDest[i] - baseline);
    }
    std::sort(&fWork[0],&fWork[digit.GetSampleCount()]);
    int inoise = 0.7*digit.GetSampleCount();
    // A threshold will be set in terms of standard deviations of the noise.
    // Peaks less than this are rejected as noise.
    double noise = fWork[inoise];

#define FILL_HISTOGRAM
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        static TH1F* gNoiseHistogram = NULL;
        if (!gNoiseHistogram) {
            gNoiseHistogram
                = new TH1F("peakSearchNoise",
                           "Peak height sigma for all wires",
                           100, 0.0, 5000);
        }
        gNoiseHistogram->Fill(noise);
#endif

        // Get the peak positions, and determine which ones are above threshold.
    float* xx = spectrum->GetPositionX();
    std::vector<float> peaks;

    int digitEndSkip = 500;
    for (int i=0; i<found; ++i)  {
        // No peaks at the ends of the digit.
        if (xx[i] < digitEndSkip) continue;
        if (xx[i] > digit.GetSampleCount()-digitEndSkip - 1) continue;
        // Check the peak size and deconvolution power.
        int index = (int) (xx[i] + 0.5);
        // Apply a cut to the overall peak size. (The digit has had the
        // baseline remove and is in units of charge.
        if (digit.GetSample(index) < fPeakMaximumCut) continue;
        // Apply a cut to the minimum value of the found peak.
        if (fDest[index] < fPeakDeconvolutionCut) continue;
        // Apply a cut based on the noise in the peak search.
        if (fDest[index] < noise*fNoiseThresholdCut + baseline) continue;
        peaks.push_back(xx[i]);
    }

    // Sort the peaks by time.
    std::sort(peaks.begin(), peaks.end());

    // This ugly bit of code is taking a found peak and making sure that if
    // it's too wide (i.e. the RMS is too big), it's split into smaller hits.
    // The heuristic is that a peak contains all charge for bins that are more
    // than "chargeThresh" times the peak value.  If more than one peak is
    // found, then the digit is split at the halfway point between the peaks.
    // If the peak is too wide, then the peak is split.
    double chargeThresh = 0.1;
    int hitCount = 0;
    for (std::vector<float>::iterator p = peaks.begin();
         p != peaks.end(); ++p) {
        int i = (int) (*p + 0.5);
        int beginIndex = 0;
        int endIndex = digit.GetSampleCount();
        double peak = digit.GetSample(i);

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

        // Search down from the peak center to find the actual integration
        // bounds based on the actual extent of the peak.  This extends beyond
        // the edge of the peak to pick up some of the background.
        int j = beginIndex;
        int endCount = -1;
        for (beginIndex = i-1; j <= beginIndex; --beginIndex) {
            double v = digit.GetSample(beginIndex);
            if (endCount > 0 && --endCount == 0) break;
            if (endCount < 0 && v < chargeThresh*peak) endCount = 10;
        }

        // Search up from the peak center to find the actual integration
        // bounds based on the actual extent of the peak.  This extends beyond
        // the edge of the peak to pick up some of the background.
        j = endIndex;
        endCount = -1;
        for (endIndex = i+1; endIndex < j; ++endIndex) {
            double v = digit.GetSample(endIndex);
            if (endCount > 0 && --endCount == 0) break;
            if (endCount < 0 && v < chargeThresh*peak) endCount = 10;
        }

        double charge = 0.0;
        double time = 0.0;
        double timeSquared = 0.0;
        double samples = 1.0;

        // Find the raw RMS for the peak.
        for (int j=beginIndex; j<endIndex; ++j) {
            double v = digit.GetSample(j);
            if (v<0.1) v = 0.1;
            charge += v;
            time += v*j;
            timeSquared += v*j*j;
            ++samples;
        }
        time /= charge;
        timeSquared /= charge;

        // Find the time per sample in the digit.
        double digitStep = digit.GetLastSample()-digit.GetFirstSample();
        digitStep /= digit.GetSampleCount();

        double rawRMS = (timeSquared - time*time);
        rawRMS = digitStep*std::sqrt(rawRMS+1.0);

        int split = (int) (rawRMS/fPeakRMSLimit + 1.0);
        double step = 1.0*(endIndex-beginIndex)/split;

        /// The loop looks a little odd since it's being done with doubles,
        /// but this means that the "digitized" split size gets distributed
        /// over the entire range instead of all at one end.
        for (double baseIndex = beginIndex; 
             baseIndex<endIndex; 
             baseIndex += step) {
            int i = (int) baseIndex;
            int j = (int) (baseIndex + step);
            CP::THandle<CP::THit> newHit = MakeHit(digit, digitStep, t0,
                                                   baselineSigma,
                                                   sampleSigma,
                                                   i, j,
                                                   (split > 1));
            if (!newHit) continue;
            hits.push_back(newHit);
            if (++hitCount > 50) return;
        }            
    }
}
