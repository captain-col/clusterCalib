#include "TWireMakeHits.hxx"

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

CP::TWireMakeHits::TWireMakeHits() {
    fSource = NULL;
    fDest = NULL;
    fNSource = 0;

    fPeakMaximumCut 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.charge");

    fPeakDeconvolutionCut 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.deconvolution");

}
CP::TWireMakeHits::~TWireMakeHits() {
    if (fSource) delete fSource;
    if (fDest) delete fDest;
}

void CP::TWireMakeHits::operator() (CP::THitSelection& hits, 
                                    const CP::TCalibPulseDigit& digit) {

    CP::TGeometryId geomId 
        = CP::TChannelInfo::Get().GetGeometry(digit.GetChannelId());

    // Make sure we have enough memory allocated for the spectrum.
    if (fNSource < (int) digit.GetSampleCount()) {
        if (fSource) delete fSource;
        fNSource = 2*digit.GetSampleCount();
        fSource = new float[fNSource];
        fDest = new float[fNSource];

    }
    
    // Fill the spectrum.
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        fSource[i] = digit.GetSample(i) + 1000.0;
        fDest[i] = 0.0;
    }

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* sourceHist 
        = new TH1F((digit.GetChannelId().AsString()+"-source").c_str(),
                   ("Source for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   0.0, 1.0*digit.GetSampleCount());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        sourceHist->SetBinContent(i+1,fSource[i]);
    }
#endif
        
    TSpectrum* spectrum = new TSpectrum(1000);

    // Find the time per sample in the digit.
    double digitStep = digit.GetLastSample()-digit.GetFirstSample();
    digitStep /= digit.GetSampleCount();

    // Set the parameters to use for the peak search.  These should be in the
    // parameters file, but the search isn't extremely sensitive to them.  The
    // most sensitive parameter is the expected peak width.  The width only
    // affects the search since the actual peak width is calculated after the
    // peaks are found.
    double sigma = 2.0;
    double threshold = 1;
    bool removeBkg = true;
    int iterations = 15;
    bool useMarkov = true;
    int window = 3;
    int found = spectrum->SearchHighRes(fSource,fDest,digit.GetSampleCount(),
                                        sigma, threshold,
                                        removeBkg, iterations, 
                                        useMarkov, window);

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* destHist 
            = new TH1F((digit.GetChannelId().AsString()+"-dest").c_str(),
                       ("Peaks found for " 
                        + digit.GetChannelId().AsString()).c_str(),
                       digit.GetSampleCount(),
                       0.0, 1.0*digit.GetSampleCount());
        for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
            destHist->SetBinContent(i+1,fDest[i]);
        }
#endif
        
    // Get the peak positions and sort by time.
    float* xx = spectrum->GetPositionX();
    std::vector<float> peaks;
    for (int i=0; i<found; ++i)  peaks.push_back(xx[i]);

    // This ugly bit of code is calculating the amount of charge in each peak
    // and the charge uncertainty; calculating the peak time RMS and time
    // uncertainty; and using that to create a hit.  The heuristic is that a
    // peak contains all charge for bins that are more than "chargeThresh"
    // times the peak value.  If there is more than one peak found, then the
    // chanrge is split at the halfway point between the peaks. 
    //
    // THIS IS A QUICK AND DIRTY HACK TO GET SOMETHING "OUT THE DOOR". 
    //
    // UGLY CODE ALERT.  This does not make me proud.
    //
    double chargeThresh = 0.1;
    for (std::vector<float>::iterator p = peaks.begin();
         p != peaks.end(); ++p) {
        // No peaks at the end of a digit.
        if (*p < 3) continue;
        if (*p > digit.GetSampleCount()-3) continue;
        double center = *p;
        int i = *p + 0.5;
        int begin = 0;
        int end = digit.GetSampleCount();
        double peak = digit.GetSample(i);
        double deconv = fDest[i];

        // Only look at peaks above the thresholds.
        if (peak < fPeakMaximumCut) continue;
        if (deconv < fPeakDeconvolutionCut) continue;

        if (peaks.size() > 1) {
            for (std::vector<float>::iterator o = peaks.begin(); 
                 o != peaks.end(); ++o) {
                if (o == p) continue;
                if (*o < *p) {
                    int t = (*p+*o)/2;
                    if (begin<t) begin = t;
                }
                else {
                    int t = (*p+*o)/2;
                    if (end>t) end = t;
                }
            }
        }
        double charge = digit.GetSample(i);
        double time = i*digit.GetSample(i);
        double timeSquared = i*i*digit.GetSample(i);
        double samples = 1.0;
        // Sum on the low side of the peak.
        for (int j = i-1; begin <= j; --j) {
            double v = digit.GetSample(j);
            if (v < chargeThresh*peak) break;
            charge += v;
            time += v*j;
            timeSquared += v*j*j;
            ++samples;
        }
        for (int j = i+1; j < end; ++j) {
            double v = digit.GetSample(j);
            if (v < chargeThresh*peak) break;
            charge += v;
            time += v*j;
            timeSquared += v*j*j;
            ++samples;
        }
        double chargeUnc = std::sqrt(charge);
        center *= digitStep;

        time /= charge;
        timeSquared /= charge;

        double rms = (timeSquared - time*time);
        time = time*digitStep;
        rms = digitStep*std::sqrt(rms + 1);
        double timeUnc = rms/std::sqrt(samples) + std::abs(center - time);
        rms = rms + std::abs(center-time);

        center += digit.GetFirstSample() + 0.5*digitStep;

        CP::TWritableDataHit hit;
        hit.SetGeomId(geomId);
        hit.SetDigit(digit.GetParent());
        hit.SetCharge(charge);
        hit.SetChargeUncertainty(chargeUnc);
        hit.SetTime(center);
        hit.SetTimeRMS(rms);
        hit.SetTimeUncertainty(timeUnc);
        hits.push_back(CP::THandle<CP::TDataHit>(new CP::TDataHit(hit)));
    }
}
