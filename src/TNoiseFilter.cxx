#include "TNoiseFilter.hxx"

#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"

#include <TChannelId.hxx>
#include <TCaptLog.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>

#include <algorithm>
#include <iostream>

CP::TNoiseFilter::~TNoiseFilter() {}

CP::TNoiseFilter::TNoiseFilter(double noisePower, double spikePower)
    : fNoisePower(noisePower), fSpikePower(spikePower), fIsNoisy(false) {}

void CP::TNoiseFilter::Calculate(CP::TChannelId id,
                                 CP::TElectronicsResponse& elecFreq,
                                 CP::TWireResponse& wireFreq,
                                 TVirtualFFT& measFreq) {

    fIsNoisy = false;

    if (fFilter.size() != elecFreq.GetSize()) {
        fFilter.resize(elecFreq.GetSize());
        fWork.resize(fFilter.size());
        fAverage.resize(fFilter.size());
        fSigma.resize(fFilter.size());
    }

    // Copy the power at every frequency in to a work area.
    double maxResponse = 0.0;
    double minResponse = std::abs(elecFreq.GetFrequency(fFilter.size()/2)
                                  *wireFreq.GetFrequency(fFilter.size()/2));
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        double rl, im;
        measFreq.GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        fWork[i] = std::abs(c);
        fFilter[i] = 1.0;
        double res = std::abs(elecFreq.GetFrequency(i)
                              *wireFreq.GetFrequency(i));
        maxResponse = std::max(maxResponse,res);
    }

    // Make a "smoothed" estimate of the power at every bin and the sigma.
    int halfWindow = 50;
    int excludeWindow = 10;
    double maxPower = 0.0;
    double minPower = 1E+20;
    for (int i = 0; i<(int) fFilter.size(); ++i) {
        // Make a first estimate of the local average.
        double rawAverage = 0.0;
        double rawAverage2 = 0.0;
        double maxBin = 0.0;
        for (int j = i - halfWindow; j < i+halfWindow+1; ++j) {
            int index = j % fFilter.size();
            if (index < 0) index += fFilter.size();
            rawAverage += fWork[index];
            rawAverage2 += fWork[index]*fWork[index];
            if (fWork[index]>maxBin) maxBin = fWork[index];
        }
        rawAverage = (rawAverage - maxBin)/(2.0*halfWindow);
        rawAverage2 = (rawAverage2 - maxBin*maxBin)/(2.0*halfWindow);
        double rawSigma = 1.0;
        if (rawAverage2 > rawAverage*rawAverage) {
            rawSigma = std::sqrt(rawAverage2 - rawAverage*rawAverage);
        }
        double truncAverage = 0.0;
        double truncSigma = 0.0;
        double truncWeight = 0.0;
        for (int j = i - halfWindow; j < i+halfWindow+1; ++j) {
            if (std::abs(j-i) < excludeWindow) continue;
            int index = j % fFilter.size();
            if (index < 0) index += fFilter.size();
            if (fWork[index] > rawAverage + 20.0*rawSigma) continue;
            truncAverage += fWork[index];
            truncSigma += fWork[index]*fWork[index];
            truncWeight += 1.0;
        }
        truncAverage /= truncWeight;
        truncSigma /= truncWeight;
        fAverage[i] = truncAverage;
        fSigma[i] = 1.0;
        if (truncSigma> truncAverage*truncAverage) {
            fSigma[i] = std::sqrt(truncSigma - truncAverage*truncAverage);
        }
        minPower = std::min(minPower, fAverage[i]);
        maxPower = std::max(maxPower, fAverage[i]);
    }

#ifdef REJECT_NOISY_CHANNELS
    // Make an estimate of the channel noise.  This should reject anything
    // where the noise is much beyond 20 ADC counts.
    double a = minPower/fFilter.size();
    double b = 2.0*a*a;
    double c = b*fFilter.size();
    double d = sqrt(c);
    if (d > 400.0) {
        fIsNoisy = true;
    }
#endif
    
    // Find the Gaussian noise estimate.
    double noise = fNoisePower;
    noise *= (maxPower*minResponse-minPower*maxResponse)/(minPower-maxPower);

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* sigHist 
        = new TH1F((id.AsString()+"-signal").c_str(),
                   ("Signal Estimate for " 
                    + id.AsString()).c_str(),
                   fFilter.size(),
                   0,fFilter.size());
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        sigHist->SetBinContent(i+1,fAverage[i]); //  + fSpikePower*fSigma[i]);
    }
#endif

    int spikeCount = 0;
    for (std::size_t i = 1; i<fFilter.size(); ++i) {
        double measPower = fWork[i];
        double respPower = std::abs(elecFreq.GetFrequency(i)
                                    *wireFreq.GetFrequency(i));
        double sigPower = fAverage[i];
        // Apply filter for any spikes in the FFT.
        if (fWork[i] > fAverage[i]+fSpikePower*fSigma[i]) {
            ++spikeCount;
            fFilter[i] *= sigPower*sigPower
                /(sigPower*sigPower+measPower*measPower);
        }
        // Apply filter for Gaussian Noise.
        fFilter[i] *= respPower*respPower/(respPower*respPower + noise*noise);
    }

}


