#include "TNoiseFilter.hxx"

#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"
#include "TChannelCalib.hxx"

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

    TChannelCalib channelCalib;
    
    fIsNoisy = false;

    if (fFilter.size() != elecFreq.GetSize()) {
        fFilter.resize(elecFreq.GetSize());
        fWork.resize(fFilter.size());
        fAverage.resize(fFilter.size());
        fSigma.resize(fFilter.size());
    }

    // Copy the power at every frequency in to a work area.
    double maxResponse = 0.0;
    double minResponse = 1E+6;
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        double rl, im;
        measFreq.GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        fWork[i] = std::abs(c);
        fFilter[i] = 1.0;
        double res = std::abs(elecFreq.GetFrequency(i)
                              *wireFreq.GetFrequency(i));
        maxResponse = std::max(maxResponse,res);
        minResponse = std::min(minResponse,res);
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
        if (truncWeight > 0.0) {
            truncAverage /= truncWeight;
            truncSigma /= truncWeight;
            fAverage[i] = truncAverage;
            fSigma[i] = 1.0;
            if (truncSigma > truncAverage*truncAverage) {
                fSigma[i] = std::sqrt(truncSigma - truncAverage*truncAverage);
            }
        }
        else {
            fAverage[i] = 0.0;
            fSigma[i] = 1.0;
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
        CaptError("Channel " << id << " rejected as too noisy");
        fIsNoisy = true;
    }
#endif
    
    // Find the Gaussian noise estimate.
    double maxSignal = 0.0;
    double maxSignalWeight = 0.0;
    double minSignal = 0.0;
    double minSignalWeight = 0.0;
    for (int i=0; i<fWork.size(); ++i) {
        double res = std::abs(elecFreq.GetFrequency(i)
                              *wireFreq.GetFrequency(i));
        double w = res/maxResponse;
        w = w*w;
        double p = fWork[i]*fWork[i]/fWork.size()/fWork.size();
        maxSignalWeight += w;
        maxSignal += w*p;
        w = 1.0-res/maxResponse;
        w *= w;
        minSignalWeight += w;
        minSignal += w*p;
    }
    maxSignal = maxSignal/maxSignalWeight;
    minSignal = minSignal/minSignalWeight;
    minSignal = std::sqrt(minSignal*fWork.size()/2.0);
    
    // Estimate the noise power relative to the averaged signal.
    double noise = 0.0;
    if (minSignal < maxSignal) {
        noise = maxSignal*minResponse-maxResponse*minSignal;
        noise /= (minSignal-maxSignal);
    }
        
    // Protect against numeric problems.
    if (noise < 0.0) noise = 0.0;

    // Add in an offset to deal with the digitization.  This limits the RMS to
    // be greater than or equal to 2 ADC counts.  It's then normalized to be a
    // fraction of the expected MIP deposited charge.
    double gain = channelCalib.GetGainConstant(id,1);
    double slope = channelCalib.GetDigitizerConstant(id,1);
    double adcNoise = 2.0/gain/slope;
    adcNoise = adcNoise/(4.0*unit::fC);
    double signalNoise = minSignal/(4.0*unit::fC);
    noise = std::sqrt(noise*noise
                      + adcNoise*adcNoise
                      + minResponse*minResponse
                      + signalNoise*signalNoise);
    noise *= fNoisePower;
    
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* sigHist 
        = new TH1F((id.AsString()+"-signal").c_str(),
                   ("Signal Estimate for " 
                    + id.AsString()).c_str(),
                   fFilter.size(),
                   0,0.5/channelCalib(id));
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        sigHist->SetBinContent(i+1,fAverage[i]);
   }
#endif

    int spikeCount = 0;
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        double measPower = fWork[i];
        double respPower = std::abs(elecFreq.GetFrequency(i)
                                    *wireFreq.GetFrequency(i));
        double sigPower = fAverage[i];
#ifdef SPIKE_FILTER
        // Apply filter for any spikes in the FFT.
        if (fWork[i] > fAverage[i]+fSpikePower*fSigma[i]) {
            ++spikeCount;
            fFilter[i] *= sigPower*sigPower
                /(sigPower*sigPower+measPower*measPower);
        }
#endif
        // Apply filter for Gaussian Noise.
        fFilter[i] *= respPower*respPower/(respPower*respPower + noise*noise);
        if (!std::isfinite(fFilter[i])) {
            CaptError("Filter not finite at " << i << " " << fFilter[i]);
            fFilter[i] = 1.0;
        }
    }
}

