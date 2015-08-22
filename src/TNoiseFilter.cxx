#include "TNoiseFilter.hxx"

#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"
#include "TChannelCalib.hxx"

#include <TChannelId.hxx>
#include <TCaptLog.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>
#include <TSpectrum.h>

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
    }
    std::fill(fAverage.begin(), fAverage.end(), 0.0);
    std::fill(fFilter.begin(), fFilter.end(), 1.0);

    // Copy the power at every frequency in to a work area.
    double maxResponse = 0.0;
    double minResponse = 1E+6;
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        double rl, im;
        measFreq.GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        fWork[i] = std::abs(c);
        double res = std::abs(elecFreq.GetFrequency(i));
        maxResponse = std::max(maxResponse,res);
        minResponse = std::min(minResponse,res);
    }

    // Make a "smoothed" estimate of the power at every frequency.  This is
    // based on the hypothesis that we have a smooth "signal" frequency (we
    // do), and that most of the noise comes at discreet frequencies.  This
    // only happens if the spike filter is enabled by specifying the
    // clusterCalib.deconvolution.spikePower in the config file to be
    // positive.
    if (fSpikePower > 1E-4) {
        // This is fairly time consuming, so don't calculate it if we aren't
        // removing the spike power.  This uses the TSpectrum implementation
        // of the Sensitive Nonlinear Iterative Peak clipping algorithm [NIM
        // B34 (1988) 396--402].
        TSpectrum spectrum;
        std::copy(fWork.begin(), fWork.end(), fAverage.begin());
        spectrum.Background(fAverage.data(),fAverage.size(), 50,
                            TSpectrum::kBackIncreasingWindow,
                            TSpectrum::kBackOrder2, // Less interpolation
                            true,
                            TSpectrum::kBackSmoothing15, // Lots of smoothing
                            false);
    }

    // Find the Gaussian noise estimate.  This uses the fact that Gaussian
    // noise has a flat frequency spectrum while the signal is not.  By using
    // a game of "averages", this solves based on the hypothesis that the
    // measured power is the sum of the power from the signal and a flat
    // noise.  The signal power spectrum is estimated using the electronics
    // response power spectrum.  It takes the total power in the measured
    // spectrum and averages it over all bins, then it averages each frequency
    // bin using the power that should be in the signal.  The noise is found
    // by looking at the excess power in the measured spectrum.
    double maxSignal = 0.0;
    double maxSignalWeight = 0.0;
    double minSignal = 0.0;
    double minSignalWeight = 0.0;
    for (std::size_t i=0; i<fWork.size(); ++i) {
        double res = std::abs(elecFreq.GetFrequency(i));
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
    maxSignal = std::sqrt(maxSignal*fWork.size()/2.0);
    minSignal = std::sqrt(minSignal*fWork.size()/2.0);

    // Estimate the noise power relative to the averaged signal.
    double noise = 0.0;
    if (minSignal < maxSignal) {
        noise = maxSignal*minResponse-maxResponse*minSignal;
        noise /= (minSignal-maxSignal);
    }

    // Protect against numeric problems.  The equations used above can give a
    // negative solution, so don't let that happen.
    if (noise < 0.0) noise = 0.0;

    // Add in the "noise" introduced by digitization.  This limits the RMS to
    // be greater than or equal to 2 ADC counts.  It's then normalized to be a
    // fraction of the expected MIP deposited charge.
    double gain = channelCalib.GetGainConstant(id,1);
    double slope = channelCalib.GetDigitizerConstant(id,1);
    double adcNoise = 2.0/gain/slope;

    /// Combine the different noise components.  The 4 fC comes the nominal
    /// MIP charge on a wire averaged over all track geometries.  Should be a
    /// parameter.  In truth, I didn't do the average, I just took a SWAG and
    /// it seems to be "good enough".
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
                   0,1.0);
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        sigHist->SetBinContent(i+1,fAverage[i]);
   }
#endif

    /// The form of filter being used below is based on an optimal filter
    /// (e.g. a Weiner Filter) where the filter value is F = S^2/(S^2 + N^2)
    /// (S being the true signal, and N being the true noise).  F will be
    /// between 0.0 and 1.0 (0.0 is more filtering, 1.0 is no filtering).
    double maxFilter = 0.0;
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        double measPower = fWork[i];
        double respPower = std::abs(elecFreq.GetFrequency(i));
        // Apply filter for any spikes in the FFT.
        if (fSpikePower > 1E-3) {
            double sigPower = fAverage[i]/fSpikePower;
            double excessPower = std::max(0.0, measPower-sigPower);
            fFilter[i] *= sigPower*sigPower
                /(sigPower*sigPower+excessPower*excessPower);            
        }
        // Apply filter for Gaussian Noise.
        fFilter[i] *= respPower*respPower/(respPower*respPower + noise*noise);
        if (!std::isfinite(fFilter[i])) {
            CaptError("Filter not finite at " << i << " " << fFilter[i]);
            fFilter[i] = 1.0;
        }
        else {
            maxFilter = std::max(maxFilter, fFilter[i]);
        }
    }

    /// Normalize the filter so that it has a maximum value of 1.0 (no
    /// filtering).  This is because we don't know the real power of the
    /// signal and noise and this helps maintain the overall normalization of
    /// the charge in the peak.
    double filterNorm = maxFilter;
    for (std::size_t i = 0; i<fFilter.size(); ++i) {
        fFilter[i] /= filterNorm;
    }

}
