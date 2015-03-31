#include "TPulseDeconvolution.hxx"
#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"
#include "TNoiseFilter.hxx"
#include "TChannelCalib.hxx"

#include <TCalibPulseDigit.hxx>
#include <TEvent.hxx>
#include <TEventFolder.hxx>
#include <TRuntimeParameters.hxx>
#include <TMCChannelId.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>

#include <memory>
#include <numeric>
#include <cmath>

CP::TPulseDeconvolution::TPulseDeconvolution(int sampleCount) {
    fSampleCount = sampleCount;
    fSmoothingWindow = CP::TRuntimeParameters::Get().GetParameterI(
        "clusterCalib.deconvolution.smoothing");
    fSmoothingWindow = std::max(1,fSmoothingWindow + 1);
    fFluctuationCut 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.deconvolution.fluctuationCut");
    fBaselineCut 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.deconvolution.baselineCut");
    fCoherenceZone 
        = CP::TRuntimeParameters::Get().GetParameterI(
            "clusterCalib.deconvolution.coherenceZone");
    fCoherenceCut 
        = CP::TRuntimeParameters::Get().GetParameterI(
            "clusterCalib.deconvolution.coherenceCut");
    fFFT = NULL;
    fInverseFFT = NULL;
    fElectronicsResponse = NULL;
    fWireResponse = NULL;
    double noisePower = CP::TRuntimeParameters::Get().GetParameterD(
        "clusterCalib.deconvolution.noisePower");
    double spikePower = CP::TRuntimeParameters::Get().GetParameterD(
        "clusterCalib.deconvolution.spikePower");
    fNoiseFilter = new CP::TNoiseFilter(noisePower, spikePower);
    fBaselineSigma = 0.0;
    fSampleSigma = 0.0;
    Initialize();
}

CP::TPulseDeconvolution::~TPulseDeconvolution() {}

void CP::TPulseDeconvolution::Initialize() {
    fSampleCount = 2*(1+fSampleCount/2);
    int nSize = fSampleCount;
    CaptLog("Initialize pulse deconvolution to " << nSize << " entries");
    
    if (fFFT) delete fFFT;
    fFFT = TVirtualFFT::FFT(1, &nSize, "R2C M K");
    if (nSize != fSampleCount) {
        CaptError("Invalid length for FFT");
        CaptError("     original length: " << fSampleCount);
        CaptError("     allocated length: " << nSize);
    }

    if (fInverseFFT) delete fInverseFFT;
    fInverseFFT = TVirtualFFT::FFT(1, &nSize, "C2R M K");
    if (nSize != fSampleCount) {
        CaptError("Invalid length for inverse FFT");
        CaptError("     original length: " << fSampleCount);
        CaptError("     allocated length: " << nSize);
    }

    if (fElectronicsResponse) delete fElectronicsResponse;
    fElectronicsResponse = new CP::TElectronicsResponse(nSize);

    if (fWireResponse) delete fWireResponse;
    fWireResponse = new CP::TWireResponse(nSize);

    fBaselineSigma = 0.0;
}

CP::TCalibPulseDigit* CP::TPulseDeconvolution::operator() 
    (const CP::TCalibPulseDigit& calib) {

    CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
    TChannelCalib channelCalib;

    if (fSampleCount < (int) calib.GetSampleCount()) {
        CaptLog("Deconvolution size needs to be increased from "
                << fSampleCount << " to " << calib.GetSampleCount());
        fSampleCount = calib.GetSampleCount();
        Initialize();
    }

    fElectronicsResponse->Calculate(ev->GetContext(),
                                    calib.GetChannelId());
    fWireResponse->Calculate(ev->GetContext(),
                             calib.GetChannelId());

    // Copy the digit into the buffer for the FFT and then apply the
    // transformation.  The FFT buffer may be longer than the number of
    // samples read for this digit.  If there aren't enough samples, then make
    // sure the FFT buffer smoothly goes to zero.
    for (int i=0; i<fSampleCount; ++i) {
        // Apply error checking for invalid calibration results.  This should
        // never happen (and after fixing bugs it doesn't seem to).
        double p = calib.GetSample(i);
        if (!std::isfinite(p)) {
            CaptError("Channel " << calib.GetChannelId() 
                      << " w/ invalid sample " << i);
        }
        // Check if the FFT buffer is being filled past the end of the digit.
        // If it is, then smoothly go to zero.
        if ((int) calib.GetSampleCount() <= i) {
            int last = calib.GetSampleCount()-1;
            double scale = 10.0;
            double delta = (i-last)/scale; delta *= delta;
            double val = 0.0;
            if (delta < 40) {
                val = fFFT->GetPointReal(last,kTRUE)*std::exp(-delta);
            }
            fFFT->SetPoint(i,val);
            continue;
        }
        // Apply smoothing to the input signal.  The signal is not smoothed if
        // fSmoothingWindow is one.
        double val = 0.0;
        double weight = 0.0;
        for (int j=0; j<fSmoothingWindow; ++j) {
            double w = fSmoothingWindow - j;
            if (j == 0) {
                val += w * calib.GetSample(i);
                weight += w;
                continue;
            }
            if (0 <= (i-j)) {
                val += w * calib.GetSample(i-j);
                weight += w;
            }
            if ((i+j) < (int) calib.GetSampleCount()) {
                val += w * calib.GetSample(i+j);
                weight += w;
            }
        }
        val /= weight;
        // Make sure that the signal is zero at the very start.  This
        // distorts the beginning of the signal, but there is a long
        // buffer of noise there anyway.
        double scale = 10.001;
        if (i<scale) val *= 1.0*i/scale;
        fFFT->SetPoint(i,val);
    }

    fFFT->Transform();

    /// Make an optimal filter based on the electronics response, and the
    /// observed signal.
    fNoiseFilter->Calculate(calib.GetChannelId(), *fElectronicsResponse,
                            *fWireResponse, *fFFT);

    // Check if this is a valid channel.
    if (fNoiseFilter->IsNoisy()) {
        CaptLog("Noisy channel: " << calib.GetChannelId());
        return NULL;
    }
    
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    {
        double timeStep = channelCalib.GetTimeConstant(calib.GetChannelId());
    
        TH1F* fftHist 
            = new TH1F((calib.GetChannelId().AsString()+"-fft").c_str(),
                       ("FFT for " 
                        + calib.GetChannelId().AsString()).c_str(),
                       fSampleCount/2,
                       0,
                       (0.5/timeStep)/unit::hertz);
        for (std::size_t i = 0; i<fSampleCount/2; ++i) {
            double rl, im;
            fFFT->GetPointComplex(i,rl,im);
            std::complex<double> c(rl,im);
        fftHist->SetBinContent(i+1, std::abs(c));
        }
    }
#endif
        
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    {
        double timeStep = channelCalib.GetTimeConstant(calib.GetChannelId());
        TH1F* respHist 
            = new TH1F((calib.GetChannelId().AsString()+"-respFFT").c_str(),
                       ("FFT of Response for " 
                        + calib.GetChannelId().AsString()).c_str(),
                       fSampleCount/2,
                       0,
                       (0.5/timeStep)/unit::hertz);
        for (std::size_t i = 0; i<fSampleCount/2; ++i) {
            std::complex<double> freq = fElectronicsResponse->GetFrequency(i);
            std::complex<double> wire = fWireResponse->GetFrequency(i);
            std::complex<double> c = freq*wire;
            respHist->SetBinContent(i+1, std::abs(c));
        }
    }
#endif

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    {
        double timeStep = channelCalib.GetTimeConstant(calib.GetChannelId());
        TH1F* filterHist 
            = new TH1F((calib.GetChannelId().AsString()+"-filter").c_str(),
                       ("FFT Filter for " 
                        + calib.GetChannelId().AsString()).c_str(),
                       fSampleCount/2,
                       0,
                       (0.5/timeStep)/unit::hertz);
        for (std::size_t i = 0; i<fSampleCount/2; ++i) {
            filterHist->SetBinContent(i+1, fNoiseFilter->GetFilter(i));
        }
    }
#endif

    // Use the transformed values to do the deconvolution.  This fills the
    // buffer for the inverse FFT.
    for (int i=0; i<fSampleCount; ++i) {
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        c /= fElectronicsResponse->GetFrequency(i);
        c /= fWireResponse->GetFrequency(i);
        c *= fNoiseFilter->GetFilter(i);
        fInverseFFT->SetPoint(i,c.real(), c.imag());
    }

    fInverseFFT->Transform();
    std::auto_ptr<CP::TCalibPulseDigit> deconv(new CP::TCalibPulseDigit(calib));

    // Set the samples into the calibrated pulse digit.
    for (std::size_t i=0; i<deconv->GetSampleCount(); ++i) {
        double v = fInverseFFT->GetPointReal(i)/fSampleCount;
        deconv->SetSample(i,v);
    }

    RemoveBaseline(*deconv);

    return deconv.release();
}

void CP::TPulseDeconvolution::RemoveBaseline(CP::TCalibPulseDigit& digit) {
    std::vector<double> diff;
    diff.resize(digit.GetSampleCount());

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* offsetHist 
        = new TH1F((digit.GetChannelId().AsString()+"-offset").c_str(),
                   ("Deconvoluted signal before baseline removal for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   digit.GetFirstSample(), digit.GetLastSample());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        offsetHist->SetBinContent(i+1,digit.GetSample(i));
    }
#endif
        
    // Find the sample median and it's "sigma".
    for (std::size_t i=1; i<digit.GetSampleCount(); ++i) {
        diff[i] = digit.GetSample(i);
    }
    std::sort(diff.begin(), diff.end());

    double baselineMedian = diff[0.5*digit.GetSampleCount()];
    double baselineSigma = diff[0.16*digit.GetSampleCount()];
    baselineSigma = std::abs(baselineSigma-baselineMedian);

    // The minimum baseline fluctuation is 1 electron charge.
    if (baselineSigma < 1) baselineSigma = 1.0;
    
    // Define a maximum separation between the sample and the median sample.
    // If it's more than this, then the sample is not baseline.  This is one
    // sided.  This is looking at the overall fluctation of the baseline and
    // prevents large plateaus from being cut.
    double baselineCut = baselineMedian + fBaselineCut*baselineSigma;

    // Find the median sample to sample difference.  Regions where the samples
    // stay withing a small difference don't have a "feature of interest".
    for (std::size_t i=1; i<digit.GetSampleCount(); ++i) {
        double delta = std::abs(digit.GetSample(i) - digit.GetSample(i-1));
        diff[i] = delta;
    }
    std::sort(diff.begin(), diff.end());

    // The 52 percentile is a of the sample-to-sample differences is a good
    // estimate of the distribution sigma.  The 52% "magic" number comes
    // looking at the differences of two samples (that gives a sqrt(2)) and
    // wanting the RMS (which is the 68%).  This gives 48% (which is
    // 68%/sqrt(2)).  Then because of the ordering after the sort, the bin we
    // look at is 1.0-52%.  The minimum allowable fluctuation is 1 electron
    // charge.
    fSampleSigma = diff[0.52*digit.GetSampleCount()];
    if (fSampleSigma < 1.0) fSampleSigma = 1.0;
    
    // Define a cut for the maximum change between two samples.  If the
    // difference is greater than this, then the samples are not "coherent".
    // If too many samples in a "coherence zone" are not considered coherent,
    // then the sample is not considered part of the baseline.
    double deltaCut = fSampleSigma*fFluctuationCut;

    // Define a cut for the maximum drift in one direction within a coherence
    // zone.  If the drift is more than this, then the region is not baseline.
    TChannelCalib channelCalib;
    double driftSigma = fSampleSigma;
    if (channelCalib.IsBipolarSignal(digit.GetChannelId())) {
        driftSigma  *= std::sqrt(fCoherenceZone);
    }
    double driftCut = fFluctuationCut * driftSigma;

    // Refill the differences... This has to be done since we are reusing an
    // existing variable (instead of allocating extra space).
    diff[0] = 0;
    for (std::size_t i=1; i < diff.size(); ++i) {
        double delta = std::abs(digit.GetSample(i) - digit.GetSample(i-1));
        diff[i] = delta;
    }


#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* diffHist 
        = new TH1F((digit.GetChannelId().AsString()+"-diff").c_str(),
                   ("Differences from previous sample for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   digit.GetFirstSample(), digit.GetLastSample());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        diffHist->SetBinContent(i+1,diff[i]/fSampleSigma);
    }
#endif
        
    // Estimate the baseline for regions where there isn't much change.
    std::vector<double> baseline;
    baseline.resize(digit.GetSampleCount());
    const double unfilledBaseline = std::numeric_limits<double>::quiet_NaN();
    std::fill(baseline.begin(), baseline.end(), unfilledBaseline);
    
    std::vector<double> drift;
    drift.resize(digit.GetSampleCount());

    // Fill the drifts for the samples.
    for (std::size_t i=0; i<drift.size(); ++i) {
        // Find the start of the coherence zone.
        int startZone = i - fCoherenceZone;

        if (startZone < 0) startZone = 0;

        double deltaSample = 0.0;
        for (std::size_t j = startZone; j<i; ++j) {
            deltaSample = std::max(deltaSample,
                                   std::abs(digit.GetSample(i)
                                            -digit.GetSample(j)));
        }
        drift[i] = deltaSample;
    }

    for (std::size_t i=drift.size()-1; 0 < i; --i) {
        // Find the start of the coherence zone.
        std::size_t startZone = i + fCoherenceZone;

        double deltaSample = 0.0;
        for (std::size_t j = i+1; j<=startZone; ++j) {
            deltaSample = std::max(deltaSample,
                                   std::abs(digit.GetSample(i)
                                            -digit.GetSample(j)));
        }
        drift[i] = std::max(deltaSample, drift[i]);
    }

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* driftHist 
        = new TH1F((digit.GetChannelId().AsString()+"-drift").c_str(),
                   ("Drift for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   digit.GetFirstSample(), digit.GetLastSample());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        driftHist->SetBinContent(i+1, drift[i]/driftSigma);
    }
#endif
        
    // First look forward through the samples.
    int coh = 0;
    for (std::size_t i=0; i < diff.size(); ++i) {
        // Find the start of the coherence zone.
        int startZone = i - fCoherenceZone;

        // Find the number of differences in the coherence zone that are less
        // than the difference cut.  This is a running sum that increments the
        // count if the new bin is less than the cut, and decrements it when
        // the bin leaves the zone.
        if (diff[i] < deltaCut) coh += 1;
        if (0 <= startZone && diff[startZone] < deltaCut) coh -= 1;

        // Handle the special cases.  This could be done more compactly, but I
        // want to be very explicit about which criteria is triggered.
        if (startZone < 0) {
            // The first several bins are always considered baseline...
            baseline[i] = digit.GetSample(i);
            continue;
        }

        // Don't consider very high signals to be baseline
        if (digit.GetSample(i) > baselineCut) {
            continue;
        }

        // If a region has changed too far in one direction, then it's not
        // part of the baseline.
        if (drift[i] > driftCut) {
            continue;
        }
        
        // Look at the number of samples since a measurement was more than a
        // statistical fluctuation.  If it's only been a few samples
        // (i.e. less than fCoherenceCut samples) then this sample isn't part
        // of the baseline.
        if (coh < fCoherenceCut) {
            continue;
        }

        // We've found a sample that should be considered part of the
        // baseline, so add it to the vector for book keeping.
        int offset = (int) (0.5*fCoherenceZone);
        baseline[i-offset] = digit.GetSample(i-offset);
    }

    coh = 0;
    // Now look backwards through the differences.
    for (std::size_t i=diff.size()-1; 0 < i; --i) {
        // Find the start of the coherence zone.
        std::size_t startZone = i + fCoherenceZone;

        // Find the number of differences in the coherence zone that are less
        // than the difference cut.
        if (diff[i] < deltaCut) {
            coh += 1;
        }
        if (startZone < diff.size() && diff[startZone] < deltaCut) {
            coh -= 1;
        }

        // Handle the special cases.  This could be done more compactly, but I
        // want to be very explicit about which criteria is triggered.
        if (diff.size() <= startZone) {
            // The first several bins are always considered baseline...
            baseline[i] = digit.GetSample(i);
            continue;
        }

        // Don't consider very high signals to be baseline
        if (digit.GetSample(i)> baselineCut) {
            continue;
        }

        // If a region has changed too far in one direction, then it's not
        // part of the baseline.
        if (drift[i] > driftCut) {
            continue;
        }
        
        // Look at the number of samples since a measurement was more than a
        // statistical fluctuation.  If it's only been a few samples
        // (i.e. less than fCoherenceCut samples) then this sample isn't part
        // of the baseline.
        if (coh < fCoherenceCut) {
            continue;
        }

        // We've found a sample that should be considered part of the
        // baseline, so add it to the vector for book keeping.
        int offset = (int) (0.5*fCoherenceZone);
        baseline[i+offset] = digit.GetSample(i+offset);

    }

    // The baseline will equal to unfilledBaseline for any point that is not
    // baseline like with gaps that correspond to where the sample is changing
    // more rapidly that would be expected from pure statistics.  The next
    // step is to tentatively interpolate into those zones and check if the
    // interpolated baseline is bigger than the sample.  If the interpolated
    // baseline is bigger, then the sample becomes the new baseline for that
    // point (but not if it's a big fluctuation from it's neighbors).  The
    // drift variable will take the new candidate baseline values.  This
    // overwrites the drift vector to save an allocation (i.e. false
    // optimization!).  After the end of this section, drift will have a value
    // for anyplace where the signal is not baseline like.
    for (std::size_t i = 1; i<baseline.size()-1; ++i) {
        drift[i] = unfilledBaseline;
        if (std::isfinite(baseline[i])) continue;
        std::size_t j = i;
        while (0<j && !std::isfinite(baseline[j])) --j;
        std::size_t k = i;
        while (k<baseline.size() && !std::isfinite(baseline[k])) ++k;
        double interp = (k-i)*baseline[j] + (i-j)*baseline[k];
        interp /= 1.0*(k-j);
        if (digit.GetSample(i) < interp
            && diff[i]/fSampleSigma < 3.0
            && diff[i+1]/fSampleSigma < 3.0) {
            drift[i] = digit.GetSample(i);
        }
    }

    // Reject any new baseline regions where its only a few bins below the
    // interpolated baseline.  This is done by resetting drift to
    // unfilledBaseline
    for (std::size_t i = 0; i<baseline.size(); ++i) {
        if (!std::isfinite(drift[i])) continue;
        std::size_t j = i;
        while (j<drift.size() && std::isfinite(drift[j])) ++j;
        if ((j-i) < 3) {
            for (std::size_t k = i; k<j; ++k) drift[k] = unfilledBaseline;
        }
        i = j;
    }

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* interpHist 
        = new TH1F((digit.GetChannelId().AsString()+"-interp").c_str(),
                   ("Interp for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   digit.GetFirstSample(), digit.GetLastSample());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        if (!std::isfinite(drift[i])) continue;
        interpHist->SetBinContent(i+1, drift[i]);
    }
#endif
        
    // Smooth the baseline to get rid of the highest frequency components, and
    // put the smoothed components into drift.  This fills in the gaps in
    // drift around where the baseline had already been estimated.
    int baselineWindow = 5;
    for (int i = 0; i<(int)baseline.size(); ++i) {
        // Skip non-baseline like regions.
        if (!std::isfinite(baseline[i])) continue;
        // Apply smoothing to the input signal.  This is controlled with
        // clusterCalib.smoothing.wire
        double val = 0.0;
        double weight = 0.0;
        for (int j=0; j<baselineWindow; ++j) {
            double w = baselineWindow - j;
            if (j == 0 && std::isfinite(baseline[i])) {
                val += w * baseline[i];
                weight += w;
                continue;
            }
            if (0 <= (i-j) && std::isfinite(baseline[i-j])) {
                val += w * baseline[i-j];
                weight += w;
            }
            if ((i+j) < (int) baseline.size() && std::isfinite(baseline[i+j])) {
                val += w * baseline[i+j];
                weight += w;
            }
        }
        if (weight < 0.001) continue;
        drift[i] = val/weight;
    }

    // copy what we know back into baseline.
    std::copy(drift.begin(), drift.end(), baseline.begin());

    // Now do a final interpolation to get the baseline everywhere.
    for (std::size_t i = 0; i<baseline.size(); ++i) {
        if (std::isfinite(drift[i])) continue;
        std::size_t j = i;
        while (0<j && !std::isfinite(drift[j])) --j;
        std::size_t k = i;
        while (k<drift.size() && !std::isfinite(drift[k])) ++k;
        double interp = (k-i)*drift[j] + (i-j)*drift[k];
        interp /= 1.0*(k-j);
        baseline[i] = interp;
    }

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* bkgHist 
        = new TH1F((digit.GetChannelId().AsString()+"-baseline").c_str(),
                   ("Estimated baseline for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   digit.GetFirstSample(), digit.GetLastSample());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        if (!std::isfinite(baseline[i])) {
            CaptError("Baseline not filled at " << i);
        }
        bkgHist->SetBinContent(i+1,baseline[i]);
    }
#endif

    // Now remove the baseline from the sample and calculate the sigma for the
    // baseline (relative to the mean baseline).
    fBaselineSigma = 0.0;
    double aveBaseline = 0.0;
    for (std::size_t i=0; i<digit.GetSampleCount(); ++i) {
        if (!std::isfinite(baseline[i])) {
            CaptError("Baseline not filled at " << i);
        }
        fBaselineSigma += baseline[i]*baseline[i];
        aveBaseline += baseline[i];
        if (channelCalib.IsBipolarSignal(digit.GetChannelId())) {
            // Only remove the baseline if it's a bipolar channel.
            double d = digit.GetSample(i) - baseline[i];
            digit.SetSample(i,d);
        }
    }
    fBaselineSigma /= digit.GetSampleCount();
    aveBaseline /= digit.GetSampleCount();
    if (fBaselineSigma > 0.0) {
        fBaselineSigma = std::sqrt(fBaselineSigma-aveBaseline*aveBaseline);
    }
    else {
        fBaselineSigma = 0.0;
    }

}
