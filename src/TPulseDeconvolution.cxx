#include "TPulseDeconvolution.hxx"
#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"
#include "TChannelCalib.hxx"

#include <TCalibPulseDigit.hxx>
#include <TEvent.hxx>
#include <TEventFolder.hxx>
#include <TRuntimeParameters.hxx>
#include <TMCChannelId.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>

#include <memory>
#include <cmath>

CP::TPulseDeconvolution::TPulseDeconvolution(int sampleCount) {
    fSampleCount = sampleCount;
    fSmoothingWindow = CP::TRuntimeParameters::Get().GetParameterI(
        "clusterCalib.smoothing.wire");
    fSmoothingWindow = std::max(1,fSmoothingWindow + 1);
    fFluctuationCut = 3.0;
    fCoherenceZone = 10;
    fCoherenceFraction = 1.1;
    fFFT = NULL;
    fInverseFFT = NULL;
    fElectronicsResponse = NULL;
    fWireResponse = NULL;
    fBaselineSigma = 0.0;
    Initialize();
}

CP::TPulseDeconvolution::~TPulseDeconvolution() {}

void CP::TPulseDeconvolution::Initialize() {
    fSampleCount = 2*(1+fSampleCount/2);
    int nSize = fSampleCount;
    CaptLog("Initialize pulse deconvolution to " << nSize << " entries");
    
    if (fFFT) delete fFFT;
    fFFT = TVirtualFFT::FFT(1, &nSize, "R2C ES K");
    if (nSize != fSampleCount) {
        CaptError("Invalid length for FFT");
        CaptError("     original length: " << fSampleCount);
        CaptError("     allocated length: " << nSize);
    }

    if (fInverseFFT) delete fInverseFFT;
    fInverseFFT = TVirtualFFT::FFT(1, &nSize, "C2R ES K");
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
    // transformation.
    for (int i=0; i<fSampleCount; ++i) {
        if (i < (int) calib.GetSampleCount()) {
            // Apply smoothing to the input signal.  This is controlled with
            // clusterCalib.smoothing.wire
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
        else {
            int last = calib.GetSampleCount()-1;
            double scale = 10.0;
            double delta = (i-last)/scale; delta *= delta;
            double val = 0.0;
            if (delta < 40) {
                val = fFFT->GetPointReal(last,kTRUE)*std::exp(-delta);
            }
            fFFT->SetPoint(i,val);
        }
    }
    fFFT->Transform();
    
    for (int i=0; i<fSampleCount; ++i) {
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        c /= fElectronicsResponse->GetFrequency(i);
        c /= fWireResponse->GetFrequency(i);
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

    // Find the median sample to sample difference.  Regions where the samples
    // stay withing a small difference don't have a "feature of interest".
    double avgDelta = 0.0;
    for (std::size_t i=1; i<digit.GetSampleCount(); ++i) {
        double delta = std::abs(digit.GetSample(i) - digit.GetSample(i-1));
        diff[i] = delta;
    }
    std::sort(diff.begin(), diff.end());

    avgDelta = diff[0.5*digit.GetSampleCount()];

    double deltaCut = avgDelta*fFluctuationCut;
    CP::TMCChannelId channel(digit.GetChannelId());

    // Refill the differences...
    diff[0] = 0;
    for (std::size_t i=1; i < diff.size(); ++i) {
        double delta = std::abs(digit.GetSample(i) - digit.GetSample(i-1));
        diff[i] = delta;
    }

    // Estimate the baseline for regions where there isn't much change.
    std::vector<double> baseline;
    baseline.resize(digit.GetSampleCount());

    // First look forward through the samples.
    int coh = 0;
    for (std::size_t i=0; i < diff.size(); ++i) {
        // Find the start of the coherence zone.
        int startZone = i - fCoherenceZone;

        // Find the number of differences in the coherence zone that are less
        // than the difference cut.
        if (diff[i] < deltaCut) coh += 1.0;
        if (0 <= startZone && diff[startZone] < deltaCut) coh -= 1.0;

        // Handle the special cases.  This could be done more compactly, but I
        // want to be very explicit about which criteria is triggered.
        if (startZone < 0) {
            // The first several bins are always considered baseline...
            baseline[i] = digit.GetSample(i);
        }
        else if (fCoherenceZone <= coh) {
            baseline[i] = digit.GetSample(i);
        }

    }

    coh = 0;
    // Now look backwards through the differences.
    for (std::size_t i=diff.size()-1; 0 < i; --i) {
        // Find the start of the coherence zone.
        std::size_t startZone = i + fCoherenceZone;

        // Find the number of differences in the coherence zone that are less
        // than the difference cut.
        if (diff[i] < deltaCut) {
            coh += 1.0;
        }
        if (startZone < diff.size() && diff[startZone] < deltaCut) {
            coh -= 1.0;
        }

        // Handle the special cases.  This could be done more compactly, but I
        // want to be very explicit about which criteria is triggered.
        if (diff.size() <= startZone) {
            // The first several bins are always considered baseline...
            baseline[i] = digit.GetSample(i);
        }
        else if (fCoherenceZone <= coh) {
            baseline[i] = digit.GetSample(i);
        }

    }

    // The baseline will be non-zero for any point that is baseline like, but
    // there will be gaps of zeros that correspond to where the sample is
    // changing more rapidly that would be expected from pure statistics.  The
    // next step is to tentatively interpolate into those zones and check if
    // the interpolated baseline is bigger than the sample.  If the
    // interpolated baseline is bigger, then the sample becomes the new
    // baseline for that point.
    for (std::size_t i = 0; i<baseline.size(); ++i) {
        if (std::abs(baseline[i]) > 1E-6) continue;
        std::size_t j = i;
        while (0<j && std::abs(baseline[j]) < 1E-6) --j;
        std::size_t k = i;
        while (k<baseline.size() && std::abs(baseline[k]) < 1E-6) ++k;
        double interp = (k-i)*baseline[j] + (i-j)*baseline[k];
        interp /= 1.0*(k-j);
        if (digit.GetSample(i) < interp) baseline[i] = digit.GetSample(i);
    }

    // Now do a final interpolation to get the baseline everywhere.
    for (std::size_t i = 0; i<baseline.size(); ++i) {
        if (std::abs(baseline[i]) > 1E-6) continue;
        std::size_t j = i;
        while (0<j && std::abs(baseline[j]) < 1E-6) --j;
        std::size_t k = i;
        while (k<baseline.size() && std::abs(baseline[k]) < 1E-6) ++k;
        double interp = (k-i)*baseline[j] + (i-j)*baseline[k];
        interp /= 1.0*k-j;
        baseline[i] = interp;
    }
    
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* bkgHist 
        = new TH1F((digit.GetChannelId().AsString()+"-bkg").c_str(),
                   ("Bkg for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   0.0, 1.0*digit.GetSampleCount());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        bkgHist->SetBinContent(i+1,baseline[i]);
    }
#endif
        
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* diffHist 
        = new TH1F((digit.GetChannelId().AsString()+"-diff").c_str(),
                   ("Diff for " 
                    + digit.GetChannelId().AsString()).c_str(),
                   digit.GetSampleCount(),
                   0.0, 1.0*digit.GetSampleCount());
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        diffHist->SetBinContent(i+1,diff[i]);
    }
#endif
        
    // Now remove the baseline from the sample and calculate the sigma for the
    // baseline (relative to the mean baseline).
    fBaselineSigma = 0.0;
    for (std::size_t i=0; i<digit.GetSampleCount(); ++i) {
        fBaselineSigma += baseline[i]*baseline[i];
        double d = digit.GetSample(i) - baseline[i];
        digit.SetSample(i,d);
    }
    fBaselineSigma /= digit.GetSampleCount();
    if (fBaselineSigma > 0.0) {
        fBaselineSigma = std::sqrt(fBaselineSigma);
    }
    else {
        fBaselineSigma = 0.0;

    }

}
