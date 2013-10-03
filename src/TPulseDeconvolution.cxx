#include "TPulseDeconvolution.hxx"
#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"

#include <TCalibPulseDigit.hxx>
#include <TEvent.hxx>
#include <TEventFolder.hxx>
#include <TRuntimeParameters.hxx>
#include <TMCChannelId.hxx>

#include <TVirtualFFT.h>

#include <memory>
#include <cmath>

CP::TPulseDeconvolution::TPulseDeconvolution(int sampleCount) {
    fSampleCount = sampleCount;
    fSmoothingWindow = CP::TRuntimeParameters::Get().GetParameterI(
        "clusterCalib.smoothing.wire");
    fSmoothingWindow = std::max(1,fSmoothingWindow + 1);
    fBaselineWindow = 20;
    fFluctuationCut = 3.0;
    fFFT = NULL;
    fInverseFFT = NULL;
    fElectronicsResponse = NULL;
    fWireResponse = NULL;
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

    // Calculate the average f all of the deltas.
    double avgDelta = 0.0;
    double deltaCount = 0.0;
    for (std::size_t i=1; i<digit.GetSampleCount(); ++i) {
        avgDelta += std::abs(digit.GetSample(i) - digit.GetSample(i-1));
        deltaCount += 1.0;
    }
    avgDelta /= deltaCount;
    
    // Now calculate the average of the deltas, but exclude the outliers
    double  truncatedDelta = 0.0;
    deltaCount = 0.0;
    for (std::size_t i=1; i<digit.GetSampleCount(); ++i) {
        double d = std::abs(digit.GetSample(i) - digit.GetSample(i-1));
        if (d>fFluctuationCut*avgDelta) continue;
        truncatedDelta += d;
        deltaCount += 1.0;
    }
    truncatedDelta /= deltaCount;
    
    std::vector<double> baseline;
    baseline.resize(digit.GetSampleCount());

    for (std::size_t i=0; i<digit.GetSampleCount(); ++i) {
        // Find a region of just random fluctuations below the current
        // element.
        std::size_t low = 0;
        if (i > 2*fBaselineWindow) {
            low = i - 2*fBaselineWindow;
            while (low>1) {
                bool found = true;
                for (std::size_t j=0; j<fBaselineWindow; ++j) {
                    double d = std::abs(digit.GetSample(low+j)
                                        -digit.GetSample(low+j+1));
                    if (d>fFluctuationCut*truncatedDelta) {
                        found = false;
                        break;
                    }
                }
                if (found) break;
                --low;
            }
        }
        
        // Find a region of just random fluctuations above the current
        // element.
        std::size_t high = digit.GetSampleCount();
        if (i + 2*fBaselineWindow < high) {
            high = i + 2*fBaselineWindow;
            while (high < digit.GetSampleCount()) {
                bool found = true;
                for (std::size_t j=0; j<fBaselineWindow; ++j) {
                    double d = std::abs(digit.GetSample(high-j)
                                        -digit.GetSample(high-j-1));
                    if (d>fFluctuationCut*truncatedDelta) {
                        found = false;
                        break;
                    }
                }
                if (found) break;
                ++high;
            }
        }
        
        double base=0.0;
        double baseWght = 0.0;
        // Find the lower average.
        for (std::size_t j = 0; j<fBaselineWindow; ++j) {
            double r = (j+1)*(fBaselineWindow-j+1);
            base += r*digit.GetSample(low+j);
            baseWght += r;
        }
        // Find the upper average.
        for (std::size_t j = 0; j<fBaselineWindow; ++j) {
            double r = (j+1)*(fBaselineWindow-j+1);
            base += r*digit.GetSample(high-j);
            baseWght += r;
        }
        baseline[i] = base/baseWght;
    }

    for (int i=0; i < (int) digit.GetSampleCount(); ++i) {
        double v = digit.GetSample(i) - baseline[i];
        digit.SetSample(i,v);
    }

}
