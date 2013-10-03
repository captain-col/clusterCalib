#include "TPulseDeconvolution.hxx"
#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"

#include <TCalibPulseDigit.hxx>
#include <TEvent.hxx>
#include <TEventFolder.hxx>
#include <TRuntimeParameters.hxx>

#include <TVirtualFFT.h>

#include <memory>
#include <cmath>

CP::TPulseDeconvolution::TPulseDeconvolution(int sampleCount) {
    fSampleCount = sampleCount;
    fSmoothingWindow = CP::TRuntimeParameters::Get().GetParameterI(
        "clusterCalib.smoothing.wire");
    fSmoothingWindow = std::max(1,fSmoothingWindow + 1);
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

#define SMOOTH
#ifdef SMOOTH
    std::vector<double> baseline;
    baseline.resize(deconv->GetSampleCount());
    // Set the samples into the calibrated pulse digit.
    int window = std::min(81,(int) deconv->GetSampleCount()/10);
    for (int i=0; i < (int) deconv->GetSampleCount(); ++i) {
        int low = std::max(0,i-window/2);
        int high = low + window;
        if (high > (int) deconv->GetSampleCount()) {
            low = deconv->GetSampleCount() - window;
            high = deconv->GetSampleCount();
        }
        double v = 0.0;
        double w = 0.0;
        for (int j = low; j < high; ++j) {
            double r = 1 - std::abs((j - 0.5*(high+low))/(0.5*(high-low+1)));
            v += r*deconv->GetSample(j);
            w += r;
        }
        v /= w;
        baseline[i]=v;
    }
    for (int i=0; i < (int) deconv->GetSampleCount(); ++i) {
        double v = deconv->GetSample(i) - baseline[i];
        deconv->SetSample(i,v);
    }
#endif

    return deconv.release();
}
