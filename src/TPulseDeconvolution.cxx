#include "TPulseDeconvolution.hxx"
#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"

#include <TCalibPulseDigit.hxx>
#include <TEvent.hxx>
#include <TEventFolder.hxx>

#include <TVirtualFFT.h>

#include <memory>

CP::TPulseDeconvolution::TPulseDeconvolution(int sampleCount) {
    fSampleCount = sampleCount;
    fSampleCount = 2*(1+fSampleCount/2);
    int nSize = fSampleCount;
    
    fFFT = TVirtualFFT::FFT(1, &nSize, "R2C EX K");
    if (nSize != fSampleCount) {
        CaptError("Invalid length for FFT");
        CaptError("     original length: " << fSampleCount);
        CaptError("     allocated length: " << nSize);
    }

    fInverseFFT = TVirtualFFT::FFT(1, &nSize, "C2R EX K");
    if (nSize != fSampleCount) {
        CaptError("Invalid length for inverse FFT");
        CaptError("     original length: " << fSampleCount);
        CaptError("     allocated length: " << nSize);
    }

    fElectronicsResponse = new CP::TElectronicsResponse(nSize);
    fWireResponse = new CP::TWireResponse(nSize);
}

CP::TPulseDeconvolution::~TPulseDeconvolution() {}

CP::TCalibPulseDigit* CP::TPulseDeconvolution::operator() 
    (const CP::TCalibPulseDigit& calib) {
    
    CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
        
    fElectronicsResponse->Calculate(ev->GetContext(),
                                    calib.GetChannelId());
    fWireResponse->Calculate(ev->GetContext(),calib.GetChannelId());

    // This should get the wrap around right so that this is periodic and
    // there isn't a step.
    for (int i=0; i<fSampleCount; ++i) {
        if (i < (int) calib.GetSampleCount()) {
            fFFT->SetPoint(i,calib.GetSample(i));
        }
        else {
            fFFT->SetPoint(i,0.0);
        }
    }
    fFFT->Transform();
    
    for (int i=0; i<fSampleCount; ++i) {
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        std::complex<double> dc = c/fElectronicsResponse->GetFrequency(i);
        fInverseFFT->SetPoint(i,dc.real(), dc.imag());
    }
    fInverseFFT->Transform();
    std::auto_ptr<CP::TCalibPulseDigit> deconv(new CP::TCalibPulseDigit(calib));
    
    for (std::size_t i=0; i<deconv->GetSampleCount(); ++i) {
        double v = fInverseFFT->GetPointReal(i)/fSampleCount;
        deconv->SetSample(i,v);
    }

    return deconv.release();
}
