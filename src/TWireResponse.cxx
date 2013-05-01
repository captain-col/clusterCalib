#include "TWireResponse.hxx"

#include <TEvent.hxx>
#include <TEventFolder.hxx>
#include <TEventContext.hxx>
#include <THandle.hxx>
#include <TRealDatum.hxx>
#include <TMCChannelId.hxx>
#include <TChannelId.hxx>

#include <TVirtualFFT.h>

#include <cmath>
#include <memory>

CP::TWireResponse::TWireResponse() 
    : fMustRecalculate(true), fWireClass(0) {}

CP::TWireResponse::TWireResponse(int size) 
    : fMustRecalculate(true), fWireClass(0) {
    fResponse.resize(size);
    fFrequency.resize(size);
}

CP::TWireResponse::~TWireResponse() { }

void CP::TWireResponse::SetSize(int s) {
    if (s == (int) fResponse.size()) return;
    fMustRecalculate = true;
    fResponse.resize(s);
    fFrequency.resize(s);
}

CP::TWireResponse::Response
CP::TWireResponse::GetResponse(int i) {
    if (i<0) return 0.0;
    if ((int) fResponse.size() <= i) return 0.0;
    if (fMustRecalculate) Calculate();
    return fResponse[i];
}

CP::TWireResponse::Frequency 
CP::TWireResponse::GetFrequency(int i) {
    if (i<0) return Frequency(0,0);
    if ((int) fFrequency.size() <= i) return Frequency(0,0);
    if (fMustRecalculate) Calculate();
    return fFrequency[i];
}

bool CP::TWireResponse::Calculate() {
    fMustRecalculate = false;

    // Fill the response function.  This explicitly normalizes.
    switch (fWireClass) {
    case 0:
        for (std::size_t i=0; i<fResponse.size(); ++i) {
            fResponse[i] = 0;
            fFrequency[i] = Frequency(1,0);
        }
        fResponse[0] = 1;
        return true;
    case 1:
        for (std::size_t i=0; i<fResponse.size(); ++i) {
            fResponse[i] = 0;
            fFrequency[i] = Frequency(1,0);
        }
        fResponse[0] = -1;
        fResponse[fResponse.size()-1] = 1;
        break;
    default:
        CaptError("Unknown wire class");
        exit(1);
    }

    // Calculate the FFT of the response 
    int size = fResponse.size();
    std::auto_ptr<TVirtualFFT> fft(TVirtualFFT::FFT(1,&size,"R2C ES K"));
    if (size != (int) fResponse.size()) {
        CaptError("Invalid length for FFT");
        CaptError("     original length: " << fResponse.size());
        CaptError("     allocated length: " << size);
    }
    for (std::size_t i=0; i<fResponse.size(); ++i) {
        fft->SetPoint(i,fResponse[i]);
    }
    fft->Transform();
    for (std::size_t i=0; i<fResponse.size(); ++i) {
        double rl,im;
        fft->GetPointComplex(i,rl,im);
        fFrequency[i] = std::complex<double>(rl,im);
    }
    fFrequency[0] = 1.0;

    return true;
}

bool CP::TWireResponse::Calculate(const CP::TEventContext& context,
                                  const CP::TChannelId channel) {
    fEventContext = context;
    fChannelId = channel;

    // Get the channel calibrations.
    if (channel.IsMCChannel()) {
        TMCChannelId mc(channel);
        int index = -1;
        if (mc.GetType() == 0) index = mc.GetSequence();
        else if (mc.GetType() == 1) index = 3;
        else {
            CaptError("Unknown channel: " << mc);
            return false;
        }
        switch (index) {
        case 0:
        case 3:
            if (fWireClass != 0) {
                fWireClass = 0;
                fMustRecalculate = true;
            }
            break;
        case 1:
        case 2:
            if (fWireClass != 1) {
                fWireClass = 1;
                fMustRecalculate = true;
            }
            break;
        default:
            CaptError("Unknown wire type");
            std::exit(1);
        }
    }
    
    return true;
}
