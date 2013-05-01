#include "TElectronicsResponse.hxx"

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

CP::TElectronicsResponse::TElectronicsResponse() 
    : fMustRecalculate(true), fPeakingTime(-1.0), fSampleTime(-1.0) {}

CP::TElectronicsResponse::TElectronicsResponse(int size) 
    : fMustRecalculate(true), fPeakingTime(-1.0), fSampleTime(-1.0) {
    fResponse.resize(size);
    fFrequency.resize(size);
}

CP::TElectronicsResponse::~TElectronicsResponse() { }

void CP::TElectronicsResponse::SetSize(int s) {
    if (s == (int) fResponse.size()) return;
    fMustRecalculate = true;
    fResponse.resize(s);
    fFrequency.resize(s);
}

void CP::TElectronicsResponse::SetStep(double s) {
    do {
        if (fSampleTime <= 0.0) {
            fMustRecalculate = true;
            break;
        }
        double change = (s-fSampleTime)/std::max(s,fSampleTime);
        if (change > 1E-8) {
            fMustRecalculate = true;
            break;
        }
    } while (false);
    if (fMustRecalculate) fSampleTime = s;
}

void CP::TElectronicsResponse::SetPeakingTime(double s) {
    do {
        if (fPeakingTime <= 0.0) {
            fMustRecalculate = true;
            break;
        }
        double change = (s-fPeakingTime)/std::max(s,fPeakingTime);
        if (change > 1E-8) {
            fMustRecalculate = true;
            break;
        }
    } while (false);
    if (fMustRecalculate) fPeakingTime = s;
}

CP::TElectronicsResponse::Response
CP::TElectronicsResponse::GetResponse(int i) {
    if (i<0) return 0.0;
    if ((int) fResponse.size() <= i) return 0.0;
    if (fMustRecalculate) Calculate(GetPeakingTime());
    return fResponse[i];
}

CP::TElectronicsResponse::Frequency 
CP::TElectronicsResponse::GetFrequency(int i) {
    if (i<0) return Frequency(0,0);
    if ((int) fFrequency.size() <= i) return Frequency(0,0);
    if (fMustRecalculate) Calculate(GetPeakingTime());
    return fFrequency[i];
}

bool CP::TElectronicsResponse::Calculate(double peakingTime) {
    fPeakingTime = peakingTime;
    fMustRecalculate = false;

    // Fill the response function.  This explicitly normalizes.
    double sum = 0.0;
    for (std::size_t i=0; i<fResponse.size(); ++i) {
        double arg = fSampleTime*i/fPeakingTime;
        double v = (arg<40)? arg*std::exp(-arg): 0.0;
        fResponse[i] = v;
        sum += v;
    }
    if (sum < 1E-20) return false;
    for (std::vector<Response>::iterator r = fResponse.begin();
         r != fResponse.end(); ++r) {
        (*r) /= sum;
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

    return true;
}

bool CP::TElectronicsResponse::Calculate(const CP::TEventContext& context,
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
        
        CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
        
        // Get the electronics shape.
        CP::THandle<CP::TRealDatum> shapeVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/shape");
        double val = (*shapeVect)[index];
        SetPeakingTime(val);

        // get the length of each sample.
        CP::THandle<CP::TRealDatum> stepVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/digitStep");
        val = (*stepVect)[index];
        SetStep(val);
    }
    
    return true;
}
