#include "TElectronicsResponse.hxx"

#include <TEvent.hxx>
#include <TEventFolder.hxx>
#include <TEventContext.hxx>
#include <THandle.hxx>
#include <TRealDatum.hxx>
#include <TMCChannelId.hxx>
#include <TChannelId.hxx>
#include <TChannelCalib.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>

#include <cmath>
#include <memory>

CP::TElectronicsResponse::TElectronicsResponse() 
    : fMustRecalculate(true) {}

CP::TElectronicsResponse::TElectronicsResponse(int size) 
    : fMustRecalculate(true) {
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

CP::TElectronicsResponse::Response
CP::TElectronicsResponse::GetResponse(int i) {
    if (i<0) return 0.0;
    if ((int) fResponse.size() <= i) return 0.0;
    Calculate();
    return fResponse[i];
}

CP::TElectronicsResponse::Frequency 
CP::TElectronicsResponse::GetFrequency(int i) {
    if (i<0) return Frequency(0,0);
    if ((int) fFrequency.size() <= i) return Frequency(0,0);
    Calculate();
    return fFrequency[i];
}

bool CP::TElectronicsResponse::Calculate() {
    if (!fMustRecalculate) return false;
    fMustRecalculate = false;

    TChannelCalib calib;
    
    // Fill the response function.  This explicitly normalizes.
    double normalization = 0.0;
    for (std::size_t i=0; i<fResponse.size(); ++i) {
        double arg = calib.GetTimeConstant(fChannelId)*(1.0*i+0.5);
        double v = calib.GetPulseShape(fChannelId,arg);
        fResponse[i] = v;
#ifdef UNIT_NORMALIZATION
        normalization += v;
#else
        normalization = std::max(normalization,v);
#endif
    }
    if (normalization < 1E-20) return false;
    for (std::vector<Response>::iterator r = fResponse.begin();
         r != fResponse.end(); ++r) {
        (*r) /= normalization;
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

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* elecResp = new TH1F("elecResp",
                              "Electronics Response",
                              fResponse.size(),
                              0.0, 1.0*fResponse.size());
    for (std::size_t i=0; i<fResponse.size(); ++i) {
        elecResp->Fill(i+0.5, std::abs(fResponse[i]));
    }
    TH1F* elecFreq = new TH1F("elecFreq",
                              "Electronics Frequency",
                              fFrequency.size(),
                              0.0, 1.0*fFrequency.size());
    for (std::size_t i=0; i<fFrequency.size(); ++i) {
        elecFreq->Fill(i+0.5, std::abs(fFrequency[i]));
    }
#endif

    return true;
}

bool CP::TElectronicsResponse::Calculate(const CP::TEventContext& context,
                                         const CP::TChannelId channel) {

    // Check for all the things that might have changed.
    if (channel.IsMCChannel() != fChannelId.IsMCChannel()) {
        fMustRecalculate = true;
    }

    if (channel.IsTPCChannel() != fChannelId.IsTPCChannel()) {
        fMustRecalculate = true;
    }

    if (context != fEventContext) {
        fMustRecalculate = true;
    }
        
    fEventContext = context;
    fChannelId = channel;
    return true;
}
