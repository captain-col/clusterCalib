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

#define USE_AVERAGE_SHAPE
#ifdef USE_AVERAGE_SHAPE
    Cache::iterator cache = fCache.begin();
#else
    Cache::iterator cache = fCache.find(fChannelId);
#endif
    if (cache != fCache.end()) {
        fResponse = cache->second.first;
        fFrequency = cache->second.second;
        return true;
    }
    
    TChannelCalib calib;

#ifdef USE_AVERAGE_SHAPE
    CaptError("Using Average Pulse Shape: "
              << " Peak: " << calib.GetAveragePulseShapePeakTime(fChannelId)
              << " Rise: " << calib.GetAveragePulseShapeRise(fChannelId)
              << " Fall: " << calib.GetAveragePulseShapeFall(fChannelId));
#endif
    
    // Fill the response function.  This explicitly normalizes so that the
    // pulse shaping is amplitude conserving (the pulse shaping for a
    // "delta-function" sample doesn't change the pulse height).
    double normalization = 0.0;
    for (std::size_t i=0; i<fResponse.size(); ++i) {
        double arg = 1.0*calib.GetTimeConstant(fChannelId)*i;
#ifdef USE_AVERAGE_SHAPE
        double v = calib.GetAveragePulseShape(fChannelId,arg);
#else
        double v = calib.GetPulseShape(fChannelId,arg);
#endif
        fResponse[i] = v;
        normalization = std::max(normalization,v);
    }
    if (normalization < 1E-20) return false;
    for (std::vector<Response>::iterator r = fResponse.begin();
         r != fResponse.end(); ++r) {
        (*r) /= normalization;
    }

    // Calculate the FFT of the response 
    int size = fResponse.size();
    std::unique_ptr<TVirtualFFT> fft(TVirtualFFT::FFT(1,&size,"R2C ES K"));
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

    fCache[fChannelId].first = fResponse;
    fCache[fChannelId].second = fFrequency;

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* elecResp = new TH1F(
        (fChannelId.AsString()+"-elec").c_str(),
        ("Electronics Response for " + fChannelId.AsString()).c_str(),
        fResponse.size(),
        0.0, fResponse.size()*calib.GetTimeConstant(fChannelId));
    for (std::size_t i=0; i<fResponse.size(); ++i) {
        elecResp->SetBinContent(i+1, std::abs(fResponse[i]));
    }
    TH1F* elecFreq = new TH1F(
        (fChannelId.AsString()+"-elecFFT").c_str(),
        ("FFT of the electronics Response for" + fChannelId.AsString()).c_str(),
        fFrequency.size(),
        0.0, (0.5/calib.GetTimeConstant(fChannelId))/unit::hertz);
    for (std::size_t i=0; i<fFrequency.size(); ++i) {
        elecFreq->SetBinContent(i+1, std::abs(fFrequency[i]));
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

    if (context.GetRun() != fEventContext.GetRun()) {
        fMustRecalculate = true;
    }
        
    if (!channel.IsMCChannel()
        && channel.IsTPCChannel()
        && channel != fChannelId) {
        fMustRecalculate = true;
    }
    fEventContext = context;
    fChannelId = channel;

    return true;
}
