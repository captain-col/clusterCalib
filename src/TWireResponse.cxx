#include "TWireResponse.hxx"

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

CP::TWireResponse::TWireResponse() 
    : fMustRecalculate(true), fWireClass(kDeltaFunction) {}

CP::TWireResponse::TWireResponse(int size) 
    : fMustRecalculate(true), fWireClass(kDeltaFunction) {
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

    CP::TChannelCalib calib;

    for (std::size_t i=0; i<fResponse.size(); ++i) {
        fResponse[i] = 0;
        fFrequency[i] = Frequency(1,0); // Initialize, will be overwritten.
    }

    // Fill the response function.  This explicitly normalizes.
    switch (fWireClass) {
    case kDeltaFunction:
        CaptLog("Calculate delta function wire response");
        // Use a delta function in time for the wire response.  This is an
        // appropriate approximation for the collection wires (i.e. the X
        // plane).
        fResponse[0] = 1;
        break;
    case kTimeDerivative:
        CaptLog("Calculate differential function wire response");
        // Use the derivative of the charge w.r.t. time for the wire response.
        // This is appropriate for the induction wires (i.e. the U and V
        // planes).
        {
            double norm = 1.0/std::sqrt(2.0);
            // Correct for the bin size.
            norm *= calib.GetTimeConstant(fChannelId)/(1*unit::microsecond);
            fResponse[0] = -1.0/norm;
            fResponse[fResponse.size()-1] = 1.0/norm;
        }
        break;
    default:
        CaptError("Unknown wire class");
        exit(1);
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

    // The wire response is symetric around zero, so there is no offset.  That
    // means that the "zero" frequency is zero power, and the deconvolution
    // blows up with a divide by zero.  Fix that by adding an artificial
    // offset.  The offset will be removed later.
    fFrequency[0] = fFrequency[1];

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    TH1F* elecResp = new TH1F(
        (fChannelId.AsString()+"-wire").c_str(),
        ("Wire Response for " + fChannelId.AsString()).c_str(),
        fResponse.size(),
        0.0, 1.0*fResponse.size());
    for (std::size_t i=0; i<fResponse.size(); ++i) {
        elecResp->Fill(i+0.5, fResponse[i]);
    }
    TH1F* elecFreq = new TH1F(
        (fChannelId.AsString()+"-wireFFT").c_str(),
        ("FFT of the Wire Response for" + fChannelId.AsString()).c_str(),
        fFrequency.size(),
        0.0, 1.0*fFrequency.size());
    for (std::size_t i=0; i<fFrequency.size(); ++i) {
        elecFreq->Fill(i+0.5, std::abs(fFrequency[i]));
    }
#endif
    return true;
}

bool CP::TWireResponse::Calculate(const CP::TEventContext& context,
                                  const CP::TChannelId channel) {
    fEventContext = context;
    fChannelId = channel;


    CP::TChannelCalib calib;
    if (calib.IsBipolarSignal(channel) && fWireClass != kTimeDerivative) {
        fWireClass = kTimeDerivative;
        fMustRecalculate = true;
    }
    if (not calib.IsBipolarSignal(channel) && fWireClass != kDeltaFunction) {
        fWireClass = kDeltaFunction;
        fMustRecalculate = true;
    }

    return true;
}
