#include "TClusterCalib.hxx"
#include "TPulseCalib.hxx"
#include "TPulseDeconvolution.hxx"

#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TRuntimeParameters.hxx>
#include <HEPUnits.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>

#include <memory>

CP::TClusterCalib::TClusterCalib() {
    fDigitStep 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.digitization.step");

    fPulseLength
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.digitization.pulse");

    fResponseLength
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.digitization.response");

    fSampleCount = (fPulseLength+fResponseLength)/fDigitStep;
    fSampleCount = 2*(1+fSampleCount/2);

    fCalibrate = new TPulseCalib();

    fDeconvolution = new TPulseDeconvolution(fSampleCount);

}

CP::TClusterCalib::~TClusterCalib() {}

bool CP::TClusterCalib::operator()(CP::TEvent& event) {
    CP::THandle<CP::TDigitContainer> pmt
        = event.Get<CP::TDigitContainer>("~/digits/pmt");
    CP::THandle<CP::TDigitContainer> drift
        = event.Get<CP::TDigitContainer>("~/digits/drift");

    CaptLog("Process " << event.GetContext());
    
    if (!pmt) {
        CaptLog("No PMT signals for this event " << event.GetContext());
        return false;
    }
    
    if (!drift) {
        CaptLog("No drift signals for this event " << event.GetContext());
        return false;
    }
    
    // Calibrate the PMT pulses.
    for (CP::TDigitContainer::const_iterator d = pmt->begin();
         d != pmt->end(); ++d) {
        const CP::TPulseDigit* pulse = dynamic_cast<const CP::TPulseDigit*>(*d);
        if (!pulse) {
            CaptError("Non-pulse in PMT digits");
            continue;
        }
        std::auto_ptr<CP::TCalibPulseDigit> calib((*fCalibrate)(pulse));
        std::auto_ptr<CP::TCalibPulseDigit> deconv((*fDeconvolution)(*calib));
        
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* calibHist 
            = new TH1F((calib->GetChannelId().AsString()+"-calib").c_str(),
                       ("Calibration for " 
                        + calib->GetChannelId().AsString()).c_str(),
                       calib->GetSampleCount(),
                       0.0, 1.*calib->GetSampleCount());
        for (std::size_t i = 0; i<calib->GetSampleCount(); ++i) {
            calibHist->SetBinContent(i+1,calib->GetSample(i));
        }
#endif
        
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* deconvHist 
            = new TH1F((deconv->GetChannelId().AsString()+"-deconv").c_str(),
                       ("Deconvolution for " 
                        + deconv->GetChannelId().AsString()).c_str(),
                       deconv->GetSampleCount(),
                       0.0, 1.*deconv->GetSampleCount());
        for (std::size_t i = 0; i<deconv->GetSampleCount(); ++i) {
            deconvHist->SetBinContent(i+1,deconv->GetSample(i));
        }
#endif
    }
        // Calibrate the drift pulses.
    for (CP::TDigitContainer::const_iterator d = drift->begin();
         d != drift->end(); ++d) {
        const CP::TPulseDigit* pulse = dynamic_cast<const CP::TPulseDigit*>(*d);
        if (!pulse) {
            CaptError("Non-pulse in drift digits");
            continue;
        }

        std::auto_ptr<CP::TCalibPulseDigit> calib((*fCalibrate)(pulse));
        std::auto_ptr<CP::TCalibPulseDigit> deconv((*fDeconvolution)(*calib));

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* calibHist 
            = new TH1F((calib->GetChannelId().AsString()+"-calib").c_str(),
                       ("Calibration for " 
                        + calib->GetChannelId().AsString()).c_str(),
                       calib->GetSampleCount(),
                       0.0, 1.*calib->GetSampleCount());
        for (std::size_t i = 0; i<calib->GetSampleCount(); ++i) {
            calibHist->SetBinContent(i+1,calib->GetSample(i));
        }
#endif
        
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* deconvHist 
            = new TH1F((deconv->GetChannelId().AsString()+"-deconv").c_str(),
                       ("Deconvolution for " 
                        + deconv->GetChannelId().AsString()).c_str(),
                       deconv->GetSampleCount(),
                       0.0, 1.*deconv->GetSampleCount());
        for (std::size_t i = 0; i<deconv->GetSampleCount(); ++i) {
            deconvHist->SetBinContent(i+1,deconv->GetSample(i));
        }
#endif
    }
    
    return true;
}
