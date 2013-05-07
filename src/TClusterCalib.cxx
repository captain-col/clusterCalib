#include "TClusterCalib.hxx"
#include "TPulseCalib.hxx"
#include "TPulseDeconvolution.hxx"
#include "TPMTMakeHits.hxx"
#include "TWireMakeHits.hxx"

#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TRuntimeParameters.hxx>
#include <HEPUnits.hxx>

#include <TChannelInfo.hxx>

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
    CP::TChannelInfo::Get().SetContext(event.GetContext());    

    if (!pmt) {
        CaptLog("No PMT signals for this event " << event.GetContext());
        return false;
    }
    
    if (!drift) {
        CaptLog("No drift signals for this event " << event.GetContext());
        return false;
    }

    std::auto_ptr<CP::THitSelection> pmtHits(new CP::THitSelection("pmt"));

    CP::TPMTMakeHits makePMTHits;

    // Calibrate the PMT pulses and turn them into hits.
    for (std::size_t d = 0; d < pmt->size(); ++d) {
        const CP::TPulseDigit* pulse 
            = dynamic_cast<const CP::TPulseDigit*>((*pmt)[d]);
        if (!pulse) {
            CaptError("Non-pulse in PMT digits");
            continue;
        }
        CP::TDigitProxy proxy(*pmt,d);
        std::auto_ptr<CP::TCalibPulseDigit> calib((*fCalibrate)(proxy));
        std::auto_ptr<CP::TCalibPulseDigit> deconv((*fDeconvolution)(*calib));
        makePMTHits(*pmtHits,*deconv);

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* calibHist 
            = new TH1F((calib->GetChannelId().AsString()+"-calib").c_str(),
                       ("Calibration for " 
                        + calib->GetChannelId().AsString()).c_str(),
                       fSampleCount,
                       0.0, 1*fSampleCount);
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
        deconvHist->ShowPeaks(2.0,"",0.1);
#endif
    }
    if (pmtHits->size() > 0) {
        // Add the pmt hits to the output.
        CP::THandle<CP::TDataVector> hits
            = event.Get<CP::TDataVector>("~/hits");
        hits->AddDatum(pmtHits.release());
        
    }

    std::auto_ptr<CP::THitSelection> driftHits(new CP::THitSelection("drift"));
    CP::TWireMakeHits makeWireHits;

    // Calibrate the drift pulses.
    for (std::size_t d = 0; d < drift->size(); ++d) {
        const CP::TPulseDigit* pulse 
            = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
        if (!pulse) {
            CaptError("Non-pulse in drift digits");
            continue;
        }

        CP::TDigitProxy proxy(*drift,d);
        std::auto_ptr<CP::TCalibPulseDigit> calib((*fCalibrate)(proxy));
        std::auto_ptr<CP::TCalibPulseDigit> deconv((*fDeconvolution)(*calib));
        makeWireHits(*driftHits,*deconv);

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* calibHist 
            = new TH1F((calib->GetChannelId().AsString()+"-calib").c_str(),
                       ("Calibration for " 
                        + calib->GetChannelId().AsString()).c_str(),
                       fSampleCount,
                       0.0, 1*fSampleCount);
        for (std::size_t i = 0; i<calib->GetSampleCount(); ++i) {
            calibHist->SetBinContent(i+1,calib->GetSample(i));
        }
#endif
        
#ifdef  FILL_HISTOGRAM
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
        deconvHist->ShowPeaks(3.0,"",0.15);
#endif
    }
    if (driftHits->size() > 0) {
        // Add the drift hits to the output.
        CP::THandle<CP::TDataVector> hits
            = event.Get<CP::TDataVector>("~/hits");
        hits->AddDatum(driftHits.release());
    }

    
    return true;
}
