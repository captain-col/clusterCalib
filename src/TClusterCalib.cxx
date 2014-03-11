#include "TClusterCalib.hxx"
#include "TPulseCalib.hxx"
#include "TPulseDeconvolution.hxx"
#include "TPMTMakeHits.hxx"
#include "TWireMakeHits.hxx"

#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TRuntimeParameters.hxx>
#include <HEPUnits.hxx>
#include <CaptGeomId.hxx>
#include <TUnitsTable.hxx>

#include <TChannelInfo.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>
#include <TSpectrum.h>

#include <memory>

namespace {
    TH1F* gClusterCalibXCharge = NULL;
    TH1F* gClusterCalibVCharge = NULL;
    TH1F* gClusterCalibUCharge = NULL;
}

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

    fSaveCalibratedPulses = true;
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
        makePMTHits(*pmtHits,*calib);

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* calibHist 
            = new TH1F((calib->GetChannelId().AsString()+"-calib").c_str(),
                       ("Calibration for " 
                        + calib->GetChannelId().AsString()).c_str(),
                       calib->GetSampleCount(),
                       calib->GetFirstSample(), calib->GetLastSample());
        for (std::size_t i = 0; i<calib->GetSampleCount(); ++i) {
            calibHist->SetBinContent(i+1,calib->GetSample(i));
        }
#endif
        
    }
    double t0 = 1E+50;
    if (pmtHits->size() > 0) {
        for (CP::THitSelection::iterator p = pmtHits->begin();
             p != pmtHits->end(); ++p) {
            t0 = std::min((*p)->GetTime(),t0);
        }
        // Add the pmt hits to the output.
        CP::THandle<CP::TDataVector> hits
            = event.Get<CP::TDataVector>("~/hits");
        hits->AddDatum(pmtHits.release());
        
    }

    std::auto_ptr<CP::THitSelection> driftHits(new CP::THitSelection("drift"));
    CP::TWireMakeHits makeWireHits;

    // Make a container to hold the deconvoluted digits.
    if (fSaveCalibratedPulses) {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-deconv");
        dv->AddDatum(dg);
        CP::THandle<CP::TDigitContainer> driftDeconv
            = event.Get<CP::TDigitContainer>("~/digits/drift-deconv");
    }
    CP::THandle<CP::TDigitContainer> driftDeconv
        = event.Get<CP::TDigitContainer>("~/digits/drift-deconv");

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
        makeWireHits(*driftHits,*deconv,t0,
                     fDeconvolution->GetBaselineSigma(),
                     fDeconvolution->GetSampleSigma());

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* calibHist 
            = new TH1F((calib->GetChannelId().AsString()+"-calib").c_str(),
                       ("Calibration for " 
                        + calib->GetChannelId().AsString()).c_str(),
                       calib->GetSampleCount(),
                       calib->GetFirstSample(), calib->GetLastSample());
        for (std::size_t i = 0; i<calib->GetSampleCount(); ++i) {
            calibHist->SetBinContent(i+1,calib->GetSample(i));
        }
#endif
        
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        TH1F* deconvHist 
            = new TH1F((deconv->GetChannelId().AsString()+"-deconv").c_str(),
                       ("Final deconvolution for " 
                        + deconv->GetChannelId().AsString()).c_str(),
                       deconv->GetSampleCount(),
                       deconv->GetFirstSample(), deconv->GetLastSample());
        for (std::size_t i = 0; i<deconv->GetSampleCount(); ++i) {
            deconvHist->SetBinContent(i+1,deconv->GetSample(i));
        }
        TSpectrum spectrum;
        TSpectrum::SetDeconIterations(30);
        spectrum.Search(deconvHist,2.0,"",0.01);
#endif

        if (driftDeconv) driftDeconv->push_back(deconv.release());
    }

#define FILL_HISTOGRAM
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    double maxCharge = 50;
    if (!gClusterCalibXCharge) {
        gClusterCalibXCharge = new TH1F("clusterCalibXCharge",
                                        "Calibrated Charge for the X wires",
                                        100, 0.0, maxCharge);
    }
    if (!gClusterCalibUCharge) {
        gClusterCalibUCharge = new TH1F("clusterCalibUCharge",
                                        "Calibrated Charge for the U wires",
                                        100, 0.0, maxCharge);
    }
    if (!gClusterCalibVCharge) {
        gClusterCalibVCharge = new TH1F("clusterCalibVCharge",
                                        "Calibrated Charge for the V wires",
                                        100, 0.0, maxCharge);
    }

    for (CP::THitSelection::iterator h = driftHits->begin();
         h != driftHits->end(); ++h) {
        TGeometryId id = (*h)->GetGeomId();
        TH1F* hist = NULL;
        switch (CP::GeomId::Captain::GetWirePlane(id)) {
        case 0: hist = gClusterCalibXCharge; break;
        case 1: hist = gClusterCalibVCharge; break;
        case 2: hist = gClusterCalibUCharge; break;
        default: std::exit(1);
        }
        hist->Fill((*h)->GetCharge()/(*h)->GetTimeRMS());
    }
#endif

    // Summarize the charges.
    double xWireCharge = 0.0;
    double xWireUnc = 0.0;
    double vWireCharge = 0.0;
    double vWireUnc = 0.0;
    double uWireCharge = 0.0;
    double uWireUnc = 0.0;
    for (CP::THitSelection::iterator h = driftHits->begin();
         h != driftHits->end(); ++h) {
        TGeometryId id = (*h)->GetGeomId();
        switch (CP::GeomId::Captain::GetWirePlane(id)) {
        case 0: {
            xWireCharge += (*h)->GetCharge(); 
            double u = (*h)->GetChargeUncertainty();
            xWireUnc += u*u;
            break;
        }
        case 1: {
            vWireCharge += (*h)->GetCharge(); 
            double u = (*h)->GetChargeUncertainty();
            vWireUnc += u*u;
            break;
        }
        case 2: {
            uWireCharge += (*h)->GetCharge(); 
            double u = (*h)->GetChargeUncertainty();
            uWireUnc += u*u;
            break;
        }
        default: std::exit(1);
        }
    }
    xWireUnc = std::sqrt(xWireUnc);
    vWireUnc = std::sqrt(vWireUnc);
    uWireUnc = std::sqrt(uWireUnc);
    double xv = (xWireCharge-vWireCharge);
    xv /= std::sqrt(xWireUnc*xWireUnc + vWireUnc*vWireUnc);
    double xu = (xWireCharge-uWireCharge);
    xu /= std::sqrt(xWireUnc*xWireUnc + uWireUnc*uWireUnc);
    double vu = (vWireCharge-uWireCharge);
    vu /= std::sqrt(vWireUnc*vWireUnc + uWireUnc*uWireUnc);
    CaptNamedInfo("TClusterCalib",
                  "X: "<< unit::AsString(xWireCharge,xWireUnc,"pe")
                  << "   V: " << unit::AsString(vWireCharge,vWireUnc,"pe")
                  << "   U: " << unit::AsString(uWireCharge,uWireUnc,"pe"));
    CaptNamedInfo("TClusterCalib",
                  "X-V/sigma: " << xv
                  << "   X-U/sigma: " << xu
                  << "   V-U/sigma: " << vu);
    
    if (driftHits->size() > 0) {
        // Add the drift hits to the output.
        CP::THandle<CP::TDataVector> hits
            = event.Get<CP::TDataVector>("~/hits");
        hits->AddDatum(driftHits.release());
    }

    return true;
}
