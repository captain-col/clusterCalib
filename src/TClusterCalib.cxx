#include "TClusterCalib.hxx"
#include "TPulseCalib.hxx"
#include "TPulseDeconvolution.hxx"
#include "TPMTMakeHits.hxx"
#include "TWireSpectrum.hxx"
#include "TWirePeaks.hxx"

#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TRuntimeParameters.hxx>
#include <HEPUnits.hxx>
#include <CaptGeomId.hxx>
#include <TUnitsTable.hxx>
#include <TChannelCalib.hxx>
#include <TChannelId.hxx>
#include <TTPCChannelId.hxx>
#include <TPDSChannelId.hxx>
#include <TDataHit.hxx>

#include <TChannelInfo.hxx>

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <memory>
#include <limits>

// The next two includes are for debugging against MC input.  The are needed
// to fill some diagnostic histograms.
#include <TPulseMCDigit.hxx>
#include <TMCChannelId.hxx>

CP::TClusterCalib::TClusterCalib() {

    double pulseLength
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.digitization.pulse");

    double responseLength
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.digitization.response");

    // Estimate the minimum number of samples that will be used.  This is only
    // needed for initialization and will be overridden by the real value
    // later.
    fSampleCount = (pulseLength+responseLength)/(500*unit::ns);
    fSampleCount = 2*(1+fSampleCount/2);

    fCalibrate = new TPulseCalib();

    fDeconvolution = new TPulseDeconvolution(fSampleCount);

    fSaveCalibratedPulses = false;
    fApplyDriftCalibration = false;
    fApplyEfficiencyCalibration = true;
}

CP::TClusterCalib::~TClusterCalib() {}

bool CP::TClusterCalib::operator()(CP::TEvent& event) {
    CaptLog("Process " << event.GetContext());
    CP::TChannelInfo::Get().SetContext(event.GetContext());

    ///////////////////////////////////////////////////////////////////////
    // PMT Calibration
    ///////////////////////////////////////////////////////////////////////
    CP::THandle<CP::TDigitContainer> pmt
        = event.Get<CP::TDigitContainer>("~/digits/pmt");

    std::auto_ptr<CP::THitSelection> pmtHits(new CP::THitSelection("pmt"));
    CP::TPMTMakeHits makePMTHits;

    if (pmt) {
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
        }
    }

#define PMT_TRIGGER_ON_TPC_CHANNEL 
#ifdef PMT_TRIGGER_ON_TPC_CHANNEL
    ///////////////////////////////////////////////////////////////////////
    // Insert fact PMT hits for the PDS trigger signals.  These are in PMT 
    // 1000 (for the geomID).
    if (event.Get<CP::TDigitContainer>("~/digits/drift")) {
        CP::THandle<CP::TDigitContainer> drift
            = event.Get<CP::TDigitContainer>("~/digits/drift");
        for (std::size_t d = 0; d < drift->size(); ++d) {
            const CP::TPulseDigit* pulse
                = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
            if (!pulse) {
                CaptError("Non-pulse in drift digits");
                continue;
            }
            if (pulse->GetChannelId() != CP::TTPCChannelId(2,7,27)) continue;
            CP::TChannelCalib calib;
            double timeOffset = calib.GetTimeConstant(pulse->GetChannelId(),0);
            double digitStep = calib.GetTimeConstant(pulse->GetChannelId(),1);
            CP::TDigitProxy proxy(*drift,d);
            for(int s = 1; s<pulse->GetSampleCount(); ++s) {
                int delta = pulse->GetSample(s-1) - pulse->GetSample(s);
                if (delta > 100) {
                    std::cout << "Trigger w/ delta " << delta <<"  at " << s
                              << std::endl;
                    CP::TWritableDataHit hit;
                    hit.SetGeomId(CP::GeomId::Captain::Photosensor(1000));
                    hit.SetDigit(proxy);
                    hit.SetCharge(1);
                    hit.SetChargeUncertainty(1);
                    hit.SetTime(digitStep*s + timeOffset);
                    hit.SetTimeUncertainty(500*unit::ns);
                    hit.SetTimeRMS(500*unit::ns);
                    pmtHits->push_back(
                        CP::THandle<CP::TDataHit>(new CP::TDataHit(hit)));
                }
            }
        }
    }
#endif
    
    if (pmtHits->size() > 0) {
        // If there are any PMT hits found, add them to the output.
        CP::THandle<CP::TDataVector> hits
            = event.Get<CP::TDataVector>("~/hits");
        hits->AddDatum(pmtHits.release());
    }

    ///////////////////////////////////////////////////////////////////////
    // Find the event time zero.
    ///////////////////////////////////////////////////////////////////////
    double pmtT0 = 0.0;
    CP::THandle<CP::THitSelection> pmtSelection = event.GetHits("pmt");
    if (pmtSelection) {
        std::vector<double> pmtTimes;
        for (CP::THitSelection::iterator p = pmtSelection->begin();
             p != pmtSelection->end(); ++p) {
            pmtTimes.push_back((*p)->GetTime());
        }
        std::sort(pmtTimes.begin(), pmtTimes.end());
        int maxHits = 0;
        for (std::vector<double>::iterator t = pmtTimes.begin();
             t != pmtTimes.end(); ++t) {
            std::vector<double>::iterator h = t;
            while (++h != pmtTimes.end()) {
                if (*h - *t > 2*unit::microsecond) break;
            }
            int dh = h - t;
            if (dh > maxHits) {
                maxHits = dh;
                pmtT0 = *t;
            }
        }
    }
#ifdef USE_PMT_TIME_ZERO
    double t0 = pmtT0;
#else
    double t0 = 0.0;
#endif
    
    ///////////////////////////////////////////////////////////////////////
    /// Drift Calibration
    ///////////////////////////////////////////////////////////////////////

    // Check to see if the deconvoluted digits are going to be saved.  If they
    // are, then create the container to hold them.
    if (fSaveCalibratedPulses) {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-deconv");
        dv->AddDatum(dg);
    }

    CP::THandle<CP::TDigitContainer> drift
        = event.Get<CP::TDigitContainer>("~/digits/drift");
    if (!drift) {
        CaptLog("No drift signals for this event " << event.GetContext());
        return false;
    }

    std::auto_ptr<CP::THitSelection> driftHits(new CP::THitSelection("drift"));
#ifdef USE_WIRE_SPECTRUM
    CP::TWireSpectrum makeWireHits(fApplyDriftCalibration,
                                   fApplyEfficiencyCalibration);
#else
    CP::TWirePeaks makeWireHits(fApplyDriftCalibration,
                                fApplyEfficiencyCalibration);
#endif

    // Calibrate the drift pulses.
    for (std::size_t d = 0; d < drift->size(); ++d) {
        const CP::TPulseDigit* pulse
            = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
        if (!pulse) {
            CaptError("Non-pulse in drift digits");
            continue;
        }

        CP::TChannelCalib channelCalib;
        if (!channelCalib.IsGoodChannel(pulse->GetChannelId())) {
            continue;
        }
        
        CP::TGeometryId pulseGeom
            = CP::TChannelInfo::Get().GetGeometry(pulse->GetChannelId());

        if (!pulseGeom.IsValid()) {
            continue;
        }

        if (d%100 == 0) {
            CaptLog("Calibrate " << pulse->GetChannelId().AsString()
                     << " " << std::setw(40) << pulseGeom << std::setw(0));
        }

        CP::TDigitProxy proxy(*drift,d);
        std::auto_ptr<CP::TCalibPulseDigit> calib((*fCalibrate)(proxy));

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        const CP::TPulseDigit* raw = proxy.As<const CP::TPulseDigit>();
        TH1F* rawHist
            = new TH1F((raw->GetChannelId().AsString()+"-raw").c_str(),
                       ("Raw ADC for "
                        + raw->GetChannelId().AsString()).c_str(),
                       raw->GetSampleCount(),
                       raw->GetFirstSample(),
                       raw->GetFirstSample()+raw->GetSampleCount());
        for (std::size_t i = 0; i<raw->GetSampleCount(); ++i) {
            rawHist->SetBinContent(i+1,raw->GetSample(i));
        }
#endif

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
        TH1F* projHist
            = new TH1F((calib->GetChannelId().AsString()+"-proj").c_str(),
                       ("Projection for "
                        + calib->GetChannelId().AsString()).c_str(),
                       100, -50000.0, 50000.0);
        for (std::size_t i = 0; i<calib->GetSampleCount(); ++i) {
            projHist->Fill(calib->GetSample(i));
        }
#endif

#define STANDARD_HISTOGRAM
#ifdef STANDARD_HISTOGRAM
#undef STANDARD_HISTOGRAM
        static TH1F* gClusterCalibXPedestal = NULL;
        static TH1F* gClusterCalibVPedestal = NULL;
        static TH1F* gClusterCalibUPedestal = NULL;
        if (!gClusterCalibXPedestal) {
            gClusterCalibXPedestal
                = new TH1F("clusterCalibXPedestal",
                           "Pedestal for the X wires",
                           5000, 0.0, 5000);
        }
        if (!gClusterCalibUPedestal) {
            gClusterCalibUPedestal
                = new TH1F("clusterCalibUPedestal",
                           "Pedestal for the U wires",
                           5000, 0.0, 5000);
        }
        if (!gClusterCalibVPedestal) {
            gClusterCalibVPedestal
                = new TH1F("clusterCalibVPedestal",
                           "Pedestal for the V wires",
                           5000, 0.0, 5000);
        }

        {
            TH1F* hist = NULL;
            switch (CP::GeomId::Captain::GetWirePlane(pulseGeom)) {
            case 0: hist = gClusterCalibXPedestal; break;
            case 1: hist = gClusterCalibVPedestal; break;
            case 2: hist = gClusterCalibUPedestal; break;
            default: std::exit(1);
            }
            hist->Fill(fCalibrate->GetPedestal());
        }
#endif

#define STANDARD_HISTOGRAM
#ifdef STANDARD_HISTOGRAM
#undef STANDARD_HISTOGRAM
        static TH1F* gClusterCalibXSigma = NULL;
        static TH1F* gClusterCalibVSigma = NULL;
        static TH1F* gClusterCalibUSigma = NULL;
        if (!gClusterCalibXSigma) {
            gClusterCalibXSigma
                = new TH1F("clusterCalibXSigma",
                           "Sigma for the X wires",
                           200, 0.0, 100);
        }
        if (!gClusterCalibUSigma) {
            gClusterCalibUSigma
                = new TH1F("clusterCalibUSigma",
                           "Sigma for the U wires",
                           200, 0.0, 100);
        }
        if (!gClusterCalibVSigma) {
            gClusterCalibVSigma
                = new TH1F("clusterCalibVSigma",
                           "Sigma for the V wires",
                           200, 0.0, 100);
        }

        {
            TH1F* hist = NULL;
            switch (CP::GeomId::Captain::GetWirePlane(pulseGeom)) {
            case 0: hist = gClusterCalibXSigma; break;
            case 1: hist = gClusterCalibVSigma; break;
            case 2: hist = gClusterCalibUSigma; break;
            default: std::exit(1);
            }
            hist->Fill(fCalibrate->GetSigma());
        }
#endif

#define STANDARD_HISTOGRAM
#ifdef STANDARD_HISTOGRAM
#undef STANDARD_HISTOGRAM
        static TH1F* gClusterCalibXGaussian = NULL;
        static TH1F* gClusterCalibVGaussian = NULL;
        static TH1F* gClusterCalibUGaussian = NULL;
        if (!gClusterCalibXGaussian) {
            gClusterCalibXGaussian
                = new TH1F("clusterCalibXGaussian",
                           "Gaussian sigma for the X wires",
                           200, 0.0, 100);
        }
        if (!gClusterCalibUGaussian) {
            gClusterCalibUGaussian
                = new TH1F("clusterCalibUGaussian",
                           "Gaussian sigma for the U wires",
                           200, 0.0, 100);
        }
        if (!gClusterCalibVGaussian) {
            gClusterCalibVGaussian
                = new TH1F("clusterCalibVGaussian",
                           "Gaussian sigma for the V wires",
                           200, 0.0, 100);
        }

        {
            TH1F* hist = NULL;
            switch (CP::GeomId::Captain::GetWirePlane(pulseGeom)) {
            case 0: hist = gClusterCalibXGaussian; break;
            case 1: hist = gClusterCalibVGaussian; break;
            case 2: hist = gClusterCalibUGaussian; break;
            default: std::exit(1);
            }
            hist->Fill(fCalibrate->GetGaussianSigma());
        }
#endif

        std::auto_ptr<CP::TCalibPulseDigit> deconv((*fDeconvolution)(*calib));
        if (!deconv.get()) continue;

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
#endif

        double wireCharge
            = makeWireHits(*driftHits,*calib,*deconv,
                           t0,
                           fDeconvolution->GetBaselineSigma(),
                           fDeconvolution->GetSampleSigma());

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
        // Check if this is an MC digit, and if it is, then check the overall
        // charge normalization.
        const CP::TPulseMCDigit* mcPulse
            = dynamic_cast<const CP::TPulseMCDigit*>((*drift)[d]);
        if (mcPulse) {
            static TH1F* gClusterCalibXTrueVFound = NULL;
            static TH1F* gClusterCalibUVTrueVFound = NULL;
            if (!gClusterCalibXTrueVFound) {
                gClusterCalibXTrueVFound = new TH1F(
                    "clusterCalibXTrueVFound",
                    "True and Calibrated Charge Fraction Diff. X wires",
                    300, -1.0, 1.0);
            }

            if (!gClusterCalibUVTrueVFound) {
                gClusterCalibUVTrueVFound = new TH1F(
                    "clusterCalibUVTrueVFound",
                    "True and Calibrated Charge Fraction Diff. UV wires",
                    300, -1.0, 1.0);
            }

            if (mcPulse->GetInformation().size()>0
                && mcPulse->GetInformation().front() > 0) {
                CP::TMCChannelId mcChan(mcPulse->GetChannelId());
                double trueCharge = mcPulse->GetInformation().front();
                double delta = (wireCharge - trueCharge)/trueCharge;
                TH1F* hist = NULL;
                if (mcChan.GetSequence() == 0) hist = gClusterCalibXTrueVFound;
                if (mcChan.GetSequence() != 0) hist = gClusterCalibUVTrueVFound;
                if (hist) hist->Fill(delta);
            }
        }
#endif
        
        // Check to see if the drift-deconv digit container exists.  If it does
        // exist, then it will be filled with the deconvoluted digits.
        CP::THandle<CP::TDigitContainer> driftDeconv
            = event.Get<CP::TDigitContainer>("~/digits/drift-deconv");
        if (driftDeconv) driftDeconv->push_back(deconv.release());
    }

#define STANDARD_HISTOGRAM
#ifdef STANDARD_HISTOGRAM
#undef STANDARD_HISTOGRAM
    static TH1F* gClusterCalibXCharge = NULL;
    static TH1F* gClusterCalibVCharge = NULL;
    static TH1F* gClusterCalibUCharge = NULL;
    double minCharge = std::log10(4000);
    double maxCharge = std::log10(50*unit::fC);
    if (!gClusterCalibXCharge) {
        gClusterCalibXCharge = new TH1F("clusterCalibXCharge",
                                        "Calibrated Charge for the X wires",
                                        100, minCharge, maxCharge);
    }
    if (!gClusterCalibUCharge) {
        gClusterCalibUCharge = new TH1F("clusterCalibUCharge",
                                        "Calibrated Charge for the U wires",
                                        100, minCharge, maxCharge);
    }
    if (!gClusterCalibVCharge) {
        gClusterCalibVCharge = new TH1F("clusterCalibVCharge",
                                        "Calibrated Charge for the V wires",
                                        100, minCharge, maxCharge);
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
        double q = std::log10((*h)->GetCharge());
        q = std::min(std::max(minCharge, q), maxCharge);
        hist->Fill(q);
    }
#endif

#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    static TH2F* gClusterCalibXTime = NULL;
    static TProfile* gClusterCalibXTimeProfile = NULL;
    if (!gClusterCalibXTime) {
        gClusterCalibXTime = new TH2F("clusterCalibXTime",
                                      "Charge of X hits vs Time",
                                      100,0.0,1*unit::ms,
                                      100,0.0,50000);
        gClusterCalibXTimeProfile = new TProfile("clusterCalibXTimeProfile",
                                      "Mean Charge of X hits vs Time",
                                      100,0.0,1*unit::ms,
                                      5000.0,50000.0);
    }
    for (CP::THitSelection::iterator h = driftHits->begin();
         h != driftHits->end(); ++h) {
        TGeometryId id = (*h)->GetGeomId();
        if (CP::GeomId::Captain::GetWirePlane(id) != 0) continue;
        int neighbors = 0;
        for (CP::THitSelection::iterator i = driftHits->begin();
             i != driftHits->end(); ++i) {
            TGeometryId wid = (*i)->GetGeomId();
            if (CP::GeomId::Captain::GetWirePlane(wid) != 0) continue;
            if (wid == id) continue;
            int dWire = CP::GeomId::Captain::GetWireNumber(id)
                - CP::GeomId::Captain::GetWireNumber(wid);
            if (std::abs(dWire) > 1) continue;
            double delta = (*h)->GetTime() - (*i)->GetTime();
            if (std::abs(delta) > 2*unit::microsecond) continue;
            ++neighbors;
        }
        if (neighbors < 1) continue;
        gClusterCalibXTime->Fill((*h)->GetTime(),(*h)->GetCharge());
        gClusterCalibXTimeProfile->Fill((*h)->GetTime(),(*h)->GetCharge());
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
    CaptLog("X: "<< unit::AsString(xWireCharge,xWireUnc,"charge")
            << "   V: " << unit::AsString(vWireCharge,vWireUnc,"charge")
            << "   U: " << unit::AsString(uWireCharge,uWireUnc,"charge"));
    CaptNamedInfo("TClusterCalib",
                  "X-V/sigma: " << xv
                  << "   X-U/sigma: " << xu
                  << "   V-U/sigma: " << vu);

#define FILL_HISTOGRAM
#ifdef FILL_HISTOGRAM
#undef FILL_HISTOGRAM
    static TH1F* gClusterCalibURatio = NULL;
    static TH1F* gClusterCalibVRatio = NULL;
    static TH1F* gClusterCalibXRatio = NULL;
    if (!gClusterCalibVRatio) {
        gClusterCalibURatio = new TH1F("clusterCalibURatio",
                                      "Ratio of U charge to X charge",
                                       100,0.0,10);
        gClusterCalibVRatio = new TH1F("clusterCalibVRatio",
                                      "Ratio of V charge to charge",
                                       100,0.0,10);
        gClusterCalibXRatio = new TH1F("clusterCalibXRatio",
                                      "Ratio of X charge to U-V average",
                                       100,0.0,10);
    }
    if (xWireCharge > 0) {
        gClusterCalibURatio->Fill(uWireCharge/xWireCharge);
        gClusterCalibVRatio->Fill(vWireCharge/xWireCharge);
    }
    if (uWireCharge + vWireCharge > 0) {
        gClusterCalibXRatio->Fill(2.0*xWireCharge/(uWireCharge+vWireCharge));
    }
#endif

    if (driftHits->size() > 0) {
        // Add the drift hits to the output.
        CP::THandle<CP::TDataVector> hits
            = event.Get<CP::TDataVector>("~/hits");
        hits->AddDatum(driftHits.release());
    }

    return true;
}
