#include "TClusterCalib.hxx"
#include "TPulseCalib.hxx"
#include "TPulseDeconvolution.hxx"
#include "TPMTMakeHits.hxx"
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
    fSaveDecorrelatedPulses = false;
    fSaveDeconvolvedPulses = false;
    fCalibrateAllChannels = false;
    fApplyEfficiencyCalibration = true;
    fRemoveCorrelatedPedestal = true;
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

    std::unique_ptr<CP::THitSelection> pmtHits(new CP::THitSelection("pmt"));
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
            std::unique_ptr<CP::TCalibPulseDigit> calib((*fCalibrate)(proxy));
            makePMTHits(*pmtHits,*calib);
        }
    }

#define PMT_TRIGGER_ON_TPC_CHANNEL 
#ifdef PMT_TRIGGER_ON_TPC_CHANNEL
    ///////////////////////////////////////////////////////////////////////
    // Insert fake PMT hits for the PDS trigger signals.  These are in PMT 
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
            for(std::size_t s = 1; s<pulse->GetSampleCount(); ++s) {
                int delta = pulse->GetSample(s-1) - pulse->GetSample(s);
                if (delta > 100) {
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

    double t0 = 0.0;
    
    ///////////////////////////////////////////////////////////////////////
    /// Drift Calibration
    ///////////////////////////////////////////////////////////////////////

    CP::THandle<CP::TDigitContainer> drift
        = event.Get<CP::TDigitContainer>("~/digits/drift");
    if (!drift) {
        CaptLog("No drift signals for this event " << event.GetContext());
        return false;
    }

    CP::THandle<CP::TDigitContainer> driftCalib;

    driftCalib = CalibrateChannels(event, drift);
    
    driftCalib = RemoveCorrelatedPedestal(event,driftCalib);

    driftCalib = DeconvolveSignals(event,driftCalib);
    
    // Create a hit selection for hits that are found.
    std::unique_ptr<CP::THitSelection> driftHits(
        new CP::THitSelection("drift"));
    
    CP::TWirePeaks wirePeaks(fApplyEfficiencyCalibration);

    // Loop over all of the calibrated pulse digits and find the hits.
    for (std::size_t d = 0; d < driftCalib->size(); ++d) {
        const CP::TCalibPulseDigit* calib
            = dynamic_cast<const CP::TCalibPulseDigit*>((*driftCalib)[d]);
        if (d%100 == 0) {
            CaptInfo("Make Hits " << calib->GetChannelId().AsString());
        }

        // Find any peaks in the deconvoluted pulse.
        wirePeaks(*driftHits,*calib,t0);
    }

    if (driftHits->size() > 0) {
        // Add the drift hits to the output.
        CP::THandle<CP::TDataVector> hits
            = event.Get<CP::TDataVector>("~/hits");
        hits->AddDatum(driftHits.release());
    }

    return true;
}

CP::THandle<CP::TDigitContainer> CP::TClusterCalib::RemoveCorrelatedPedestal(
    CP::TEvent& event, CP::THandle<CP::TDigitContainer> driftCalib) {
    if (!fRemoveCorrelatedPedestal) return driftCalib;
    
    CaptInfo("Remove Correlations " << event.GetContext());

    // Find the RMS for each wire.  The pedestal is already removed so the
    // expected (mean) value is zero.
    static std::vector<double> sigmas;
    if (sigmas.size() != driftCalib->size()) {
        sigmas.resize(driftCalib->size());
    }
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        CP::TCalibPulseDigit* calib1
            = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        double sigma = 0.0;
        for (std::size_t i=0; i<(std::size_t) calib1->GetSampleCount(); ++i) {
            double p = calib1->GetSample(i);
            sigma += p*p;
        }
        sigma /= calib1->GetSampleCount();
        sigmas[d1] = std::sqrt(sigma);
    }

    // Fill a "2d" array of wire to wire correlations.  The array indexing is
    // done by hand so I can use a vector. 
    static std::vector<double> correlations;
    if (correlations.size() != (driftCalib->size()*driftCalib->size())) {
        correlations.resize(driftCalib->size()*driftCalib->size());
    }
    std::size_t pulseSamples = 0;
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        if (sigmas[d1]<1.0) continue;
        CP::TCalibPulseDigit* calib1
            = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        if (d1%100 == 0) {
            CaptInfo("Find Correlations " << calib1->GetChannelId().AsString());
        }
        pulseSamples = std::max(pulseSamples,calib1->GetSampleCount());
        for (std::size_t d2 = d1+1; d2 < driftCalib->size(); ++d2) {
            if (sigmas[d2]<1.0) continue;
            if (d1 == d2) continue;
            CP::TCalibPulseDigit* calib2
                = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d2]);
            double corr12 = 0.0;
            int steps = 0;
            for (std::size_t i=0; i<calib2->GetSampleCount(); i += 4) {
                double p1 = calib1->GetSample(i);
                double p2 = calib2->GetSample(i);
                corr12 += p1*p2;
                ++steps;
            }
            corr12 /= steps;
            corr12 /= (sigmas[d1]*sigmas[d2]);
            correlations[d1*driftCalib->size()+d2] = corr12;
            correlations[d2*driftCalib->size()+d1] = corr12;
        }
    }

    // Find the average "correlated" pedestal for each wire.  For a wire, this
    // is the correlation weighted average of the samples in all of the other
    // wires.  This needs to be calculated separately so when the new pedestal
    // subtraction is done, the subtraction does affect the pedestal
    // calculation.
    static std::vector<double> pedestals;
    if (pedestals.size() != driftCalib->size()*pulseSamples) {
        pedestals.resize(driftCalib->size()*pulseSamples);
    }
    static std::vector<double> weights;
    if (weights.size() != driftCalib->size()) {
        weights.resize(driftCalib->size());
    }
    static std::vector<double> wires;
    if (wires.size() != driftCalib->size()) {
        wires.resize(driftCalib->size());
    }
    int correlatedWires = 0;
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        if (sigmas[d1]<1.0) continue;
        CP::TCalibPulseDigit* calib1
            = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        if (d1%100 == 0) {
            CaptInfo("Estimate correlated pedestal "
                    << calib1->GetChannelId().AsString());
        }
        for (std::size_t d2 = 0; d2 < driftCalib->size(); ++d2) {
            if (sigmas[d2]<1.0) continue;
            if (d1 == d2) continue;
            CP::TCalibPulseDigit* calib2
                = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d2]);
            double w = correlations[d1*driftCalib->size()+d2];
            if (w < 0.6) continue;
            w = w*w*w; // Favor channels with high correlations.
            for (std::size_t i = 0; i< calib2->GetSampleCount(); ++i) {
                // this is going to take the average (weighted by the
                // correlation) of the scaled baseline fluctuation.
                double v = sigmas[d1]*calib2->GetSample(i)/sigmas[d2];
                pedestals[d1*pulseSamples+i] += w*v;
            }
            weights[d1] += std::abs(w);
            wires[d1] += 1.0;
        } 
        if (weights[d1] > 0.0) ++correlatedWires;
    }
    CaptInfo("Wires with strong correlations: " << correlatedWires);

    // Find the averaged pedestals for each wire.
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        const CP::TCalibPulseDigit* calib
            = dynamic_cast<const CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        // Strongly deweight pedestals for wires that have few correlated
        // wires.
        double weight = (wires[d1]-100.0)/10.0;
        if (weight < -10.0) weight = 0.0;
        else if (weight > 10.0) weight = 1.0;
        else weight = 1.0-1.0/(1.0 + std::exp(weight));
        for (std::size_t i = 0; i< calib->GetSampleCount(); ++i) {
            if (weights[d1] < 1.0) {
                pedestals[d1*pulseSamples+i] = 0.0;
                continue;
            }
            pedestals[d1*pulseSamples+i] *= weight/weights[d1];
        }
    }

    // Create a handle for the current calibrated digits being worked on.
    // This will keep track of the current stage of the calibration as it
    // progresses.
    if (fSaveDecorrelatedPulses) {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-correl");
        dv->AddDatum(dg);
    }
    else {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-correl");
        dv->AddTemporary(dg);
    }
    CP::THandle<CP::TDigitContainer> driftCorrel
        = event.Get<CP::TDigitContainer>("~/digits/drift-correl");

    // Subtract the average correlated pedestal from each channel
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        CP::TCalibPulseDigit* calib
            = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        std::unique_ptr<CP::TCalibPulseDigit> correl(
            new CP::TCalibPulseDigit(*calib));
        // Remove the correlated pedestal from the channel, but only for the
        // induction planes.  The collection plane is just copied.
        for (std::size_t i = 0; i< calib->GetSampleCount(); ++i) {
            double v = calib->GetSample(i);
            double p = pedestals[d1*pulseSamples+i];
            CP::TGeometryId pulseGeom
                = CP::TChannelInfo::Get().GetGeometry(calib->GetChannelId());
            if (CP::GeomId::Captain::IsXWire(pulseGeom)) {
                p = 0.0;
            }
            correl->SetSample(i,v-p);
        }
        driftCorrel->push_back(correl.release());
    }

    return driftCorrel;
}

CP::THandle<CP::TDigitContainer> CP::TClusterCalib::DeconvolveSignals(
    CP::TEvent& event, CP::THandle<CP::TDigitContainer> driftCalib) {

    CaptInfo("Deconvolve Event " << event.GetContext());

    // Create a container for the deconvoluted TCalibPulseDigits and check to
    // see if the deconvoluted digits are going to be saved.  If they are,
    // then create a permanent container to hold them.  Otherwise, the
    // container is temporary.
    if (fSaveDeconvolvedPulses) {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-deconv");
        dv->AddDatum(dg);
    }
    else {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-deconv");
        dv->AddTemporary(dg);
    }
    CP::THandle<CP::TDigitContainer> driftDeconv
        = event.Get<CP::TDigitContainer>("~/digits/drift-deconv");

    // Loop over all of the calibrated pulse digits and deconvolve.  The
    // deconvolution is going to apply a Weiner filter.
    for (std::size_t d = 0; d < driftCalib->size(); ++d) {
        const CP::TCalibPulseDigit* calib
            = dynamic_cast<const CP::TCalibPulseDigit*>((*driftCalib)[d]);
        std::unique_ptr<CP::TCalibPulseDigit> deconv((*fDeconvolution)(*calib));
        if (!deconv.get()) continue;
        if (d%100 == 0) {
            CaptInfo("Deconvolve " << calib->GetChannelId().AsString());
        }

        // Add the deconvoluted digits to the event (remember, they might be
        // temporary).
        driftDeconv->push_back(deconv.release());
    }

    return driftDeconv;
}

CP::THandle<CP::TDigitContainer> CP::TClusterCalib::CalibrateChannels(
    CP::TEvent& event,
    CP::THandle<CP::TDigitContainer> drift) {
    
    CaptInfo("Calibrate Event " << event.GetContext());
    
    // Create a handle for the current calibrated digits being worked on.
    // This will keep track of the current stage of the calibration as it
    // progresses.
    if (fSaveCalibratedPulses) {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-calib");
        dv->AddDatum(dg);
    }
    else {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-calib");
        dv->AddTemporary(dg);
    }
    CP::THandle<CP::TDigitContainer> driftCalib
        = event.Get<CP::TDigitContainer>("~/digits/drift-calib");

    // Calibrate the drift pulses.  This loop should only apply the
    // calibration.  The resulting calibrated digits are pointed to in the
    // handle driftCalib.  This fills driftCalib for the first time.
    for (std::size_t d = 0; d < drift->size(); ++d) {
        const CP::TPulseDigit* pulse
            = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
        if (!pulse) {
            CaptError("Non-pulse in drift digits");
            continue;
        }

        CP::TChannelCalib channelCalib;

        // Skip channels that have bad calibrations.
        if (!channelCalib.IsGoodChannel(pulse->GetChannelId())) {
            continue;
        }

        
        CP::TGeometryId pulseGeom
            = CP::TChannelInfo::Get().GetGeometry(pulse->GetChannelId());

        if (!pulseGeom.IsValid() && !fCalibrateAllChannels) {
            continue;
        }

        if (d%100 == 0) {
            CaptInfo("Calibrate " << pulse->GetChannelId().AsString()
                     << " " << std::setw(40) << pulseGeom << std::setw(0));
        }

        CP::TDigitProxy proxy(*drift,d);
        std::unique_ptr<CP::TCalibPulseDigit> calib((*fCalibrate)(proxy));
        if (!calib.get()) {
            CaptError("Channel not calibrated ");
            continue;
        }

        // Add the calibrated digits to the event (remember, they are probably
        // temporary).
        driftCalib->push_back(calib.release());
    }

    return driftCalib;
}


