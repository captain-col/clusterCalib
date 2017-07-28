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
    fCalibrateAllChannels = false;
    fApplyDriftCalibration = false;
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

    // Create a container for the calibrated, but not deconvoluted
    // TCalibPulseDigits.
    {
        CP::THandle<CP::TDataVector> dv
            = event.Get<CP::TDataVector>("~/digits");
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift-calib");
        dv->AddTemporary(dg);
    }
    CP::THandle<CP::TDigitContainer> driftCalib
        = event.Get<CP::TDigitContainer>("~/digits/drift-calib");
    
    // Create a container for the deconvoluted TCalibPulseDigits and check to
    // see if the deconvoluted digits are going to be saved.  If they are,
    // then create a permanent container to hold them.  Otherwise, the
    // container is temporary.
    if (fSaveCalibratedPulses) {
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

    // Create a hit selection for hits that are found.
    std::unique_ptr<CP::THitSelection> driftHits(
        new CP::THitSelection("drift"));

    CP::TWirePeaks wirePeaks(fApplyDriftCalibration,
                                fApplyEfficiencyCalibration);

    // Calibrate the drift pulses.  This loop should only apply the
    // calibration, and fill a bunch of histograms.
    for (std::size_t d = 0; d < drift->size(); ++d) {
        const CP::TPulseDigit* pulse
            = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
        if (!pulse) {
            CaptError("Non-pulse in drift digits");
            continue;
        }

        CP::TChannelCalib channelCalib;
        if (!channelCalib.IsGoodChannel(pulse->GetChannelId())) {
            CaptLog("Bad Channel " <<pulse->GetChannelId().AsString());
            continue;
        }
        
        CP::TGeometryId pulseGeom
            = CP::TChannelInfo::Get().GetGeometry(pulse->GetChannelId());

        if (!pulseGeom.IsValid() && !fCalibrateAllChannels) {
            continue;
        }

        if (d%100 == 0) {
            CaptLog("Calibrate " << pulse->GetChannelId().AsString()
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

    if (fRemoveCorrelatedPedestal) RemoveCorrelatedPedestal(event);
    
    // Loop over all of the calibrated pulse digits and deconvolve.  The
    // deconvolution is going to apply a Weiner filter.
    for (std::size_t d = 0; d < driftCalib->size(); ++d) {
        const CP::TCalibPulseDigit* calib
            = dynamic_cast<const CP::TCalibPulseDigit*>((*driftCalib)[d]);
        std::unique_ptr<CP::TCalibPulseDigit> deconv((*fDeconvolution)(*calib));
        if (!deconv.get()) continue;
        if (d%100 == 0) {
            CaptLog("Deconvolve " << calib->GetChannelId().AsString());
        }

        // Find any peaks in the deconvoluted pulse.
        wirePeaks(*driftHits,*deconv,t0,fDeconvolution);
        
        // Add the deconvoluted digits to the event (remember, they might be
        // temporary).
        driftDeconv->push_back(deconv.release());
    }

    if (driftHits->size() > 0) {
        // Add the drift hits to the output.
        CP::THandle<CP::TDataVector> hits
            = event.Get<CP::TDataVector>("~/hits");
        hits->AddDatum(driftHits.release());
    }

    return true;
}

void CP::TClusterCalib::RemoveCorrelatedPedestal(CP::TEvent& event) {
    CP::THandle<CP::TDigitContainer> driftCalib
        = event.Get<CP::TDigitContainer>("~/digits/drift-calib");

    // Find the RMS for each wire.  The pedestal is already removed so the
    // expected (mean) value is zero.
    std::vector<double> sigmas(driftCalib->size());
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
    std::vector<double> correlations(driftCalib->size()*driftCalib->size());
    std::size_t pulseSamples = 0;
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        if (sigmas[d1]<1.0) continue;
        CP::TCalibPulseDigit* calib1
            = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        if (d1%100 == 0) {
            CaptLog("Find Correlations " << calib1->GetChannelId().AsString());
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
    std::vector<double> pedestals(driftCalib->size()*pulseSamples);
    std::vector<double> weights(driftCalib->size());
    int correlatedWires = 0;
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        if (sigmas[d1]<1.0) continue;
        CP::TCalibPulseDigit* calib1
            = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        if (d1%100 == 0) {
            CaptLog("Estimate correlated pedestal "
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
        } 
        if (weights[d1] > 0.0) ++correlatedWires;
    }
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        const CP::TCalibPulseDigit* calib1
            = dynamic_cast<const CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        for (std::size_t i = 0; i< calib1->GetSampleCount(); ++i) {
            if (weights[d1] < 0.1) pedestals[d1*pulseSamples+i] = 0.0;
            else pedestals[d1*pulseSamples+i] /= weights[d1];
        }
    }
    CaptLog("Wires with strong correlations: " << correlatedWires);

    // Subtract the average correlated pedestal from each channel
    for (std::size_t d1 = 0; d1 < driftCalib->size(); ++d1) {
        CP::TCalibPulseDigit* calib1
            = dynamic_cast<CP::TCalibPulseDigit*>((*driftCalib)[d1]);
        for (std::size_t i = 0; i< calib1->GetSampleCount(); ++i) {
            double v = calib1->GetSample(i);
            double p = pedestals[d1*pulseSamples+i];
            calib1->SetSample(i,v-p);
        }
    }

}
