#include "TMakeWireHit.hxx"
#include "TPulseDeconvolution.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TFADCHit.hxx>
#include <HEPUnits.hxx>
#include <THandle.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>
#include <CaptGeomId.hxx>
#include <TChannelCalib.hxx>

#include <TSpectrum.h>
#include <TH1F.h>
#include <TH2F.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <sstream>
#include <iostream>

CP::TMakeWireHit::TMakeWireHit(bool correctLifetime,
                               bool correctEfficiency) {
    fCorrectElectronLifetime = correctLifetime;
    fCorrectCollectionEfficiency = correctEfficiency;
}
CP::TMakeWireHit::~TMakeWireHit() { }

CP::THandle<CP::THit>
CP::TMakeWireHit::operator()(const CP::TCalibPulseDigit& digit,
                             double step, double t0,
                             const CP::TPulseDeconvolution* pulseDeconvolution,
                             std::size_t beginIndex, std::size_t endIndex,
                             bool split) {
    double charge = 0.0;
    double sample = 0.0;
    double sampleSquared = 0.0;
    double samples = 1.0;

    /// Protect against out of bounds problems.
    if (digit.GetSampleCount() <= endIndex) endIndex = digit.GetSampleCount()-1;
    if (endIndex < beginIndex) return CP::THandle<THit>();
    
    // Find the start and stop time.
    double startTime = beginIndex*step + digit.GetFirstSample();
    double stopTime = (endIndex+1)*step + digit.GetFirstSample();

    // Copy the samples into the hit.
    std::vector<double> sampleCharge(endIndex-beginIndex+1);
    for (std::size_t i = beginIndex; i<=endIndex; ++i) {
        sampleCharge[i-beginIndex] = digit.GetSample(i);
    }

    // Find the peak charge, average sample index and sample index squared for
    // the peak.  The average is charge weighted.
    for (std::size_t j=beginIndex; j<=endIndex; ++j) {
        double v = digit.GetSample(j);
        if (v <= 0) continue;
        charge += v;
        sample += v*j;
        sampleSquared += v*j*j;
        ++samples;
    }

    // Protect against empty regions.
    if (charge<1) return CP::THandle<CP::THit>();

    sample /= charge;
    sampleSquared /= charge;

    // Convert the sample index into a time.
    double time = (sample + 0.5)*step + digit.GetFirstSample();
    
    // Find the sample RMS, and then convert into a time RMS.
    double rms = step*std::sqrt(sampleSquared - sample*sample + 1.0);

#define FWHM_OVERRIDE
#ifdef FWHM_OVERRIDE
    {
        // Estimate the RMS using the FWHM.
        std::vector<double>::iterator maxBin = sampleCharge.end();
        for (std::vector<double>::iterator s = sampleCharge.begin();
             s != sampleCharge.end(); ++s) {
            if (maxBin==sampleCharge.end() || *maxBin < *s) maxBin = s;
        }
        if (maxBin != sampleCharge.end()) {
            std::vector<double>::iterator lowBin = maxBin;
            while (lowBin != sampleCharge.begin()) {
                --lowBin;
                if (*lowBin < 0.5*(*maxBin)) break;
            }
            std::vector<double>::iterator hiBin = maxBin;
            while (hiBin != sampleCharge.end()) {
                if (*hiBin < 0.5*(*maxBin)) break;
                if (hiBin+1 == sampleCharge.end()) break;
                ++hiBin;
            }
            double lowVal=(0.5*(*maxBin)-(*lowBin))/(*(lowBin+1)-(*lowBin));
            double hiVal = (0.5*(*maxBin)-(*hiBin))/(*(hiBin-1)-(*hiBin));
            int diff = hiBin-lowBin;
            // Convert FWHM into an RMS.  The 2.36 factor is the ratio between
            // the rms and the FWHM for a Gaussian peak.
            double bins = 1.0*(diff - lowVal - hiVal)/2.36;
            if (bins > 1) rms = bins*step;
        }
    }
#endif
    
#define TIME_OVERRIDE
#ifdef TIME_OVERRIDE
    {
        // Find the time by looking at samples near the peak.
        std::vector<double>::iterator maxBin = sampleCharge.end();
        for (std::vector<double>::iterator s = sampleCharge.begin();
             s != sampleCharge.end(); ++s) {
            if (maxBin==sampleCharge.end() || *maxBin < *s) maxBin = s;
        }
        if (maxBin != sampleCharge.end()) {
            std::vector<double>::iterator lowBin = maxBin;
            while (lowBin != sampleCharge.begin()) {
                --lowBin;
                if (2.0*rms < step*(maxBin-lowBin)) break;
            }
            std::vector<double>::iterator hiBin = maxBin;
            while (hiBin != sampleCharge.end()) {
                if (2.0*rms < step*(hiBin-maxBin)) break;
                ++hiBin;
            }
            time = 0.0;
            double w = 0.0;
            int bin = 0;
            for (std::vector<double>::iterator b = lowBin; b != hiBin; ++b) {
                time += (*b)*bin;
                w += (*b);
                ++bin;
            }
            time *= step/w;
            time += step*(lowBin-sampleCharge.begin()) + startTime + step/2.0;
        }
    }
#endif
    
    // Base the uncertainty in the time on the number of samples used to find
    // the RMS.
    double timeUnc = rms/std::sqrt(samples);

    // Calculate the charge uncertainty (the calculate starts with chargeUnc
    // being the variance, and then takes the sqrt).  The fundamental
    // uncertainty is in the number of electrons (Poisson distributed).
    double chargeUnc = charge/unit::eplus;

    // Add uncertainty for the sample to sample RMS.  These are fluctuations
    // introduced by the gaussian fluctuations of the electronics (including
    // shot noise).
    double sig = pulseDeconvolution->GetSampleSigma()*std::sqrt(samples);
    double sigC = pulseDeconvolution->GetSampleSigma(samples);
    
    chargeUnc += sigC*sigC;

    // Channels that have had baseline subtraction, will have a positive
    // baslineSigma (this only affects the induction planes).  In this case
    // add a correlated uncertainty for the "wandering of the baseline".  This
    // is analogous to the deviation of the track from a "bestfit" straight
    // line (the s_plane parameter in the PDG notation).  The RMS is
    // correlated between all the channels, but is uncorrelated with the other
    // uncertainties (number of electrons, and electronics noise).  The
    // uncertainty for the "baseline wander" is approximated as the
    // sample sigma.  This should be calculated separately, but since the
    // wander is from the same physics as the sample sigma, it should be a
    // good approximation.
    double sigB = 0.0;
    if (pulseDeconvolution->GetBaselineSigma() > 0) {
        sigB = sigC/4.0;
        chargeUnc += sigB*sigB - sigC*sigC;
    }

    if (chargeUnc < 0.0) chargeUnc = sigC*sigC;
    
    // Now take the sqrt of the variance to get the uncertainty.
    chargeUnc = std::sqrt(chargeUnc);

    CaptNamedInfo("TMakeWireHit", digit.GetChannelId() << " Sigma " << charge
                  << " " << samples
                  << " " << pulseDeconvolution->GetSampleSigma()
                  << " " << sig
                  << " " << sigC
                  << " " << sigB
                  << " " << chargeUnc);
    
    CP::TGeometryId geomId
        = CP::TChannelInfo::Get().GetGeometry(digit.GetChannelId());

    if (!geomId.IsValid()) return CP::THandle<CP::THit>();

    if (!std::isfinite(rms) || rms <= 0.0) {
        CaptError("Time RMS for " << digit.GetChannelId()
                  << " is not positive and finite ("
                  << rms << " out of " << sampleCharge.size() << ")");
        timeUnc = 1*unit::second;
    }

    if (!std::isfinite(timeUnc) || timeUnc <= 0.0) {
        CaptError("Time uncertainty for " << digit.GetChannelId()
                  << " is not positive and finite (" << timeUnc << ")");
        timeUnc = 1*unit::second;
    }

        if (!std::isfinite(chargeUnc) || chargeUnc <= 0.0) {
        CaptError("Charge uncertainty for " << digit.GetChannelId()
                  << " is not positive and finite " << samples);
        chargeUnc = 1*unit::coulomb;
    }

    TChannelCalib calib;

    // Correct for the wire collection efficiency.  Before this correction,
    // the wire is calibrated in terms of the charged collected (or induced)
    // on the wire.  After this correction, the wire is calibrated in terms of
    // "sensed" electrons which is the number of electrons passing in the
    // vicinity of the wire.  For the collection wires, the measured electrons
    // and sensed electrons are the same thing.
    if (fCorrectCollectionEfficiency) {
        charge /= calib.GetCollectionEfficiency(digit.GetChannelId());
        chargeUnc /= calib.GetCollectionEfficiency(digit.GetChannelId());
    }

    // Correct for the electron drift lifetime.
    double deltaT = time - t0;
    if (fCorrectElectronLifetime && deltaT > 0.0) {
        charge *= std::exp(deltaT/calib.GetElectronLifetime());
        chargeUnc *= std::exp(deltaT/calib.GetElectronLifetime());
    }

    // Check the hit validity.
    bool isValid = true;
    if (!std::isfinite(charge)) {
        CaptError("Invalid hit charge for " << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(chargeUnc) || chargeUnc <= 0.0) {
        CaptError("Invalid hit charge uncertainty for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(time)) {
        CaptError("Invalid hit time for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(timeUnc)) {
        CaptError("Invalid hit time uncertainty for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(rms)) {
        CaptError("Invalid hit time RMS for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(startTime)) {
        CaptError("Invalid hit start time for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (!std::isfinite(stopTime)) {
        CaptError("Invalid hit stop time for "
                  << digit.GetChannelId());
        isValid = false;
    }
    if (sampleCharge.size()<1) {
        CaptError("Invalid hit charge samples "
                  << digit.GetChannelId()
                  << " (samples == "<< sampleCharge.size() << ")");
        isValid = false;
    }

    if (!isValid) return CP::THandle<CP::TFADCHit>();

    // Build the hit.
    CP::TWritableFADCHit hit;
    hit.SetGeomId(geomId);
    hit.SetDigit(digit.GetParent());
    hit.SetCharge(charge);
    hit.SetChargeUncertainty(chargeUnc);
    hit.SetTime(time);
    hit.SetTimeRMS(rms);
    hit.SetTimeStart(startTime);
    hit.SetTimeStop(stopTime);
    hit.SetTimeSamples(sampleCharge.begin(), sampleCharge.end());
    hit.SetTimeUncertainty(timeUnc);

    return CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit));
}
