#include "TMakeWireHit.hxx"
#include "TChannelCalib.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TFADCHit.hxx>
#include <HEPUnits.hxx>
#include <THandle.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>
#include <CaptGeomId.hxx>

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
                             double baselineSigma, double sampleSigma,
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

    // Base the uncertainty in the time on the number of samples used to find
    // the RMS.
    double timeUnc = rms/std::sqrt(samples);

    // The charge uncertainty is calculated assuming the limit of Poisson
    // statistics, but assumes that the background uncertainties are
    // correlated.
    double sig = sampleSigma*std::sqrt(samples);
    double chargeUnc = std::sqrt(charge + sig*sig) + baselineSigma*samples;

    CP::TGeometryId geomId
        = CP::TChannelInfo::Get().GetGeometry(digit.GetChannelId());

    if (!geomId.IsValid()) return CP::THandle<CP::THit>();

    if (!std::isfinite(timeUnc) || timeUnc <= 0.0) {
        CaptError("Time uncertainty for " << digit.GetChannelId()
                  << " is not positive or and finite ");
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
