#ifndef TClusterCalib_hxx_seen
#define TClusterCalib_hxx_seen

#include <TEvent.hxx>
#include <TChannelId.hxx>
#include <TCaptLog.hxx>

#include <memory>

namespace CP {
    class TClusterCalib;
    class TPulseCalib;
    class TPulseDeconvolution;
};

/// Apply the calibration to an event and find pulse to make hits.
class CP::TClusterCalib {
public:
    typedef std::vector<double> DoubleVector;
    TClusterCalib();
    virtual ~TClusterCalib();

    bool operator()(CP::TEvent& event);

    /// Set a flag to save the calibrated pulses.
    void SaveCalibratedPulses(bool value = true) {
        CaptLog("Save Calibrated Pulses set to " << value);
        fSaveCalibratedPulses = value;
    }

    /// Set a flag to save the decorrelated pulses.
    void SaveDecorrelatedPulses(bool value = true) {
        CaptLog("Save Decorrelated Pulses set to " << value);
        fSaveDecorrelatedPulses = value;
    }

    /// Set a flag to save the deconvolved pulses.
    void SaveDeconvolvedPulses(bool value = true) {
        CaptLog("Save Deconvolved Pulses set to " << value);
        fSaveDeconvolvedPulses = value;
    }

    /// Set a flag to calibrate all channels (even channels not attached to a
    /// wire).  This defaults to false.
    void CalibrateAllChannels(bool value = true) {
        fCalibrateAllChannels = value;
    }

    /// Set a flag to correct the wire collection efficiency.  If this is
    /// true, then the normal calibration is applied, but if it's false, then
    /// the collection efficiency correction is not applied.  The normal
    /// setting is true.  This is false when running the collection efficiency
    /// calibration is being calculated.
    void ApplyEfficiencyCalibration(bool value = true) {
        fApplyEfficiencyCalibration = value;
    }

    /// Set a flag to trigger removing the correlated pedestal.
    void RemoveCorrelations(bool value = true) {
        fRemoveCorrelatedPedestal = value;
    }
private:

    /// Apply the channel calibrations to all channels.  This takes a
    /// container of raw digits and returns a container of calibrated digits.
    CP::THandle<CP::TDigitContainer>
    CalibrateChannels(CP::TEvent& event,
                     CP::THandle<CP::TDigitContainer> drift);

    /// Removes baseline fluctuations shared by channels.  It calculates the
    /// correlations, and then uses the correlations to calculate a pedestal
    /// based on correlated channels.  It's pretty slow.
    CP::THandle<CP::TDigitContainer>
    RemoveCorrelatedPedestal(CP::TEvent& event,
                             CP::THandle<CP::TDigitContainer> driftCalib);

    /// Applies the deconvolution to all channels.  It's pretty slow, but must
    /// be done before the peaks are found.
    CP::THandle<CP::TDigitContainer>
    DeconvolveSignals(CP::TEvent& event,
                      CP::THandle<CP::TDigitContainer> driftCalib);

    /// The maximum pulse length (in time).  The fft (time sequence) length is
    /// the pulse length plus the response length.  The exact value isn't
    /// important as long as it is bigger, or equal to, the actual pulse
    /// length. The number of samples is the fft length divided by the digit
    /// step.
    double fPulseLength;

    /// The maximum length of the response function.  The exact value of this
    /// variable isn't important as long it is bigger than the actual response
    /// function length.  See fPulseLength for more documentation.
    double fResponseLength;

    /// A convenient holder for the number of samples used in the FFT.  This
    /// must equivalent to fFFT->GetN()[0], or there is a bug.  The
    /// equivalence is guarranteed in the constructor.
    int fSampleCount;

    /// A class to calibrate the digital pulses into calibrated charge pulses. 
    CP::TPulseCalib* fCalibrate;

    /// A class to deconvolute the electronics shape 
    CP::TPulseDeconvolution* fDeconvolution;

    /// A flag that the calibrated pulse digits should be saved on the output.
    bool fSaveCalibratedPulses;

    /// A flag that the calibrated pulse digits after subtracting correlated
    /// pedestals should be saved on the output.
    bool fSaveDecorrelatedPulses;
    
    /// A flag that the calibrated pulse digits after applying the
    /// deconvolution should be saved on the output.
    bool fSaveDeconvolvedPulses;
    
    /// A flag to trigger calibrate channels not attached to a wire.
    bool fCalibrateAllChannels;

    /// A flag to determine if the electron lifetime calibration is applied.
    /// This is normally true, and the false value is reserved for calculating
    /// the drift calibration constant.
    bool fApplyDriftCalibration;

    /// A flag to determine if the efficiency calibration is applied.
    /// This is normally true, and the false value is reserved for calculating
    /// calibration constant.
    bool fApplyEfficiencyCalibration;

    /// A flag to remove baseline fluctuations shared by channels.  It
    /// calculates the correlations, and then uses the correlations to
    /// calculate a pedestal based on correlated channels.  It's pretty slow.
    bool fRemoveCorrelatedPedestal;
};
#endif
