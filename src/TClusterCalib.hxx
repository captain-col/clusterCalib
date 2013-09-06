#ifndef TClusterCalib_hxx_seen
#define TClusterCalib_hxx_seen

#include <TEvent.hxx>
#include <TMCChannelId.hxx>

#include <memory>

namespace CP {
    class TClusterCalib;
    class TPulseCalib;
    class TPulseDeconvolution;
};

/// This is a very simplistic electronics simulation.  It is not intended for
/// doing physic, but does capture enough of the behavior to develop software.
/// On output, the hits for the PMT are saved in the "pmt" hit selection, and
/// the wire hits are contained in "drift".
class CP::TClusterCalib {
public:
    typedef std::vector<double> DoubleVector;
    TClusterCalib();
    virtual ~TClusterCalib();

    bool operator()(CP::TEvent& event);

    /// Set a flag to save the calibrated pulses.
    void SaveCalibratedPulses(bool value = true) {
        fSaveCalibratedPulses = value;
    }
    
private:
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

    /// The minimum time per sample in the ellectronics.  This is used to
    /// calculate the length of the FFT, and the exact value isn't important
    /// (as long as it is equal to, or shorter, than the actual time per
    /// sample
    double fDigitStep;

    /// A convenient holder for the number of samples used in the FFT.  This
    /// must equivalent to fFFT->GetN()[0], or there is a bug.  The
    /// equivalence is guarranteed in the constructor.
    int fSampleCount;

    /// A class to calibrate the digital pulses into calibrated charge pulses. 
    CP::TPulseCalib* fCalibrate;

    /// A class to deconvolute the electronics shape 
    CP::TPulseDeconvolution* fDeconvolution;

    /// A flag that the calibrated pulse digits should be saved on the output.
    /// This is the input into the peak finding.
    bool fSaveCalibratedPulses;

};
#endif
