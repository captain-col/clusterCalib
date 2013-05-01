#ifndef TClusterCalib_hxx_seen
#define TClusterCalib_hxx_seen

#include <TEvent.hxx>
#include <TMCChannelId.hxx>

namespace CP {
    class TElectronicsResponse;
    class TWireResponse;
    class TClusterCalib;
};

class TVirtualFFT;

/// This is a very simplistic electronics simulation.  It is not intended for
/// doing physic, but does capture enough of the behavior to develop software.
class CP::TClusterCalib {
public:
    typedef std::vector<double> DoubleVector;
    TClusterCalib();
    virtual ~TClusterCalib();

    bool operator()(CP::TEvent& event);

private:
    /// The maximum pulse length (in time).  The fft (time sequence) length is
    /// the pulse length plus the response length.  The number of samples is
    /// the fft lenght divided by the digit step.
    double fPulseLength;

    /// The maximum length of the response function.  See fPulseLength for
    /// more documentation.
    double fResponseLength;

    /// The shortest time per sample in the ellectronics.  This is used to
    /// calculate the length of the FFT.
    double fDigitStep;

    /// A convenient holder for the number of samples used in the FFT.  This
    /// must equivalent to fFFT->GetN()[0], or there is a bug.  The
    /// equivalence is guarranteed in the constructor.
    int fSampleCount;

    /// The fft class to take the fourier transform.
    TVirtualFFT *fFFT;

    /// The fft class to take the inverse fourier transform.
    TVirtualFFT *fInverseFFT;

    /// The electronics response.  For now, there is only one object, but if
    /// the difference classes of electronics end up with different responses,
    /// then it may pay to have more than one.
    TElectronicsResponse* fElectronicsResponse;

    /// The wire response.  For now, there is only one object, but if
    /// the difference types of wires end up with different responses,
    /// then it may pay to have more than one.
    TWireResponse* fWireResponse;

};
#endif
