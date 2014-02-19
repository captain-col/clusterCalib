#ifndef TPulseDeconvolution_hxx_seen
#define TPulseDeconvolution_hxx_seen

#include <TCalibPulseDigit.hxx>
#include <TChannelId.hxx>

namespace CP {
    class TPulseDeconvolution;
    class TElectronicsResponse;
    class TWireResponse;
};

class TVirtualFFT;

/// Take a digit, that must be derived from a TCalibPulseDigit, and
/// deconvolute the distribution to get an inferred current. The resulting
/// digit will have the time, and pedestal subtracted charge in the standard
/// HEPUnits definitions.  That's ns and the positron charge (0.00016 fC).
/// The time is measured relative to the start of the digitization period, and
/// not the time of the particle.  The charge is not corrected for the
/// electron lifetime.  The caller owns the resulting
/// CP::TCalibPulseDigit pointer and is responsible for deleting it.
class CP::TPulseDeconvolution {
public:
    explicit TPulseDeconvolution(int sampleCount);
    virtual ~TPulseDeconvolution();

    /// Do the actual calibration.
    CP::TCalibPulseDigit* operator()(const CP::TCalibPulseDigit& digit);

    /// Get the number of samples in the FFT.
    int GetSampleCount() {return fSampleCount;}

    /// Get the baseline uncertainty.  This is calculated relative to the
    /// average baseline and represents the per sample uncertainty in the
    /// deconvolution.
    double GetBaselineSigma() {return fBaselineSigma;}

    /// Get the sample to sample sigma.  Represents the uncorrelated
    /// uncertainty between samples.
    double GetSampleSigma() {return fSampleSigma;}

private:

    /// Initialize the class.
    void Initialize();

    /// Remove the baseline drift that is introduced by random for the
    /// induction wires.  This is applied after the deconvolution, and doesn't
    /// affect the collection wires.
    void RemoveBaseline(CP::TCalibPulseDigit& digit);

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

    /// Hold the number of side band samples to use during the smoothing.
    /// This is set using clusterCalib.smoothing.wire, and is calculated to be
    /// (clusterCalib.smoothing.wire + 1).
    int fSmoothingWindow;

    /// Hold the cut for random fluctuations.
    double fFluctuationCut;

    /// When determining where the peaks probably are so they can be excluded
    /// from the baseline estimate, this is the number of samples that have to
    /// stay inside the fFluctuationCut.
    int fCoherenceZone;

    /// The number of samples that have to be within fFluctuationCut for a
    /// region to be considered "coherent".
    double fCoherenceFraction;

    /// The baseline sigma relative to the average baseline.
    double fBaselineSigma;
    
    /// The sample to sample sigma.
    double fSampleSigma;
};
#endif
