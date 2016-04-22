#ifndef TPulseDeconvolution_hxx_seen
#define TPulseDeconvolution_hxx_seen

#include <TCalibPulseDigit.hxx>
#include <TChannelId.hxx>

namespace CP {
    class TPulseDeconvolution;
    class TElectronicsResponse;
    class TWireResponse;
    class TNoiseFilter;
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
    int GetSampleCount() const {return fSampleCount;}

    /// Get the baseline uncertainty.  This is calculated relative to the
    /// average baseline and represents the per sample uncertainty in the
    /// deconvolution.
    double GetBaselineSigma() const {return fBaselineSigma;}

    /// Get the sample to sample sigma.  Represents the uncorrelated
    /// uncertainty between samples.
    double GetSampleSigma(int i=0) const {
        if (i<1) return fSampleSigma[0];
        if (i<50) return fSampleSigma[i];
        return fSampleSigma[49];
    }

private:

    /// Initialize the class.
    void Initialize();

    /// Remove the baseline drift from the deconvolution.  This looks at the
    /// sample to sample fluctuations to estimate the background.
    void RemoveBaseline(CP::TCalibPulseDigit& digit,
                        const CP::TCalibPulseDigit& calib);

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

    /// An optimal filter can be applied to reduce the effect of high
    /// frequency noise.  This sets the frequency cut-off as a fraction of the
    /// Nyquist frequency.  A "cut" of 1.0 means no cut-off.  Since the
    /// Nyquist frequency is nominally 1 MHz and our signal is about 500 kHz
    /// (i.e. half the Nyquist frequency), the nominal cut value will probably
    /// be somewhere between 0.8 and 0.95.  It is set using the parameter
    /// clusterCalib.deconvolution.nyquistFraction.
    TNoiseFilter* fNoiseFilter;
    
    /// An optimal filter can be applied to reduce the effect of high
    /// frequency noise.  This sets the frequency cut-off as a fraction of the
    /// Nyquist frequency.  A "cut" of 1.0 means no cut-off.  Since the
    /// Nyquist frequency is nominally 1 MHz and our signal is about 500 kHz
    /// (i.e. half the Nyquist frequency), the nominal cut value will probably
    /// be somewhere between 0.8 and 0.95.  It is set using the parameter
    /// clusterCalib.deconvolution.nyquistFraction.
    double fNyquistFraction;

    /// Hold the range of the loop used to smooth a sample. This is set using
    /// clusterCalib.deconvolution.smoothing which sets the number of
    /// side-band samples to be used.  The loop range is zero to max(1,
    /// clusterCalib.deconvolution.smoothing+1), so smoothing can be disabled
    /// by using a value of zero (that means no side-band samples are used).
    int fSmoothingWindow;

    /// Set the minimum allowed RMS for the calibrated samples around the
    /// baseline.  It's used when the sample RMS is zero (usually because the
    /// MC is being run without any noise.  This sets the minimum value for
    /// fSamlpeSigma.  If this has a negative value, then it is directly used
    /// (or rather abs(fMinimumSigma) is directly used.
    double fMinimumSigma;

    /// Hold the cut for random fluctuations.  This is in units of the RMS of
    /// the sample fluctuations.  This is set using
    /// clusterCalib.deconvolution.fluctuationCut.
    double fFluctuationCut;

    /// Hold the cut on how high the baseline can get.  This is in units of
    /// the baseline standard deviation relative to the median baseline.  This
    /// is set using clusterCalib.deconvolution.baselineCut.
    double fBaselineCut;

    /// When determining where the peaks probably are so they can be excluded
    /// from the baseline estimate, this is the number of samples that have to
    /// stay inside the fFluctuationCut.  This is set using the parameter
    /// clusterCalib.deconvolution.coherenceZone.
    int fCoherenceZone;

    /// The number of samples that have to be within fFluctuationCut for a
    /// region to be considered "coherent".  This is set using the parameter
    /// clusterCalib.deconvolutioncoherenceCut.
    double fCoherenceCut;

    /// When determining where the peaks probably are so they can be excluded
    /// from the baseline estimate, this is the number of samples over which
    /// the drift is estimated.  This is set using the parameter
    /// clusterCalib.deconvolution.driftZone.
    int fDriftZone;

    /// The maximum drift in terms of sigma that will be considered a baseline
    /// fluctuation.
    double fDriftCut;

    /// The baseline sigma relative to the average baseline.
    double fBaselineSigma;
    
    /// The sample to sample sigma.
    double fSampleSigma[50];
};
#endif
