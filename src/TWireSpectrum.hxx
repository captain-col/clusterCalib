#ifndef TWireSpectrum_hxx_seen
#define TWireSpectrum_hxx_seen

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>

namespace CP {
    class TWireSpectrum;
};

/// This takes a TCalibPulseDigit and turns it into one or more TDataHit
/// object.  The hit will contain the time, and charge of the hit.  This has
/// the option of applying the electron lifetime correction, or to just return
/// the raw integrated charge.
class CP::TWireSpectrum {
public:

    /// If the parameter is true, then the electron lifetime is corrected.
    /// This is the normal setting and equalizes the response across the
    /// detector.  If the first parameter is false, then this runs in
    /// "calibration" mode and the drift correction is not applied.  If the
    /// second parameter is false, then this runs in calibration mode, and the
    /// collection efficiency is not applied.
    explicit TWireSpectrum(bool correctDrift=true,
                           bool correctEfficiency=true);
    ~TWireSpectrum();
    
    /// Take a TCalibPulseDigit object with the calibrated charges, and the
    /// deconvolution of the calibrated charges, to find any hits in the
    /// pulse.  Any hits that are found are appended to output THitSelection.
    /// This takes a hit selection of PMT hits that will be used as a hint for
    /// the expected dispersion of the wire signals.  If the sample
    /// uncertainty is provided, this is the uncertainty of each measurement
    /// in the TCalibPulseDigit.  This returns the total charge on the wire.
    double operator () (CP::THitSelection& hits, 
                      const CP::TCalibPulseDigit& calib,
                      const CP::TCalibPulseDigit& deconv,
                      double t0,
                      double baseSigma = 0.0,
                      double sampleSigma = 0.0);

private:

    /// Determine the bounds of the hit.
    std::pair<int, int> HitExtent(int peakIndex,
                                  const CP::TCalibPulseDigit& deconv,
                                  double baselineSigma, double sampleSigma);

    /// If this is true, then the electron lifetime correction is applied.
    /// The correction should normally be applied, but for certain
    /// calibrations it needs to be turned off.  This is controlled in the
    /// constructor.
    bool fCorrectElectronLifetime;

    /// If this is true, then the collection efficiency correction is applied.
    /// The correction should normally be applied, but for certain
    /// calibrations it needs to be turned off.  This is controlled in the
    /// constructor.
    bool fCorrectCollectionEfficiency;

    /// The size of the buffer for the spectrum
    int fNSource;

    /// A buffer for the spectrum.
    float* fSource;

    /// A buffer for the found peaks.
    float* fDest;

    /// A buffer for local work.
    float* fWork;

    /// The required charge in the sample at the peak required for it to be
    /// considered valid.  This is pedestal subtracted and after the response
    /// function has been deconvoluted.
    double fPeakMaximumCol;
    double fPeakMaximumInd;

    /// The required total charge in the peak.
    double fPeakAreaCol;
    double fPeakAreaInd;

    /// The expected width for a drifting point charge.  The deconvolution
    /// removes the electronics response, but does not remove the width
    /// introduced by the drift physics.
    double fPeakWidthCol;
    double fPeakWidthInd;
    
    /// The noise threshold in "RMS" of the measured noise for the channel
    /// (i.e. "sigma" fluctuation before it's not considered noise).  The
    /// noise is not Gaussian, so a 3 sigma cut isn't 99.8 percent.
    double fNoiseThresholdCut;
    
    /// Set the limit on how wide a peak can be in drift time before it's
    /// split into multiple hits.  Peaks wider than this are split up into
    /// multiple hits by drift time.
    double fPeakRMSLimit;

    /// The maximum number of hits allowed per wire.
    int fMaxPeaks;

    /// The number of samples to skip at the beginning and ending of the
    /// digit.  This is needed since the first and last run of samples are
    /// contaminated by FFT "wrap around".  There should not be any signal in
    /// that part of the event anyway.
    int fDigitEndSkip;
    
    /// Once a peak has been found, the charge is calculated by summing "out
    /// from the peak" until the sample charges are below the charge
    /// threshold.  This sets a threshold as a fraction of the peak height.
    double fIntegrationChargeThreshold;

    /// Once a peak has been found, the charge is calculated by summing "out
    /// from the peak" until the sample charges are below the charge
    /// threshold.  This sets a threshold as a sigma  of the baseline and
    /// sample noise (calculated in PulseDeconvolution).
    double fIntegrationNoiseThreshold;

};

#endif
