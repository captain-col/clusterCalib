#ifndef TWireMakeHits_hxx_seen
#define TWireMakeHits_hxx_seen

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>

namespace CP {
    class TWireMakeHits;
};

/// This takes a TCalibPulseDigit and turns it into one or more TDataHit
/// object.  The hit will contain the time, and charge of the hit.  This has
/// the option of applying the electron lifetime correction, or to just return
/// the raw integrated charge.
class CP::TWireMakeHits {
public:

    /// If the parameter is true, then the electron lifetime is corrected.
    /// This is the normal setting and equalizes the response across the
    /// detector.  If the parameter is false, then this runs in "calibration"
    /// mode and the raw integrated charge is returned.
    explicit TWireMakeHits(bool correctDrift=true);
    ~TWireMakeHits();
    
    /// Take a TCalibPulseDigit reference and find any hits in the pulse.  Any
    /// hits that are found are appended to output THitSelection.  This takes
    /// a hit selection of PMT hits that will be used as a hint for the
    /// expected dispersion of the wire signals.  If the sample uncertainty is
    /// provided, this is the uncertainty of each measurement in the
    /// TCalibPulseDigit.
    void operator () (CP::THitSelection& hits, 
                      const CP::TCalibPulseDigit& digits,
                      double t0,
                      double baseSigma = 0.0,
                      double sampleSigma = 0.0);

private:

    /// Build a hit out of the digit samples between beginIndex and endIndex.
    /// The digit step size is provided as an input.
    CP::THandle<CP::THit> 
    MakeHit(const CP::TCalibPulseDigit& digit, 
            double digitStep, double t0,
            double baselineSigma, double sampleSigma,
            int beginIndex, int endIndex, bool split);

    /// If this is true, then the electron lifetime correction is applied.
    /// The correction should normally be applied, but for certain
    /// calibrations it needs to be turned off.  This is controlled in the
    /// constructor.
    bool fCorrectElectronLifetime;

    /// The size of the buffer for the spectrum
    int fNSource;

    /// A buffer for the spectrum.
    float* fSource;

    /// A buffer for the found peaks.
    float* fDest;

    /// The charge required at the peak for it to be considered valid.
    double fPeakMaximumCut;

    /// The power in the TSpectrum deconvolution for the peak for it to be
    /// considered valid.
    double fPeakDeconvolutionCut;

    /// Set the limit on how wide a peak can be in drift time before it's
    /// split into multiple hits.  Peaks wider than this are split up into
    /// multiple hits by drift time.
    double fPeakRMSLimit;
};

#endif
