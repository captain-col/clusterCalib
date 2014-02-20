#ifndef TWireMakeHits_hxx_seen
#define TWireMakeHits_hxx_seen

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>

namespace CP {
    class TWireMakeHits;
};

/// This takes a TCalibPulseDigit and turns it into one or more TDataHit
/// object.  The hit will contain the time, and charge of the hit.
class CP::TWireMakeHits {
public:
    TWireMakeHits();
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
