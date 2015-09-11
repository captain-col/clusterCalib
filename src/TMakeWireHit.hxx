#ifndef TMakeWireHit_hxx_seen
#define TMakeWireHit_hxx_seen

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>

namespace CP {
    class TMakeWireHit;
};

/// Build a THit out of a sequence of samples in a TCalibPulseDigit.  The hit
/// will contain the time, and charge of the hit.  This has the option of
/// applying the electron lifetime correction, or to just return the raw
/// integrated charge.
class CP::TMakeWireHit {
public:

    /// If the parameter is true, then the electron lifetime is corrected.
    /// This is the normal setting and equalizes the response across the
    /// detector.  If the first parameter is false, then this runs in
    /// "calibration" mode and the drift correction is not applied.  If the
    /// second parameter is false, then this runs in calibration mode, and the
    /// collection efficiency is not applied.
    explicit TMakeWireHit(bool correctDrift=true,
                          bool correctEfficiency=true);
    ~TMakeWireHit();
    
    /// Build a hit out of the digit samples between beginIndex and endIndex
    /// (inclusive).  The digit step size is provided as an input.
    CP::THandle<CP::THit> 
    operator ()(const CP::TCalibPulseDigit& digit, 
                double digitStep, double t0,
                double baselineSigma, double sampleSigma,
                std::size_t beginIndex, std::size_t endIndex, bool split);
    
private:

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

};
#endif