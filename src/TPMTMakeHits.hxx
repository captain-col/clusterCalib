#ifndef TPMTMakeHits_hxx_seen
#define TPMTMakeHits_hxx_seen

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>

namespace CP {
    class TPMTMakeHits;
};

/// This takes a TCalibPulseDigit and turns it into one or more TDataHit
/// object.  The hit will contain the time, and charge of the hit.
class CP::TPMTMakeHits {
public:
    TPMTMakeHits();
    ~TPMTMakeHits();
    
    /// Take a TCalibPulseDigit reference and find any hits in the pulse.  Any
    /// hits that are found are appended to output THitSelection.
    void operator () (CP::THitSelection& output, 
                      const CP::TCalibPulseDigit& input);

private:
    
    /// The threshold for defining a hit.
    double fThreshold;

};

#endif
