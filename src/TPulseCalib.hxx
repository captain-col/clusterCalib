#ifndef TPulseCalib_hxx_seen
#define TPulseCalib_hxx_seen

#include <TCalibPulseDigit.hxx>
#include <TChannelId.hxx>

namespace CP {
    class TPulseCalib;
};

/// Take a digit, that must be derived from a TPulseDigit, and apply
/// calibrations.  The resulting digit will have the time, and pedestal
/// subtracted charge in the standard HEPUnits definitions.  That's ns and
/// the positron charge (0.00016 fC).  The time is measured relative to
/// the start of the digitization period, and not the time of the
/// particle.  The charge is not corrected for the electron lifetime.  The
/// caller owns the resulting CP::TCalibPulseDigit pointer and is
/// responsible for deleting it.
class CP::TPulseCalib {
public:
    TPulseCalib();
    virtual ~TPulseCalib();

    /// Do the actual calibration.
    CP::TCalibPulseDigit* operator()(const CP::TDigitProxy& digit);

};
#endif
