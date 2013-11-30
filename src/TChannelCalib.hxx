#ifndef TChannelCalib_hxx_seen
#define TChannelCalib_hxx_seen

#include "EClusterCalib.hxx"

namespace CP {
    class TChannelCalib;
    class TChannelId;
};

namespace CP {
    EXCEPTION(EChannelCalib,ECore);
    EXCEPTION(EChannelCalibUnknownType, EChannelCalib);
};

/// Provide generic an interface to get the calibration coefficients.  This
/// works for calibration constants provided in an MC file, and for constants
/// for the data.
class CP::TChannelCalib {
public:
    TChannelCalib();
    ~TChannelCalib();

    /// Get the gain constants for a channel.  The second parameter is the
    /// order of the constant.  Normally, order 0 is the pedestal, order 1 is
    /// linear, order 2 is quadratic, etc.
    double GetGainConstant(CP::TChannelId id, int order=1);

    /// Get the time constants for a channel.  The second parameter is the
    /// order of the constant.  Normally, order 0 is the pedestal, order 1 is
    /// linear, order 2 is quadratic, etc.
    double GetTimeConstant(CP::TChannelId id, int order=1);

};

#endif
