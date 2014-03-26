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

    /// This returns true if the signal is a bipolar signal.  The collection
    /// wires and PMTs are unipolar.  The induction wires are bipolar.
    bool IsBipolarSignal(CP::TChannelId id);

    /// Get the amplifier gain constants for a channel.  The second parameter
    /// is the order of the constant.  Normally, order 0 is the pedestal,
    /// order 1 is linear, order 2 is quadratic, etc.
    double GetGainConstant(CP::TChannelId id, int order=1);

    /// Get the time constants for a channel.  The second parameter is the
    /// order of the constant.  Normally, order 0 is the pedestal, order 1 is
    /// linear, order 2 is quadratic, etc.
    double GetTimeConstant(CP::TChannelId id, int order=1);

    /// Get the digitizer slope of ADC/(input voltage).  The second parameter
    /// is the order of the constant.  Normally, order 0 is the offset (always
    /// zero), order 1 is linear, order 2 is quadratic, etc.
    double GetDigitizerConstant(CP::TChannelId id, int order=1);

    /// Get the electron lifetime.
    double GetElectronLifetime();

    /// Get the electron drift velocity.  This can be calculated by disabling
    /// the drift calibration for CLUSTERCALIB.exe using the -O no-drift flag.
    /// The average X charge vs time can then be fit to determine the drift
    /// lifetime.  
    double GetElectronDriftVelocity();

    /// Get the wire collection efficiency.  This can be calculated by
    /// disabling the efficiency calibration for CLUSTERCALIB.exe using the -O
    /// no-efficiency flag.  The collection efficiency is then the ratio
    /// between the U/X and V/X ratios.  This can also include wire-to-wire
    /// differences, but those are not calculated with CLUSTERCALIB.exe.
    double GetCollectionEfficiency(CP::TChannelId id);

};

#endif
