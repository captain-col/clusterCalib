#ifndef TActivityFilter_hxx_seen
#define TActivityFilter_hxx_seen

#include <TEvent.hxx>
#include <TChannelId.hxx>

#include <memory>

namespace CP {
    class TActivityFilter;
};

/// This is a very simplistic activity filter.  This returns true if there is
/// sufficient activity in the detector to make the event "interesting".  It
/// only uses the X planes and requires at least a segment of "track" to be
/// found.  There is a very loose definition of "track".
class CP::TActivityFilter {
public:
    TActivityFilter();
    virtual ~TActivityFilter();

    bool operator()(CP::TEvent& event);

    /// Set the required significance of a channel above the baseline for the
    /// channel to be counted as active.  This is a cut developed by Yujing
    /// Sun and is a way to find the channels that are significantly above the
    /// noise in an event.
    void SetRequiredSignificance(double v) {
        std::cout << "TActivityFilter:: Set the required sigificance: "
                  << v
                  << std::endl;
        fRequiredSignificance = v;
    }

    /// Set the minimum number of ADC counts above the baseline required for a
    /// channel to be considered significant.  This sets an absolute charge
    /// requirement on the channel.
    void SetMinimumSignal(double v) {
        std::cout << "TActivityFilter:: Set the minimum signal: "
                  << v
                  << std::endl;
        fMinimumSignal = v;
    }

    /// Set the number of hits found in a cluster before an event is
    /// considered to contain activity.
    void SetRequiredHits(int i) {
        std::cout << "TActivityFilter:: Set the required hits: "
                  << i
                  << std::endl;
        fRequiredHits = i;
    }
    
    /// Set the maximum number of hits in an event.
    void SetMaximumAllowedHits(int i) {
        std::cout << "TActivityFilter:: Set maximum allowed hits: "
                  << i
                  << std::endl;
        fMaximumAllowedHits = i;
    }

private:

    /// The minimum number of ADC counts above the baseline for a channel to be
    /// counted as active.
    double fMinimumSignal;

    /// The required significance of the fluctuation of a channel above the
    /// baseline for the channel to be counted as active.  This is in "sigma".
    double fRequiredSignificance;

    /// The number of hits that must be in a cluster for the
    /// event to be considered to have activity.
    int fRequiredHits;

    /// The maximum number of hits that are allowed to register activity
    int fMaximumAllowedHits;
    
};
#endif
