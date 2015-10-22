#ifndef TActivityFilter_hxx_seen
#define TActivityFilter_hxx_seen

#include <TEvent.hxx>
#include <TChannelId.hxx>

#include <memory>

namespace CP {
    class TActivityFilter;
};

/// This is a very simplistic activity filter.  This returns true if there is
/// sufficient activity in the detector to make the event "interesting".
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

    /// Set the number of active channels required before an event is
    /// considered to contain activity.
    void SetRequiredChannels(int i) {
        std::cout << "TActivityFilter:: Set the required channels: "
                  << i
                  << std::endl;
        fRequiredChannels = i;
    }
    
    /// Set the maximum number of active channels allowed in an event
    void SetMaximumRequiredChannels(int i) {
        std::cout << "TActivityFilter:: Set maximum required channels: "
                  << i
                  << std::endl;
        fMaximumRequiredChannels = i;
    }

private:

    /// The minimum number of ADC counts above the baseline for a channel to be
    /// counted as active.
    double fMinimumSignal;

    /// The required significance of the fluctuation of a channel above the
    /// baseline for the channel to be counted as active.  This is in "sigma".
    double fRequiredSignificance;

    /// The number of channels that must be active in the detector for the
    /// event to be considered to have activity.
    int fRequiredChannels;

    /// The maximum number of channels that are allowed to register activity
    int fMaximumRequiredChannels;
    
};
#endif
