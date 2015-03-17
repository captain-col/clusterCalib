#ifndef TActivityFilter_hxx_seen
#define TActivityFilter_hxx_seen

#include <TEvent.hxx>
#include <TChannelId.hxx>

#include <memory>

namespace CP {
    class TActivityFilter;
};

/// This is a very simplistic activity filter.  This returns true if any of
/// the collection plane wires show signs of activity.
class CP::TActivityFilter {
public:
    TActivityFilter();
    virtual ~TActivityFilter();

    bool operator()(CP::TEvent& event);

private:
};
#endif
