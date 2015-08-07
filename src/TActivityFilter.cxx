#include "TActivityFilter.hxx"

#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>

#include <CaptGeomId.hxx>
#include <HEPUnits.hxx>

#include <vector>

CP::TActivityFilter::TActivityFilter() {}

CP::TActivityFilter::~TActivityFilter() {}

bool CP::TActivityFilter::operator() (CP::TEvent& event) {
    CaptLog("Filter event " << event.GetContext());
    CP::TChannelInfo::Get().SetContext(event.GetContext());

#ifdef PMT_ACTIVITY_CHECK
    /// Check if there is any PMT activity
    CP::THandle<CP::TDigitContainer> pmt
        = event.Get<CP::TDigitContainer>("~/digits/pmt");


#endif
    
    CP::THandle<CP::TDigitContainer> drift
        = event.Get<CP::TDigitContainer>("~/digits/drift");
    std::vector<double> work;
    
    double minimumSignal = 30; // about a 1/3 mip.
    double significance = 7.0;
    int significantChannels = 0;
    for (std::size_t d = 0; d < drift->size(); ++d) {
        const CP::TPulseDigit* digit
            = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
        if (!digit) {
            CaptError("Non-pulse in drift digits");
            continue;
        }

        // Check that we have a wire;
        CP::TChannelId cid = digit->GetChannelId();
        CP::TGeometryId geometryId = CP::TChannelInfo::Get().GetGeometry(cid);
        if (!CP::GeomId::Captain::IsWire(geometryId)) continue;
        
        if (work.size() != digit->GetSampleCount()) {
            work.resize(digit->GetSampleCount());
        }

        // Calculate the pedestal for the channel.
        for (std::size_t i=1; i< digit->GetSampleCount(); ++i) {
            work[i] = digit->GetSample(i);
        }
        std::sort(work.begin(), work.end());
        double pedestal = work[0.5*digit->GetSampleCount()];
        if (pedestal < 200) continue;
        
        // Calculate the raw channel sigma.
        double channelSigma = 0.0;
        for (std::size_t i=0; i< digit->GetSampleCount(); ++i) {
            double v = digit->GetSample(i) - pedestal;
            channelSigma += v*v;
        }
        channelSigma /= digit->GetSampleCount();
        channelSigma = std::sqrt(channelSigma);
        if (channelSigma < 1) continue;
        
        // Calculate the Gaussian sigma for the channel.
        for (std::size_t i=1; i< digit->GetSampleCount(); ++i) {
            work[i] = std::abs(digit->GetSample(i)-digit->GetSample(i-1));
        }
        work[0] = 0.0;
        std::sort(work.begin(), work.end());
        int iMedian = 0.52*digit->GetSampleCount();
        int iLow = iMedian;
        while (iLow > 0) {
            if (work[iLow] != work[iMedian]) break;
            --iLow;
        }
        int iHigh = iMedian;
        while ((iHigh+1)<work.size()) {
            if (work[iHigh] != work[iMedian]) break;
            ++iHigh;
        }
        double lowDiff = iMedian-iLow;
        double diff = iHigh - iLow;
        double gaussianSigma = work[iMedian] + lowDiff/diff;

        // Check if there is activity on this channel.
        for (std::size_t i=0.05*digit->GetSampleCount();
             i < 0.95*digit->GetSampleCount(); ++i) {
            double signal = digit->GetSample(i)-pedestal;
            if (significance*channelSigma <= signal
                && minimumSignal <= signal) {
                ++significantChannels;
                break;
            }
        }
    }

    std::cout << significantChannels << " " << significance << std::endl;
    if (significantChannels > 6) return true;
    return false;
}
