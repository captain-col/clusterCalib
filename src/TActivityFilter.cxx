#include "TActivityFilter.hxx"

#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>

#include <TRuntimeParameters.hxx>
#include <CaptGeomId.hxx>
#include <HEPUnits.hxx>

#include <vector>

CP::TActivityFilter::TActivityFilter() {
    fRequiredChannels
        = CP::TRuntimeParameters::Get().GetParameterI(
            "clusterCalib.filter.requiredChannels");
    fRequiredSignificance
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.filter.requiredSignificance");
    fMinimumSignal 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.filter.minimumSignal");
    fMaximumRequiredChannels
        = CP::TRuntimeParameters::Get().GetParameterI(
            "clusterCalib.filter.maximumRequiredChannels");
}

CP::TActivityFilter::~TActivityFilter() {}

bool CP::TActivityFilter::operator() (CP::TEvent& event) {
    CP::TChannelInfo::Get().SetContext(event.GetContext());

#ifdef PMT_ACTIVITY_CHECK
    /// Check if there is any PMT activity
    CP::THandle<CP::TDigitContainer> pmt
        = event.Get<CP::TDigitContainer>("~/digits/pmt");
#endif


    /// There wasn't PMT activity, so check for wire activity.
    CP::THandle<CP::TDigitContainer> drift
        = event.Get<CP::TDigitContainer>("~/digits/drift");
    if (!drift) return false;

    std::vector<double> work;
    
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

        // Reject channels that have a strange pedestal
        if (pedestal < 200) continue;
        if (pedestal > 3000) continue;
        
        // Calculate the raw channel sigma.
        double channelSigma = 0.0;
        for (std::size_t i=0; i< digit->GetSampleCount(); ++i) {
            double v = digit->GetSample(i) - pedestal;
            channelSigma += v*v;
        }
        channelSigma /= digit->GetSampleCount();
        channelSigma = std::sqrt(channelSigma);

        
        // Reject channels that don't show any activity.
        if (channelSigma < 0.5) continue;

#ifdef CALCULATE_GAUSSIAN_SIGMA
        // Calculate the Gaussian sigma for the channel.
        for (std::size_t i=1; i< digit->GetSampleCount(); ++i) {
            work[i] = std::abs(digit->GetSample(i)-digit->GetSample(i-1));
        }
        work[0] = 0.0;
        std::sort(work.begin(), work.end());
        int iMedian = 0.52*digit->GetSampleCount();

        // Count the number of channels that have the same value as the median
        // ADC difference.
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

        // Interpolate between the integer values for the ADC counts to find
        // the interpolated median.
        double gaussianSigma = work[iMedian] + lowDiff/diff;
#endif        

        // Check if there is activity on this channel.  This ignores the very
        // beginning and very end of the time window.
        for (std::size_t i=0.05*digit->GetSampleCount();
             i < 0.95*digit->GetSampleCount(); ++i) {
            double signal = digit->GetSample(i)-pedestal;
            if (fRequiredSignificance*channelSigma <= signal
                && fMinimumSignal <= signal) {
                ++significantChannels;
                break;
            }
        }
    }

    if ((significantChannels >= fRequiredChannels) && (significantChannels <= fMaximumRequiredChannels)) {
        CaptLog("Save event: " << event.GetContext());
        return true;
    }

    return false;
}
