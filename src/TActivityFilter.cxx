#include "TActivityFilter.hxx"
#include "TTmplDensityCluster.hxx"

#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TChannelCalib.hxx>

#include <TRuntimeParameters.hxx>
#include <CaptGeomId.hxx>
#include <HEPUnits.hxx>

#include <vector>
#include <utility>

namespace {
    class HitDistanceMetric {
    public:
        double operator() (const std::pair<double,double>& lhs,
                           const std::pair<double,double>& rhs) {
            double x = std::abs(lhs.first-rhs.first);
            double z = std::abs(lhs.second-rhs.second);
            return std::sqrt(x*x+z*z);
        }
    };
    typedef CP::TTmplDensityCluster<std::pair<double,double>,
                                    HitDistanceMetric> HitCluster;
}

CP::TActivityFilter::TActivityFilter() {
    fRequiredHits
        = CP::TRuntimeParameters::Get().GetParameterI(
            "clusterCalib.filter.requiredHits");
    fMaximumAllowedHits
        = CP::TRuntimeParameters::Get().GetParameterI(
            "clusterCalib.filter.maximumAllowedHits");
    fRequiredSignificance
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.filter.requiredSignificance");
    fMinimumSignal 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.filter.minimumSignal");
}

CP::TActivityFilter::~TActivityFilter() {}

bool CP::TActivityFilter::operator() (CP::TEvent& event) {
    CP::TChannelInfo::Get().SetContext(event.GetContext());

    CP::THandle<CP::TDigitContainer> drift
        = event.Get<CP::TDigitContainer>("~/digits/drift");
    if (!drift) return false;

    std::vector<double> work;
    std::vector< std::pair<double,double> > collectionWireHits;

    CP::TChannelCalib calib;
    for (std::size_t d = 0; d < drift->size(); ++d) {
        const CP::TPulseDigit* digit
            = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
        if (!digit) {
            CaptError("Non-pulse in drift digits");
            continue;
        }

        CP::TChannelId cid = digit->GetChannelId();

        // Check if this is a bad channel.
        if (!calib.IsGoodChannel(cid)) continue;
        
        // Check that we have a collection wire;
        CP::TGeometryId geometryId = CP::TChannelInfo::Get().GetGeometry(cid);
        if (!CP::GeomId::Captain::IsWire(geometryId)) continue;
        if (calib.IsBipolarSignal(cid)) continue;

        if (work.size() != digit->GetSampleCount()) {
            work.resize(digit->GetSampleCount());
        }

        // Calculate the pedestal for the channel.
        for (std::size_t i=1; i< digit->GetSampleCount(); ++i) {
            work[i] = digit->GetSample(i);
        }
        std::sort(work.begin(), work.end());
        double pedestal = work[0.5*digit->GetSampleCount()];

        // Reject channels that have a strange pedestals
        if (pedestal < 200) continue;
        if (pedestal > 3000) continue;
        
        // Calculate the raw and Gaussian sigma for the channel.
        double channelSigma = 0.0;
        for (std::size_t i=1; i< digit->GetSampleCount(); ++i) {
            double v = digit->GetSample(i) - pedestal;
            channelSigma += v*v;
            work[i] = std::abs(digit->GetSample(i)-digit->GetSample(i-1));
        }
        channelSigma /= digit->GetSampleCount();
        channelSigma = std::sqrt(channelSigma);

        // Reject channels that don't show any activity.
        if (channelSigma < 0.5) continue;

        // Find the "RMS" based on sample to sample fluctuations.
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
        while ((std::size_t) (iHigh+1)<work.size()) {
            if (work[iHigh] != work[iMedian]) break;
            ++iHigh;
        }
        double lowDiff = iMedian-iLow;
        double diff = iHigh - iLow;

        // Interpolate between the integer values for the ADC counts to find
        // the interpolated median.
        double gaussianSigma = work[iMedian] + lowDiff/diff;

        // Check if there is activity on this channel.  This ignores the very
        // beginning and very end of the time window.
        std::pair<double,double> candidatePair(-1,0);
        double candidatePeak = 0.0;
        for (std::size_t i=0.05*digit->GetSampleCount();
             i < 0.95*digit->GetSampleCount(); ++i) {
            double signal = digit->GetSample(i)-pedestal;
            double rise = digit->GetSample(i-1)-pedestal;
            if (rise >= signal) continue;
            rise = digit->GetSample(i-2)-pedestal;
            if (rise >= signal) continue;
            double fall = digit->GetSample(i+1)-pedestal;
            if (fall >= signal) continue;
            fall = digit->GetSample(i+2)-pedestal;
            if (fall >= signal) continue;
            if (fRequiredSignificance*gaussianSigma <= signal
                && fMinimumSignal <= signal) {
                double x = 3.0*CP::GeomId::Captain::GetWireNumber(geometryId);
                double z = 1.6*0.5*i; // 1.6 mm/us * 0.5 us/sample * samples
                if (candidatePair.first < 0) {
                    candidatePair = std::make_pair(x,z);
                    candidatePeak = signal;
                    continue;
                }
                if (std::abs(candidatePair.first-x) > 1) {
                    collectionWireHits.push_back(candidatePair);
                    candidatePair = std::make_pair(x,z);
                    candidatePeak = signal;
                    continue;
                }
                if (std::abs(candidatePair.second-z) > 10) {
                    collectionWireHits.push_back(candidatePair);
                    candidatePair = std::make_pair(x,z);
                    candidatePeak = signal;
                    continue;
                }
                if (candidatePeak < signal) {
                    candidatePair = std::make_pair(x,z);
                    candidatePeak = signal;
                }
            }
        }
        if (candidatePair.first >= 0) {
            collectionWireHits.push_back(candidatePair);
        }
    }

    CaptLog("Filter with " << collectionWireHits.size() << " hits");
    if (collectionWireHits.size() > (std::size_t) fMaximumAllowedHits) {
        return false;
    }

    // Form a cluster requiring two neighbors, and that the hits be in
    // nearby wires (wire spacing is 3mm).
    HitCluster hitCluster(2,8.5);
    hitCluster.Cluster(collectionWireHits.begin(), collectionWireHits.end());
    CaptLog("Clusters found " << hitCluster.GetClusterCount());
    if (hitCluster.GetClusterCount() < 1) return false;
    
    CaptLog("Size of biggest cluster " << hitCluster.GetCluster(0).size());
    if (hitCluster.GetCluster(0).size()<(std::size_t)fRequiredHits) {
        return false;
    }

    return true;
}
