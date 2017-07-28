#include "TWirePeaks.hxx"
#include "TMakeWireHit.hxx"
#include "TChannelCalib.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TFADCHit.hxx>
#include <HEPUnits.hxx>
#include <THandle.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>
#include <CaptGeomId.hxx>

#include <TH1F.h>
#include <TH2F.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <sstream>
#include <iostream>

CP::TWirePeaks::TWirePeaks(bool correctLifetime,
                           bool correctEfficiency) {
    fCorrectElectronLifetime = correctLifetime;
    fCorrectCollectionEfficiency = correctEfficiency;
    fMaxPeaks = 50;
    
    fPeakMaximumCol
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.charge.collection");

    fPeakMaximumInd
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.charge.induction");

    fPeakAreaCol
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.area.collection");

    fPeakAreaInd
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.area.induction");

    fPeakWidthCol
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.width.collection");

    fPeakWidthInd
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.width.induction");

    fNoiseThresholdCut
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.noise");

    fPeakRMSLimit
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.peakSearch.rmsLimit");

    fDigitEndSkip
        = CP::TRuntimeParameters::Get().GetParameterI(
            "clusterCalib.peakSearch.endSkip");
    if (fDigitEndSkip < 1) fDigitEndSkip = 1;
    
    fIntegrationChargeThreshold = 0.001; // effectively zero.

    fIntegrationShoulderThreshold = 0.5;

    fIntegrationValleyThreshold = 0.1;

}
CP::TWirePeaks::~TWirePeaks() {}

std::pair<int, int>
CP::TWirePeaks::PeakExtent(int peakIndex,
                           const CP::TCalibPulseDigit& deconv) {

    double peak = deconv.GetSample(peakIndex);
    // The lowest point found so far.
    double valley = peak;
    // A boolean that is true when the shoulder threshold has been passed.
    bool passedShoulder = false;
    // A threshold to say that the main body of the peak has been passed.
    double shoulderThreshold = peak*fIntegrationShoulderThreshold;
    // A threshold to say that we are climbing the next peak.
    double valleyThreshold = peak*fIntegrationValleyThreshold;
    
    std::size_t beginIndex = peakIndex;
    for (std::size_t i=peakIndex; i>0; --i) {
        if (0.5 < fWork[i]) break;
        double v = deconv.GetSample(i);
        if (v < valley) {
            valley = v;
            beginIndex = i;
        }
        if (v < shoulderThreshold) passedShoulder = true;
        if (v < fIntegrationChargeThreshold*peak) break;
        if (passedShoulder && v > valley+valleyThreshold) break;
    }

    valley = peak;
    passedShoulder = false;

    std::size_t endIndex = peakIndex;
    for (std::size_t i=peakIndex; i<deconv.GetSampleCount(); ++i) {
        if (0.5 < fWork[i]) break;
        double v = deconv.GetSample(i);
        if (v < valley) {
            valley = v;
            endIndex = i;
        }
        if (v < shoulderThreshold) passedShoulder = true;
        if (v < fIntegrationChargeThreshold*peak) break;
        if (passedShoulder && v > valley+valleyThreshold) break;
    }
    return std::make_pair(beginIndex,endIndex);
}

double CP::TWirePeaks::PeakFWHM(int peakIndex,
                                const std::pair<int,int>& extent,
                                const CP::TCalibPulseDigit& deconv) {
    if (peakIndex < extent.first) return -1;
    if (peakIndex > extent.second) return -1;
    double peakHeight = deconv.GetSample(peakIndex);
    
    // Find the interpolated full width half max.
    int lowBound = peakIndex;
    for (int j=peakIndex-1; extent.first<=j; --j) {
        if (deconv.GetSample(j) < 0.5*peakHeight) break;
        lowBound = j;
    }
    int highBound = peakIndex;
    for (int j = peakIndex+1; j<=extent.second; ++j) {
        if (deconv.GetSample(j) < 0.5*peakHeight) break;
        highBound = j;
    }
    double lowValue = lowBound -
        (0.5*peakHeight-deconv.GetSample(lowBound-1))
        /(deconv.GetSample(lowBound)-deconv.GetSample(lowBound-1));
    double highValue = highBound +
        (deconv.GetSample(highBound)-0.5*peakHeight)
        /(deconv.GetSample(highBound)-deconv.GetSample(highBound+1));
    double fwhm = highValue - lowValue;
    return fwhm;
}

double CP::TWirePeaks::operator() (CP::THitSelection& hits,
                                   const CP::TCalibPulseDigit& deconv,
                                   double t0) {
    double wireCharge = 0.0;

    // Find the time per sample in the digit.
    double digitStep = deconv.GetLastSample()-deconv.GetFirstSample();
    digitStep /= deconv.GetSampleCount();

    if (fWork.capacity() < deconv.GetSampleCount()) {
        fWork.resize(deconv.GetSampleCount());
    }
    
    CP::TChannelCalib channelCalib;
    bool wireIsBipolar = channelCalib.IsBipolarSignal(deconv.GetChannelId());

    double peakMaximumCut = fPeakMaximumCol;
    double peakAreaCut = fPeakAreaCol;
    double peakWidthCut = fPeakWidthCol;
    if (wireIsBipolar) {
        peakMaximumCut = fPeakMaximumInd;
        peakAreaCut = fPeakAreaInd;
        peakWidthCut = fPeakWidthInd;
    }

#ifdef TPC_WIRE    
    std::ostringstream histName;
    histName << "wire-"
             << std::setw(3)
             << std::setfill('0')
             << CP::TChannelInfo::Get().GetWireNumber(deconv.GetChannelId());

    std::ostringstream histTitle;
    histTitle << "Wire "
             << CP::TChannelInfo::Get().GetWireNumber(deconv.GetChannelId())
             << " (" << deconv.GetChannelId().AsString() << ")";
#else
    std::ostringstream histName;
    CP::TGeometryId id
        = CP::TChannelInfo::Get().GetGeometry(deconv.GetChannelId());

    if (CP::GeomId::Captain::IsUWire(id)) histName << "wire-u";
    if (CP::GeomId::Captain::IsVWire(id)) histName << "wire-v";
    if (CP::GeomId::Captain::IsXWire(id)) histName << "wire-x";
    histName << "-" << std::setw(3) << std::setfill('0')
             << CP::GeomId::Captain::GetWireNumber(id);

    std::ostringstream histTitle;
    if (CP::GeomId::Captain::IsUWire(id)) histTitle << "Wire U";
    if (CP::GeomId::Captain::IsVWire(id)) histTitle << "Wire V";
    if (CP::GeomId::Captain::IsXWire(id)) histTitle << "Wire X";
    histTitle << "-" << CP::GeomId::Captain::GetWireNumber(id)
              << " (" << deconv.GetChannelId().AsString() << ")";
#endif

    // Find the magnitude of noise for this channel.
    for (std::size_t i = 0; i<deconv.GetSampleCount(); ++i) {
        fWork[i] = std::abs(deconv.GetSample(i));
    }
    std::sort(&fWork[fDigitEndSkip],
              &fWork[deconv.GetSampleCount()-fDigitEndSkip]);
    int inoise = fDigitEndSkip + 0.7*(deconv.GetSampleCount()-2*fDigitEndSkip);
    // A threshold will be set in terms of standard deviations of the noise.
    // Peaks less than this are rejected as noise.
    double noise = fWork[inoise];

    // Protect against a "zero" channel.
    if (noise < 3) {
        CaptLog("Wire with no signal: " << deconv.GetChannelId()
                << " noise: " << noise
                << " max: " << fWork[deconv.GetSampleCount()-1]);
        // return wireCharge;
    }
    
    // Don't bother with channels that have crazy big noise.  This says if the
    // noise is bigger than the peak selection threshold, don't save the wire.
    if (noise > peakAreaCut) {
        CaptInfo("Wire with large peak noise: " << deconv.GetChannelId()
                 << "  noise: " << noise);
        // return wireCharge;
    }

    // Now use fWork to mask out peaks that are part of another peak.
    std::fill(&fWork[0],&fWork[deconv.GetSampleCount()], 0.0);
    
    typedef std::pair<double, std::size_t> PeakPosition;

    // Find all the possible peak positions and sort them so the biggest are
    // last.
    std::vector<PeakPosition> peakCandidates;
    // Find the candidate peaks.
    for (std::size_t i = fDigitEndSkip;
         i<deconv.GetSampleCount()-fDigitEndSkip;
         ++i) {
        if (deconv.GetSample(i) < deconv.GetSample(i-1)) continue;
        if (deconv.GetSample(i) < deconv.GetSample(i+1)) continue;
        if (deconv.GetSample(i) < peakMaximumCut) continue;
        peakCandidates.push_back(std::make_pair(deconv.GetSample(i), i));
    }
    std::sort(peakCandidates.begin(), peakCandidates.end());
    
    // The vector of peak positions that need to be made into hits.
    std::vector< std::pair<int,int> > peaks;
    
    // Look at all of the peak candidates and find the "real" hits.
    while (!peakCandidates.empty()) {
        // The height of the current candidate
        int candidateHeight = peakCandidates.back().first;
        int candidateIndex = peakCandidates.back().second;
        peakCandidates.pop_back();
        // Skip peaks that close to other peaks.
        if (fWork[candidateIndex] > 0.5) {
            continue;
        }
        // Cut peaks that are smaller than NoiseThresholdCut times the RMS of
        // the channel.
        if (candidateHeight < noise*fNoiseThresholdCut) continue;
        // Find the extent of the peak.
        std::pair<int, int> extent = PeakExtent(candidateIndex, deconv);
        // Find the charge in the peak candidate.
        double charge = 0.0;
        for (int j=extent.first; j<=extent.second; ++j) {
            double r = deconv.GetSample(j); 
            charge += r;            
        }
        // Find the fwhm.
        double fwhm = digitStep*PeakFWHM(candidateIndex, extent, deconv);

        // Apply a cut to the overall peak size. (The digit has had the
        // baseline remove and is in units of charge).
        if (candidateHeight<peakMaximumCut) {
            continue;
        }
        // Apply a cut on the charge in the peak.
        if (charge < peakAreaCut) {
            continue;
        }
        // Apply a cut on the overall peak width.
        if (fwhm < peakWidthCut) {
            continue;
        }
        // This looks like a real hit, so mark the samples as used, and save
        // the extent for later.
        for (int j = extent.first; j <= extent.second; ++j) {
            fWork[j] = 1.0;
        }
        peaks.push_back(extent);
        // Stop if we are getting to many peaks.  The peaks are built in order
        // of height, so this tends to keep the biggest hits.  This only cuts
        // hits on really noisy wires.
        if (fMaxPeaks > 0 && (std::size_t) fMaxPeaks <= peaks.size()) {
            CaptError("Found more than " << fMaxPeaks
                      << " hits for channel "
                      << deconv.GetChannelId());
            break;
        }
    }

    // Sort the peaks by sample number.
    if (!peaks.empty()) std::sort(peaks.begin(), peaks.end());

#ifdef CHECK_FOR_PEAK_OVERLAP
    std::vector< std::pair<int,int> >::iterator last = peaks.end();
    for(std::vector< std::pair<int,int> >::iterator p = peaks.begin();
        p != peaks.end(); ++p) {
        if (last != peaks.end()) {
            if (p->first <= last->second) {
                CaptError("Overlapping hits for "
                          << deconv.GetChannelId()
                          << " " << last->first << "--" << last->second
                          << " overlaps "
                          << " " << p->first << "--" << p->second);
            }
        }
        last = p;
    }
#endif

    // Make all the hits.  
    CP::TMakeWireHit makeHit(fCorrectElectronLifetime,
                             fCorrectCollectionEfficiency);
    for(std::vector< std::pair<int,int> >::iterator p = peaks.begin();
        p != peaks.end(); ++p) {
        CP::THandle<CP::THit> newHit
            = makeHit(deconv, digitStep, t0,
                      p->first, p->second,
                      false);
        if (!newHit) continue;
        wireCharge += newHit->GetCharge();
        hits.push_back(newHit);
    }

#ifdef  CHECK_FOR_OVERLAPPED_HITS
    for (CP::THitSelection::iterator h = hits.begin(); h != hits.end(); ++h) {
        for (CP::THitSelection::iterator i = h+1; i != hits.end(); ++i) {
            // The hits aren't the same channel, so they can't overlap.
            if ((*i)->GetChannelId() != (*h)->GetChannelId()) continue;
            // "i" starts after "h" ends, so they can't overlap.
            if ((*i)->GetTimeStart() >= (*h)->GetTimeStop()) continue;
            // "i" ends before "h" begins, so they can't overlap.
            if ((*i)->GetTimeStop() <= (*h)->GetTimeStart()) continue;
            CaptError("Overlapping hit on " << (*h)->GetChannelId());
            CaptError("   hit " << (*h)->GetTimeStart()
                      << " " << (*h)->GetTimeStop());
            CaptError("   hit " << (*i)->GetTimeStart()
                      << " " << (*i)->GetTimeStop());
        }
    }
#endif

    return wireCharge;
}
