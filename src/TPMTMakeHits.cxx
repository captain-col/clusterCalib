#include "TPMTMakeHits.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TDataHit.hxx>
#include <HEPUnits.hxx>
#include <THandle.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>

#include <cmath>

CP::TPMTMakeHits::TPMTMakeHits() {
    fThreshold 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "clusterCalib.pmt.threshold");
}
CP::TPMTMakeHits::~TPMTMakeHits() {}

void CP::TPMTMakeHits::operator() (CP::THitSelection& hits, 
                                   const CP::TCalibPulseDigit& digit) {
    CP::TGeometryId geomId 
        = CP::TChannelInfo::Get().GetGeometry(digit.GetChannelId());

    double digitStep = digit.GetLastSample()-digit.GetFirstSample();
    digitStep /= digit.GetSampleCount();

    // This should do a peak search, but the PMT peaks are really (really,
    // really) sharp, and extremely narrow.  This is probably good enough for
    // now.
    for (std::size_t i=0; i<digit.GetSampleCount(); ++i) {
        if (digit.GetSample(i) < fThreshold)  continue;
        double hitTime = digitStep*i + digit.GetFirstSample();
        std::size_t hitStart = i;
        if (i>2) hitStart -= 2;
        else hitStart = 0;
        std::size_t belowCount = 0;
        while (i < digit.GetSampleCount() && belowCount < 4) {
            if (digit.GetSample(i) > fThreshold) {
                belowCount = 0;
                ++i;
                continue;
            }
            ++belowCount;
            ++i;
        }
        std::size_t hitStop = i;
        double hitCharge = 0.0;
        for (std::size_t j = hitStart; j<hitStop; ++j) {
            hitCharge += digit.GetSample(j);
        }
        double hitLength = (hitStop-hitStart)*digitStep;
        CP::TWritableDataHit hit;
        hit.SetGeomId(geomId);
        hit.SetDigit(digit.GetParent());
        hit.SetCharge(hitCharge);
        hit.SetChargeUncertainty(std::sqrt(hitCharge));
        hit.SetTime(hitTime);
        hit.SetTimeUncertainty(digitStep/std::sqrt(12.0));
        hit.SetTimeRMS(hitLength/2.0);
        hits.push_back(CP::THandle<CP::TDataHit>(new CP::TDataHit(hit)));
    }

}

