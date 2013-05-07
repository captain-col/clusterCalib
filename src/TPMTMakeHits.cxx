#include "TPMTMakeHits.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TDataHit.hxx>
#include <HEPUnits.hxx>
#include <THandle.hxx>
#include <TCaptLog.hxx>

#include <cmath>

CP::TPMTMakeHits::TPMTMakeHits() {}
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
        double s = digit.GetSample(i);
        if (s < 1E+5)  continue;
        CP::TWritableDataHit hit;
        hit.SetGeomId(geomId);
        hit.SetDigit(digit.GetParent());
        hit.SetCharge(s);
        hit.SetChargeUncertainty(std::sqrt(s));
        hit.SetTime(digitStep*i);
        hit.SetTimeUncertainty(digitStep/std::sqrt(12.0));
        hit.SetTimeRMS(digitStep/std::sqrt(12.0));
        hits.push_back(CP::THandle<CP::TDataHit>(new CP::TDataHit(hit)));
    }

}

