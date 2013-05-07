#include "TWireMakeHits.hxx"

#include <THitSelection.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TDataHit.hxx>
#include <HEPUnits.hxx>
#include <THandle.hxx>
#include <TCaptLog.hxx>

#include <TSpectrum.h>

#include <cmath>
#include <vector>

CP::TWireMakeHits::TWireMakeHits() {
    fSource = NULL;
    fDest = NULL;
    fNSource = 0;
}
CP::TWireMakeHits::~TWireMakeHits() {
    if (fSource) delete fSource;
    if (fDest) delete fDest;
}

void CP::TWireMakeHits::operator() (CP::THitSelection& hits, 
                                    const CP::TCalibPulseDigit& digit) {

    CP::TGeometryId geomId 
        = CP::TChannelInfo::Get().GetGeometry(digit.GetChannelId());

    // Make sure we have enough memory allocated for the spectrum.
    if (fNSource < (int) digit.GetSampleCount()) {
        if (fSource) delete fSource;
        fNSource = 2*digit.GetSampleCount();
        fSource = new float[fNSource];
        fDest = new float[fNSource];

    }
    
    // Fill the spectrum.
    for (std::size_t i = 0; i<digit.GetSampleCount(); ++i) {
        fSource[i] = digit.GetSample(i);
        fDest[i] = 0.0;
    }

    TSpectrum* spectrum = new TSpectrum;
    double sigma = 4.0;
    double threshold = 15;
    bool removeBkg = false;
    int iterations = 3;
    bool useMarkov = false;
    int window = 3;
    
    int found = spectrum->SearchHighRes(fSource,fDest,digit.GetSampleCount(),
                                       sigma, threshold,
                                       removeBkg, iterations, 
                                       useMarkov, window);

    // Get the peak positions and sort by time.
    float* xx = spectrum->GetPositionX();
    std::vector<float> peaks;
    for (int i=0; i<found; ++i)  peaks.push_back(xx[i]);

    // UGLY CODE ALERT.  This does not make me proud.
    for (std::vector<float>::iterator p = peaks.begin();
         p != peaks.end(); ++p) {
        // No peaks at the end of a digit.
        if (*p < 3) continue;
        if (*p > digit.GetSampleCount()-3) continue;
        double center = *p;
        int i = *p + 0.5;
        int b = 0;
        int e = digit.GetSampleCount();
        double peak = digit.GetSample(i);
        peak += digit.GetSample(i-1);
        peak += digit.GetSample(i+1);
        peak /= 3.0;
        if (peaks.size() > 1) {
            for (std::vector<float>::iterator o = peaks.begin(); 
                 o != peaks.end(); ++o) {
                if (o == p) continue;
                if (*o < *p) {
                    int t = (*p+*o)/2;
                    if (b<t) b = t;
                }
                else {
                    int t = (*p+*o)/2;
                    if (e>t) e = t;
                }
            }
        }
        double q = 0.0;
        double t = 0.0;
        double tt = 0.0;
        for (int j = b; j < e; ++j) {
            double v = digit.GetSample(j);
            if (v < 0.1*peak) continue;
            q += v;
            t += v*j;
            tt += v*j*j;
        }
        double qUnc = std::sqrt(q);
        center *= 250.0*unit::ns;

        t /= q;
        tt /= q;

        double rms = (tt - t*t);
        rms = 250.0*unit::ns*std::sqrt(rms + 1);
        double tUnc = rms + std::abs(center - 250.0*unit::ns*t);
        tUnc /= qUnc;
        center += digit.GetFirstSample();

        CP::TWritableDataHit hit;
        hit.SetGeomId(geomId);
        hit.SetDigit(digit.GetParent());
        hit.SetCharge(q);
        hit.SetChargeUncertainty(qUnc);
        hit.SetTime(center);
        hit.SetTimeRMS(rms);
        hit.SetTimeUncertainty(tUnc);
        hits.push_back(CP::THandle<CP::TDataHit>(new CP::TDataHit(hit)));
    }

}
