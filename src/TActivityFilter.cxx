#include "TActivityFilter.hxx"

#include <THitSelection.hxx>
#include <CaptGeomId.hxx>
#include <HEPUnits.hxx>

#include <vector>

CP::TActivityFilter::TActivityFilter() {

}

CP::TActivityFilter::~TActivityFilter() {

}

bool CP::TActivityFilter::operator() (CP::TEvent& event) {
    CaptLog("Filter event " << event.GetContext());

    CP::THandle<CP::THitSelection> drift(event.GetHits("drift"));   
    CP::THitSelection uHits("uHits");
    CP::THitSelection vHits("vHits");
    CP::THitSelection xHits("xHits");
    
    for (CP::THitSelection::iterator h = drift->begin(); h != drift->end();
         ++h) {
        int plane = CP::GeomId::Captain::GetWirePlane((*h)->GetGeomId());
        if (plane == CP::GeomId::Captain::kUPlane) {
            uHits.push_back(*h);
        }
        if (plane == CP::GeomId::Captain::kVPlane) {
            vHits.push_back(*h);
        }
        if (plane == CP::GeomId::Captain::kXPlane) {
            xHits.push_back(*h);
        }
    }

    CaptLog("Hits found: " << uHits.size()
            << " " << vHits.size()
            << " " << xHits.size()
            << " " << drift->size());

    std::vector<double> vTimes;
    int vNeighbors = 0;
    for (CP::THitSelection::iterator h = vHits.begin();
         h != vHits.end(); ++h) {
        int w1 = CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
        int timeNeighbors = 0;
        int wireNeighbors = 0;
        for (CP::THitSelection::iterator i = h;
             i != vHits.end(); ++i) {
            int w2 = CP::GeomId::Captain::GetWireNumber((*i)->GetGeomId());
            int d = w2-w1;
            if (d == 0) continue;
            double dT = std::abs((*h)->GetTime() - (*i)->GetTime());
            if (dT > 10.0*unit::microsecond) continue;
            ++timeNeighbors;
            if (1 < std::abs(d)) continue;
            ++wireNeighbors;
        }
        if (timeNeighbors > 3) continue;
        if (wireNeighbors < 1) continue;
        ++vNeighbors;
        vTimes.push_back((*h)->GetTime());
    }

    std::vector<double> uTimes;
    int uNeighbors = 0;
    for (CP::THitSelection::iterator h = uHits.begin();
         h != uHits.end(); ++h) {
        int w1 = CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
        int timeNeighbors = 0;
        int wireNeighbors = 0;
        for (CP::THitSelection::iterator i = h;
             i != uHits.end(); ++i) {
            int w2 = CP::GeomId::Captain::GetWireNumber((*i)->GetGeomId());
            int d = w2-w1;
            if (d == 0) continue;
            double dT = std::abs((*h)->GetTime() - (*i)->GetTime());
            if (dT > 10.0*unit::microsecond) continue;
            ++timeNeighbors;
            if (1 < std::abs(d)) continue;
            ++wireNeighbors;
        }
        if (timeNeighbors > 4) continue;
        if (wireNeighbors < 1) continue;
        ++uNeighbors;
        uTimes.push_back((*h)->GetTime());
    }

    CaptLog("Neighbors " << uNeighbors << " " << vNeighbors);
    int planeNeighbors = 0;
    double avgT = 0.0;
    double avgW = 0.0;
    double minT = 1E+20;
    for (std::vector<double>::iterator uT = uTimes.begin();
         uT != uTimes.end(); ++uT) {
        for (std::vector<double>::iterator vT = vTimes.begin();
             vT != vTimes.end(); ++vT) {
            double dT = std::abs((*uT)-*(vT));
            minT = std::min(minT,dT);
            if (dT > 20*unit::microsecond) continue;
            ++ planeNeighbors;
            avgT += (*vT);
            avgW += 1.0;
        }
    }
    if (avgW > 0) avgT /= avgW;
    else avgT = 0;
    CaptLog("Plane Neighbors " << planeNeighbors
            << " " << avgT/unit::microsecond
            << " " << minT/unit::microsecond);

    if (planeNeighbors > 0) return true;
    return false;
}
