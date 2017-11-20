#include "GaussianNoise.hxx"
#include "FindPedestal.hxx"

#include <eventLoop.hxx>
#include <HEPUnits.hxx>
#include <CaptGeomId.hxx>
#include <THandle.hxx>
#include <THitSelection.hxx>
#include <TTPCChannelId.hxx>

#include <TChannelInfo.hxx>

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <TCanvas.h>
#include <TStyle.h>

#include <memory>
#include <sstream>
#include <algorithm>

namespace CP {class TRootOutput;};

namespace {
    const TVector3 PositionXY(const CP::THandle<CP::THit>& hit1,
                              const CP::THandle<CP::THit>& hit2) {
        double x1 = hit1->GetPosition().X();
        double y1 = hit1->GetPosition().Y();
        double dx1 = hit1->GetYAxis().X();
        double dy1 = hit1->GetYAxis().Y();
        
        double x2 = hit2->GetPosition().X();
        double y2 = hit2->GetPosition().Y();
        double dx2 = hit2->GetYAxis().X();
        double dy2 = hit2->GetYAxis().Y();
        
        // Solve 
        //      x1 + s1*dx1 = x2 + s2*dx2;
        //      y1 + s1*dy1 = y2 + s2*dy2; 
        // for s1, s2
        
        // Only the first shift is used, but here is the full solution (for s1
        // and s2):
        double s1 = -(dx2*(y1-y2)+dy2*x2-dy2*x1)/(dx2*dy1-dx1*dy2);
        // double s2 = -(dx1*(y1-y2)+dy1*x2-dy1*x1)/(dx2*dy1-dx1*dy2);
        
        return TVector3( (x1+s1*dx1), (y1+s1*dy1), 0.0);
    }
}

/// Validate the X wire channel map by comparing against the expected X
/// position of a wire based on the U and V wire crossing point. This uses the
/// idea that the U and V wire positions were validated against laser data in
/// previous runs while the X plane was not working.  This is run as a
/// standard event loop program, and needs to read a calibrated event file.
/// It produces a root file with histograms summarizing the information.
class TCalibHistsLoop: public CP::TEventLoopFunction {
public:
    TCalibHistsLoop() {
        fPrintFile = "";
    }

    virtual ~TCalibHistsLoop() {};

    void Usage(void) {
        std::cout << "   -O draw=<base-name> Output histograms to pdf file."
                  << std::endl;
    }

    void Initialize() {
        fXYPositions = new TH2F("xyPositions","X-Y position of U-V cross",
                                140,-700,700,
                                140,-700,700);
        fXYCrosses = new TH2F("xyCrosses","X-Y position of U-V-X cross",
                                140,-700,700,
                                140,-700,700);
        fUVXCrosses = new TH2F("uvxCrosses","X-Y position of U-V-X cross",
                                140,-700,700,
                                140,-700,700);
        fHitTimes = new TH1F("hitTimes", "hit times",
                             100, -1*unit::ms, 3*unit::ms);
        fXPositions = new TH2F("xPositions","X Position vs U&V Crossing Point",
                               332,0,332,
                               332,0,332);
        fXPositions->SetXTitle("Current X Wire Position");
        fXPositions->SetYTitle("Correct X Wire Position from U&V Wire Crossing Point");
        fXCorrections
            = new TH2F("xCorrections",
                       "Correction to X Based On U&V Crossing Point",
                       332,0,332,
                       101,-50,51
                );
        fXCorrections->SetXTitle("Current X Wire Position");
        fXCorrections->SetYTitle("Correction to X Wire Position");

        fXChannels
            = new TH2F("xChannels",
                       "Wire for Channel Based On U&V Crossing Point",
                       1,0,1,
                       332,0,332
                );
        fXChannels->SetYTitle("Wire");

}
    
    virtual bool SetOption(std::string option,std::string value="") {
        if (option=="draw") {
            fPrintFile=value;
            return true;
        }
        return false;
    }

    void Finalize(CP::TRootOutput * const file) {
        fXChannels->LabelsDeflate("x");
        fXChannels->LabelsOption("a","x");
    }

    bool operator () (CP::TEvent& event) {
        CP::THandle<CP::THitSelection> drift = event.GetHits("drift");
        if (!drift) {
            CaptError("No drift signals for this event " << event.GetContext());
            return false;
        }

        CP::TChannelInfo& chanInfo = CP::TChannelInfo::Get();
        chanInfo.SetContext(event.GetContext());

        std::ostringstream titlePrefix;
        titlePrefix << "Event " << event.GetContext().GetRun()
                    << "." << event.GetContext().GetEvent() << ": ";

        std::unique_ptr<CP::THitSelection> uHits(new CP::THitSelection());
        std::unique_ptr<CP::THitSelection> vHits(new CP::THitSelection());
        std::unique_ptr<CP::THitSelection> xHits(new CP::THitSelection());
        for (CP::THitSelection::iterator h = drift->begin();
             h != drift->end(); ++h) {
            CP::TGeometryId id = (*h)->GetGeomId();
            int w = CP::GeomId::Captain::GetWireNumber(id);
            switch(CP::GeomId::Captain::GetWirePlane(id)) {
            case CP::GeomId::Captain::kUPlane: 
                uHits->push_back(*h);
                break;
            case CP::GeomId::Captain::kVPlane:
                vHits->push_back(*h);
                break;
            case CP::GeomId::Captain::kXPlane:
                xHits->push_back(*h);
                break;
            default:
                std::cout << "Unknown wire plane" << std::endl;
            }
        }

        for (CP::THitSelection::iterator hu = uHits->begin();
             hu != uHits->end(); ++hu) {
            double tu = (*hu)->GetTime()
                - (*hu)->GetPosition().Z()/(1.6*unit::mm/unit::microsecond);
            if (tu < 0.0 || tu > 700*unit::microsecond) continue;
            CP::THitSelection::iterator cv = vHits->begin();
            double tv = (*cv)->GetTime()
                - (*cv)->GetPosition().Z()/(1.6*unit::mm/unit::microsecond);
            double ct = std::abs(tu-tv);
            for (CP::THitSelection::iterator hv = vHits->begin();
                 hv != vHits->end(); ++hv) {
                tv = (*hv)->GetTime()
                    - (*hv)->GetPosition().Z()/(1.6*unit::mm/unit::microsecond);
                if (tv < 0.0 || tv > 700*unit::microsecond) continue;
                double dt = std::abs(tu-tv);
                if (dt < ct) {
                    cv = hv;
                    ct = dt;
                }
            }
            if (ct > 1.5*unit::microsecond) continue;
            tv = (*cv)->GetTime()
                - (*cv)->GetPosition().Z()/(1.6*unit::mm/unit::microsecond);
            fHitTimes->Fill(tu);
            TVector3 uvPos = PositionXY((*hu),(*cv));
            fXYPositions->Fill(uvPos.X(),uvPos.Y());

            CP::THitSelection::iterator cx = xHits->begin();
            double tx = (*cx)->GetTime()
                - (*cx)->GetPosition().Z()/(1.6*unit::mm/unit::microsecond);
            ct = std::max(std::abs(tu-tx),std::abs(tv-tx));
            for (CP::THitSelection::iterator hx = xHits->begin();
                 hx != xHits->end(); ++hx) {
                tx = (*hx)->GetTime()
                    - (*hx)->GetPosition().Z()/(1.6*unit::mm/unit::microsecond);
                double dt = std::max(std::abs(tu-tx),std::abs(tv-tx));
                if (dt < ct) {
                    cx = hx;
                    ct = dt;
                }
            }
            if (ct > 1*unit::microsecond) continue;
            CP::TGeometryId id = (*cx)->GetGeomId();
            int x = CP::GeomId::Captain::GetWireNumber(id);
            int w = (uvPos.X()+497)/3.0;
            fXPositions->Fill(1.0*x+0.1,1.0*w + 0.1);
            CP::TTPCChannelId cid = (*cx)->GetChannelId();
            fXChannels->Fill(cid.AsString().c_str(),1.0*w + 0.1,1.0);
            fXCorrections->Fill(1.0*x+0.1,1.0*(w-x)+0.1);
            fXYCrosses->Fill(uvPos.X(),uvPos.Y());
            TVector3 uxPos = PositionXY((*hu),(*cx));
            TVector3 vxPos = PositionXY((*cv),(*cx));
            if ((uxPos-vxPos).Mag() < 10*unit::mm) {
                TVector3 avgPos = 0.33333333*(uvPos + uxPos + vxPos);
                fUVXCrosses->Fill(avgPos.X(),avgPos.Y());
            }
        }

        return false;
    }

private:
    /// The prefix for file names draw the wire histograms in
    std::string fPrintFile;

    TH2F *fXYPositions;

    TH2F *fXYCrosses;

    TH1F *fHitTimes;

    TH2F *fXPositions;

    TH2F *fXCorrections;

    TH2F *fXChannels;

    TH2F *fUVXCrosses;

};

int main(int argc, char **argv) {
    TCalibHistsLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
