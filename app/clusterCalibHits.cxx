#include "GaussianNoise.hxx"
#include "FindPedestal.hxx"

#include <eventLoop.hxx>
#include <HEPUnits.hxx>
#include <CaptGeomId.hxx>
#include <THitSelection.hxx>
#include <TTPCChannelId.hxx>

#include <TChannelInfo.hxx>

#include <TH1F.h>
#include <TH2F.h>

#include <TCanvas.h>
#include <TStyle.h>

#include <memory>
#include <sstream>
#include <algorithm>

namespace CP {class TRootOutput;};

/// Make diagnostic histograms for a file.  This summarizes the hits that are
/// found on each channel.
class TCalibHistsLoop: public CP::TEventLoopFunction {
public:
    TCalibHistsLoop() {
        fPrintFile = "";
        fGeomUHits = NULL;
        fGeomVHits = NULL;
        fGeomXHits = NULL;
    }

    virtual ~TCalibHistsLoop() {};

    void Usage(void) {
        std::cout << "   -O draw=<base-name> Output histograms to pdf file."
                  << std::endl;
    }

    void Initialize() {
    }
    
    virtual bool SetOption(std::string option,std::string value="") {
        if (option=="draw") {
            fPrintFile=value;
            return true;
        }
        return false;
    }

    void Finalize(CP::TRootOutput * const file) {}

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
        
        if (!fGeomUHits) {
            fGeomUHits = new TH1F("geomUHits",
                                  "Hits on the U plane",
                                  338, 0, 338);
            fGeomVHits = new TH1F("geomVHits",
                                  "Hits on the V plane",
                                  338, 0, 338);
            fGeomXHits = new TH1F("geomXHits",
                                  "Hits on the X plane",
                                  338, 0, 338);
            const int channels = 4096;
            fUnconnectedHits = new TH1F("unconnectedHits",
                                        "Hits on Unconnected Channels",
                                        1, 0, 1);

            fGeomUCharge = new TH2F("geomUCharge",
                                    "Charge on the U plane",
                                    338, 0, 338,
                                    50, 0, 50000);
            fGeomVCharge = new TH2F("geomVCharge",
                                    "Charge on the V plane",
                                    338, 0, 338,
                                    50, 0, 50000);
            fGeomXCharge = new TH2F("geomXCharge",
                                    "Charge on the X plane",
                                    338, 0, 338,
                                    50, 0, 50000);
            fUnconnectedCharge = new TH2F("unconnectedCharge",
                                          "Charge on Unconnected Channels",
                                          1, 0, 1,
                                          50, 0, 50000);
        }

        for (CP::THitSelection::iterator h = drift->begin();
             h != drift->end(); ++h) {
            CP::TGeometryId id = (*h)->GetGeomId();
            int w = CP::GeomId::Captain::GetWireNumber(id);
            double q = (*h)->GetCharge();
            if (q > 49999) q = 49999;
            if (CP::GeomId::Captain::GetWirePlane(id)
                == CP::GeomId::Captain::kUPlane) {
                fGeomUHits->Fill(w);
                fGeomUCharge->Fill(w,q);
            }
            else if (CP::GeomId::Captain::GetWirePlane(id)
                == CP::GeomId::Captain::kVPlane) {
                fGeomVHits->Fill(w);
                fGeomVCharge->Fill(w,q);
            }
            else if (CP::GeomId::Captain::GetWirePlane(id)
                == CP::GeomId::Captain::kXPlane) {
                int w = CP::GeomId::Captain::GetWireNumber(id);
                fGeomXHits->Fill(w);
                fGeomXCharge->Fill(w,q);
            }
            else {
                CP::TTPCChannelId cid = (*h)->GetChannelId();
                std::cout << cid.AsString()
                          << std::endl;
                fUnconnectedHits->Fill(cid.AsString().c_str(),1.0);
                fUnconnectedCharge->Fill(cid.AsString().c_str(),q,1.0);
            }
        }
        
        return true;
    }

private:
    /// The prefix for file names draw the wire histograms in
    std::string fPrintFile;

    /// A histograms of the number of hits by geometry wire.
    TH1F* fGeomUHits;
    TH1F* fGeomVHits;
    TH1F* fGeomXHits;
    TH1F* fUnconnectedHits;
    
    TH2F* fGeomUCharge;
    TH2F* fGeomVCharge;
    TH2F* fGeomXCharge;
    TH2F* fUnconnectedCharge;

};

int main(int argc, char **argv) {
    TCalibHistsLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
