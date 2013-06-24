#include "TClusterCalib.hxx"

#include <eventLoop.hxx>

/// Run the cluster calibration on a file.
class TClusterCalibLoop: public CP::TEventLoopFunction {
public:
    TClusterCalibLoop() {
        fClusterCalib = NULL;
    }

    virtual ~TClusterCalibLoop() {};

    void Usage(void) {     }

    virtual bool SetOption(std::string option,std::string value="") {
        return true;
    }

    bool operator () (CP::TEvent& event) {
        // Make sure the electronics simulated is created.
        if (!fClusterCalib) fClusterCalib = new CP::TClusterCalib();

        // Run the simulation on the event.
        (*fClusterCalib)(event);

        // Save everything.
        return true;
    }

private:
    CP::TClusterCalib* fClusterCalib;
};

int main(int argc, char **argv) {
    TClusterCalibLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
