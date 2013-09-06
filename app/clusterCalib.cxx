#include "TClusterCalib.hxx"

#include <eventLoop.hxx>

/// Run the cluster calibration on a file.
class TClusterCalibLoop: public CP::TEventLoopFunction {
public:
    TClusterCalibLoop() {
        fClusterCalib = NULL;
        fSavePulses = true;
    }

    virtual ~TClusterCalibLoop() {};

    void Usage(void) {
        std::cout << "    -O pulse     Save the calibrated pulses"
                  << std::endl;
    }

    virtual bool SetOption(std::string option,std::string value="") {
        if (value != "") return false;
        if (option == "nopulse") fSavePulses = false;
        else if (option == "pulse") fSavePulses = true;
        return true;
    }

    bool operator () (CP::TEvent& event) {
        // Make sure the electronics simulated is created.
        if (!fClusterCalib) fClusterCalib = new CP::TClusterCalib();

        // Set the action for the calibrated pulse digits.
        fClusterCalib->SaveCalibratedPulses(fSavePulses);

        // Run the simulation on the event.
        (*fClusterCalib)(event);

        // Save everything.
        return true;
    }

private:
    CP::TClusterCalib* fClusterCalib;

    bool fSavePulses;
};

int main(int argc, char **argv) {
    TClusterCalibLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
