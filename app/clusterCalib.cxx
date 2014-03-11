#include "TClusterCalib.hxx"

#include <eventLoop.hxx>

/// Run the cluster calibration on a file.
class TClusterCalibLoop: public CP::TEventLoopFunction {
public:
    TClusterCalibLoop() {
        fClusterCalib = NULL;
        fSavePulses = true;
        fApplyDriftCalibration = true;
    }

    virtual ~TClusterCalibLoop() {};

    void Usage(void) {
        std::cout << "    -O pulse     Save the calibrated pulses"
                  << std::endl;
        std::cout << "    -O nopulse   Don't save the calibrated pulses"
                  << std::endl;
        std::cout << "    -O drift     Apply the drift calibration (default)"
                  << std::endl;
        std::cout << "    -O nodrift   Don't apply the drift calibration"
                  << std::endl;
    }

    virtual bool SetOption(std::string option,std::string value="") {
        if (value != "") return false;
        if (option == "nopulse") fSavePulses = false;
        else if (option == "pulse") fSavePulses = true;
        else if (option == "nodrift") fApplyDriftCalibration = false;
        else if (option == "drift") fApplyDriftCalibration = true;
        else return false;
        return true;
    }

    bool operator () (CP::TEvent& event) {
        // Make sure the electronics simulated is created.
        if (!fClusterCalib) fClusterCalib = new CP::TClusterCalib();

        // Set the action for the calibrated pulse digits.
        fClusterCalib->SaveCalibratedPulses(fSavePulses);
        fClusterCalib->ApplyDriftCalibration(fApplyDriftCalibration);

        // Run the simulation on the event.
        (*fClusterCalib)(event);

        // Save everything.
        return true;
    }

private:
    CP::TClusterCalib* fClusterCalib;

    bool fSavePulses;
    bool fApplyDriftCalibration;
};

int main(int argc, char **argv) {
    TClusterCalibLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
