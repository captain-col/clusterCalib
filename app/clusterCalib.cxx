#include "TClusterCalib.hxx"
#include "TActivityFilter.hxx"

#include <eventLoop.hxx>

#include <sstream>

/// Run the cluster calibration on a file.
class TClusterCalibLoop: public CP::TEventLoopFunction {
public:
    TClusterCalibLoop() {
        fClusterCalib = NULL;
        fActivityFilter = NULL;
        fSavePulses = false;
        fApplyDriftCalibration = false;
        fApplyEfficiencyCalibration = true;
    }

    virtual ~TClusterCalibLoop() {};

    void Usage(void) {
        std::cout << "   -O filter[=c:S:s] Apply a simple event activity filter"
                  << std::endl
                  << "        c: Number of required channels"
                  << std::endl
                  << "        S: Required ADC above baseline"
                  << std::endl
                  << "        s: Required significance above baseline"
                  << std::endl;
        std::cout << "   -O pulse      "
                  << "Save the calibrated (after deconvolution) pulses"
                  << std::endl;
        std::cout << "   -O no-pulse   Don't save the calibrated pulses"
                  << std::endl;
        std::cout << "   -O drift      Apply the drift calibration (default)"
                  << std::endl;
        std::cout << "   -O no-drift   Don't apply the drift calibration"
                  << std::endl;
        std::cout << "   -O efficiency     Apply the efficiency (default)"
                  << std::endl;
        std::cout << "   -O no-efficiency  Don't apply the efficiency"
                  << std::endl;
    }

    virtual bool SetOption(std::string option,std::string value="") {
        if (option == "no-pulse") fSavePulses = false;
        else if (option == "pulse") fSavePulses = true;
        else if (option == "no-drift") {
            CaptLog("Do not apply the electron lifetime calibration");
            fApplyDriftCalibration = false;
        }
        else if (option == "drift") fApplyDriftCalibration = true;
        else if (option == "no-efficiency") {
            CaptLog("Do not apply the collection efficiency calibration");
            fApplyEfficiencyCalibration = false;
        }
        else if (option == "efficiency") fApplyEfficiencyCalibration = true;
        else if (option == "filter") {
            fActivityFilter = new CP::TActivityFilter();
            if (value!="") {
                std::istringstream vStr(value);
                int chan;
                vStr >> chan;
                if (chan>0) fActivityFilter->SetRequiredChannels(chan);
                char colon;
                vStr >> colon;
                if (colon != ':') return false;
                double v;
                vStr >> v;
                if (v>0) fActivityFilter->SetMinimumSignal(v);
                vStr >> colon;
                if (colon != ':') return false;
                vStr >> v;
                if (v>0) fActivityFilter->SetRequiredSignificance(v);
            }
        }
        else return false;
        return true;
    }

    bool operator () (CP::TEvent& event) {
        if (!fClusterCalib) {
            fClusterCalib = new CP::TClusterCalib();
            // Set the action for the calibrated pulse digits.
            fClusterCalib->SaveCalibratedPulses(fSavePulses);
            fClusterCalib->ApplyDriftCalibration(fApplyDriftCalibration);
            fClusterCalib->ApplyEfficiencyCalibration(
                fApplyEfficiencyCalibration);
        }

        // Possibly run a filter to reject noise events using uncalibrated
        // data.
        if (fActivityFilter) {
            bool result = (*fActivityFilter)(event);
            CaptLog("Reject " << event.GetContext());
            if (!result) return false;
        }
        
        // Run the calibration on the event.
        (*fClusterCalib)(event);

        // Save everything that gets to here.
        return true;
    }

private:
    CP::TClusterCalib* fClusterCalib;
    CP::TActivityFilter* fActivityFilter;
    
    bool fSavePulses;
    bool fApplyDriftCalibration;
    bool fApplyEfficiencyCalibration;

};

int main(int argc, char **argv) {
    TClusterCalibLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
