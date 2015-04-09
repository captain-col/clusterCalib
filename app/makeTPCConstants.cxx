#include <eventLoop.hxx>
#include <TPulseDigit.hxx>

#include <sstream>

/// Run the cluster calibration on a file.
class TMakeTPCConstantsLoop: public CP::TEventLoopFunction {
public:
    TMakeTPCConstantsLoop() {
        fSaveEvents = false;
        fPulseVoltage = -1.0;
        fGainConfig = 2;
        fShapeConfig = 1;
    }

    virtual ~TMakeTPCConstantsLoop() {};

    void Usage(void) {
        std::cout << "    -O pulse=<mV>    Set the input pulse (required)"
                  << std::endl;
        std::cout << "    -O rise=<0..3>   Set the risetime config. (def: 1)"
                  << std::endl;
        std::cout << "    -O gain=<0..3>   Set the gain config. (def: 2)"
                  << std::endl;
        std::cout << "    -O save          Save the events (def: no-save)"
                  << std::endl;
        std::cout << "    -O no-save       Do not save the events."
                  << std::endl;

    }

    virtual bool SetOption(std::string option,std::string value="") {
        std::cout << "OPTION: " << option << " " << value << std::endl;
        if (option == "rise") {
            std::istringstream val(value);
            val >> fShapeConfig;
            return true;
        }
        else if (option == "gain") {
            std::istringstream val(value);
            val >> fGainConfig;
            return true;
        }
        else if (option == "pulse") {
            std::istringstream val(value);
            val >> fPulseVoltage;
            return true;
        }
        else if (option == "save") {
            fSaveEvents = true;
        }
        else if (option == "no-save") {
            fSaveEvents = false;
        }
        return false;
    }

    bool operator () (CP::TEvent& event) {
        CP::THandle<CP::TDigitContainer> drift
            = event.Get<CP::TDigitContainer>("~/digits/drift");

        if (!drift) {
            CaptError("No drift signals for this event " << event.GetContext());
            return false;
        }

        for (std::size_t d = 0; d < drift->size(); ++d) {
            const CP::TPulseDigit* pulse 
                = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
            if (!pulse) {
                CaptError("Non-pulse digit in drift digits");
                continue;
            }
        }
        return true;
    }

private:

    

    /// Flag to determine if the events are saved to an output file
    bool fSaveEvents;

    /// The calibration pulse voltage in millivolts
    double fPulseVoltage;

    /// The ASIC pulse shape configuration.  The configuration is: 0) 500 ns, 1)
    /// 1000ns, 2) ??? 3) ???
    int fShapeConfig;

    /// The ASIC gain configuration.  The configuration is: 0) 4.7 mV/fC 1)
    /// 7.8 mV/fC 2) 14 mV/fC 3) 25 mV/fC.
    int fGainConfig;
};

int main(int argc, char **argv) {
    TMakeTPCConstantsLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
