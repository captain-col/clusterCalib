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
        fSaveCalib = false;
        fSaveDecorrel = false;
        fSaveDeconv = false;
        fCalibrateAllChannels = false;
        fApplyEfficiencyCalibration = true;
        fRemoveCorrelatedPedestal = true;
    }

    virtual ~TClusterCalibLoop() {};

    void Usage(void) {
        std::cout << "   -O all       Calibrate all channels,"
                  << " including disconected"
                  << std::endl;
        std::cout << "   -O filter[=S]:h]:s:]H] Apply an event activity filter"
                  << std::endl
                  << "        S: Override number of ADC above baseline"
                  << std::endl
                  << "        h: Override number of hits in cluster"
                  << std::endl
                  << "        s: Override significance above baseline"
                  << std::endl
                  << "        H: Override maximum allowed number of hits"
                  << std::endl;
        std::cout << "   -O save-deconv     "
                  << "Save the calibrated (after deconvolution) pulses"
                  << std::endl;
        std::cout << "   -O save-decorrel   "
                  << "Save the calibrated (after decorrelation) pulses"
                  << std::endl;
        std::cout << "   -O save-calib      "
                  << "Save the calibrated pulses"
                  << std::endl;
        std::cout << "   -O efficiency     Apply the efficiency (default)"
                  << std::endl;
        std::cout << "   -O no-efficiency  Don't apply the efficiency"
                  << std::endl;
        std::cout << "   -O no-correlation  No extra pedestal handling"
                  << std::endl;
        std::cout << "   -O correlation     "
                  << "Remove wire to wire correlations"
                  << " in the pedestal (slow)"
                  << std::endl;
    }

    virtual bool SetOption(std::string option,std::string value="") {
        if (option == "save-calib") fSaveCalib = true;
        else if (option == "save-decorrel") fSaveDecorrel = true;
        else if (option == "save-deconv") fSaveDeconv = true;
        else if (option.find("no-corr") == 0) {
            fRemoveCorrelatedPedestal = false;
        }
        else if (option.find("corr") == 0) {
            fRemoveCorrelatedPedestal = true;
        }
        else if (option == "no-efficiency") {
            CaptLog("Do not apply the collection efficiency calibration");
            fApplyEfficiencyCalibration = false;
        }
        else if (option == "efficiency") fApplyEfficiencyCalibration = true;
        else if (option == "all") fCalibrateAllChannels = true;
        else if (option == "filter") {
            fActivityFilter = new CP::TActivityFilter();
            if (value!="") {
                std::istringstream vStr(value);
                double v;
                vStr >> v;
                if (v>0) fActivityFilter->SetMinimumSignal(v);
                char colon;
                vStr >> colon;
                if (colon != ':') return true;

                int chan;
                vStr >> chan;
                if (chan>0) fActivityFilter->SetRequiredHits(chan);
                vStr >> colon;
                if (colon != ':') return true;

                vStr >> v;
                if (v>0) fActivityFilter->SetRequiredSignificance(v);
                vStr >> colon;
                if (colon != ':') return true;
                
                vStr >> chan;
                if (chan>0) fActivityFilter->SetMaximumAllowedHits(chan);
            }
        }
        else return false;
        return true;
    }

    bool operator () (CP::TEvent& event) {
        if (!fClusterCalib) {
            fClusterCalib = new CP::TClusterCalib();
            // Set the action for the calibrated pulse digits.
            fClusterCalib->SaveCalibratedPulses(fSaveCalib);
            fClusterCalib->SaveDecorrelatedPulses(fSaveDecorrel);
            fClusterCalib->SaveDeconvolvedPulses(fSaveDeconv);
            fClusterCalib->CalibrateAllChannels(fCalibrateAllChannels);
            fClusterCalib->ApplyEfficiencyCalibration(
                fApplyEfficiencyCalibration);
            fClusterCalib->RemoveCorrelations(fRemoveCorrelatedPedestal);
        }

        // Possibly run a filter to reject noise events using uncalibrated
        // data.
        if (fActivityFilter) {
            bool result = (*fActivityFilter)(event);
            if (!result) {
                CaptLog("Reject " << event.GetContext());
                return false;
            }
        }
        
        // Run the calibration on the event.
        (*fClusterCalib)(event);

        // Save everything that gets to here.
        return true;
    }

private:
    CP::TClusterCalib* fClusterCalib;
    CP::TActivityFilter* fActivityFilter;
    
    bool fSaveDeconv;
    bool fSaveDecorrel;
    bool fSaveCalib;
    bool fCalibrateAllChannels;
    bool fApplyEfficiencyCalibration;
    bool fRemoveCorrelatedPedestal;
};

int main(int argc, char **argv) {
    TClusterCalibLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
