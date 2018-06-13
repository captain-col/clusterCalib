#include "GaussianNoise.hxx"
#include "FindPedestal.hxx"

#include <eventLoop.hxx>
#include <TPulseCalib.hxx>
#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TChannelCalib.hxx>
#include <TPulseCalib.hxx>
#include <TDigitProxy.hxx>
#include <HEPUnits.hxx>
#include <CaptGeomId.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <memory>
#include <complex>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace CP {class TRootOutput;};

/// Make diagnostic histograms for a file.
class TCalibHistsLoop: public CP::TEventLoopFunction {
public:
    TCalibHistsLoop() {
        fFFT = NULL;
        fSampleCount = 0;
        fPrintFile = "";
        fSampleTime = 500*unit::ns;
        fCalibrate = new CP::TPulseCalib();
        fAvgFFTHist = NULL;
        fAvgPowerHist = NULL;
        fAvgPhaseHist = NULL;
    }

    virtual ~TCalibHistsLoop() {};

    void Usage(void) {
        std::cout << "   -O draw=<base-name> Output histograms to pdf file."
                  << std::endl;
    }

    void Initialize() {
    }
    
    virtual bool SetOption(std::string option,std::string value="") {
        std::cout << "OPTION: " << option << " " << value << std::endl;
        if (option=="draw") {
            fPrintFile=value;
            return true;
        }
        return false;
    }

    void Finalize(CP::TRootOutput * const file) {
        std::cout << "Finalize the run " << fAvgEntries << std::endl;

        for (int i=0; i<fAvgFFT.size()/2; ++i) {
            fAvgFFT[i] /= fAvgEntries;
            fAvgFFTHist->SetBinContent(i,std::abs(fAvgFFT[i]));
        }
        
        if (fPrintFile.size() < 1) return;

        TCanvas c("clusterCalibFFT");

        gPad->Print((fPrintFile + ".pdf[").c_str());
        gPad->Print(((fPrintFile+".pdf") + "]").c_str());
    }

    bool operator () (CP::TEvent& event) {
        CP::THandle<CP::TDigitContainer> drift
            = event.Get<CP::TDigitContainer>("~/digits/drift");

        if (!drift) {
            CaptError("No drift signals for this event " << event.GetContext());
            return false;
        }

        CP::TChannelInfo& chanInfo = CP::TChannelInfo::Get();
        chanInfo.SetContext(event.GetContext());

        std::ostringstream titlePrefix;
        titlePrefix << "Event " << event.GetContext().GetRun()
                    << "." << event.GetContext().GetEvent() << ": ";
        
        for (std::size_t d = 0; d < drift->size(); ++d) {
            const CP::TPulseDigit* pulse 
                = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
            if (!pulse) {
                CaptError("Non-pulse in drift digits");
                continue;
            }

            // Check that this is a good channel attached to a wire.
            CP::TChannelCalib channelCalib;
            if (!channelCalib.IsGoodChannel(pulse->GetChannelId())) continue;
            CP::TGeometryId pulseGeom
                = CP::TChannelInfo::Get().GetGeometry(pulse->GetChannelId());
            if (!pulseGeom.IsValid()) continue;

            CP::TDigitProxy proxy(*drift,d);
            std::unique_ptr<CP::TCalibPulseDigit> calib((*fCalibrate)(proxy));
            if (!calib.get()) {
                CaptError("Channel not calibrated ");
                continue;
            }

            // Check that the FFT is properly created.
            int nSize = 2*(1+pulse->GetSampleCount()/2);
            fNyquistFreq = (0.5/fSampleTime)/unit::hertz;
            fDeltaFreq = 2.0*fNyquistFreq/nSize;
            if (!fFFT || nSize != fSampleCount) {
                if (!fFFT) delete fFFT;
                CaptLog("Nyquist frequency: " << fNyquistFreq);
                CaptLog("Frequency bin size: " << fDeltaFreq);
                fFFT = TVirtualFFT::FFT(1, &nSize, "R2C M K");
                fSampleCount = nSize;
                fAvgFFT.resize(nSize);
                fAvgEntries = 0;
                fAvgFFTHist
                    = new TH1F("averageFFT", "Average FFT",
                               nSize/2, 0.0, fNyquistFreq);
            }

            int wire = chanInfo.GetWireNumber(pulse->GetChannelId());

            // Find the pedestal for the wire
            double pedestal = CP::FindPedestal(pulse->begin(), pulse->end());

            // Fill the FFT and decompose.
            for (std::size_t i=0; i<(std::size_t) nSize; ++i) {
                double p = 0.0;
                if (i<pulse->GetSampleCount()) {
                    p = pulse->GetSample(i)-pedestal;
                }
                fFFT->SetPoint(i,p);
            }
            fFFT->Transform();

            std::ostringstream histBaseName;
            if (0 <= wire) {
                histBaseName << "W:"
                             << std::setw(4) 
                             << std::setfill('0') << wire<<"-";
            }
            std::string chanName;
            if (!pulse->GetChannelId().IsMCChannel()) {
                chanName = pulse->GetChannelId().AsString().substr(4);
            }
            else {
                chanName = pulse->GetChannelId().AsString();
            }            
            histBaseName << chanName;

            std::ostringstream histBaseTitle;
            if (0 <= wire) histBaseTitle << "wire number " << wire << " ";
            histBaseTitle << chanName;
            // std::cout << histBaseTitle.str() << " " << nSize << std::endl;

            if (!fAvgPowerHist) {
                fAvgPowerHist
                    = new TProfile((histBaseName.str()+"-power").c_str(),
                                   (titlePrefix.str() 
                                    + "Average Power for " 
                                    + histBaseTitle.str()).c_str(),
                                   nSize/2, 0.0, fNyquistFreq);
            }

            if (!fAvgPhaseHist) {
                fAvgPhaseHist
                    = new TProfile((histBaseName.str()+"-phase").c_str(),
                                   (titlePrefix.str() 
                                    + "Average Phase for " 
                                    + histBaseTitle.str()).c_str(),
                                   nSize/2, 0.0, fNyquistFreq);
            }
            
            

            for (std::size_t i = 0; i<nSize; ++i) {
                double rl, im;
                fFFT->GetPointComplex(i,rl,im);
                std::complex<double> fft(rl,im);
                fAvgPowerHist->Fill(fDeltaFreq*(i+0.5),std::abs(fft));
                fAvgPhaseHist->Fill(fDeltaFreq*(i+0.5),std::arg(fft),
                                    std::abs(fft));
                fAvgFFT[i] += fft;
            }
            ++fAvgEntries;
        }

        return false;
    }

private:
    /// An fft to get frequencies.
    TVirtualFFT* fFFT;

    int fSampleCount;

    /// The sample time.
    double fSampleTime;

    /// The nyquist frequency
    double fNyquistFreq;

    /// The frequency step between FFT samples
    double fDeltaFreq;
    
    /// A file to draw the wire histogram in (usually png).
    std::string fPrintFile;

    /// The class to calibration the digit.
    CP::TPulseCalib* fCalibrate;

    /// The average FFT
    std::vector< std::complex<double> > fAvgFFT;

    /// The number of entries added to the average.
    int fAvgEntries;
    
    TH1F* fAvgFFTHist;
    
    /// The average power at each frequency
    TProfile* fAvgPowerHist;

    /// The average phase at each frequency
    TProfile* fAvgPhaseHist;
};

int main(int argc, char **argv) {
    TCalibHistsLoop userCode;
    CP::eventLoop(argc,argv,userCode,1);
}
