#include <TPulseDigit.hxx>
#include <HEPUnits.hxx>
#include <CaptGeomId.hxx>

#include <TChannelInfo.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <memory>
#include <complex>

#include <eventLoop.hxx>

namespace CP {class TRootOutput;};

/// Make diagnostic histograms for a file.
class TCalibHistsLoop: public CP::TEventLoopFunction {
public:
    TCalibHistsLoop() {
        fFFT = NULL;
        fNoiseHist = NULL;
        fWireHist = NULL;
        fSampleCount = 0;
        fMaxPower = -1000;
        fPrintFile = "";
        fSampleTime = 500*unit::ns;
    }

    virtual ~TCalibHistsLoop() {};

    void Usage(void) {
    }

    void Initialize() {
        fNoiseHist = new TH1F("noiseHist",
                              "Noise for each Channel",
                              100.0, 0, 50.0);
        fNoiseHist->SetXTitle("ADC");
        fGaussHist = new TH1F("gaussHist",
                              "Gaussian Noise for each Channel",
                              100.0, 0, 50.0);
        fGaussHist->SetXTitle("ADC");
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
        std::cout << "Finalize the run " << std::endl;
        if (fPrintFile.size() < 1) return;

        gStyle->SetOptStat(0);
        if (fMaxPower>0.0) fWireHist->SetMaximum(0.0);
        fWireHist->Draw("colz");
        gPad->Update();
        gPad->Print(fPrintFile.c_str());
        
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
                CaptError("Non-pulse in drift digits");
                continue;
            }

            // Check that the FFT is properly created.
            int nSize = 2*(1+pulse->GetSampleCount()/2);
            double nyquistFreq = 0.5/fSampleTime;
            double deltaFreq = 2.0*nyquistFreq/nSize;
            if (!fFFT || nSize != fSampleCount) {
                if (!fFFT) delete fFFT;
                fFFT = TVirtualFFT::FFT(1, &nSize, "R2C M K");
                fSampleCount = nSize;
                fWireHist = new TH2F("wireHist",
                                       "Log10(Power) for all wires",
                                      drift->size(), 0.0, drift->size(),
                                      nSize/2-1, deltaFreq, nyquistFreq);
                fWireHist->SetXTitle("Channel");
                fWireHist->SetYTitle("Frequency (Hz)");
            }

            // Find the average pedestal and subtract it.
            double pedestal = 0.0;
            for (std::size_t i=0; i<pulse->GetSampleCount(); ++i) {
                double p = pulse->GetSample(i);
                pedestal += p;
            }
            pedestal /= pulse->GetSampleCount();
            
            for (std::size_t i=0; i<(std::size_t) nSize; ++i) {
                double p = 0.0;
                if (i<pulse->GetSampleCount()) {
                    p = pulse->GetSample(i)-pedestal;
                }
                fFFT->SetPoint(i,p);
            }
            
            fFFT->Transform();

            TH1F* fftHist 
                = new TH1F((pulse->GetChannelId().AsString()+"-fft").c_str(),
                           ("FFT for " 
                            + pulse->GetChannelId().AsString()).c_str(),
                           nSize/2-1, deltaFreq, nyquistFreq);
            fftHist->SetYTitle("Power (#frac{ADC^{2}}{#Delta#it{f}})");
            fftHist->SetXTitle("Frequency (Hz)");

            std::vector<double> power(nSize/2);
            double amp = 0.0;
            for (std::size_t i = 1; i<(std::size_t) nSize/2; ++i) {
                double rl, im;
                fFFT->GetPointComplex(i,rl,im);
                std::complex<double> c(rl,im);
                double p = 2.0*std::abs(c*c)/nSize/nSize;
                amp += p;
                fMaxPower = std::max(fMaxPower,std::log10(p));
                power[i] = p;
                fftHist->SetBinContent(i, p);
                if (p<0.01/nSize) p = 0.01/nSize;
                fWireHist->SetBinContent(d+1,i,std::log10(p));
            }

            fNoiseHist->Fill(std::sqrt(amp));
            std::sort(power.begin(), power.end());
            int index = 0.68*power.size();
            amp = power[index]*nSize;
            fGaussHist->Fill(std::sqrt(amp));
            
        }

        return true;
    }

private:
    /// An fft to get frequencies.
    TVirtualFFT* fFFT;

    /// The number of samples in the FFT.
    int fSampleCount;

    /// The sample time.
    double fSampleTime;

    /// The power in the maximum bin, used to set a histogram max for display
    /// purposes.
    double fMaxPower;

    /// A histogram of the noise RMS on each wire.
    TH1F* fNoiseHist;

    /// A histogram of the noise on each wire assuming it's Gaussian
    TH1F* fGaussHist;

    /// The FFT for all wires in an event.
    TH2F* fWireHist;

    /// A file to draw the wire histogram in (usually png).
    std::string fPrintFile;
    
};

int main(int argc, char **argv) {
    TCalibHistsLoop userCode;
    CP::eventLoop(argc,argv,userCode,1);
}
