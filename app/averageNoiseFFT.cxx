#include "FindPedestal.hxx"

#include <eventLoop.hxx>

#include <TPulseDigit.hxx>
#include <TChannelInfo.hxx>
#include <TChannelCalib.hxx>
#include <CaptGeomId.hxx>
#include <HEPUnits.hxx>

#include <TVirtualFFT.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TSpectrum.h>

#include <memory>
#include <complex>
#include <sstream>
#include <fstream>
#include <iostream>

namespace CP {class TRootOutput;};

/// Make diagnostic histograms for a file.
class TCalibHistsLoop: public CP::TEventLoopFunction {
public:
    TCalibHistsLoop() {
        fFFT = NULL;
        fSampleCount = 0;
        fSampleTime = 500*unit::ns;
        fAvgFFTHist = NULL;
        fAvgPowerHist = NULL;
        fMeanPowerHist = NULL;
        fAvgPhaseHist = NULL;
        fAvgNoiseHist = NULL;
        fAvgUNoiseHist = NULL;
        fAvgVNoiseHist = NULL;
        fAvgXNoiseHist = NULL;
        fMaxNoiseHist = NULL;
        fAsicNoiseHist = NULL;
        fAsicMaxHist = NULL;
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
        return false;
    }

    void Finalize(CP::TRootOutput * const file) {
        std::cout << "Finalize the run " << fAvgEntries << std::endl;

        for (std::size_t i=0; i<fAvgFFT.size()/2; ++i) {
            fAvgFFT[i] /= fAvgEntries;
            fAvgFFTHist->SetBinContent(i,std::abs(fAvgFFT[i]));
        }
        
        TSpectrum *spectrum = new TSpectrum(300);

        int bins = fAvgPowerHist->GetNbinsX();
        TH1* peakHist
            = new TH1F("peakHist",
                       "Peaks of Discreet Frequency Noise",
                       bins,
                       fAvgPowerHist->GetXaxis()->GetBinLowEdge(1),
                       fAvgPowerHist->GetXaxis()->GetBinUpEdge(bins+1));

        for (int i=0; i<bins; ++i) {
            peakHist->SetBinContent(i+1,fAvgPowerHist->GetBinContent(i+1));
        }
                                  
        TH1* background = spectrum->Background(peakHist,100);
        background->SetName("noiseContinuum");
        background->SetTitle("Stocastic Noise Spectrum");

        for (int i=0; i<bins; ++i) {
            double r = peakHist->GetBinContent(i+1)
                - background->GetBinContent(i+1);
            peakHist->SetBinContent(i+1,r);
        }


       fDeltaFreq = 0.5/fSampleTime/bins/unit::hertz;

       std::ofstream peaksOutput("peaks.txt");
       peaksOutput << "# Index"
                   << " Frequency(Hz)"
                   << " Peak(arb)"
                   << " RMS(arb)"
                   << " Gamma(Hz)"
                   << " Bkg(arb)"
                   << std::endl;
        int peak = 0;
        for (int i=2; i<bins; ++i) {
            double top = peakHist->GetBinContent(i);
            double low = peakHist->GetBinContent(i-1);
            double high = peakHist->GetBinContent(i+1);
            if (top < 10.0) continue;
            if (low >= top) continue;
            if (high >= top) continue;
            int iLow = i-1;
            while (iLow>1) {
                if (peakHist->GetBinContent(iLow-1)
                    >= peakHist->GetBinContent(iLow)) break;
                --iLow;
            }
            int iHigh = i+1;
            while (iHigh < peakHist->GetNbinsX()) {
                if (peakHist->GetBinContent(iHigh+1)
                    >= peakHist->GetBinContent(iHigh)) break;
                ++iHigh;
            }
            low = peakHist->GetBinContent(iLow);
            high = peakHist->GetBinContent(iHigh);
            double delta = 0.5*(top - std::min(low,high));
            if (delta < 20.0) continue;
            int iHalfLow = i;
            while (iHalfLow>1) {
                if (peakHist->GetBinContent(iHalfLow-1)
                    <= top-delta) break;
                --iHalfLow;
            }
            double r1 = top-peakHist->GetBinContent(iHalfLow);
            double r2 = top-peakHist->GetBinContent(iHalfLow-1);
            if (r2 <= r1) continue;
            double d = (delta-r1)/(r2-r1);
            double halfLow = (iHalfLow-d);
            int iHalfHigh = i;
            while (iHalfHigh < peakHist->GetNbinsX()-3) {
                if (peakHist->GetBinContent(iHalfHigh+1)
                    <= top-delta) break;
                ++iHalfHigh;
            }
            r1 = top-peakHist->GetBinContent(iHalfHigh);
            r2 = top-peakHist->GetBinContent(iHalfHigh+1);
            if (r2 <= r1) continue;
            d = (delta-r1)/(r2-r1);
            double halfHigh = (iHalfHigh)+d;
            ++peak;
            peaksOutput << peak << " " << fDeltaFreq*i
                        << " " << top
                        << " " << fAvgPowerHist->GetBinError(i)
                        << " " << fDeltaFreq*std::max(i-halfLow,halfHigh-i)
                        << " " << background->GetBinContent(i)
                        << std::endl;
        }

        std::ofstream noiseOutput("noise.txt");
        noiseOutput << "# Plane"
                    << " Wire"
                    << " Noise(ADC)"
                    << " RMS(ADC)"
                    << std::endl;
        for (int i=1; i<fAvgUNoiseHist->GetNbinsX(); ++i) {
            double avg = fAvgUNoiseHist->GetBinContent(i);
            double err = fAvgUNoiseHist->GetBinError(i);
            noiseOutput << CP::GeomId::Captain::kUPlane
                        << " " << i-1
                        << " " << avg
                        << " " << err
                        << std::endl;
        }
        for (int i=1; i<fAvgVNoiseHist->GetNbinsX(); ++i) {
            double avg = fAvgVNoiseHist->GetBinContent(i);
            double err = fAvgVNoiseHist->GetBinError(i);
            std::cout << CP::GeomId::Captain::kVPlane
                      << " " << i-1
                      << " " << avg
                      << " " << err
                      << std::endl;
        }
        for (int i=1; i<fAvgXNoiseHist->GetNbinsX(); ++i) {
            double avg = fAvgXNoiseHist->GetBinContent(i);
            double err = fAvgXNoiseHist->GetBinError(i);
            std::cout << CP::GeomId::Captain::kXPlane
                      << " " << i-1
                      << " " << avg
                      << " " << err
                      << std::endl;
        }
        
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
        
        std::ostringstream runName;
        runName << "Run " << event.GetContext().GetRun();
        
        int maxASIC = 6*12*16;

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
            
            // Check that the FFT is properly created.
            int nSize = 2*(1+pulse->GetSampleCount()/2);
            fNyquistFreq = (0.5/fSampleTime)/unit::hertz;
            fDeltaFreq = 2.0*fNyquistFreq/nSize;
            if (!fFFT || nSize != fSampleCount) {
                if (fFFT) delete fFFT;
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
             
            if (!fAvgPowerHist) {
                fAvgPowerHist
                    = new TProfile("averagePower",
                                   ("Average Power for " +
                                    runName.str()).c_str(),
                                   nSize/2, 0.0, fNyquistFreq, "S");
                fMeanPowerHist
                    = new TProfile("meanPower",
                                   ("Mean Power for " +
                                    runName.str()).c_str(),
                                   nSize/2, 0.0, fNyquistFreq);
            }

            if (!fAvgPhaseHist) {
                fAvgPhaseHist
                    = new TProfile("averagePhase",
                                   ("Average Phase for " 
                                    + runName.str()).c_str(),
                                   nSize/2, 0.0, fNyquistFreq);
            }

            double maxAmp = 0.0;
            double power = 0.0;
            for (std::size_t i = 0; i<(std::size_t) nSize; ++i) {
                double rl, im;
                fFFT->GetPointComplex(i,rl,im);
                std::complex<double> fft(rl,im);
                fAvgPowerHist->Fill(fDeltaFreq*(i+0.5),std::abs(fft));
                fMeanPowerHist->Fill(fDeltaFreq*(i+0.5),std::abs(fft));
                fAvgPhaseHist->Fill(fDeltaFreq*(i+0.5),std::arg(fft),
                                    std::abs(fft));
                // Apply the missing normalization factor for the FFT.
                fft /= std::sqrt(1.0*nSize);
                fAvgFFT[i] += fft;
                power += std::norm(fft);
                maxAmp = std::max(maxAmp,std::abs(fft));
            }
            
            if (!fAvgNoiseHist) {
                fAvgNoiseHist
                    = new TProfile("averageNoise",
                                   ("Average RMS for " 
                                    + runName.str()).c_str(),
                                   3*333, 0.0, 3*333);
                fAvgUNoiseHist
                    = new TProfile("averageUNoise",
                                   ("Average U RMS for " 
                                    + runName.str()).c_str(),
                                   333, 0.0, 333, "S");
                fAvgVNoiseHist
                    = new TProfile("averageVNoise",
                                   ("Average V RMS for " 
                                    + runName.str()).c_str(),
                                   333, 0.0, 333, "S");
                fAvgXNoiseHist
                    = new TProfile("averageXNoise",
                                   ("Average X RMS for " 
                                    + runName.str()).c_str(),
                                   333, 0.0, 333, "S");
            }

            if (!fMaxNoiseHist) {
                fMaxNoiseHist
                    = new TProfile("maxNoise",
                                   ("Max Amplitude for " 
                                    + runName.str()).c_str(),
                                   3*333, 0.0, 3*333);
            }

            if (!fAsicNoiseHist) {
                fAsicNoiseHist
                    = new TProfile("asicNoise",
                                   ("Average ASIC RMS for " 
                                    + runName.str()).c_str(),
                                   maxASIC, 0.0, maxASIC);
            }

            if (!fAsicMaxHist) {
                fAsicMaxHist
                    = new TProfile("asicMax",
                                   ("Max ASIC for " 
                                    + runName.str()).c_str(),
                                   maxASIC, 0.0, maxASIC);
            }

            double noise = std::sqrt(power/nSize);
            double wire =  333*CP::GeomId::Captain::GetWirePlane(pulseGeom)
                + CP::GeomId::Captain::GetWireNumber(pulseGeom);
            fAvgNoiseHist->Fill(wire+0.5,noise);
            if (CP::GeomId::Captain::IsUWire(pulseGeom)) {
                fAvgUNoiseHist->Fill(
                    CP::GeomId::Captain::GetWireNumber(pulseGeom)+0.5,noise);
            }
            if (CP::GeomId::Captain::IsVWire(pulseGeom)) {
                fAvgVNoiseHist->Fill(
                    CP::GeomId::Captain::GetWireNumber(pulseGeom)+0.5,noise);
            }
            if (CP::GeomId::Captain::IsXWire(pulseGeom)) {
                fAvgXNoiseHist->Fill(
                    CP::GeomId::Captain::GetWireNumber(pulseGeom)+0.5,noise);
            }
            fMaxNoiseHist->Fill(wire+0.5,maxAmp);
            
            int asic = chanInfo.GetMotherboard(pulse->GetChannelId())-1;
            asic = 12*asic+chanInfo.GetASIC(pulse->GetChannelId())-1;
            asic = 16*asic+chanInfo.GetASICChannel(pulse->GetChannelId())-1;
            fAsicNoiseHist->Fill(asic+0.5,noise);
            fAsicMaxHist->Fill(asic+0.5,maxAmp);
            
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
    
    /// The average FFT
    std::vector< std::complex<double> > fAvgFFT;

    /// The number of entries added to the average.
    int fAvgEntries;
    
    TH1F* fAvgFFTHist;
    
    /// The average power at each frequency
    TProfile* fAvgPowerHist;

    /// The average power at each frequency
    TProfile* fMeanPowerHist;

    /// The average phase at each frequency
    TProfile* fAvgPhaseHist;

    /// The average RMS for each wire.
    TProfile* fAvgNoiseHist;
    TProfile* fAvgUNoiseHist;
    TProfile* fAvgVNoiseHist;
    TProfile* fAvgXNoiseHist;

    /// The average RMS for each wire.
    TProfile* fMaxNoiseHist;

    TProfile* fAsicNoiseHist;
    TProfile* fAsicMaxHist;
    
};

int main(int argc, char **argv) {
    TCalibHistsLoop userCode;
    CP::eventLoop(argc,argv,userCode,1);
}
