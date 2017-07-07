#include "GaussianNoise.hxx"
#include "FindPedestal.hxx"

#include <eventLoop.hxx>
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
#include <sstream>
#include <algorithm>

namespace CP {class TRootOutput;};

/// Make diagnostic histograms for a file.
class TCalibHistsLoop: public CP::TEventLoopFunction {
public:
    TCalibHistsLoop() {
        fFFT = NULL;
        fNoiseHist = NULL;
        fGaussHist = NULL;
        fPowerHist = NULL;
        fExtremaHist = NULL;
	fCountPeaksHist = NULL;
        fMaxFFTHist = NULL;
        fChanFFTHist = NULL;
        fWireFFTHist = NULL;
        fASICFFTHist = NULL;
        fCorrChanHist = NULL;
        fCorrASICHist = NULL;
        fCorrWireHist = NULL;
        fPeakChanHist = NULL;
        fPeakASICHist = NULL;
        fPeakWireHist = NULL;
        fSampleCount = 0;
        fPrintFile = "";
        fSampleTime = 500*unit::ns;
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
        std::cout << "Finalize the run " << std::endl;
        if (fPrintFile.size() < 1) return;

         TCanvas c("clusterCalibFFT");

        gPad->Print((fPrintFile + ".pdf[").c_str());

        fNoiseHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-noiseHist.png").c_str());

        fPowerHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-gaussHist.png").c_str());

        fExtremaHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-extremaHist.png").c_str());

        gPad->SetLogy(true);
        fCountPeaksHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-countPeaksHist.png").c_str());
        gPad->SetLogy(false);

        fPeakChanHist->SetMinimum(fNoiseHist->GetMean());
        gPad->SetLogy(true);
        fPeakChanHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-peakChanHist.png").c_str());
        gPad->SetLogy(false);

        fPeakWireHist->SetMinimum(fNoiseHist->GetMean());
        gPad->SetLogy(true);
        fPeakWireHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-peakWireHist.png").c_str());
        gPad->SetLogy(false);

        fPeakASICHist->SetMinimum(fNoiseHist->GetMean());
        gPad->SetLogy(true);
        fPeakASICHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-peakASICHist.png").c_str());
        gPad->SetLogy(false);

        fMaxCorrChanHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-maxCorrChanHist.png").c_str());

        fMaxCorrWireHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-maxCorrWireHist.png").c_str());

        fMaxCorrASICHist->Draw();
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-maxCorrASICHist.png").c_str());

        fMaxFFTHist->Draw();
        gPad->SetLogy(true);
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-maxFFTHist.png").c_str());
        gPad->SetLogy(false);

        fChanFFTHist->Draw("colz");
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-chanFFTHist.png").c_str());

        fWireFFTHist->Draw("colz");
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-wireFFTHist.png").c_str());

        fASICFFTHist->Draw("colz");
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-asicFFTHist.png").c_str());

        fCorrChanHist->Draw("colz");
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-corrChanHist.png").c_str());

        fCorrWireHist->Draw("colz");
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-corrWireHist.png").c_str());

        fCorrASICHist->Draw("colz");
        gPad->Print((fPrintFile+".pdf").c_str());
        gPad->Print((fPrintFile+"-corrASICHist.png").c_str());

        gPad->Print(((fPrintFile+".pdf") + "]").c_str());
    }

    bool operator () (CP::TEvent& event) {
        int maxASIC = 6*12*16;
        int maxWire = 996;
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
        
        if (!fNoiseHist) {
            fNoiseHist = new TH1F("noiseHist",
                                  (titlePrefix.str() +
                                   "RMS of ADC values").c_str(),
                                  100.0, 0, 50.0);
            fNoiseHist->SetXTitle("ADC");
        }

        if (!fGaussHist) {
            fGaussHist = new TH1F("gaussHist",
                                  (titlePrefix.str() +
                                   "Estimate of ADC Gaussian noise"
                                   + " from differences").c_str(),
                                  100.0, 0, 50.0);
            fGaussHist->SetXTitle("ADC");
        }

        if (!fPowerHist) {
            fPowerHist = new TH1F("powerHist",
                                  (titlePrefix.str() +
                                   "Estimate of ADC Gaussian noise"
                                   + " from power spectrum").c_str(),
                                  100.0, 0, 50.0);
            fPowerHist->SetXTitle("ADC");
        }

        if (!fExtremaHist) {
            fExtremaHist = new TH1F("extremaHist",
                                    (titlePrefix.str() +
                                     "RMS of ADC extrema values").c_str(),
                                    100, 0.0, 200.0);
            fExtremaHist->SetXTitle("ADC");
        }
        
        if (!fCountPeaksHist) {
            fCountPeaksHist = new TH1F("countPeaksHist",
                                    (titlePrefix.str() +
                                     "FFT samples great than 0.5").c_str(),
                                    maxWire, 1, maxWire+1);
            fCountPeaksHist->SetXTitle("Wire Number");
            fCountPeaksHist->SetYTitle("Number of Events");
	    fCountPeaksHist->SetStats(false);
        }

        if (!fCorrChanHist) {
            fCorrChanHist = new TH2F("corrChanHist",
                                     (titlePrefix.str() +
                                      "Correlations for all channels").c_str(),
                                     drift->size(), 0.0, drift->size(),
                                     drift->size(), 0.0, drift->size());
            fCorrChanHist->GetXaxis()->SetNdivisions(drift->size()/64,false);
            fCorrChanHist->GetYaxis()->SetNdivisions(drift->size()/64,false);
            fCorrChanHist->SetStats(false);
            fCorrChanHist->SetMinimum(-1);
            fCorrChanHist->SetMaximum(1);
            fCorrChanHist->SetXTitle("Channel");
            fCorrChanHist->SetYTitle("Channel");

            fMaxCorrChanHist
                = new TH1F("maxCorrChanHist",
                           (titlePrefix.str() +
                            "Maximum correlation for all channels").c_str(),
                           drift->size(), 0.0, drift->size());
            fMaxCorrChanHist->SetStats(false);
            fMaxCorrChanHist->GetXaxis()->SetNdivisions(drift->size()/64,false);
            fMaxCorrChanHist->SetXTitle("Channel");
            fMaxCorrChanHist->SetYTitle("Max Correlation Coefficient");
        }

        if (!fCorrASICHist) {
            fCorrASICHist = new TH2F("corrASICHist",
                                     (titlePrefix.str() +
                                      "Correlations for all ASICs").c_str(),
                                     maxASIC, 1.0, maxASIC+1,
                                     maxASIC, 1.0, maxASIC+1);
            fCorrASICHist->GetXaxis()->SetNdivisions(6 + 12*100,false);
            fCorrASICHist->GetYaxis()->SetNdivisions(6 + 12*100,false);
            fCorrASICHist->SetStats(false);
            fCorrASICHist->SetMinimum(-1);
            fCorrASICHist->SetMaximum(1);
            fCorrASICHist->SetXTitle("ASIC");
            fCorrASICHist->SetYTitle("ASIC");

            fMaxCorrASICHist
                = new TH1F("maxCorrASICHist",
                           (titlePrefix.str() +
                            "Max correlation for all ASICs").c_str(),
                           maxASIC, 1.0, maxASIC+1);
            fMaxCorrASICHist->SetStats(false);
            fMaxCorrASICHist->GetXaxis()->SetNdivisions(6 + 12*100,false);
            fMaxCorrASICHist->SetXTitle("ASIC");
            fMaxCorrASICHist->SetYTitle("Max Correlation Coefficient");
        }

        if (!fCorrWireHist) {
            fCorrWireHist
                = new TH2F("corrWireHist",
                           (titlePrefix.str() 
                            + "Correlations for all wires").c_str(),
                           maxWire, 1, maxWire+1,
                           maxWire, 1, maxWire+1);
            fCorrWireHist->SetStats(false);
            fCorrWireHist->SetMinimum(-1);
            fCorrWireHist->SetMaximum(1);
            fCorrWireHist->SetXTitle("Wire");
            fCorrWireHist->SetYTitle("Wire");

            fMaxCorrWireHist
                = new TH1F("maxCorrWireHist",
                           (titlePrefix.str() +
                            + "Max correlation for all wires").c_str(),
                           maxWire, 1.0, maxWire+1);
            fMaxCorrWireHist->SetXTitle("Wire");
            fMaxCorrWireHist->SetYTitle("Max Correlation Coefficient");
            fMaxCorrWireHist->SetStats(false);
        }

        if (!fPeakWireHist) {
            fPeakWireHist = new TH1F("peakWireHist",
                                     (titlePrefix.str()
                                      + "RMS peak height (top quintile)"
                                      + " on each wire").c_str(),
                                     maxWire, 1, maxWire+1);
            fPeakWireHist->SetXTitle("Wire Number");
            fPeakWireHist->SetYTitle("|ADC|");
            fPeakWireHist->SetStats(false);
        }
        if (!fPeakASICHist) {
            fPeakASICHist
                = new TH1F("peakASICHist",
                           (titlePrefix.str()
                            + "RMS peak height (top quintile)"
                            + " on each ASIC").c_str(),
                           maxASIC, 1, maxASIC+1);
            fPeakASICHist->SetXTitle("ASIC Number");
            fPeakASICHist->SetYTitle("|ADC|");
            fPeakASICHist->GetXaxis()->SetNdivisions(6 + 12*100,false);
            fPeakASICHist->SetStats(false);
        }
        if (!fPeakChanHist) {
            fPeakChanHist = new TH1F("peakChanHist",
                                     (titlePrefix.str()
                                      + "RMS peak height (top quintile)"
                                      + " on each channel").c_str(),
                                     drift->size(), 0, drift->size());
            fPeakChanHist->SetXTitle("Channel Number");
            fPeakChanHist->SetYTitle("|ADC|");
            fPeakChanHist->SetStats(false);
            fPeakChanHist->GetXaxis()->SetNdivisions(drift->size()/64,false);
        }
        
        std::vector<double> pedestals(drift->size());
        std::vector<double> sigmas(drift->size());
        std::vector<double> differences;
        
        std::vector<float> powerRange;

        for (std::size_t d = 0; d < drift->size(); ++d) {
            const CP::TPulseDigit* pulse 
                = dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
            if (!pulse) {
                CaptError("Non-pulse in drift digits");
                continue;
            }
            // Check that the FFT is properly created.
            int nSize = 2*(1+pulse->GetSampleCount()/2);
            double nyquistFreq = (0.5/fSampleTime)/unit::hertz;
            double deltaFreq = 2.0*nyquistFreq/nSize;

            int wire = chanInfo.GetWireNumber(pulse->GetChannelId());

            if (differences.size() != pulse->GetSampleCount()) {
                differences.resize(pulse->GetSampleCount());
            }
            
            int asic = chanInfo.GetMotherboard(pulse->GetChannelId())-1;
            asic = 12*asic+chanInfo.GetASIC(pulse->GetChannelId())-1;
            asic = 16*asic+chanInfo.GetASICChannel(pulse->GetChannelId())-1;
            
            if (!fFFT || nSize != fSampleCount) {
                if (!fFFT) delete fFFT;
                CaptLog("Nyquist frequency: " << nyquistFreq);
                CaptLog("Frequency bin size: " << deltaFreq);
                fFFT = TVirtualFFT::FFT(1, &nSize, "R2C M K");
                fSampleCount = nSize;
                int overSample = 20;
                
                fMaxFFTHist = new TH1F("maxFFTHist",
                                        (titlePrefix.str() 
                                         + "Max FFT power vs Frequency"
                                         + " in any channel").c_str(),
                                       (nSize/2-1)/overSample,
                                       deltaFreq, nyquistFreq);
                fMaxFFTHist->SetXTitle("Frequency (Hz)");
                fMaxFFTHist->SetYTitle(
                    "FFT Power (#frac{ADC^{2}}{#Delta#it{f}})");
                fMaxFFTHist->SetStats(false);
                
                fChanFFTHist = new TH2F("chanFFTHist",
                                        (titlePrefix.str() 
                                         + "Max FFT power in window"
                                         + " for channels").c_str(),
                                        drift->size(), 0.0, drift->size(),
                                        (nSize/2-1)/overSample,
                                        deltaFreq, nyquistFreq);
                fChanFFTHist->SetXTitle("Channel");
                fChanFFTHist->SetYTitle("Frequency (Hz)");
                fChanFFTHist->GetXaxis()->SetNdivisions(drift->size()/64,false);
                fChanFFTHist->SetStats(false);
                fWireFFTHist
                    = new TH2F("wireFFTHist",
                               (titlePrefix.str() 
                                + "Max FFT power in window"
                                + " for wires").c_str(),
                               maxWire, 1.0, maxWire+1,
                               (nSize/2-1)/overSample,
                               deltaFreq, nyquistFreq);
                fWireFFTHist->SetXTitle("Wire");
                fWireFFTHist->SetYTitle("Frequency (Hz)");
                fWireFFTHist->SetStats(false);
                fASICFFTHist = new TH2F("asicFFTHist",
                                        (titlePrefix.str()
                                         + "Max FFT power in window"
                                         + " for asics").c_str(),
                                        maxASIC, 0.0, maxASIC,
                                        (nSize/2-1)/overSample,
                                        deltaFreq, nyquistFreq);
                fASICFFTHist->SetXTitle("ASIC");
                fASICFFTHist->SetYTitle("Frequency (Hz)");
                fASICFFTHist->GetXaxis()->SetNdivisions(6 + 12*100,false);
                fASICFFTHist->SetStats(false);
                powerRange.reserve(drift->size()*nSize/2);
            }

            // Find the average pedestal and subtract it.
            double pedestal = CP::FindPedestal(pulse->begin(), pulse->end());
            pedestals[d] = pedestal;

            double gaussNoise = CP::GaussianNoise(pulse->begin(), pulse->end());
            fGaussHist->Fill(std::max(gaussNoise,0.0));

            // Find the minima and maxima and use them to tabulate the peak
            // heights.
            std::vector<double> localExtrema;
            for (std::size_t i=1; i<pulse->GetSampleCount(); ++i) {
                if (pulse->GetSample(i) < pulse->GetSample(i-1)
                    && pulse->GetSample(i) < pulse->GetSample(i+1)
                    && pulse->GetSample(i) < pedestal) {
                    localExtrema.push_back(pedestal - pulse->GetSample(i));
                }
                else if (pulse->GetSample(i) > pulse->GetSample(i-1)
                         && pulse->GetSample(i) > pulse->GetSample(i+1)
                    && pulse->GetSample(i) > pedestal) {
                    localExtrema.push_back(pulse->GetSample(i)-pedestal);
                }
            }
            std::sort(localExtrema.begin(), localExtrema.end());
            double rmsExtrema = 0.0;
            double weightExtrema = 0.0;
            for (int i = 0.8*localExtrema.size();
                 i<0.99*localExtrema.size(); ++i) {
                double p = localExtrema[i];
                rmsExtrema += p*p;
                weightExtrema += 1.0;
            }
            if (weightExtrema>0.5) {
                rmsExtrema /= weightExtrema;
            }
            rmsExtrema = sqrt(rmsExtrema);
            fExtremaHist->Fill(rmsExtrema);
            if (wire>=0) {
                fPeakWireHist->SetBinContent(wire+0.1,rmsExtrema);
            }
            fPeakChanHist->SetBinContent(asic+0.1,rmsExtrema);
            fPeakASICHist->SetBinContent(d+0.1,rmsExtrema);
            
            double sigma = 0.0;
            for (std::size_t i=0; i<(std::size_t) nSize; ++i) {
                double p = 0.0;
                if (i<pulse->GetSampleCount()) {
                    p = pulse->GetSample(i)-pedestal;
                    sigma += p*p;
                }
                fFFT->SetPoint(i,p);
            }
            sigma /= pulse->GetSampleCount();
            sigmas[d] = std::sqrt(sigma);
            
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
            if (0 <= wire) histBaseTitle << "wire " << wire << " ";
            histBaseTitle << chanName;
            
            TH1F* fftHist 
                = new TH1F((histBaseName.str()+"-fft").c_str(),
                           (titlePrefix.str() 
                            + "FFT for " 
                            + histBaseTitle.str()).c_str(),
                           nSize/2-1, deltaFreq, nyquistFreq);
            fftHist->SetYTitle("FFT Power (#frac{ADC^{2}}{#Delta#it{f}})");
            fftHist->SetXTitle("Frequency (Hz)");

	    int peakCounter = 0;
            std::vector<double> power(nSize/2);
            double amp = 0.0;
            for (std::size_t i = 1; i<(std::size_t) nSize/2; ++i) {
                double rl, im;
                fFFT->GetPointComplex(i,rl,im);
                std::complex<double> c(rl,im);
                double p = 2.0*std::abs(c*c)/nSize/nSize;
                amp += p;
                power[i] = p;
                fftHist->SetBinContent(i, p);
		
		// Count the number of peaks in the power spectrum.
		if (p > 0.5) peakCounter++;

                // Fill the powerRange vector, which will be sorted, and then
                // its values used for the scale of the 2D histos
		if (p>0.01/nSize) powerRange.push_back(p);
                else p = 0.01/nSize;

		// Fill 2D histograms.  Be careful with the "Y" binning
                int b = 1+2*fChanFFTHist->GetNbinsY()*i/nSize;
		double m = fChanFFTHist->GetBinContent(d+1,b);
                if (p > m || m == 0.0) {
                    fChanFFTHist->SetBinContent(d+1,b,p);
                }
                b = 1+2*fWireFFTHist->GetNbinsY()*i/nSize;
                if (0 <= wire) {
                    m = fWireFFTHist->GetBinContent(wire,b);
                    if (p > m || m == 0.0) {
                        fWireFFTHist->SetBinContent(wire,b,p);
                    }
                }
                b = 1+2*fASICFFTHist->GetNbinsY()*i/nSize;
                m = fASICFFTHist->GetBinContent(asic,b);
                if (p > m || m == 0.0) {
                    fASICFFTHist->SetBinContent(asic,b,p);
                }
            }

	    if (0 <= wire) {
                fCountPeaksHist->SetBinContent(wire+0.1,peakCounter);
	    }

            fNoiseHist->Fill(std::sqrt(amp));
            std::sort(power.begin(), power.end());
            int index = 0.68*power.size();
            amp = power[index]*nSize;
            fPowerHist->Fill(std::sqrt(amp));
        }

        std::sort(powerRange.begin(), powerRange.end());
        double minPower = powerRange[0.70*powerRange.size()];
        double maxPower = powerRange[0.98*powerRange.size()];
        if (minPower < 1E-6) minPower = 1E-6;
        if (maxPower < 1E-4) maxPower = 1E-4;
        
        fChanFFTHist->SetMinimum(minPower);
        fChanFFTHist->SetMaximum(maxPower);

        fWireFFTHist->SetMinimum(minPower);
        fWireFFTHist->SetMaximum(maxPower);

        fASICFFTHist->SetMinimum(minPower);
        fASICFFTHist->SetMaximum(maxPower);

        // Fill the "maximum" power vs frequency bin.  This is calculated from
        // the fChanFFTHist which has the maximum power in each frequency bin
        // for each channel.  Notice that this isn't the absolute maximum
        // since it throws out the top few outlying channels for each
        // frequency bin.
        for (int f1=1; f1<=fChanFFTHist->GetNbinsY(); ++f1) {
            // Loop over the channels.
            std::vector<double> cPower(fChanFFTHist->GetNbinsX());
            cPower.clear();
            for (int c1=1; c1<=fChanFFTHist->GetNbinsX(); ++c1) {
                cPower.push_back(fChanFFTHist->GetBinContent(c1,f1));
            }
            if (cPower.empty()) continue;
            std::sort(cPower.begin(), cPower.end());
            fMaxFFTHist->SetBinContent(f1,cPower[0.95*cPower.size()]);
        }
        
        // Calculate the correlations.
        for (std::size_t d1 = 0; d1 < drift->size(); ++d1) {
            const CP::TPulseDigit* pulse1
                = dynamic_cast<const CP::TPulseDigit*>((*drift)[d1]);
            if (!pulse1) {
                CaptError("Non-pulse in drift digits");
                continue;
            }
            
            CP::TChannelId id1 = pulse1->GetChannelId();

            int asic1 = chanInfo.GetMotherboard(id1)-1;
            asic1 = 12*asic1 + chanInfo.GetASIC(id1)-1;
            asic1 = 16*asic1+ chanInfo.GetASICChannel(id1)-1;

            int wire1 = chanInfo.GetWireNumber(id1);
            double pedestal1 = pedestals[d1];
            double sigma1 = sigmas[d1];

            double maxCorr = 0.0;
            for (std::size_t d2 = 0; d2 < drift->size(); ++d2) {
                if (d1 == d2) continue;
                const CP::TPulseDigit* pulse2
                    = dynamic_cast<const CP::TPulseDigit*>((*drift)[d2]);
                if (!pulse2) {
                    CaptError("Non-pulse in drift digits");
                    continue;
                }
        
                CP::TChannelId id2 = pulse2->GetChannelId();
                
                int asic2 = chanInfo.GetMotherboard(id2)-1;
                asic2 = 12*asic2 + chanInfo.GetASIC(id2)-1;
                asic2 = 16*asic2+ chanInfo.GetASICChannel(id2)-1;

                int wire2 = chanInfo.GetWireNumber(id2);
                
                double pedestal2 = pedestals[d2];
                double sigma2 = sigmas[d2];

                double corr12 = 0.0;
                for (std::size_t i=0; i<pulse2->GetSampleCount(); ++i) {
                    double p1 = pulse1->GetSample(i)-pedestal1;
                    double p2 = pulse2->GetSample(i)-pedestal2;
                    corr12 += p1*p2;
                }
                corr12 /= pulse2->GetSampleCount();
                corr12 /= (sigma1*sigma2);

                fCorrChanHist->SetBinContent(d1,d2,corr12);
                fCorrASICHist->SetBinContent(asic1,asic2,corr12);
                fCorrWireHist->SetBinContent(wire1,wire2,corr12);
                if (std::abs(maxCorr) < std::abs(corr12)) maxCorr = corr12;
            }
            fMaxCorrChanHist->SetBinContent(d1,std::abs(maxCorr));
            fMaxCorrASICHist->SetBinContent(asic1,std::abs(maxCorr));
            fMaxCorrWireHist->SetBinContent(wire1,std::abs(maxCorr));
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

    /// A histogram of the noise RMS on each wire.
    TH1F* fNoiseHist;

    /// A histogram of the noise on each wire assuming it's Gaussian
    TH1F* fGaussHist;

    /// A histogram of the noise on each wire assuming it's Gaussian
    TH1F* fPowerHist;

    /// A histogram of the extrema value fluctuations.
    TH1F* fExtremaHist;

    /// A histogram of the number of peaks above a certain power spectra value
    /// Hardcoded as above 0.5 for now
    TH1F* fCountPeaksHist;

    /// The maximum power at each frequency.  This is the max for all channels
    /// and can be used to pick out the important frequencies.
    TH1F* fMaxFFTHist;
    
    /// The FFT for all wires in an event.
    TH2F* fChanFFTHist;
    TH2F* fWireFFTHist;
    TH2F* fASICFFTHist;

    /// The correlations for all wires in an event.
    TH2F* fCorrChanHist;

    /// The correlations for all wires in an event.
    TH2F* fCorrASICHist;

    /// The correlations for all wires in an event.
    TH2F* fCorrWireHist;

    /// The noise vs channe.
    TH1F* fPeakWireHist;
    TH1F* fPeakChanHist;
    TH1F* fPeakASICHist;

    /// The maximum correlation.
    TH1F* fMaxCorrChanHist;
    TH1F* fMaxCorrASICHist;
    TH1F* fMaxCorrWireHist;
    
    /// A file to draw the wire histogram in (usually png).
    std::string fPrintFile;
    
};

int main(int argc, char **argv) {
    TCalibHistsLoop userCode;
    CP::eventLoop(argc,argv,userCode,1);
}
