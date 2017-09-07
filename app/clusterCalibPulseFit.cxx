///////////////////////////////////////////////////////////////////
// Take a file with pulses injected into the calibration capacitor and
// determine the gain and pulse shape.  This produces a file suitable for
// loading into the DB (i.e. a SQL file of calibration constants).
#include "GaussianNoise.hxx"
#include "FindPedestal.hxx"

#include <eventLoop.hxx>
#include <TPulseDigit.hxx>
#include <TChannelId.hxx>
#include <TTPCChannelId.hxx>
#include <HEPUnits.hxx>
#include <CaptGeomId.hxx>

#include <TChannelInfo.hxx>
#include <TTPC_Channel_Calib_Table.hxx>

#include <TF1.h>
#include <TPad.h>
#include <TString.h>
#include <TProfile.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <ctime>

namespace CP {class TRootOutput;};

double gain(double p) {
  return 400.0*std::exp(p);
}

double peakingTime(double p) {
  double v = 2.0+0.05*std::atan(p)/1.6;
  return v;
}

double fallShape(double p) {
  return std::exp(std::atan(p));
}

double riseShape(double p) {
  return std::exp(0.3*std::atan(p));
}

double PulseShape(double* x, double* par) {
  double ped = par[0];
  double t = x[0] - 3200.0 - 2.0*std::atan(par[1])/1.6;
  double arg = t/peakingTime(par[5]);
  if (arg < 0.0) arg =0.0;
  else if (arg < 1.0) arg = std::pow(arg,riseShape(par[3]));
  else arg = std::pow(arg,fallShape(par[4]));
  double v = (arg<40)? arg*std::exp(-arg)*std::exp(1.0): 0.0;
  v *= gain(par[2]);
  return v + ped;
}

int getStatus(std::vector<int> status_vector) {
  int status = 0;
  int statusBitmask = 0;
  for (unsigned int i=0; i<status_vector.size(); i++) {
    int s = status_vector[i];
    status += (s!=0);
    statusBitmask |= s;
  }
  if (status/status_vector.size() > 0.5)
    return statusBitmask;
  else
    return 0;
}

float getAvg(std::vector<double> parameter_vector) {
  double parAvg = 0.0;
  for (unsigned int i=0; i<parameter_vector.size();i++) {
    double p = parameter_vector[i];
    parAvg += p;
  }
  return parAvg/parameter_vector.size();
}

float getRMS(std::vector<double> parameter_vector, double mean) {
  double parRMS = 0.0;
  for (unsigned int i=0; i<parameter_vector.size();i++) {
    double p = parameter_vector[i];
    parRMS += ((p - mean)*(p - mean));
  }
  return sqrt(parRMS/parameter_vector.size());
}

/// Make diagnostic histograms for a file.
class TCalibPulseFitLoop: public CP::TEventLoopFunction {
public:
  TCalibPulseFitLoop() {
    // Chuck "fixed" the input voltage.  I think it was 0.8 V, but this
    // needs to be checked.
    fInjectedPulse = 0.8*unit::volt;
    // The ASICs have a calibrated capacitor.  I don't recall the exact
    // capacitance.  This needs to be checked.
    fCapacitance = 17.5*unit::femto*unit::farad;
    // Turn the pulse voltage and capacitance into the injected charge.
    fInjectedCharge = fCapacitance*fInjectedPulse;
    // This is fixed, and could be made a constant. 
    fDigitizerSlope = 2.5/unit::millivolt;
  }

  virtual ~TCalibPulseFitLoop() {};

  void Usage(void) {
  }

  void Initialize() {
  }
    
  virtual bool SetOption(std::string option, std::string value="") {
    return false;
  }

  void Finalize(CP::TRootOutput * const file) {
    // Open the output file (with the date in the name), and generate a
    // string with the creation date.
    std::time_t t = std::time(0);
    char buf[256];
    std::strftime(buf,sizeof(buf),"%Y-%m-%d %H:%M:%S",std::gmtime(&t));
    std::string creationDate(buf);
    std::ostringstream outputName;
    outputName << "TPC_CHANNEL_CALIB_TABLE"
	       << "_" << fTimeStamp.substr(0,fTimeStamp.find(" "))
	       << "_" << fRunNumber
	       << "_" << creationDate.substr(0,creationDate.find(" "))
	       << ".update";
            
    std::ofstream output(outputName.str());

    // Write the header for the calibration table.
    output << "BEGIN_TABLE TPC_CHANNEL_CALIB_TABLE"
	   << " '" << fTimeStamp << "'"
	   << " " << "'2037-09-01 00:00:00'" // end-date
	   << " " << 0                       // Aggregate Number
	   << " '" << creationDate << "'"
	   << " " << 0                      // task ???
	   << " " << "DETECTORMASK=mCAPTAIN"
	   << " " << "SIMMASK=Data"
	   << " " << "EPOCH=0"
	   << std::endl;

    // Fit the calibration constants for every channel and write the
    // constants to the calibration table.
    //int h_i = 0;

    int status = 0;
    double par0 = 0.0;
    double gainValue = 0.0;
    double peakTime = 0.0;
    double riseShp = 0.0;
    double fallShp = 0.0;
    
    for (std::map<CP::TChannelId,std::vector<int>>::iterator
	   h = fChanStatus.begin();
	 h != fChanStatus.end(); ++h) {

      CP::TTPCChannelId chanId(h->first);
      
      status = getStatus((h->second));
      par0 = getAvg(fChanPar0[chanId]);
      gainValue = getAvg(fChanGainValue[chanId]);
      peakTime = getAvg(fChanPeakingTime[chanId]);
      riseShp = getAvg(fChanRiseShape[chanId]);
      fallShp = getAvg(fChanFallShape[chanId]);
      
      output << chanId.GetCrate()
	     << "," << chanId.GetFEM()
	     << "," << chanId.GetChannel()
	     << "," << status
	     << "," << par0
	     << "," << gainValue
	     << ","
	     << peakTime
	     << "," << riseShp
	     << "," << fallShp
	     << std::endl;
      


      // std::cout << chanId.GetCrate()
      // 		<< ", " << chanId.GetFEM()
      // 		<< ", " << chanId.GetChannel()
      // 		<< ", " << status
      // 		<< ", " << par0 << "+/-" << getRMS(fChanPar0[chanId],par0)
      // 		<< ", " << gainValue << "+/-" << getRMS(fChanGainValue[chanId],gainValue)
      // 		<< ", " << peakTime << "+/-" << getRMS(fChanPeakingTime[chanId],peakTime)
      // 		<< ", " << riseShp <<"+/-" << getRMS(fChanRiseShape[chanId],riseShp)
      // 		<< ", " << fallShp <<"+/-" << getRMS(fChanFallShape[chanId],fallShp)
      // 		<< std::endl;      
      
    }



  }

  bool operator () (CP::TEvent& event) {
    CP::THandle<CP::TDigitContainer> drift
      = event.Get<CP::TDigitContainer>("~/digits/drift");

    fTimeStamp = event.GetContext().GetTimeStampString();
    fTimeStamp = fTimeStamp.substr(1,fTimeStamp.size()-2);
    fRunNumber = event.GetContext().GetRun();
	
    if (!drift) {
      CaptError("No drift signals for this event " << event.GetContext());
      return false;
    }

    CP::TChannelInfo& chanInfo = CP::TChannelInfo::Get();
    chanInfo.SetContext(event.GetContext());

    /// For each digit... (but only the first one for debugging).   
    for (std::size_t d = 0; d < drift->size(); ++d) {
      const CP::TPulseDigit* pulse
	= dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
      if (!pulse) {
	CaptError("Non-pulse in drift digits");
	continue;
      }

      CP::TChannelId chanId = pulse->GetChannelId();
      TProfile *chanProfile = new TProfile(chanId.AsString().c_str(),
					   chanId.AsString().c_str(),
					   500, 3190, 3220.0,"i");
            
	    
      for (std::size_t i=3180; i<3230; ++i) {
	double p = pulse->GetSample(i);
	chanProfile->Fill(i+0.5, p);
      }

      TF1* pulseShape = new TF1("pulseShape",PulseShape,3190.,3220.,6);
      
      double pedestal = chanProfile->GetMean(2) - 0.5*chanProfile->GetRMS(2);
      pulseShape->SetParameter(0,pedestal);
      pulseShape->SetParameter(1,0.0);
      pulseShape->SetParameter(2,0.0);
      pulseShape->SetParameter(3,0.0);
      pulseShape->SetParameter(4,0.0);
      pulseShape->SetParameter(5,0.0);
      chanProfile->Fit("pulseShape","q0","",3190.,3220.);
      chanProfile->SetMinimum(chanProfile->GetMean(2)-chanProfile->GetRMS(2));
      int status = 0;
      double gainValue =  gain(pulseShape->GetParameter(2));
      gainValue /= fDigitizerSlope;  // V
      gainValue /= fInjectedCharge; // V/C
      gainValue /= (unit::mV/unit::femto/unit::coulomb);
      if (gainValue < 1.0) {
	gainValue = 0.0;
	status |= CP::TTPC_Channel_Calib_Table::kNoSignal;
      }
      else if (gainValue < 10) {
	status |= CP::TTPC_Channel_Calib_Table::kLowGain;
      }
      else if (gainValue > 20) {
	status |= CP::TTPC_Channel_Calib_Table::kHighGain;
      }
      double peakTime = peakingTime(pulseShape->GetParameter(5));
      peakTime *= 500.0*unit::ns;
      if (peakTime < 700*unit::ns) {
	status |= CP::TTPC_Channel_Calib_Table::kBadPeak;
      }


      float goodnessOfFit = 0.0;

      for (int i=1; i < chanProfile->GetXaxis()->GetNbins(); i++) {
	float xWire = chanProfile->GetXaxis()->GetBinCenter(i);
	if (chanProfile->GetBinContent(i) != 0) {
	  //std::cout<<xWire<< "," << float(chanProfile->GetBinContent(i))<< "," << pulseShape->Eval(xWire) <<std::endl;
	  goodnessOfFit += fabs((float(chanProfile->GetBinContent(i)) - pulseShape->Eval(xWire)));
	}
      }
      
      if (goodnessOfFit > 400) {
	status |= CP::TTPC_Channel_Calib_Table::kBadFit;
      }
      
      CP::TTPCChannelId TPCchanId(chanId);
      
      // std::cout << TPCchanId.GetCrate()
      // << "," << TPCchanId.GetFEM()
      // << "," << TPCchanId.GetChannel()
      // << "," << status
      // << "," << pulseShape->GetParameter(0)
      // << "," << gainValue
      // << ","
      // << peakingTime(pulseShape->GetParameter(5))*500*unit::ns
      // << "," << riseShape(pulseShape->GetParameter(3))
      // << "," << fallShape(pulseShape->GetParameter(4))
      // << std::endl;	    
      

      std::vector<int> p0 = fChanStatus[chanId];
      std::vector<double> p1 = fChanPar0[chanId];
      std::vector<double> p2 = fChanGainValue[chanId];
      std::vector<double> p3 = fChanPeakingTime[chanId];
      std::vector<double> p4 = fChanRiseShape[chanId];
      std::vector<double> p5 = fChanFallShape[chanId];

      p0.push_back(status);
      p1.push_back(pulseShape->GetParameter(0));
      p2.push_back(gainValue);
      p3.push_back(peakingTime(pulseShape->GetParameter(5))*500*unit::ns);
      p4.push_back(riseShape(pulseShape->GetParameter(3)));
      p5.push_back(fallShape(pulseShape->GetParameter(4)));
      
      fChanStatus[chanId] = p0;
      fChanPar0[chanId] = p1;
      fChanGainValue[chanId] = p2;
      fChanPeakingTime[chanId] = p3;
      fChanRiseShape[chanId] = p4;
      fChanFallShape[chanId] = p5;

      chanProfile->Delete();
      
    } 

    return true;
  }

private:
  // The time stamp for the data generating this calibration.
  std::string fTimeStamp;

  // The run number of the calibration file.
  int fRunNumber;
    
  // The pulse into the capacitor
  double fInjectedPulse;

  // The capacitance
  double fCapacitance;

  // The injected charge;
  double fInjectedCharge;

  // The digitizer slope
  double fDigitizerSlope;

  std::map<CP::TChannelId,std::vector<int>> fChanStatus;
  std::map<CP::TChannelId,std::vector<double>> fChanPar0;
  std::map<CP::TChannelId,std::vector<double>> fChanGainValue;
  std::map<CP::TChannelId,std::vector<double>> fChanPeakingTime;
  std::map<CP::TChannelId,std::vector<double>> fChanRiseShape;
  std::map<CP::TChannelId,std::vector<double>> fChanFallShape;


};

int main(int argc, char **argv) {
  TCalibPulseFitLoop userCode;
  CP::eventLoop(argc,argv,userCode);
}
