#include "TPulseCalib.hxx"
#include "TChannelCalib.hxx"
#include "FindPedestal.hxx"
#include "GaussianNoise.hxx"
#include "TruncatedRMS.hxx"

#include <TCalibPulseDigit.hxx>
#include <TPulseDigit.hxx>
#include <TChannelId.hxx>
#include <TMCChannelId.hxx>
#include <TEventFolder.hxx>
#include <TRealDatum.hxx>
#include <HEPUnits.hxx>
#include <TEvent.hxx>

CP::TPulseCalib::TPulseCalib() 
    : fPedestal(0), fAverage(0), fSigma(0), fGaussianSigma(0) {}
    
CP::TPulseCalib::~TPulseCalib() {}

CP::TCalibPulseDigit* 
CP::TPulseCalib::operator()(const CP::TDigitProxy& digit) {
    const CP::TPulseDigit* pulse = digit.As<const CP::TPulseDigit>();
    if (!pulse) return NULL;

    TChannelId chan(pulse->GetChannelId());

    CP::TChannelCalib calib;

    double timeOffset = calib.GetTimeConstant(chan,0);
    double digitStep = calib.GetTimeConstant(chan,1);
    fPedestal = calib.GetDigitizerConstant(chan,0);
    double slope = calib.GetDigitizerConstant(chan,1);
    double gain = calib.GetGainConstant(chan,1);

    double startTime = digitStep*pulse->GetFirstSample() + timeOffset;
    double stopTime 
        = digitStep*(pulse->GetFirstSample()+pulse->GetSampleCount()) 
        + timeOffset;

    // Use the median ADC value as the pedestal
    fPedestal = CP::FindPedestal(pulse->begin(), pulse->end());;

    std::pair<double,double> avgRMS
        = CP::TruncatedRMS(pulse->begin(), pulse->end());

    fAverage = avgRMS.first;
    fSigma = avgRMS.second;
    fGaussianSigma = CP::GaussianNoise(pulse->begin(), pulse->end());
    
    // Actually apply the calibration.
    CP::TCalibPulseDigit::Vector samples(pulse->GetSampleCount());
    for (std::size_t i=0; i< pulse->GetSampleCount(); ++i) {
        double p = 1.0*(pulse->GetSample(i)-fPedestal)/gain/slope;
        samples[i] = p;
        if (!std::isfinite(p)) {
            CaptError("Channel " << pulse->GetChannelId() 
                      << " w/ invalid sample " << i);
        }
    }
    return new CP::TCalibPulseDigit(digit,startTime,stopTime,samples);
}
