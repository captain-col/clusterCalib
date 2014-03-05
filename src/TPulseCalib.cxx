#include "TPulseCalib.hxx"
#include "TChannelCalib.hxx"

#include <TCalibPulseDigit.hxx>
#include <TPulseDigit.hxx>
#include <TChannelId.hxx>
#include <TMCChannelId.hxx>
#include <TEventFolder.hxx>
#include <TRealDatum.hxx>
#include <HEPUnits.hxx>
#include <TEvent.hxx>

CP::TPulseCalib::TPulseCalib() {}
    
CP::TPulseCalib::~TPulseCalib() {}

CP::TCalibPulseDigit* 
CP::TPulseCalib::operator()(const CP::TDigitProxy& digit) {
    const CP::TPulseDigit* pulse = digit.As<const CP::TPulseDigit>();
    if (!pulse) return NULL;

    TChannelId chan(pulse->GetChannelId());

    CP::TChannelCalib calib;

    double timeOffset = calib.GetTimeConstant(chan,0);
    double digitStep = calib.GetTimeConstant(chan,1);
    double pedestal = calib.GetGainConstant(chan,0);
    double gain = calib.GetGainConstant(chan,1);
    double slope = calib.GetDigitizerConstant(chan,1);

    double startTime = digitStep*pulse->GetFirstSample() + timeOffset;
    double stopTime 
        = digitStep*(pulse->GetFirstSample()+pulse->GetSampleCount()) 
        + timeOffset;
    CP::TCalibPulseDigit::Vector samples(pulse->GetSampleCount());
    for (std::size_t i=0; i< pulse->GetSampleCount(); ++i) {
        samples[i] = 1.0*(pulse->GetSample(i)-pedestal)/gain/slope;
    }
    return new CP::TCalibPulseDigit(digit,startTime,stopTime,samples);
}
