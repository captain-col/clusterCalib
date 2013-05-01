#include "TPulseCalib.hxx"

#include <TCalibPulseDigit.hxx>
#include <TPulseDigit.hxx>
#include <TChannelId.hxx>
#include <TMCChannelId.hxx>
#include <TEventFolder.hxx>
#include <TRealDatum.hxx>

#include <TEvent.hxx>

CP::TPulseCalib::TPulseCalib() {}
    
CP::TPulseCalib::~TPulseCalib() {}

CP::TCalibPulseDigit* CP::TPulseCalib::operator()(const CP::TDigit* digit) {
    const CP::TPulseDigit* pulse = dynamic_cast<const CP::TPulseDigit*>(digit);
    if (!pulse) return NULL;

    TChannelId chan(digit->GetChannelId());

    double digitStep;
    double pedestal;
    double gain;

    // Get the channel calibrations.
    if (chan.IsMCChannel()) {
        TMCChannelId mc(chan);
        int index = -1;
        if (mc.GetType() == 0) index = mc.GetSequence();
        else if (mc.GetType() == 1) index = 3;
        else {
            CaptError("Unknown channel: " << mc);
            return NULL;
        }
            
        CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
        
        // Get the digitization step.
        CP::THandle<CP::TRealDatum> stepVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/digitStep");
        digitStep = (*stepVect)[index];

        // Get the pedestal
        CP::THandle<CP::TRealDatum> pedVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/pedestal");
        pedestal = (*pedVect)[index];
        
        // Get the gain
        CP::THandle<CP::TRealDatum> gainVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/gain");
        gain = (*gainVect)[index];
    }
    else {
        // We don't have any calibrations yet!
        return NULL;
    }

    double startTime = digitStep*pulse->GetFirstSample();
    double stopTime 
        = digitStep*(pulse->GetFirstSample()+pulse->GetSampleCount());
    CP::TCalibPulseDigit::Vector samples(pulse->GetSampleCount());
    for (std::size_t i=0; i< pulse->GetSampleCount(); ++i) {
        samples[i] = 1.0*(pulse->GetSample(i)-pedestal)/gain;
    }
    return new CP::TCalibPulseDigit(pulse,startTime,stopTime,samples);
}
