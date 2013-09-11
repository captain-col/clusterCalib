#include "TPulseCalib.hxx"

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

    double digitStep;
    double timeOffset;
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

        // The offset for the digitization time for each type of MC channel.
        // This was determined "empirically", and depends on the details of
        // the simulation.  It should be in the parameter file.  The exact
        // value depends on the details of the time clustering.
        if (index == 1) timeOffset = -digitStep + 23*unit::ns;
        else if (index == 2) timeOffset = -digitStep + 52*unit::ns;
        else timeOffset = 0.0;

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

    double startTime = digitStep*pulse->GetFirstSample() + timeOffset;
    double stopTime 
        = digitStep*(pulse->GetFirstSample()+pulse->GetSampleCount()) 
        + timeOffset;
    CP::TCalibPulseDigit::Vector samples(pulse->GetSampleCount());
    for (std::size_t i=0; i< pulse->GetSampleCount(); ++i) {
        samples[i] = 1.0*(pulse->GetSample(i)-pedestal)/gain;
    }
    return new CP::TCalibPulseDigit(digit,startTime,stopTime,samples);
}
