#include "TChannelCalib.hxx"

#include <TEvent.hxx>
#include <THandle.hxx>
#include <TRealDatum.hxx>
#include <TEventFolder.hxx>
#include <TCaptLog.hxx>
#include <TChannelId.hxx>
#include <TMCChannelId.hxx>
#include <HEPUnits.hxx>
#include <TRuntimeParameters.hxx>

#define SKIP_DATA_CALIBRATION

CP::TChannelCalib::TChannelCalib() { }

CP::TChannelCalib::~TChannelCalib() { } 

bool CP::TChannelCalib::IsBipolarSignal(CP::TChannelId id) {
    if (id.IsMCChannel()) {
        TMCChannelId mc(id);

        int index = -1;
        if (mc.GetType() == 0) index = mc.GetSequence();
        else if (mc.GetType() == 1) index = 3;
        else {
            CaptError("Unknown channel: " << id);
            throw CP::EChannelCalibUnknownType();
        }
            
        switch (index) {
        case 1: case 2: return true;
        default: return false;
        }
    }

#ifdef SKIP_DATA_CALIBRATION
    return false;
#endif

    CaptError("Unknown channel: " << id);
    throw EChannelCalibUnknownType();
    return false;
}

double CP::TChannelCalib::GetGainConstant(CP::TChannelId id, int order) {
    if (id.IsMCChannel()) {
        TMCChannelId mc(id);

        int index = -1;
        if (mc.GetType() == 0) index = mc.GetSequence();
        else if (mc.GetType() == 1) index = 3;
        else {
            CaptError("Unknown channel: " << id);
            throw CP::EChannelCalibUnknownType();
        }
            
        CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
        
        if (order == 0) {
            // Get the pedestal
            CP::THandle<CP::TRealDatum> pedVect
                = ev->Get<CP::TRealDatum>("~/truth/elecSimple/pedestal");
            return (*pedVect)[index];
        }
        else if (order == 1) {
            // Get the gain
            CP::THandle<CP::TRealDatum> gainVect
                = ev->Get<CP::TRealDatum>("~/truth/elecSimple/gain");
            return (*gainVect)[index];
        }

        return 0.0;
    }

#ifdef SKIP_DATA_CALIBRATION
    if (order == 1) return 4.7*unit::mV/unit::fC;
    return 2048.0;
#endif

    CaptError("Unknown channel: " << id);
    throw EChannelCalibUnknownType();
    return 0.0;
}

double CP::TChannelCalib::GetPulseShape(CP::TChannelId id, double t) {
    if (id.IsMCChannel()) {
        if (t < 0.0) return 0.0;

        TMCChannelId mc(id);
        int index = -1;
        if (mc.GetType() == 0) index = mc.GetSequence();
        else if (mc.GetType() == 1) index = 3;
        else {
            CaptError("Unknown channel: " << id);
            throw CP::EChannelCalibUnknownType();
        }
            
        CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
        
        // Get the rise time.
        CP::THandle<CP::TRealDatum> shapeVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/shape");
        double peakingTime = (*shapeVect)[index];

        
        double arg = t/peakingTime;
        double v = (arg<40)? arg*std::exp(-arg): 0.0;
        return v;
    }

#ifdef SKIP_DATA_CALIBRATION
    if (t < 0.0) return 0.0;
    double arg = std::pow(t/(1.0*unit::microsecond), 2.0);
    double v = (arg<40)? arg*std::exp(-arg): 0.0;
    return v;
#endif
    
    CaptError("Unknown channel: " << id);
    throw EChannelCalibUnknownType();
    return 0.0;
}

double CP::TChannelCalib::GetTimeConstant(CP::TChannelId id, int order) {
    if (id.IsMCChannel()) {
        TMCChannelId mc(id);

        int index = -1;
        if (mc.GetType() == 0) index = mc.GetSequence();
        else if (mc.GetType() == 1) index = 3;
        else {
            CaptError("Unknown channel: " << id);
            throw CP::EChannelCalibUnknownType();
        }
            
        CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
        
        // Get the digitization step.
        CP::THandle<CP::TRealDatum> stepVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/digitStep");
        double digitStep = (*stepVect)[index];

        if (order == 0) {
            // The offset for the digitization time for each type of MC channel.
            // This was determined "empirically", and depends on the details of
            // the simulation.  It should be in the parameter file.  The exact
            // value depends on the details of the time clustering.
            double timeOffset = 0.0;
            if (index == 1) timeOffset = -digitStep + 23*unit::ns;
            else if (index == 2) timeOffset = -digitStep + 52*unit::ns;
            return timeOffset;
        }
        else if (order == 1) {
            return digitStep;
        }

        return 0.0;
    }

#ifdef SKIP_DATA_CALIBRATION
    if (order == 1) return 500.0*unit::ns;
    return 0.0;
#endif

    CaptError("Unknown channel: " << id);
    throw EChannelCalibUnknownType();
    return 0.0;
}

double CP::TChannelCalib::GetDigitizerConstant(CP::TChannelId id, int order) {
    if (id.IsMCChannel()) {
        TMCChannelId mc(id);

        int index = -1;
        if (mc.GetType() == 0) index = mc.GetSequence();
        else if (mc.GetType() == 1) index = 3;
        else {
            CaptError("Unknown channel: " << id);
            throw CP::EChannelCalibUnknownType();
        }
            
        CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
        
        if (order == 0) {
            return 0.0;
        }
        else if (order == 1) {
            // Get the digitizer slope
            CP::THandle<CP::TRealDatum> slopeVect
                = ev->Get<CP::TRealDatum>("~/truth/elecSimple/slope");
            return (*slopeVect)[index];
        }

        return 0.0;
    }

#ifdef SKIP_DATA_CALIBRATION
    if (order == 1) return 2.5/unit::mV;
    return 0.0;
#endif

    CaptError("Unknown channel: " << id);
    throw EChannelCalibUnknownType();
    return 0.0;
}

double CP::TChannelCalib::GetElectronLifetime() {
    CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
    
    // Get the electron lifetime.
    CP::THandle<CP::TRealDatum> stepVect
        = ev->Get<CP::TRealDatum>("~/truth/elecSimple/argon");
    if (!stepVect) return 3.14E+8*unit::second;
    return (*stepVect)[1];
}

double CP::TChannelCalib::GetElectronDriftVelocity() {
    CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
    
    // Get the electron lifetime.
    CP::THandle<CP::TRealDatum> stepVect
        = ev->Get<CP::TRealDatum>("~/truth/elecSimple/argon");
    if (!stepVect) return 1.6*unit::millimeter/unit::microsecond;
    return (*stepVect)[0];
}

double CP::TChannelCalib::GetCollectionEfficiency(CP::TChannelId id) {
    if (id.IsMCChannel()) {
        TMCChannelId mc(id);

        int index = -1;
        if (mc.GetType() == 0) index = mc.GetSequence();
        else if (mc.GetType() == 1) return 1.0;
        else {
            CaptError("Unknown channel: " << id);
            throw CP::EChannelCalibUnknownType();
        }
            
        if (index == 0) {
            // A X channel.
            double eff = CP::TRuntimeParameters::Get().GetParameterD(
                "clusterCalib.mc.wire.collection.x");
            return eff;
        }
        else if (index == 1) {
            // A V channel.
            double eff = CP::TRuntimeParameters::Get().GetParameterD(
                "clusterCalib.mc.wire.collection.v");
            return eff;
        }
        else if (index == 2) {
            // A U channel.
            double eff = CP::TRuntimeParameters::Get().GetParameterD(
                "clusterCalib.mc.wire.collection.u");
            return eff;
        }
    }

#ifdef SKIP_DATA_CALIBRATION
    return 1.0;
#endif


    CaptError("Unknown channel: " << id);
    throw EChannelCalibUnknownType();
    return 1.0;
}

