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
#include <CaptGeomId.hxx>
#include <TChannelInfo.hxx>
#include <TTPCChannelId.hxx>

#define SKIP_DATA_CALIBRATION

CP::TChannelCalib::TChannelCalib() { }

CP::TChannelCalib::~TChannelCalib() { } 

bool CP::TChannelCalib::IsGoodChannel(CP::TChannelId id) {
    if (id == TTPCChannelId(1,13,16)) {
        CaptLog("Bad Channel " << id);
        return false;
    }
    return true;
}

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

    /// Get the event context.
    CP::TEvent* ev = CP::TEventFolder::GetCurrentEvent();
    if (!ev) {
        CaptError("No event is loaded so context cannot be set.");
        throw EChannelCalibUnknownType();
    }
    CP::TEventContext context = ev->GetContext();

    /// Get the geometry id for the current wire.
    CP::TGeometryId geomId = CP::TChannelInfo::Get().GetGeometry(id);
    if (!geomId.IsValid()) {
        CaptError("Channel is not associated with a geometry object: " << id);
        return false;
    }
    if (!CP::GeomId::Captain::IsWire(geomId)) {
        CaptError("Channel " << id << " is not a wire: " << geomId);
        return false;
    }

    // Check for special runs.  This should really be done with a table, but
    // for now, it is hardcoded.
    if (context.IsMiniCAPTAIN()
        && 4090 <= context.GetRun()
        && context.GetRun() < 10000000) {
        // The X wires are disconnected for most of these runs, but return
        // "unipolar" anyway since it makes other studies easier.
        if (CP::GeomId::Captain::IsXWire(geomId)) return false;
        if (CP::GeomId::Captain::IsVWire(geomId)) return false;
        return true;
    }

    // Normally, the X wire is the collection, and the other wires (U, V) are
    // induction (i.e. bipolar).
    if (CP::GeomId::Captain::IsXWire(geomId)) return false;
    return true;

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
    if (order == 1) return (14.0*unit::mV/unit::fC);
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

        // Get the rising edge shape.
        CP::THandle<CP::TRealDatum> riseVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/shapeRise");
        double riseShape = (*riseVect)[index];
        
        // Get the rising edge shape.
        CP::THandle<CP::TRealDatum> fallVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/shapeFall");
        double fallShape = (*fallVect)[index];
        
        double arg = t/peakingTime;
        if (arg < 1.0) arg = std::pow(arg,riseShape);
        else arg = std::pow(arg,fallShape);
        
        double v = (arg<40)? arg*std::exp(-arg): 0.0;
        return v;
    }

#ifdef SKIP_DATA_CALIBRATION
    if (t < 0.0) return 0.0;

    double peakingTime = 1.0*unit::microsecond;
    double riseShape = 2.0;
    double fallShape = 2.0;

    double arg = t/peakingTime;
    if (arg < 1.0) arg = std::pow(arg,riseShape);
    else arg = std::pow(arg,fallShape);
    
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

        // Get the trigger offset.
        CP::THandle<CP::TRealDatum> offsetVect
            = ev->Get<CP::TRealDatum>("~/truth/elecSimple/triggerOffset");
        double triggerOffset = (*offsetVect)[index];

        if (order == 0) {
            // The offset for the digitization time for each type of MC channel.
            // This was determined "empirically", and depends on the details of
            // the simulation.  It should be in the parameter file.  The exact
            // value depends on the details of the time clustering.
            double timeOffset = -triggerOffset;
            if (index == 1) timeOffset += -digitStep + 23*unit::ns;
            else if (index == 2) timeOffset += -digitStep + 52*unit::ns;
            return timeOffset;
        }
        else if (order == 1) {
            return digitStep;
        }

        return 0.0;
    }

#ifdef SKIP_DATA_CALIBRATION
    if (order == 0) return -1.600*unit::ms;
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

