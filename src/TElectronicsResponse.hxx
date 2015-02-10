#ifndef TElectronicsResponse_hxx_seen
#define TElectronicsResponse_hxx_seen

#include <TEventContext.hxx>
#include <TChannelId.hxx>
#include <HEPUnits.hxx>

#include <complex>

namespace CP {
    class TElectronicsResponse;
};

/// This class manages the shape of the delta function response for the
/// electronics.  It is used to fill an array which the response function
/// which can then be used in an FFT based deconvolution.
class CP::TElectronicsResponse {
public: 
    /// The real response.
    typedef double Response;

    /// The complex frequency.  This comes from a FFT of the response.
    typedef std::complex<double> Frequency;

    virtual ~TElectronicsResponse();

    /// Connstruct a response function without setting the number of samples,
    /// or the time per sample.  The must be set using the SetSize(), and
    /// SetStep() methods.
    TElectronicsResponse();

    /// Construct a response with a time per sample of "step", and "size"
    /// samples.  This means that the response function is for a time of
    /// size*step.
    TElectronicsResponse(int size);

    /// @{ Set (or Get) the number of samples.
    void SetSize(int s);
    std::size_t GetSize() const {return fResponse.size();}
    /// @}

    /// Get the response for a time bin.  This returns 0.0 if the bin is out
    /// of range (either negative of larger than the response distribution
    /// vector size).
    double GetResponse(int i);

    /// Get the FFT for a frequency bin.  This returns 0+0*i if the bin is out
    /// of the frequency range.
    Frequency GetFrequency(int i);

    /// Get the event context of the last event for which the response
    /// function was calculated.
    const CP::TEventContext& GetEventContext() const {return fEventContext;}

    /// Get the channel id of the last channel for which the response function
    /// was calculated.
    const CP::TChannelId& GetChannelId() const {return fChannelId;}

    /// Calculate the delta function response for a particular
    /// CP::TEventContext and CP::TChannelId.  This also calculates the FFT
    /// for use in a convolution or deconvolution.  This returns true if the
    /// response function will be recalculated, but doesn't do the actual
    /// recalculation.  The recalculation will be done on the first call to
    /// GetResponse, or GetFrequency.
    bool Calculate(const CP::TEventContext& context, const CP::TChannelId id);
    
    /// Calculate the delta function response for a particular peaking time.
    /// This version of the Calculate method always returns true, and does the
    /// actual recalculation.
    bool Calculate();
    
private:
    /// Return true if the delta function response has changed since the last
    /// time this class calculated a shape.  This assumes that there will be
    /// one TPulseShaping object per category of channel.  It is very likely
    /// that all of the collection wires can share one TElectronicsResponse
    /// object, all of the V induction wire use another, and all of the U
    /// induction wires use a third.  The another TElectronicsResponse object
    /// is required for the PMTs.
    // bool ShapeChanged(const CP::TEventContext& context, 
    // const CP::TChannelId id);

    /// Flag that something forces a recalculation.
    bool fMustRecalculate;

    /// The event context of the last calculation.
    CP::TEventContext fEventContext;

    /// The channel id of the last calculation.
    CP::TChannelId fChannelId;

    /// The response distribution
    std::vector<Response> fResponse;

    /// The frequency distribution
    std::vector<Frequency> fFrequency;
};
#endif
