#ifndef TWireResponse_hxx_seen
#define TWireResponse_hxx_seen

#include <TEventContext.hxx>
#include <TChannelId.hxx>
#include <HEPUnits.hxx>

#include <complex>

namespace CP {
    class TWireResponse;
};

/// This class manages the shape of the delta function response for the
/// electronics.  It is used to fill an array which the response function
/// which can then be used in an FFT based deconvolution.
class CP::TWireResponse {
public: 
    /// The real response.
    typedef double Response;

    /// The complex frequency.  This comes from a FFT of the response.
    typedef std::complex<double> Frequency;

    virtual ~TWireResponse();

    /// Connstruct a response function without setting the number of samples.
    /// The must be set using the SetSize() method.
    TWireResponse();

    /// Construct a response with "size" samples.
    TWireResponse(int size);

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

    /// Actually fill the frequency from the response.
    bool Calculate();
    
private:

    /// Flag that something forces a recalculation.
    bool fMustRecalculate;

    /// The event context of the last calculation.
    CP::TEventContext fEventContext;

    /// The channel id of the last calculation.
    CP::TChannelId fChannelId;

    /// The class of the wire response function shape.  0) delta function, 1)
    /// derivative.
    enum {kDeltaFunction, kTimeDerivative};
    int fWireClass;

    /// The response distribution
    std::vector<Response> fResponse;

    /// The frequency distribution
    std::vector<Frequency> fFrequency;
};
#endif
