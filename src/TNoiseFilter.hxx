#ifndef TNoiseFilter_hxx_seen
#define TNoiseFilter_hxx_seen

#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"

#include <TChannelId.hxx>
#include <HEPUnits.hxx>

#include <TVirtualFFT.h>
#include <RVersion.h>

namespace CP {
    class TNoiseFilter;
};


/// A noise filter for a wire based on the assumed electronics response, and
/// the power spectrum measured on the wire.  The filter is based on an
/// estimate of the Gaussian noise seen in the wire power spectrum, and a
/// "notch" filter to remove any fixed frequency noise present on the wire.
class CP::TNoiseFilter {
public:

    /// Create a noise filter for a particular measurement and response
    /// function.  The noisePower is a multiplicitive factor applied to the
    /// estimated Gaussian noise.  A value of zero means no filtering, a large
    /// value means more filtering.  A value of 1.0 means that the raw
    /// estimated noise is used in the filter (a good first guess).  The
    /// spikePower is a multiplicitive factor applied to the "notch" filter.
    /// Zero means that the notch filter is not applied.  A value of 1.0 means
    /// that the esimated power for the fixed frequency noise is directly used
    /// in the filter (a good first guess).
    TNoiseFilter(double noisePower, double spikePower);
    virtual ~TNoiseFilter();

    /// Figure out the optimal filter.
    void Calculate(CP::TChannelId id,
                   CP::TElectronicsResponse& elecFreq,
                   CP::TWireResponse& wireFreq,
                   TVirtualFFT& measFreq);
    
    /// Return the filter value for a particular frequency bin.
    double GetFilter(int i) const {return fFilter[i];}
    
    /// Return a boolean whether this channel should be flagged as "noisy"
    bool IsNoisy() const {return fIsNoisy;} 

private:
    /// The cached filter values.
    std::vector<double> fFilter;

    /// A work area.
    std::vector<double> fWork;
    
    /// A work area.
#if ROOT_VERSION(6,0,0) < ROOT_VERSION_CODE
    std::vector<double> fAverage;
#else
    std::vector<float> fAverage;
#endif
    
    /// The power of the gaussian noise.
    double fNoisePower;

    /// The maximum deviation of any bin in the FFT from the local average (in
    /// terms of sigma).
    double fSpikePower;

    /// A flag for if this is a noisy channel
    bool fIsNoisy;
    
};
    
#endif
