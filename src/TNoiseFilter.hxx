#ifndef TNoiseFilter_hxx_seen
#define TNoiseFilter_hxx_seen

#include "TElectronicsResponse.hxx"
#include "TWireResponse.hxx"

#include <TChannelId.hxx>
#include <HEPUnits.hxx>

#include <TVirtualFFT.h>

namespace CP {
    class TNoiseFilter;
};


class CP::TNoiseFilter {
public:

    /// Create a noise filter for a particular measurement and response
    /// function.
    TNoiseFilter(double noisePower, double spikePower);
    virtual ~TNoiseFilter();

    /// Figure out the optimal filter.
    void Calculate(CP::TChannelId id,
                   CP::TElectronicsResponse& elecFreq,
                   CP::TWireResponse& wireFreq,
                   TVirtualFFT& measFreq);
    
    /// Return the filter value for a particular frequency bin.
    double GetFilter(int i) {return fFilter[i];}
    
private:
    /// The cached filter values.
    std::vector<double> fFilter;

    /// A work area.
    std::vector<double> fWork;
    
    /// A work area.
    std::vector<double> fAverage;
    
    /// A work area.
    std::vector<double> fSigma;
    
    /// The power of the gaussian noise.
    double fNoisePower;

    /// The maximum deviation of any bin in the FFT from the local average (in
    /// terms of sigma).
    double fSpikePower;
};
    
#endif
