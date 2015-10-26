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
    fPedestal = calib.GetGainConstant(chan,0);
    double gain = calib.GetGainConstant(chan,1);
    double slope = calib.GetDigitizerConstant(chan,1);

    double startTime = digitStep*pulse->GetFirstSample() + timeOffset;
    double stopTime 
        = digitStep*(pulse->GetFirstSample()+pulse->GetSampleCount()) 
        + timeOffset;
    CP::TCalibPulseDigit::Vector samples(pulse->GetSampleCount());

    // Use the median ADC value as the pedestal
    for (std::size_t i=0; i< pulse->GetSampleCount(); ++i) {
        double p = pulse->GetSample(i);
        if (!std::isfinite(p)) continue;
        samples[i] = p;
    }
    std::sort(samples.begin(), samples.end());
    int iMedian = pulse->GetSampleCount()/2;
    std::pair<CP::TCalibPulseDigit::Vector::iterator,
              CP::TCalibPulseDigit::Vector::iterator> bounds;
    bounds = std::equal_range(samples.begin(), samples.end(),
                              samples[iMedian]);
    fPedestal = (samples.begin()+iMedian)-bounds.first;
    fPedestal /= bounds.second-bounds.first;
    fPedestal += std::min(samples[iMedian]-0.5,0.0);

    // Calculate the RMS around the average with the tails thrown out.
    fAverage = 0.0;
    fSigma = 0.0;
    double norm = 0.0;
    int lowBound = 0.05*pulse->GetSampleCount();
    int highBound = 0.95*pulse->GetSampleCount();
    for (int i=lowBound; i<highBound; ++i) {
        double p = samples[i];
        if (!std::isfinite(p)) continue;
        fAverage += p;
        fSigma += p*p;
        norm += 1.0;
    }
    fAverage = fAverage/norm;
    fSigma = fSigma/norm;
    fSigma = (fSigma - fAverage*fAverage);
    if (fSigma>0) fSigma = std::sqrt(fSigma);
    else fSigma = - std::sqrt(-fSigma);

    // Calculate the channel-to-channel Gaussian sigma.  The 52 percentile is
    // a of the sample-to-sample differences is a good estimate of the
    // distribution sigma.  The 52% "magic" number comes looking at the
    // differences of two samples (that gives a sqrt(2)) and wanting the RMS
    // (which is the 68%).  This gives 48% (which is 68%/sqrt(2)).  Then
    // because of the ordering after the sort, the bin we look at is 1.0-52%.
    for (std::size_t i=1; i< pulse->GetSampleCount(); ++i) {
        samples[i] = std::abs(pulse->GetSample(i)-pulse->GetSample(i-1));
    }
    samples[0] = 0.0;
    std::sort(samples.begin(), samples.end());
    iMedian = 0.52*pulse->GetSampleCount();
    std::pair<CP::TCalibPulseDigit::Vector::iterator,
              CP::TCalibPulseDigit::Vector::iterator> gaussianBounds;
    gaussianBounds = std::equal_range(samples.begin(), samples.end(),
                              samples[iMedian]);
    fGaussianSigma = (samples.begin()+iMedian)-gaussianBounds.first;
    fGaussianSigma /= gaussianBounds.second-gaussianBounds.first;
    fGaussianSigma += std::min(samples[iMedian]-0.5,0.0);

    // Actually apply the calibration.
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
