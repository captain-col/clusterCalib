#ifndef TruncatedRMS_hxx_seen
#define TruncatedRMS_hxx_seen
#include <algorithm>
#include <vector>
#include <utility>
#include <cmath>

namespace CP {

    /// Calculated the truncated RMS for the samples between begin and end.
    /// \code
    /// const CP::TPulseDigit* pulse = digit.As<const CP::TPulseDigit>();
    /// std::pair<double,double> avgRMS
    ///          = CP::TruncatedRMS(pulse->begin(), pulse->end());
    /// std::cout << "Average " << avgRMS.first << "  RMS " << avgRMS.second
    ///           << std::endl;
    /// \endcode
    template <typename iter>
    std::pair<double,double> TruncatedRMS(iter begin, iter end,
                                          double trunc=0.05) {
        std::size_t length = end - begin;
        std::size_t lowBound = trunc*length;
        std::size_t highBound = (1.0-trunc)*length;
        if (highBound < lowBound+1) return std::pair<double,double>(0.0,0.0);
        static std::vector<typename iter::value_type> work;
        if (work.capacity() < length) work.reserve(length);
        work.clear();
        std::copy(begin,end,std::back_inserter(work));
        std::sort(work.begin(),work.end());
        double average = 0.0;
        double rms = 0.0;
        double norm = 0.0;
        typename std::vector<typename iter::value_type>::iterator first
            = work.begin() + lowBound;
        typename std::vector<typename iter::value_type>::iterator last
            = work.begin() + highBound;
        while (first != last) {
            average += *first;
            rms += (*first)*(*first);
            norm += 1.0;
            ++first;
        }
        average /= norm;
        rms /= norm;
        rms = (rms - average*average);
        if (rms > 0.0) rms = std::sqrt(rms);
        else rms = - std::sqrt(-rms);
        return std::pair<double,double>(average,rms);
    }
}
#endif
