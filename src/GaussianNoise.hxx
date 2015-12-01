#ifndef GaussianNoise_hxx_seen
#define GaussianNoise_hxx_seen
#include <algorithm>
#include <vector>
#include <utility>

namespace CP {
    /// Calculate the channel-to-channel Gaussian sigma.  The 52 percentile is
    /// a of the sample-to-sample differences is a good estimate of the
    /// distribution sigma.  The 52% "magic" number comes looking at the
    /// differences of two samples (that gives a sqrt(2)) and wanting the RMS
    /// (which is the 68%).  This gives 48% (which is 68%/sqrt(2)).  Then
    /// because of the ordering after the sort, the bin we look at is 1.0-52%.
    template <typename iter>
    double GaussianNoise(iter begin, iter end) {
        if (begin == end) return 0.0;
        std::size_t length = end - begin;
        static std::vector<typename iter::value_type> work;
        if (work.capacity() < length) work.reserve(length);
        work.clear();
        while (begin+1 != end) {
            work.push_back(std::abs(*begin - *(begin+1)));
            ++begin;
        }
        std::sort(work.begin(),work.end());
        typename std::vector<typename iter::value_type>::iterator middle
            = work.begin() + 0.52*work.size();
        if (middle == work.begin()) return *middle;
        std::pair<iter,iter> bounds
            = std::equal_range(work.begin(), work.end(), *middle);
        if (bounds.second == work.end()) {
            return *middle;
        }
        double diff = 0.5*(*bounds.second - *middle);
        double pedestal = 1.0*(middle-bounds.first);
        pedestal /= 1.0*(bounds.second-bounds.first);
        pedestal += *middle - diff;
        return pedestal;
    }
}
#endif
