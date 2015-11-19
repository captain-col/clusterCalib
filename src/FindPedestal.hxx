#ifndef FindPedestal_hxx_seen
#define FindPedestal_hxx_seen
#include <algorithm>
#include <vector>
#include <utility>

namespace CP {
    template <typename iter>
        double FindPedestal(iter begin, iter end) {
        std::size_t length = end - begin;
        std::vector<typename iter::value_type> work;
        work.reserve(length);
        std::copy(begin,end,std::back_inserter(work));
        std::sort(work.begin(),work.end());
        typename std::vector<typename iter::value_type>::iterator middle
            = work.begin() + work.size()/2;
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
