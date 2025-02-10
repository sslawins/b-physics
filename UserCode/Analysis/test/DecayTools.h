#ifndef DECAYTOOLS_H
#define DECAYTOOLS_H

#include <vector>
#include <set>

namespace DecayTools {
    
    const std::vector<int> BsStarG = {22, 533};
    const std::vector<int> Bs = {531};
    const std::vector<int> PhiG = {333, 22};
    const std::vector<int> MuMu = {13, -13};

    inline bool isSameChannel(const std::vector<int>& dec1, const std::vector<int>& dec2) {
        return std::is_permutation(dec1.begin(), dec1.end(), dec2.begin());
    }

}
#endif