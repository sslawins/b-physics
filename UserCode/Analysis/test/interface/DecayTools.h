#ifndef DECAYTOOLS_H
#define DECAYTOOLS_H

#include <vector>
#include <set>
#include "TMath.h"

namespace DecayTools {
    
    const std::vector<int> BsStarG = {22, 533};
    const std::vector<int> Bs = {531};
    const std::vector<int> PhiG = {333, 22};
    const std::vector<int> MuMu = {13, -13};
    const std::vector<int> KK = {321, -321};

    inline bool isSameChannel(const std::vector<int>& dec1, const std::vector<int>& dec2) {
        if (dec1.size() != dec2.size()) return false;
        return std::is_permutation(dec1.begin(), dec1.end(), dec2.begin());
    }

    inline math::XYZPoint pca(math::XYZPoint pv, math::XYZPoint sv, math::XYZVectorD pMuMu){
        //s = (PV - SV) * (pMuMu) / |pMuMu|^2 
        double s = ((pv - sv).Dot(pMuMu)) / pMuMu.Mag2();
        //pca = sv + s*pMuMu
        math::XYZPoint PCA = sv + s*pMuMu;

        return PCA;
    }

}
#endif