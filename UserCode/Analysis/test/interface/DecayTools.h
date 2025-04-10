#ifndef DECAYTOOLS_H
#define DECAYTOOLS_H

#include <vector>
#include <set>
#include "TMath.h"
#include "Math/Vector4D.h"

namespace DecayTools {
    
    const std::vector<int> BsStarG = {22, 533};
    const std::vector<int> Bs = {531};
    const std::vector<int> PhiG = {333, 22};
    const std::vector<int> MuMu = {13, -13};
    const std::vector<int> KK = {321, -321};

    const std::vector<double> MuMuGmasses = {0.105658, 0.105658, 0.0};
    const std::vector<double> MuMuMasses = {0.105658, 0.105658};
    const std::vector<double> KKMasses = {0.493677, 0.493677};
    const std::vector<double> KKGmasses = {0.493677, 0.493677, 0.0};

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

    template <typename T>
    double invariantMass(const std::vector<T>& particles, const std::vector<double>& masses) {

        if (particles.size() != masses.size()) {
            throw std::invalid_argument("Size mismatch: particles and masses vectors must have the same length");
        }

        ROOT::Math::PxPyPzEVector totalP4(0, 0, 0, 0);

        for (size_t i = 0; i < particles.size(); ++i) {
            ROOT::Math::PxPyPzEVector p4(particles[i]->px(), particles[i]->py(), particles[i]->pz(),
                                        std::sqrt(masses[i] * masses[i] + particles[i]->p() * particles[i]->p()));
            totalP4 += p4;
        }

        return totalP4.M();
    }

    template <typename T>
    ROOT::Math::PxPyPzEVector fourMomenta(const std::vector<T>& particles, const std::vector<double>& masses) {
        if (particles.size() != masses.size()) {
            throw std::invalid_argument("Size mismatch: particles and masses vectors must have the same length");
        }

        ROOT::Math::PxPyPzEVector totalP4(0, 0, 0, 0);

        for (size_t i = 0; i < particles.size(); ++i) {
            ROOT::Math::PxPyPzEVector p4(
                particles[i]->px(),
                particles[i]->py(),
                particles[i]->pz(),
                std::sqrt(masses[i] * masses[i] + particles[i]->p() * particles[i]->p())
            );
            totalP4 += p4;
        }

        return totalP4;
    }

    template <typename T>
    ROOT::Math::XYZVector scaledP(const T& length, const T& direction) {
        ROOT::Math::XYZVector vecDir(direction->px(), direction->py(), direction->pz());
        double magnitude = length->p();  // długość pędu

        return vecDir.Unit() * magnitude;
    }


    template <typename T, typename MagT, typename DirT>
    double scaledInvariant(
        const std::vector<T>& particles,
        const MagT& magnitudeFrom,                  
        const DirT& directionFrom,                  
        const std::vector<double>& masses
    ) {
        if (particles.size() + 1 != masses.size()) {
            throw std::invalid_argument("Mismatch between particles and mass vector sizes");
        }
    
        ROOT::Math::PxPyPzEVector totalP4(0, 0, 0, 0);
    
        // Dodajemy cząstki z listy (miony, kaony itp.)
        for (size_t i = 0; i < particles.size(); ++i) {
            ROOT::Math::PxPyPzEVector p4(
                particles[i]->px(),
                particles[i]->py(),
                particles[i]->pz(),
                std::sqrt(masses[i] * masses[i] + particles[i]->p() * particles[i]->p())
            );
            totalP4 += p4;
        }
    
        // Photon z przeskalowanym pędem
        ROOT::Math::XYZVector scaledPhotonVec = scaledP(magnitudeFrom, directionFrom);
    
        ROOT::Math::PxPyPzEVector photonP4(
            scaledPhotonVec.X(),
            scaledPhotonVec.Y(),
            scaledPhotonVec.Z(),
            std::sqrt(masses.back() * masses.back() + scaledPhotonVec.Mag2())
        );
    
        totalP4 += photonP4;
    
        return totalP4.M();
    }
    
}

#endif