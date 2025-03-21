#include <iomanip>
#include "ParticleMatcher.h"

template <typename T>
ParticleMatcher<T>::ParticleMatcher(const std::vector<T>& recoParticles_, 
                         std::vector<const reco::Candidate*>& genParticles_,
                         double deltaRcut_)
    : recoParticles(recoParticles_), genParticles(genParticles_), recoMatched(std::vector< const reco::Candidate*>(genParticles.size(), nullptr)), deltaRcut(deltaRcut_) {}

template <typename T>
void ParticleMatcher::matchRecoToGen() {
    std::cout << "Number of all reconstructed muons/kaons from this event: " << recoParticles.size() << std::endl;
    std::cout << "Number of all generated muons/kaons from the decay: " << genParticles.size() << std::endl;

    for  (unsigned int genIter = 0; genIter < genParticles.size(); ++genIter ){
        double lowestRpart = 100.0; //high initial value
        std::cout << "Gen particle pt: " << genParticles[genIter]->pt() << ", reco particle pt: "  << std::endl;
        
        for (const auto& recoPart : recoParticles){
            std::cout << "          " << recoPart.pt();
            double deltaRVal = reco::deltaR(recoPart, *genParticles[genIter]);
            std::cout << "          " <<deltaRVal<<"/"<< lowestRpart ;
            
            if (deltaRVal < lowestRpart && deltaRVal < deltaRcut) {

                std::cout << "      pt of matched GEN particle: " << genParticles[genIter]->pt() ; 
                lowestRpart = deltaRVal;
                recoMatched[genIter] = &recoPart;
            }

            std::cout << std::endl;
        }

    }
}

template <typename T>
void ParticleMatcher::printMatchedParticles() {
    std::cout << std::left << std::setw(20) << "Gen Particle pt" 
              << std::setw(20) << "Reco Particle pt" << std::endl;
    std::cout << std::string(30, '-') << std::endl;

    for (unsigned int genIter = 0; genIter < genParticles.size(); ++genIter) {
        const auto& genPart = genParticles[genIter];
        const auto& matched = recoMatched[genIter];

        std::cout << std::left << std::setw(20) << genPart->pt()
                  << std::setw(20) << (matched ? matched->pt() : 0) << std::endl;
    }
}


