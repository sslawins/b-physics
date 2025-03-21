#ifndef PARTICLEMATCHER_H
#define PARTICLEMATCHER_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include "TH1D.h"

template <typename T>
class ParticleMatcher {
public:
    ParticleMatcher(const std::vector<T>& recoParticles, 
                std::vector< const reco::Candidate*>& genParticles,
                double deltaRcut,
                TH1D* hKK_deltaR);

    ParticleMatcher(const std::vector<T>& recoParticles, 
                std::vector< const reco::Candidate*>& genParticles,
                double deltaRcut);

    void matchRecoToGen(bool print = false);

    const std::vector<const reco::Candidate*>& getMatched() const { return recoMatched; }
    void printMatchedParticles();
    bool isSuccessful() { return    (recoMatched.size() == genParticles.size() && 
                                    (std::find(recoMatched.begin(), recoMatched.end(), nullptr) == recoMatched.end()) &&
                                    (recoMatched.size() != 0)) ; }

    ~ParticleMatcher() {};

private:
    const std::vector<T>& recoParticles;
    std::vector<const reco::Candidate*>& genParticles;
    std::vector<const reco::Candidate*> recoMatched;

    double deltaRcut;
    TH1D* hKK_deltaR;
};

template <typename T>
ParticleMatcher<T>::ParticleMatcher(const std::vector<T>& recoParticles_, 
                         std::vector<const reco::Candidate*>& genParticles_,
                         double deltaRcut_,
                         TH1D* hKK_deltaR_)
    : recoParticles(recoParticles_), genParticles(genParticles_), recoMatched(std::vector< const reco::Candidate*>(genParticles.size(), nullptr)), deltaRcut(deltaRcut_), hKK_deltaR(hKK_deltaR_) {}

template <typename T>
ParticleMatcher<T>::ParticleMatcher(const std::vector<T>& recoParticles_, 
                         std::vector<const reco::Candidate*>& genParticles_,
                         double deltaRcut_)
    : recoParticles(recoParticles_), genParticles(genParticles_), recoMatched(std::vector< const reco::Candidate*>(genParticles.size(), nullptr)), deltaRcut(deltaRcut_), hKK_deltaR(nullptr) {}

template <typename T>
void ParticleMatcher<T>::matchRecoToGen(bool print) {
    std::cout << "Number of all reconstructed muons/kaons from this event: " << recoParticles.size() << std::endl;
    std::cout << "Number of all generated muons/kaons from the decay: " << genParticles.size() << std::endl;

    for (unsigned int genIter = 0; genIter < genParticles.size(); ++genIter ){
        double lowestRpart = 100.0; //high initial value
        if(print) std::cout << "Gen particle pt: " << genParticles[genIter]->pt() << ", reco particle pt: "  << std::endl;
        
        for (const auto& recoPart : recoParticles){
            if(print) std::cout << "          " << recoPart.pt();
            //std::cout << recoPart.pdgId() << std::endl;
            double deltaRVal = reco::deltaR(recoPart, *genParticles[genIter]);
            if (hKK_deltaR) hKK_deltaR->Fill(deltaRVal);
            if(print) std::cout << "          " <<deltaRVal<<"/"<< lowestRpart ;
            
            if (deltaRVal < lowestRpart && deltaRVal < deltaRcut && recoPart.charge() == genParticles[genIter]->charge()) {

                if(print) std::cout << "      pt of matched GEN particle: " << genParticles[genIter]->pt() ; 
                lowestRpart = deltaRVal;
                recoMatched[genIter] = &recoPart;
            }

            if(print) std::cout << std::endl;
        }

    }
}

template <typename T>
void ParticleMatcher<T>::printMatchedParticles() {
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

#endif // PARTICLEMATCHER_H
