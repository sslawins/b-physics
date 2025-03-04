#include <iomanip>
#include "MuonMatcher.h"

MuonMatcher::MuonMatcher(const std::vector<reco::Muon>& recoMuons_, 
                         std::vector<const reco::Candidate*>& genMuons_,
                         double deltaRcut_)
    : recoMuons(recoMuons_), genMuons(genMuons_), recoMatchedMu(std::vector< const reco::Candidate*>(genMuons.size(), nullptr)), deltaRcut(deltaRcut_) {}
    
void MuonMatcher::matchRecoToGen() {
    std::cout << "Number of all reconstructed muons from this event: " << recoMuons.size() << std::endl;
    std::cout << "Number of all generated muons from the decay: " << genMuons.size() << std::endl;

    for  (unsigned int genIter = 0; genIter < genMuons.size(); ++genIter ){
        double lowestRmuon = 100.0; //high initial value
        std::cout << "Gen muon pt: " << genMuons[genIter]->pt() << ", reco muons pt: "  << std::endl;
        
        for (const auto& recoMu : recoMuons){
            std::cout << "          " << recoMu.pt();
            double deltaRVal = reco::deltaR(recoMu, *genMuons[genIter]);
            std::cout << "          " <<deltaRVal<<"/"<< lowestRmuon ;
            
            if (deltaRVal < lowestRmuon && deltaRVal < deltaRcut) {

                std::cout << "      pt of matched GEN muon: " << genMuons[genIter]->pt() ; 
                lowestRmuon = deltaRVal;
                recoMatchedMu[genIter] = &recoMu;
            }

            std::cout << std::endl;
        }

    }
}

void MuonMatcher::printMatchedMuons() {
    std::cout << std::left << std::setw(15) << "Gen Muon pt" 
              << std::setw(15) << "Reco Muon pt" << std::endl;
    std::cout << std::string(30, '-') << std::endl;

    for (unsigned int genIter = 0; genIter < genMuons.size(); ++genIter) {
        const auto& genMu = genMuons[genIter];
        const auto& recoMatched = recoMatchedMu[genIter];

        std::cout << std::left << std::setw(15) << genMu->pt()
                  << std::setw(15) << (recoMatched ? recoMatched->pt() : 0) << std::endl;
    }
}


