#ifndef BDecayAnalyzer_H
#define BDecayAnalyzer_H

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

class BDecayAnalyzer {
public:
    BDecayAnalyzer();
    std::vector<const reco::Candidate*> getPhiProducts();
    std::vector<const reco::Candidate*> getPhotons();
    std::vector<std::vector<const reco::Candidate*>> analyzeBDecays(const std::vector<reco::GenParticle>&, const std::vector<int>);

    void printTheTree( std::vector<std::vector<const reco::Candidate*>> );
    bool analyzeEvent( edm::EventID );

private:

    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken;
    std::vector<const reco::Candidate*> phiProducts;
    std::vector<const reco::Candidate*> genPhotons;
};

#endif