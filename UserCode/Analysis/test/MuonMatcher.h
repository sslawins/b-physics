#ifndef MUONMATCHER_H
#define MUONMATCHER_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <vector>
#include <iostream>

class MuonMatcher {
public:
    MuonMatcher(const std::vector<reco::Muon>& recoMuons, 
                std::vector< const reco::Candidate*>& genMuons,
                double deltaRcut);

    void matchRecoToGen();

    const std::vector<const reco::Candidate*>& getMatched() const { return recoMatchedMu; }
    void printMatchedMuons();
    bool isSuccessful() { return    (recoMatchedMu.size() == genMuons.size() && 
                                    (std::find(recoMatchedMu.begin(), recoMatchedMu.end(), nullptr) == recoMatchedMu.end())); }

    ~MuonMatcher() {};

private:
    const std::vector<reco::Muon>& recoMuons;
    std::vector<const reco::Candidate*>& genMuons;

    std::vector< const reco::Candidate*> recoMatchedMu;

    double deltaRcut;
};

#endif // MUONMATCHER_H
