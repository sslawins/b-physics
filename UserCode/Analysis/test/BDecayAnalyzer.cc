#include "BDecayAnalyzer.h"
#include "DecayTools.h"

BDecayAnalyzer::BDecayAnalyzer()
    : genMuons(), genPhotons()  {}

std::vector<const reco::Candidate*> BDecayAnalyzer::getMuons() {
  return genMuons;
}

std::vector<const reco::Candidate*> BDecayAnalyzer::getPhotons() {
  return genPhotons;
}

std::vector<std::vector<const reco::Candidate*>> BDecayAnalyzer::analyzeBDecays(const std::vector<reco::GenParticle>& genParticles) {
    genMuons.clear();
    genPhotons.clear();
    std::vector<std::vector<const reco::Candidate*>> bTree; //a vector to contain the tree
    bool initial = true;
    std::vector<const reco::Candidate*> lineage;

    for (const auto& part : genParticles) {
        if (abs(part.pdgId()) / 100 == 5) { // Check if the particle is a B meson
        
          for(size_t momIter = 0; momIter < part.numberOfMothers(); ++momIter){
            const reco::Candidate* mother = part.mother(momIter);
            if ( std::abs(mother->pdgId()) / 100 == 5 ) {
              // if it's an initial B, create a vector for its descendants and change the bool value
              initial = false;
              
              break; //no need to continue
            }
          }

          //if it is an initial B, upload the descending B into the vector
          if( initial ){
            lineage.clear();
            const reco::Candidate* ancestor = static_cast<const reco::Candidate*>(&part);
            lineage.push_back(ancestor); //we have a vector for the evolution of one b quark and the initial particle a the beginning

            while (ancestor != nullptr){
              std::vector<int> decay; //daughters, with sign
              bool last = true;

              for (unsigned int dauIter = 0; dauIter < ancestor->numberOfDaughters(); ++dauIter) {
                const reco::Candidate* daughter = ancestor -> daughter(dauIter);
                decay.push_back( daughter->pdgId() );

                if (last == true && (abs(daughter->pdgId()) / 100) == 5) { //mother B found

                    lineage.push_back(daughter);
                    ancestor = daughter; // new ancestor
                    
                    last = false;
                    
                    break; // stop, and enter the while loop once again
                }
              }
              
              if (last){
                //std::cout << "Pdg ID " << abs(part.pdgId()) <<std::endl;
                //if( abs(lineage.back()->pdgId()) == 531 ) hBssPt -> Fill( lineage.back()->pt() ) ;  
                //std::cout << "Pdg ID " << abs(lineage.back()->pdgId()) <<std::endl;

                if( DecayTools::isSameChannel(decay, DecayTools::PhiG)  ) {

                  for( size_t i = 0; i < lineage.back()->numberOfDaughters(); ++i ){
                    
                    if(lineage.back()-> daughter(i) ->pdgId()  == 333){ //phi
                      std::vector<int> phiDecay;
                      const reco::Candidate* phi = lineage.back()-> daughter(i);
                      
                      std::cout << " Phi decay, number of daughters: " <<  phi -> numberOfDaughters() << std::endl;
                      std::cout << "Phi daughters: " ;
                      for( size_t dauPhi = 0; dauPhi < phi -> numberOfDaughters(); ++dauPhi ){
                        std::cout << phi -> daughter(dauPhi)->pdgId() << " " ;
                        phiDecay.push_back( phi ->daughter(dauPhi)->pdgId() );
                    
                      }
                      std::cout << std::endl;
                      
                      if( DecayTools::isSameChannel(phiDecay, DecayTools::MuMu)){
                        
                        genMuons.push_back(phi->daughter(0));
                        genMuons.push_back(phi->daughter(1));
                      }

                    }

                    if(lineage.back()-> daughter(i) ->pdgId()  == 22){
                      
                      genPhotons.push_back( lineage.back()-> daughter(i));
                    }
                    
                  }

                  if( bTree.size() != 0){
                    bTree.insert(bTree.begin(), lineage);
                  }else{
                    bTree.push_back(lineage);
                  }
                }else{
                  bTree.push_back(lineage);
                }
                break;

              }
            }
          }
        }
      }  
    return bTree;
}


