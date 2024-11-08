#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <sstream>
#include <iomanip> 
#include <utility>
#include <numeric>
#include <vector>
#include <stack>


using namespace std;


//object definition
class ControlChannelPhi : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit ControlChannelPhi(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~ControlChannelPhi();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool isSameChannel(const std::vector<int>&, const std::vector<int>&);
  std::vector<std::vector<const reco::Candidate*>> bFamilyTree(const std::vector<reco::GenParticle>&);

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  //Bs and phi decay
  vector<const reco::Candidate*> genMuons;
  vector<const reco::Muon*> recoMatchedMuons;

  vector<const reco::Candidate*> genPhotons;
  vector<const reco::Photon*> recoMatchedPhotons;

  TH1D *hBPt;
  TH1D *hBssPt;

  TH1D *hBsToPhiG;
  TH1D *hPhiToMuMu;
  TH1D *hBsToB;

  TH1D *hNinit;
  TH1D *hInitial;
  TH1D *hHad; 
  TH1D *hBsParents ;
  TH1D *hBsAncestor ; 
  //TH1D *hBsProduct ; 

  TH1D *hCascade;
  TH1D *hMuMuDeltaR;
  TH1D *hMuMuDeltaPhi;
  TH1D *hMuMuDeltaEta;
  TH2D *hMuEta_MuEta ;
  TH1D *hMuPt ;
  TH1D *hMuPtAll ;

  TH1D *hPhiPt;
  TH1D *hPhiEta;
  TH2D *hGEta_GPt;
  TH2D *hPhiEta_PhiPt;
  TH1D *hGammaPt;
  TH1D *hGammaEta;

  TH1D *hGGen_GReco_DeltaR;
  TH1D *hRecoGammaPt;


  TH1D* hGen_MuPt;
  TH1D* hRecoMu_Pt;
  TH2D* hGen_MuEta_Pt;
  TH1D* hMuGen_MuReco_DeltaR;

  TH2D* hGenGammaPtVsEta;

  TH2D* hRecoMu_Eta_Pt;

  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;
  edm::EDGetTokenT < vector<reco::Muon> > theMuonToken;
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;

  std::vector<int> MuMuG = {22, 13, -13};
  std::vector<int> BsStarG= {22, 533};
  std::vector<int> Bs= {531};
  std::vector<int> PhiG= {333, 22};

};


ControlChannelPhi::ControlChannelPhi(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles" ));

}

ControlChannelPhi::~ControlChannelPhi()
{
  cout <<" DTOR" << endl;
}

bool ControlChannelPhi::isSameChannel(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}

std::vector<std::vector<const reco::Candidate*>> ControlChannelPhi::bFamilyTree(const std::vector<reco::GenParticle>& genParticles) {
    
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
            hInitial -> Fill( abs(part.pdgId()) );
            lineage.clear();
            const reco::Candidate* ancestor = static_cast<const reco::Candidate*>(&part);
            lineage.push_back(ancestor);

            while (ancestor != nullptr){
              std::vector<int> decay;
              bool last = true;

              for (unsigned int dauIter = 0; dauIter < ancestor->numberOfDaughters(); ++dauIter) {
                const reco::Candidate* daughter = ancestor -> daughter(dauIter);
                decay.push_back( daughter->pdgId() );

                if (last== true && (abs(daughter->pdgId()) / 100) == 5) { //mother B found

                    lineage.push_back(daughter);
                    ancestor = daughter; // new ancestor
                    
                    last = false;
                    
                    break; // stop, and enter the while loop once again
                }
              }
              
              if (last){
                hBPt -> Fill ( part.pt() );
                std::cout << "Pdg ID " << abs(part.pdgId()) <<std::endl;
                if( abs(lineage.back()->pdgId()) == 531 ) hBssPt -> Fill( lineage.back()->pt() ) ;  
                std::cout << "Pdg ID " << abs(lineage.back()->pdgId()) <<std::endl;

                if( isSameChannel(decay, PhiG)  ) {
                  
                  for( size_t i = 0; i < lineage.back()->numberOfDaughters(); ++i ){
                    if(lineage.back()-> daughter(i) ->pdgId()  == 333){

                      const reco::Candidate* phi = lineage.back()-> daughter(i);
                      
                      hPhiEta -> Fill(phi->eta());
                      hPhiPt -> Fill(phi->pt());
                      hPhiEta_PhiPt ->Fill(phi->eta() , lineage.back()-> daughter(i)->pt());

                      if( phi -> numberOfDaughters() == 2 && abs(phi->daughter(0) ->pdgId()) == 13 && abs(phi->daughter(1) ->pdgId()) == 13){
                        hPhiToMuMu -> Fill(1.);

                      }

                    }
                    if(lineage.back()-> daughter(i) ->pdgId()  == 22){
                      hGammaEta -> Fill(lineage.back()-> daughter(i)->eta());
                      hGammaPt -> Fill(lineage.back()-> daughter(i)->pt());
                      hGEta_GPt ->Fill(lineage.back()-> daughter(i)->eta() , lineage.back()-> daughter(i)->pt());
                    
                      genPhotons.push_back( lineage.back()-> daughter(i));
                    }
                    
                  }
                  
                  hBsToPhiG -> Fill ( 1. );

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



void ControlChannelPhi::beginJob()
{
  //create a histogram
  hBPt = new TH1D("hBPt","Transverse momentum of B mesons ; p_{T} [GeV]; Counts",1000, 0., 100.);
  hBssPt = new TH1D("hBssPt","Transverse momentum of B_{s}^{0}, #bar{B}_{s}^{0} ; p_{T} [GeV]; Counts",1000, 0., 100.);
  
  hBsToPhiG = new TH1D("hBsToPhiG", "Number of B_{s}^{0}/#bar{B}_{s}^{0} #rightarrow #Phi#gamma decays; Decays per Event; Events", 6, -0.5, 5.5);
  hPhiToMuMu = new TH1D("hPhiToMuMu", "Number of #Phi #rightarrow #mu#mu decays; Decays per Event; Events", 6, -0.5, 5.5);
  hBsToB =  new TH1D("hBsToB","Ratio of B_{s} to all B mesons ; Ratio; Events",120, -0.15, 1.05); 
  
  hNinit = new TH1D("hNinit","Number of b quarks that hadronized into B measons; per event; Events",10, -1.5, 8.5);
  hInitial =  new TH1D("hInitial","First particle formed by b quarks; Particle ID; Events",95, 505.5, 600.5);
  hHad =  new TH1D("hHad"," Hadronization of b quarks; Particle ID; Events",95, 505.5, 600.5); 
  hBsParents =  new TH1D("hBsParents","Parents of B_{s}^{0} (#Phi#gamma); Particle ID; Events",45, 505.5, 550.5);
  hBsAncestor =  new TH1D("hBsAncestor","First produced B meson that decayed into B_{s}^{0} #rightarrow #Phi#gamma; Particle ID; Events",45, 505.5, 550.5); 
  //hBsProduct =  new TH1D("hBsProduct","Production of B_{s}^{0} ; Decay Channel; Events",7, -1.5, 5.5); 

  hCascade = new TH1D("hCascade","Muon cascade; Counts per event; Events",12, -1.5 , 10.5);
  hMuMuDeltaR = new TH1D("hMuMuDeltaR","Angular separation between two muons; #Delta R; Events",1000, 0., 20.);
  hMuMuDeltaEta = new TH1D("hMuMuDeltaEta","Angular separation between two muons - #Delta #eta; #Delta #eta; Events",1000, 0., 20.);
  hMuMuDeltaPhi = new TH1D("hMuMuDeltaPhi","Angular separation between two muons - #Delta #phi; #Delta #phi; Events",1000, 0., 20.);
  hMuEta_MuEta = new TH2D("hMuEta_MuEta","Pseudorapidity of #gamma; #eta_{1}; #eta_{2}; Events",1000, -10., 10., 1000, -10., 10.);
  hMuPtAll = new TH1D("hMuPtAll","Transverse momentum of muons; p_{T} [GeV]; Events",1000, 0., 100.);
  hMuPt = new TH1D("hMuPt","Transverse momentum of muons (#mu^{+}#mu^{-}); p_{T} [GeV]; Events",1000, 0., 100.);


  //Products of the Bs decay
  hPhiPt = new TH1D("hPhiPt","Transverse momentum of #Phi; p_{T} [GeV]; Events",1000, 0., 100.);
  hPhiEta = new TH1D("hPhiEta","Pseudorapidity of #Phi; #eta; Events",1000, -10., 10.);
  hGEta_GPt = new TH2D("hGEta_GPt","Pseudorapidity and transverse momentum of #gamma; #eta_{#gamma}; p_{T}_{#gamma} [GeV]; Events",200, -10., 10., 200, 0., 50.);
  hPhiEta_PhiPt = new TH2D("hPhiEta_PhiPt","Pseudorapidity and transverse momentum of #Phi; #eta_{#Phi}; p_{T}_{#Phi} [GeV]; Events",200, -10., 10., 200, 0., 50.);
  hGammaPt = new TH1D("hGammaPt","Transverse momentum of #gamma; p_{T} [GeV]; Events",1000, 0., 100.);
  hGammaEta = new TH1D("hGammaEta","Pseudorapidity of #gamma; #eta; Events",1000, -10., 10.);
  
  hGGen_GReco_DeltaR = new TH1D("hGGen_GReco_DeltaR","Angular separation between reconstructed and generated #gamma; #Delta R; Events",1000, 0., 20.);
  hRecoGammaPt = new TH1D("hRecoGammaPt", "Reconstructed photons from B_{s}^{0}/#bar{B}_{s}^{0} #rightarrow #Phi#gamma; pT [GeV]; counts", 300, 0, 30);
  
  // products of phi decay
  hGen_MuPt = new TH1D("hRecoMu_Pt", "reco muon pT; pT [GeV]; counts", 300, 0, 30);
  hGen_MuEta_Pt = new TH2D("hGenMu_Eta_Pt", "gen muon pT vs eta; p_{T} [GeV]; #eta", 300, 0, 30, 1000, -10, 10);

  hMuGen_MuReco_DeltaR = new TH1D("hMuGen_MuReco_DeltaR","Reconstructed muons from #Phi #rightarrow #mu #mu; #Delta R; Events",1000, 0., 20.);
  hRecoMu_Pt = new TH1D("hRecoMu_Pt", "Reconstructed muons from #Phi #rightarrow #mu #mu; pT [GeV]; counts", 300, 0, 30);
  hRecoMu_Eta_Pt = new TH2D("hRecoMu_Eta_Pt", "Reconstructed muons from #Phi #rightarrow #mu #mu; p_{T} [GeV]; #eta", 300, 0, 30, 1000, -10, 10);

  cout << "HERE ControlChannelPhi::beginJob()" << endl;
}

void ControlChannelPhi::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hBPt-> Write();
  hBssPt-> Write();
  
  hBsToPhiG-> Write();
  hPhiToMuMu-> Write();
  hBsToB-> Write();
  
  hNinit ->Write();
  hInitial -> Write();
  hHad-> Write();
  hBsParents -> Write();
  hBsAncestor -> Write();
  //hBsProduct -> Write();

  hCascade -> Write();
  hMuMuDeltaR-> Write();
  hMuMuDeltaEta-> Write();
  hMuMuDeltaPhi-> Write();
  hMuEta_MuEta -> Write();
  hMuPt -> Write();
  hMuPtAll -> Write();

  hPhiPt-> Write();
  hPhiEta-> Write();
  hGEta_GPt-> Write();
  hPhiEta_PhiPt-> Write();
  hGammaPt-> Write();
  hGammaEta-> Write();

  hRecoGammaPt ->Write();
  hGGen_GReco_DeltaR -> Write();

  //phi decay
  hGen_MuPt -> Write();
  hGen_MuEta_Pt -> Write();

  hMuGen_MuReco_DeltaR -> Write();

  hRecoMu_Pt -> Write();
  hRecoMu_Eta_Pt -> Write();

  myRootFile.Close();

  delete hBPt;
  delete hBssPt;

  delete hBsToPhiG;
  delete hPhiToMuMu;
  delete hBsToB;
  
  delete hInitial;
  delete hHad;
  delete hBsParents ;
  delete hBsAncestor ;
  
  //delete hBsProduct
  delete hNinit;
  delete hCascade;
  delete hMuMuDeltaR;
  delete hMuMuDeltaPhi;
  delete hMuMuDeltaEta;
  delete hMuEta_MuEta ;
  delete hMuPt ;
  delete hMuPtAll ;

  delete hPhiPt;
  delete hPhiEta;
  delete hPhiEta_PhiPt;
  delete hGEta_GPt;
  delete hGEta_GPt;
  delete hGammaPt;
  delete hGammaEta;

  delete hGGen_GReco_DeltaR;
  delete hRecoGammaPt;

  //phi decay
  delete hGen_MuPt;
  delete hGen_MuEta_Pt;

  delete hMuGen_MuReco_DeltaR;

  delete hRecoMu_Pt;
  delete hRecoMu_Eta_Pt;

  cout << "HERE Mgr::endJob()" << endl;
}


void ControlChannelPhi::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE ControlChannelPhi::analyze "<< std::endl;

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);
  
  /*vector<const reco::Candidate*> genMuons;
  vector<const reco::Muon*> recoMatchedMuons;

  vector<const reco::Candidate*> genPhotons;
  vector<const reco::Photon*> recoMatchedPhotons;
  */
  std::vector<const reco::Candidate*> cascadeMuons;
  int nBs = 0;

  std::vector<std::vector<const reco::Candidate*>> tree = bFamilyTree(genPar);

  std::cout << "gen photons" << genPhotons.size() << std::endl;

  std::cout <<"Family tree: " << std::endl;
  for (const auto& lineage : tree) { 
    if( abs(lineage.back()->pdgId()) == 531 ) nBs++;
    for (const auto* particle : lineage) { 
        std::cout << particle->pdgId() << " "; 
        if(lineage == tree.at(0)){
          //Phi Gamma decay
        }
    }
    std::cout << std::endl; 
  }

  hBsToB -> Fill (nBs*1./tree.size());

  hBsAncestor -> Fill(tree.at(0).at(0)->pdgId() );
  if (tree.at(0).size() > 1){
    hBsParents -> Fill ( tree.at(0).at( tree.at(0).size() -2 ) ->pdgId() );
  }else{
    hBsParents -> Fill( 506.);
  }

  hNinit -> Fill(tree.size());

  std::cout << "Search for muons: " << std::endl;
  if(tree.size() == 2){ //only two B mesons

    hHad ->Fill( abs(tree.at(1).at(0)->pdgId()) );
    
    //for(const auto* particle : tree.at(1)){
    const reco::Candidate* particle = tree.at(1).at(0);
    std::cout << std::endl;

    std::cout <<"Tree component:" << particle->pdgId() << std::endl;
    std::stack<const reco::Candidate*> stack;
    stack.push(particle); //put particle onto the stack to be checked 
      
    while (!stack.empty()) {
      const reco::Candidate* current = stack.top();
      stack.pop();
        
        
      if (abs(current->pdgId()) == 13) {
          cascadeMuons.push_back(current);
          hMuPtAll -> Fill(current->pt());
          cout << "Muon pt (currently checked): " << current->pt() << endl;
      }

      std::cout << "Daughters of " <<current->pdgId() << ": " << std::endl;
      for (size_t i = 0; i < current->numberOfDaughters(); ++i) {

          stack.push(current->daughter(i));
          std::cout << current->daughter(i)->pdgId()<< "  " ;
      }
      std::cout << std::endl;
    }
    
    
  }

  hCascade -> Fill(cascadeMuons.size());
  if(cascadeMuons.size() == 2){
    const reco::Candidate* gamma1 = cascadeMuons[0];
    const reco::Candidate* gamma2 = cascadeMuons[1];
    if(gamma1->charge() *gamma2->charge() <0 ){
      hMuMuDeltaR -> Fill(std::sqrt( reco::deltaR2(gamma1->eta(), gamma1->phi(), gamma2->eta(), gamma2->phi()) ));
      hMuEta_MuEta -> Fill(gamma1->eta(), gamma2->eta());
      hMuPt -> Fill(gamma1->pt());
      hMuPt -> Fill(gamma2->pt());
    }
  }


  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(ControlChannelPhi);

