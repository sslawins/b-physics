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
//#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

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
  void findPhotonFromBs(const std::vector<reco::Photon>&);
  void findMuonsFromBs(const std::vector<reco::Muon>&);
  double energy     ( const reco::Candidate* , double );
  double invariantMass  ( const reco::Candidate* , const reco::Candidate* , const reco::Candidate* = nullptr );

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  //Bs and phi decay
  vector<const reco::Candidate*> genMuons;
  vector<const reco::Muon*> recoMatchedMuons;

  //vector<const reco::Candidate*> genMuonsOtherSide;
  vector<const reco::Muon*> recoMatchedMuonsOtherSide;

  vector<const reco::Candidate*> genPhotons;
  vector<const reco::Photon*> recoMatchedPhotons;

  //B mesons momentum
  TH1D *hBPt;
  TH1D *hBssPt;

  //check the number of decays
  TH1D *hBsToPhiG;
  TH1D *hPhiToMuMu;
  TH1D *hPhiDecayProducts;
  TH1D *hBsMass;
  TH1D *hPhiMass;

  //check how often a Bs is formed
  TH1D *hBsToB;

  //history of B meson decaying into phi gamma (5)
  TH1D *hNinit;
  TH1D *hInitial;
  TH1D *hHad; 
  TH1D *hBsParents ;
  TH1D *hBsAncestor ; 

  // muons from the other side - expected to give a trigger
  TH1D *hCascade;

  TH2D *hCascadeGen_Eta_Pt ;
  TH2D *hCascadeGen_Eta_Eta ;
  TH1D *hCascadeGen_PtAll ;
  TH1D *hCascadeGen_Pt ;

  TH1D *hCascadeMuMuDeltaR;
  TH1D *hCascadeMuMuDeltaPhi;
  TH1D *hCascadeMuMuDeltaEta;

  TH1D *hCascadeReco_Pt ;
  TH2D *hCascadeReco_Eta_Pt;

  TH1D *hCascade_RecoGenDeltaR;
  TH2D *hCascade_RecoPt_GenPt;
  
  //Products of the Bs decay
  TH1D *hPhiPt;
  TH1D *hPhiEta;
  TH2D *hPhi_Eta_Pt;

  TH2D *hDecayGamma_Eta_Pt;
  TH1D *hDecayGamma_Pt;
  TH1D *hDecayGamma_Eta;

  TH1D *hDecayGammaReco_Pt;
  TH2D *hDecayGammaReco_Eta_Pt;

  TH1D *hDecayGamma_RecoGenDeltaR;
  TH2D *hDecayGamma_RecoPt_GenPt;

  //with additional cuts
  TH1D *hCutsDecayGammaReco_Pt;
  TH2D *hCutsDecayGammaReco_Eta_Pt;
  TH1D *hCutsDecayGamma_RecoGenDeltaR;
  TH2D *hCutsDecayGamma_RecoPt_GenPt;

  TH1D *hCutsDecayGammaGen_Pt;
  TH2D *hCutsDecayGammaGen_Eta_Pt;
  TH1D *hCutsDecayGammaGenMatched_Pt;
  TH2D *hCutsDecayGammaGenMatched_Eta_Pt;

  //products of phi decay -> two muons
  TH1D* hPhiDecayGen_Pt;
  TH2D* hPhiDecayGen_Eta_Pt;
  
  TH1D* hPhiDecayReco_Pt;
  TH2D* hPhiDecayReco_Eta_Pt;

  TH1D* hPhiDecay_RecoGenDeltaR;
  TH2D *hPhiDecay_RecoPt_GenPt;

  //with additional cuts
  TH1D* hCutsPhiDecayReco_Pt;
  TH2D* hCutsPhiDecayReco_Eta_Pt;
  TH1D* hCutsPhiDecay_RecoGenDeltaR;
  TH2D *hCutsPhiDecay_RecoPt_GenPt;

  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;
  edm::EDGetTokenT < vector<reco::Muon> > theMuonToken;
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;

  std::vector<int> MuMuG = {22, 13, -13};
  std::vector<int> BsStarG= {22, 533};
  std::vector<int> Bs= {531};
  std::vector<int> PhiG= {333, 22};

  double mMu =  0.1056583755; // GeV

};


ControlChannelPhi::ControlChannelPhi(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles" ));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));
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
                //std::cout << "Pdg ID " << abs(part.pdgId()) <<std::endl;
                //if( abs(lineage.back()->pdgId()) == 531 ) hBssPt -> Fill( lineage.back()->pt() ) ;  
                //std::cout << "Pdg ID " << abs(lineage.back()->pdgId()) <<std::endl;

                if( isSameChannel(decay, PhiG)  ) {
                  hBssPt -> Fill (lineage.back()->pt());
                  for( size_t i = 0; i < lineage.back()->numberOfDaughters(); ++i ){
                    if(lineage.back()-> daughter(i) ->pdgId()  == 333){

                      const reco::Candidate* phi = lineage.back()-> daughter(i);
                      
                      hPhiEta -> Fill(phi->eta());
                      hPhiPt -> Fill(phi->pt());
                      hPhi_Eta_Pt ->Fill(phi->eta() , lineage.back()-> daughter(i)->pt());
                      cout << " Phi decay, number of daughters: " <<  phi -> numberOfDaughters() << endl;
                      if( phi -> numberOfDaughters() == 2 && abs(phi->daughter(0) ->pdgId()) == 13 && abs(phi->daughter(1) ->pdgId()) == 13){
                        hPhiToMuMu -> Fill(1.);
                        hPhiDecayGen_Pt -> Fill(phi->daughter(0)->pt());
                        hPhiDecayGen_Pt -> Fill(phi->daughter(1)->pt());
                        hPhiDecayGen_Eta_Pt -> Fill(phi->daughter(0)->eta(), phi->daughter(0)->pt());
                        hPhiDecayGen_Eta_Pt -> Fill(phi->daughter(1)->eta(), phi->daughter(1)->pt());
                        genMuons.push_back(phi->daughter(0));
                        genMuons.push_back(phi->daughter(1));
                      }else{
                        //cout << "Phi daughters: " ;
                        for( size_t dauPhi = 0; dauPhi < phi -> numberOfDaughters(); ++dauPhi ){
                          //cout << phi -> daughter(dauPhi)->pdgId() << " " ;
                          hPhiDecayProducts -> Fill( abs(phi -> daughter(dauPhi)->pdgId()) );
                        }
                        //cout << endl;
                      }

                    }
                    if(lineage.back()-> daughter(i) ->pdgId()  == 22){
                      hDecayGamma_Eta -> Fill(lineage.back()-> daughter(i)->eta());
                      hDecayGamma_Pt -> Fill(lineage.back()-> daughter(i)->pt());
                      hDecayGamma_Eta_Pt ->Fill(lineage.back()-> daughter(i)->eta() , lineage.back()-> daughter(i)->pt());
                    
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

void ControlChannelPhi::findMuonsFromBs(const std::vector<reco::Muon>& recoGamma) {
  
}
void ControlChannelPhi::findPhotonFromBs(const std::vector<reco::Photon>& recoMuon) {
}

double ControlChannelPhi::energy(const reco::Candidate* par, double mass){
  double energy = sqrt( pow(par->p() , 2.) + mass * mass);
  return energy;
}

double ControlChannelPhi::invariantMass(const reco::Candidate* c1, const reco::Candidate* c2, const reco::Candidate* c3) {
    TLorentzVector p4_1(c1->px(), c1->py(), c1->pz(), energy(c1, mMu));
    TLorentzVector p4_2(c2->px(), c2->py(), c2->pz(), energy(c2, mMu));

    TLorentzVector sum = p4_1 + p4_2;
    cout << "c1 and c2 ok" << endl;
    if (c3 != nullptr) {
      cout << "c3 " ;
        TLorentzVector p4_3(c3->px(), c3->py(), c3->pz(), energy(c3, 0.));
        cout << "ok" << endl;
        sum += p4_3;
    }
    return sum.M();
}

void ControlChannelPhi::beginJob()
{
  //B mesons momentum
  hBPt = new TH1D("hBPt","Transverse momentum of B mesons ; p_{T} [GeV]; Counts",1000, 0., 100.);
  hBssPt = new TH1D("hBssPt","Transverse momentum of B_{s}^{0}, #bar{B}_{s}^{0} ; p_{T} [GeV]; Counts",1000, 0., 100.);
  
  //check the number of decays
  hBsToPhiG = new TH1D("hBsToPhiG", "Number of B_{s}^{0}/#bar{B}_{s}^{0} #rightarrow #Phi#gamma decays; Decays per Event; Events", 6, -0.5, 5.5);
  hPhiToMuMu = new TH1D("hPhiToMuMu", "Number of #Phi #rightarrow #mu#mu decays; Decays per Event; Events", 6, -0.5, 5.5);
  hPhiDecayProducts = new TH1D("hPhiDecayProducts", "Products of #Phi decays; Particle ID; Events",350 , 0.5, 350.5);
  hBsMass = new TH1D("hBsMass", "Reconstruction of B_{s}^{0}; M_{inv} [GeV]; Counts",48 , 5., 6.2);
  hPhiMass = new TH1D("hPhiMass", "Reconstruction of #Phi ; M_{inv} [GeV]; Counts", 48 , 0.9, 1.2);

  //reconstruction

  //check how often a Bs is formed
  hBsToB =  new TH1D("hBsToB","Ratio of B_{s} to all B mesons ; Ratio; Events",120, -0.15, 1.05); 
  
  //history of B meson decaying into phi gamma (5)
  hNinit = new TH1D("hNinit","Number of b quarks that hadronized into B measons; per event; Events",10, -1.5, 8.5);
  hInitial =  new TH1D("hInitial","First particle formed by b quarks; Particle ID; Events",95, 505.5, 600.5);
  hHad =  new TH1D("hHad"," Hadronization of b quarks; Particle ID; Events",95, 505.5, 600.5); //biased
  hBsParents =  new TH1D("hBsParents","Parents of B_{s}^{0} (#Phi#gamma); Particle ID; Events",45, 505.5, 550.5);
  hBsAncestor =  new TH1D("hBsAncestor","First produced B meson (B_{s}^{0} #rightarrow #Phi#gamma); Particle ID; Events",45, 505.5, 550.5); 
  
  // muons from the other side - expected to give a trigger (11)
  hCascade = new TH1D("hCascade","Muon cascade; Counts per event; Events",12, -1.5 , 10.5);
  //gen part
  hCascadeGen_Eta_Pt = new TH2D("hCascadeGen_Eta_Pt", "Muons ; #eta ; p_{T} [GeV]", 1000, -8, 8, 500, 0, 30);
  hCascadeGen_Eta_Eta = new TH2D("hCascadeGen_Eta_Eta","Muons ; #eta_{1}; #eta_{2}; Events",1000, -8., 8., 1000, -8., 8.);
  hCascadeGen_PtAll = new TH1D("hCascadeGen_PtAll","Muons ; p_{T} [GeV]; Events",1500, 0., 15.);
  hCascadeGen_Pt = new TH1D("hCascadeGen_Pt","Muons - only #mu^{+}#mu^{-} pairs ; p_{T} [GeV]; Events",1000, 0., 10.);
  //to explain the shape of delta R
  hCascadeMuMuDeltaR = new TH1D("hCascadeMuMuDeltaR","Angular separation between two muons (#mu^{+}#mu^{-}); #Delta R; Events",1000, 0., 20.);
  hCascadeMuMuDeltaEta = new TH1D("hCascadeMuMuDeltaEta","Angular separation between two muons (#mu^{+}#mu^{-}) - #Delta #eta; #Delta #eta; Events",1000, 0. , 10.);
  hCascadeMuMuDeltaPhi = new TH1D("hCascadeMuMuDeltaPhi","Angular separation between two muons (#mu^{+}#mu^{-}) - #Delta #phi; #Delta #phi; Events",1000, 0., 6.5);
  //reco part
  hCascadeReco_Pt = new TH1D("hCascadeReco_Pt","Reconstructed muons (#mu^{+}#mu^{-}) ; p_{T} [GeV]; Events",400, 0., 20.);
  hCascadeReco_Eta_Pt = new TH2D("hCascadeReco_Eta_Pt", "Reconstructed muons (#mu^{+}#mu^{-}); #eta ; p_{T} [GeV]", 600, -3, 3, 200, 0, 20);
  //reco vs gen - to adjust the cut
  hCascade_RecoGenDeltaR = new TH1D("hCascade_RecoGenDeltaR","Reconstructed muons (#mu^{+}#mu^{-}) ;#Delta R; Events",500, 0., 2.);
  hCascade_RecoPt_GenPt = new TH2D("hCascade_RecoPt_GenPt"," Reconstructed vs generated muons (#mu^{+}#mu^{-}); p_{T}_{reco} [GeV]; p_{T}_{gen} [GeV]; Counts",500, 0., 30., 500, 0., 30.);

  //Products of the Bs decay (9)
  //gen phi
  hPhiPt = new TH1D("hPhiPt","Transverse momentum of #Phi (B_{s}^{0} #rightarrow #Phi#gamma) ; p_{T} [GeV]; Events",1000, 0., 100.);
  hPhiEta = new TH1D("hPhiEta","Pseudorapidity of #Phi (B_{s}^{0} #rightarrow #Phi#gamma); #eta; Events",200, -10., 10.);
  hPhi_Eta_Pt = new TH2D("hPhi_Eta_Pt","Pseudorapidity and transverse momentum of #Phi (B_{s}^{0} #rightarrow #Phi#gamma); #eta_{#Phi}; p_{T}_{#Phi} [GeV]; Events",200, -10., 10., 200, 0., 50.);
  //gen gamma
  hDecayGamma_Eta_Pt = new TH2D("hDecayGamma_Eta_Pt","Pseudorapidity and transverse momentum of photons (B_{s}^{0} #rightarrow #Phi#gamma); #eta_{#gamma}; p_{T}_{#gamma} [GeV]; Events",1000, -10., 10., 500, 0., 30.);
  hDecayGamma_Pt = new TH1D("hDecayGamma_Pt","Transverse momentum of photons (B_{s}^{0} #rightarrow #Phi#gamma); p_{T} [GeV]; Events", 500., 0., 30.);
  hDecayGamma_Eta = new TH1D("hDecayGamma_Eta","Pseudorapidity of photons (B_{s}^{0} #rightarrow #Phi#gamma); #eta; Events",1000, -10., 10.);
  //reco gamma
  hDecayGammaReco_Pt = new TH1D("hDecayGammaReco_Pt", "Reconstructed photons (B_{s}^{0} #rightarrow #Phi#gamma); pT [GeV]; Counts", 500, 0, 30);
  hDecayGammaReco_Eta_Pt = new TH2D("hDecayGammaReco_Eta_Pt","Pseudorapidity and transverse momentum of photons (B_{s}^{0} #rightarrow #Phi#gamma); #eta_{#gamma}; p_{T}_{#gamma} [GeV]; Events",700, -3.5, 3.5, 500, 0., 50.);
  //comparison of reco and gen gamma - to ajust the cut
  hDecayGamma_RecoGenDeltaR = new TH1D("hDecayGamma_RecoGenDeltaR","Angular separation between reconstructed and generated photons (B_{s}^{0} #rightarrow #Phi#gamma); #Delta R; Events",1000, 0., 0.1);
  hDecayGamma_RecoPt_GenPt = new TH2D("hDecayGamma_RecoPt_GenPt"," Reconstructed vs generated photons (B_{s}^{0} #rightarrow #Phi#gamma); p_{T,reco} [GeV]; p_{T,gen} [GeV]; Counts",500, 0., 50., 500, 0., 50.);

  //with additional cuts
  hCutsDecayGammaReco_Pt = new TH1D("hCutsDecayGammaReco_Pt", "Reconstructed photons, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV ; pT [GeV]; counts", 500, 0, 30);
  hCutsDecayGammaReco_Eta_Pt = new TH2D("hCutsDecayGammaReco_Eta_Pt","Reconstructed photons, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV; #eta_{#gamma}; p_{T}_{#gamma} [GeV]; Events",700, -3.5, 3.5, 500, 0., 50.);
  hCutsDecayGamma_RecoGenDeltaR = new TH1D("hCutsDecayGamma_RecoGenDeltaR","Photons, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV; #Delta R; Events",100, 0., 0.1);
  hCutsDecayGamma_RecoPt_GenPt = new TH2D("hCutsDecayGamma_RecoPt_GenPt"," Photons, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV; p_{T,reco} [GeV]; p_{T,gen} [GeV]; Counts",500, 0., 30., 500, 0., 30.);

  hCutsDecayGammaGen_Eta_Pt = new TH2D("hCutsDecayGammaGen_Eta_Pt","Generated photons, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV; #eta_{#gamma}; p_{T}_{#gamma} [GeV]; Events",700, -3.5, 3.5, 500, 0., 50.);
  hCutsDecayGammaGen_Pt = new TH1D("hCutsDecayGammaGen_Pt", "Generated photons, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV ; pT [GeV]; counts", 500, 0, 30);
  hCutsDecayGammaGenMatched_Eta_Pt = new TH2D("hCutsDecayGammaGenMatched_Eta_Pt","Generated and matched photons, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV; #eta_{#gamma}; p_{T}_{#gamma} [GeV]; Events",700, -3.5, 3.5, 500, 0., 50.);
  hCutsDecayGammaGenMatched_Pt = new TH1D("hCutsDecayGammaGenMatched_Pt", "Generated and matched photons, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV ; pT [GeV]; counts", 500, 0, 30);
  
  // products of phi decay -> two muons (5)
  //gen part
  hPhiDecayGen_Pt = new TH1D("hPhiDecayGen_Pt", "Muons (#Phi #rightarrow #mu#mu); p_{T} [GeV]; counts", 500, 0, 30);
  hPhiDecayGen_Eta_Pt = new TH2D("hPhiDecayGen_Eta_Pt", "Muons (#Phi #rightarrow #mu#mu); #eta ; p_{T} [GeV]", 1000, -10, 10, 300, 0, 30 );
  //reco part
  hPhiDecayReco_Pt= new TH1D("hPhiDecayReco_Pt", "Reconstructed muons (#Phi #rightarrow #mu#mu); p_{T} [GeV]; counts", 500, 0, 30);
  hPhiDecayReco_Eta_Pt = new TH2D("hPhiDecayReco_Eta_Pt", "Reconstructed muons (#Phi #rightarrow #mu#mu); #eta ; p_{T} [GeV]",600, -3., 3., 300, 0, 30 );
  //comparison - to adjust the cut
  hPhiDecay_RecoGenDeltaR= new TH1D("hPhiDecay_RecoGenDeltaR","Angular separation between reconstructed and generated muons (#Phi #rightarrow #mu#mu); #Delta R; Events",1000, 0., 0.1);
  hPhiDecay_RecoPt_GenPt = new TH2D("hPhiDecay_RecoPt_GenPt"," Reconstructed vs generated muons (#Phi#rightarrow #mu#mu); p_{T,reco} [GeV]; p_{T,gen} [GeV]; Counts",200, 0., 50., 200, 0., 50.);

  //with additional cuts
  hCutsPhiDecayReco_Pt = new TH1D("hCutsPhiDecayReco_Pt", "Reconstructed muons (#Phi #rightarrow #mu#mu, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV ); p_{T} [GeV]; Counts", 500, 0, 30);
  hCutsPhiDecayReco_Eta_Pt = new TH2D("hCutsPhiDecayReco_Eta_Pt", "Reconstructed muons (#Phi #rightarrow #mu#mu, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV); #eta ; p_{T} [GeV]",600, -3., 3., 300, 0, 30);
  hCutsPhiDecay_RecoGenDeltaR = new TH1D("hCutsPhiDecay_RecoGenDeltaR","Muons (#Phi #rightarrow #mu#mu, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV); #Delta R; Events",5000, 0., 5.);
  hCutsPhiDecay_RecoPt_GenPt = new TH2D("hCutsPhiDecay_RecoPt_GenPt"," Muons (#Phi #rightarrow #mu#mu, |#eta_{#mu}| < 2.4, p_{T,#mu} > 4 GeV); p_{T,reco} [GeV]; p_{T,gen} [GeV]; Counts",500, 0., 30., 500, 0., 30.);

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
  hPhiDecayProducts -> Write();
  hBsMass -> Write();
  hPhiMass -> Write();

  hBsToB-> Write();
  
  hNinit ->Write();
  hInitial -> Write();
  hHad-> Write();
  hBsParents -> Write();
  hBsAncestor -> Write();

  hCascade -> Write();

  hCascadeGen_Eta_Pt -> Write();
  hCascadeGen_Eta_Eta -> Write();
  hCascadeGen_PtAll -> Write();
  hCascadeGen_Pt -> Write();

  hCascadeMuMuDeltaR-> Write();
  hCascadeMuMuDeltaEta-> Write();
  hCascadeMuMuDeltaPhi-> Write();

  hCascadeReco_Pt -> Write();
  hCascadeReco_Eta_Pt -> Write();

  hCascade_RecoGenDeltaR -> Write();
  hCascade_RecoPt_GenPt -> Write();

  hPhiPt-> Write();
  hPhiEta-> Write();
  hPhi_Eta_Pt-> Write();

  hDecayGamma_Eta_Pt-> Write();
  hDecayGamma_Pt-> Write();
  hDecayGamma_Eta-> Write();

  hDecayGammaReco_Pt ->Write();
  hDecayGammaReco_Eta_Pt ->Write();

  hDecayGamma_RecoGenDeltaR -> Write();
  hDecayGamma_RecoPt_GenPt -> Write();

  //with additional cuts
  hCutsDecayGammaReco_Pt ->Write();
  hCutsDecayGammaReco_Eta_Pt ->Write();
  hCutsDecayGamma_RecoGenDeltaR -> Write();
  hCutsDecayGamma_RecoPt_GenPt -> Write();

  hCutsDecayGammaGen_Pt ->Write();
  hCutsDecayGammaGen_Eta_Pt ->Write();
  hCutsDecayGammaGenMatched_Pt ->Write();
  hCutsDecayGammaGenMatched_Eta_Pt ->Write();

  hPhiDecayGen_Pt -> Write();
  hPhiDecayGen_Eta_Pt -> Write();

  hPhiDecayReco_Pt-> Write();
  hPhiDecayReco_Eta_Pt -> Write();
  
  hPhiDecay_RecoGenDeltaR-> Write();
  hPhiDecay_RecoPt_GenPt -> Write();

  //with additional cuts 
  hCutsPhiDecayReco_Pt-> Write();
  hCutsPhiDecayReco_Eta_Pt -> Write();
  hCutsPhiDecay_RecoGenDeltaR-> Write();
  hCutsPhiDecay_RecoPt_GenPt -> Write();

  myRootFile.Close();

  delete hBPt;
  delete hBssPt;

  delete hBsToPhiG;
  delete hPhiToMuMu;
  delete hPhiDecayProducts;
  delete hBsMass;
  delete hPhiMass;

  delete hBsToB;
  
  //delete hBsProduct
  delete hNinit;
  delete hInitial;
  delete hHad;
  delete hBsParents;
  delete hBsAncestor;

  delete hCascade;

  delete hCascadeGen_Eta_Pt ;
  delete hCascadeGen_Eta_Eta ;
  delete hCascadeGen_PtAll ;
  delete hCascadeGen_Pt ;
  
  delete hCascadeMuMuDeltaR;
  delete hCascadeMuMuDeltaPhi;
  delete hCascadeMuMuDeltaEta;

  delete hCascadeReco_Pt ;
  delete hCascadeReco_Eta_Pt ;

  delete hCascade_RecoGenDeltaR;
  delete hCascade_RecoPt_GenPt;

  delete hPhiPt;
  delete hPhiEta;
  delete hPhi_Eta_Pt;

  delete hDecayGamma_Eta_Pt;
  delete hDecayGamma_Pt;
  delete hDecayGamma_Eta;

  delete hDecayGammaReco_Pt;
  delete hDecayGammaReco_Eta_Pt;

  delete hDecayGamma_RecoGenDeltaR;
  delete hDecayGamma_RecoPt_GenPt;

  //with additional cuts
  delete hCutsDecayGammaReco_Pt;
  delete hCutsDecayGammaReco_Eta_Pt;
  delete hCutsDecayGamma_RecoGenDeltaR;
  delete hCutsDecayGamma_RecoPt_GenPt;

  delete hCutsDecayGammaGen_Pt ;
  delete hCutsDecayGammaGen_Eta_Pt ;
  delete hCutsDecayGammaGenMatched_Pt ;
  delete hCutsDecayGammaGenMatched_Eta_Pt;

  //phi decay
  delete hPhiDecayGen_Pt;
  delete hPhiDecayGen_Eta_Pt;

  delete hPhiDecayReco_Pt;
  delete hPhiDecayReco_Eta_Pt;

  delete hPhiDecay_RecoGenDeltaR;
  delete hPhiDecay_RecoPt_GenPt;

  //with additional cuts
  delete hCutsPhiDecayReco_Pt;
  delete hCutsPhiDecayReco_Eta_Pt;
  delete hCutsPhiDecay_RecoGenDeltaR;
  delete hCutsPhiDecay_RecoPt_GenPt;

  cout << "HERE ControlChannelPhi::endJob()" << endl;
}


void ControlChannelPhi::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE ControlChannelPhi::analyze "<< std::endl;

  genMuons.clear();
  recoMatchedMuons.clear();
  //genMuonsOtherSide.clear(); ----> this is exactly what cascadeMuons is
  recoMatchedMuonsOtherSide.clear();
  genPhotons.clear();
  recoMatchedPhotons.clear();

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);

  std::vector<const reco::Candidate*> cascadeMuons;
  int nBs = 0;

  //find all the B mesons and put them in order from mother to daughter 
  std::vector<std::vector<const reco::Candidate*>> tree = bFamilyTree(genPar);

  // print the family tree for future reference
  std::cout <<"Family tree: " << std::endl;
  for (const auto& lineage : tree) { 
    if( abs(lineage.back()->pdgId()) == 531 ) nBs++;
    for (const auto* particle : lineage) { 
        std::cout << particle->pdgId() << " "; 
    }
    std::cout << std::endl; 
  }

  //number of DECAYING (to somthing else) Bs compared to number of b quarks hadronizing into B mesons
  hBsToB -> Fill (nBs*1./tree.size());
  //the first B meson that was produced with the b quark, that would later be a part o Bs0 decaying into phi gamma 
  hBsAncestor -> Fill(tree.at(0).at(0)->pdgId() );

  if (tree.at(0).size() > 1){
    // parent (another B meson) of the Bs0 ( phi gamma ), if there is any
    hBsParents -> Fill ( tree.at(0).at( tree.at(0).size() -2 ) ->pdgId() );
  }else{
    hBsParents -> Fill( 506.); // quarks and gluons, if there is no B parent
  }

  //number of b quarks that hadronized into B mesons
  hNinit -> Fill(tree.size());

  // TAG MUONS (cascade decays) //////////////////////////
  
  //std::cout << "Search for muons: " << std::endl;
  if(tree.size() == 2){ //only two B mesons

    hHad ->Fill( abs(tree.at(1).at(0)->pdgId()) ); // what type of B meson has the other b quark 
                                                   //(if there were only 2) hadronized into?
    const reco::Candidate* particle = tree.at(1).at(0);
    //std::cout << std::endl;
    //std::cout <<"Tree component:" << particle->pdgId() << std::endl;

    // Find all the muons originating from the second b quark (all the decays of all the daughters)
    std::stack<const reco::Candidate*> stack;
    stack.push(particle); 
    
    while (!stack.empty()) {
      const reco::Candidate* current = stack.top();
      stack.pop();
        
      if (abs(current->pdgId()) == 13) {
          cascadeMuons.push_back(current);
          hCascadeGen_PtAll -> Fill(current->pt());
          cout << "Muon pt (currently processed): " << current->pt() << endl;
      }

      //std::cout << "Daughters of " <<current->pdgId() << ": " << std::endl;
      for (size_t i = 0; i < current->numberOfDaughters(); ++i) {

          stack.push(current->daughter(i));
          //std::cout << current->daughter(i)->pdgId()<< "  " ;
      }
      //std::cout << std::endl;
    }
    
  }

  //fill histograms for the cascade decays
  hCascade -> Fill(cascadeMuons.size());

  if(cascadeMuons.size() == 2){
    const reco::Candidate* gamma1 = cascadeMuons[0];
    const reco::Candidate* gamma2 = cascadeMuons[1];
    if(gamma1->charge() *gamma2->charge() <0 ){
      hCascadeMuMuDeltaR -> Fill(std::sqrt( reco::deltaR2(*gamma1, *gamma2) ));
      hCascadeMuMuDeltaPhi -> Fill(abs( gamma1->phi() - gamma2->phi() ));
      hCascadeMuMuDeltaEta -> Fill(abs( gamma1->eta() - gamma2->eta() ));
      
      hCascadeGen_Eta_Eta -> Fill(gamma1->eta(), gamma2->eta());
      hCascadeGen_Eta_Pt -> Fill(gamma1->eta(), gamma1->pt());
      hCascadeGen_Eta_Pt -> Fill(gamma2->eta(), gamma2->pt());
      hCascadeGen_Pt -> Fill(gamma1->pt());
      hCascadeGen_Pt -> Fill(gamma2->pt());
    }
  }

  if (genMuons.size() == 2){
    //MUON MATCHING ////////////////////////////////////////////

    //to store; only for Phi decay for now
    std::vector<const reco::Candidate*> recoMatchedMu(genMuons.size(), nullptr);

    //print some info for debuging
    cout << "All reconstructed muons from this event: " << recoMuons.size() << endl;
    cout << "Generated muons from the decay: " << genMuons.size() << endl;

    // reco muon matching
    for (const auto& recoMu : recoMuons){
      
      double lowestRmuon = 100. ; //for each reco photon new threshold

      //muons from phi decay (there are only two in the genMuons vector, but only for phi -> mumu)
      for (unsigned int genIter = 0; genIter < genMuons.size(); ++genIter ){

        hPhiDecay_RecoGenDeltaR-> Fill(reco::deltaR(recoMu, *genMuons[genIter]));
        
        if (reco::deltaR(recoMu, *genMuons[genIter]) < lowestRmuon){
          
          cout << "Momentum of Matched reco muon  " << recoMu.pt() << endl;
          lowestRmuon = reco::deltaR(recoMu, *genMuons[genIter]);
          recoMatchedMu[genIter] = &recoMu;
          
        }
      }

      //muons from the other side matching, but only 
      if(cascadeMuons.size() == 2 && (cascadeMuons[0]->charge() *cascadeMuons[1]->charge()) <0){
        
        for (const reco::Candidate* genMu : cascadeMuons ){
          hCascade_RecoGenDeltaR -> Fill(reco::deltaR(recoMu, *genMu));
          
          if (reco::deltaR(recoMu, *genMu) < 0.01){
            //cout << "Reconstructed muons pT: " << recoMu.pt() << endl;
            //recoMatchedMuonsOtherSide.push_back(&recoMu);
            hCascadeReco_Eta_Pt -> Fill(recoMu.eta(), recoMu.pt());
            hCascadeReco_Pt -> Fill(recoMu.pt());

            hCascade_RecoPt_GenPt -> Fill (recoMu.pt(), genMu->pt());
            
          }
        }
      }
    }

    // PHOTON MATCHING////////////////////////////////////
    double lowestRphoton = 100. ;
    const reco::Candidate* recoPhoton = nullptr;
    
    // find photon with the lowest delta R with the gen photon
    for (const auto& recoPh : recoPhotons)
    {
      //cout << "Photon" << endl;
      for (const reco::Candidate* genPh : genPhotons) // only one photon
      {
        hDecayGamma_RecoGenDeltaR -> Fill(reco::deltaR(recoPh, *genPh)); // to find a suitable  cut

        if (reco::deltaR(recoPh, *genPh) < lowestRphoton) 
        {
          lowestRphoton = reco::deltaR(recoPh, *genPh);
          recoPhoton = &recoPh;
        }
      }
    }

    
    //recoMatchedMu[0] != nullptr && recoMatchedMu[1] != nullptr &&
    bool reconstructWithCuts = true;
    //reconstruct photon, Bs and Phi with additional cuts only if BOTH muons are reconstructed

    if ( recoMatchedMu[0] != recoMatchedMu[1]){ //the matched reco muon cannot be the same for both gen muons

      for(unsigned int k =0 ; k < recoMatchedMu.size(); ++k){
        
        //muon with the lowest delta R with respect to the gen muon is reconstructed IF AND ONLY IF delta R < 0.01
        //look at the hPhiDecay_RecoGenDeltaR
        if(recoMatchedMu[k] != nullptr && reco::deltaR(*genMuons[k], *recoMatchedMu[k]) < 0.01){
          
          //fill histograms for all recostructed matched muons
          hPhiDecayReco_Eta_Pt -> Fill(recoMatchedMu[k]->eta(), recoMatchedMu[k]->pt());
          hPhiDecayReco_Pt-> Fill(recoMatchedMu[k]->pt());
          hPhiDecay_RecoPt_GenPt -> Fill (recoMatchedMu[k]->pt(), genMuons[k]->pt());
        
          if( recoMatchedMu[k]-> pt() <= 4 ||  abs(recoMatchedMu[k]-> eta()) >= 2.4) {
            reconstructWithCuts = false;
          }

        }else{ reconstructWithCuts = false; }
      }
    
    }else{
      reconstructWithCuts = false;
    }


    // for BOTH muons satisfying the conditions
    if(reconstructWithCuts == true && recoMatchedMu.size() == 2){

      for(unsigned int k =0 ; k < recoMatchedMu.size(); ++k){
        hCutsPhiDecay_RecoGenDeltaR-> Fill(reco::deltaR(*genMuons[k], *recoMatchedMu[k]));
        hCutsPhiDecayReco_Eta_Pt -> Fill(recoMatchedMu[k]->eta(), recoMatchedMu[k]->pt());
        hCutsPhiDecayReco_Pt-> Fill(recoMatchedMu[k]->pt());
        hCutsPhiDecay_RecoPt_GenPt -> Fill (recoMatchedMu[k]->pt(), genMuons[k]->pt());

      }
      
      hCutsDecayGammaGen_Eta_Pt -> Fill( genPhotons[0]->eta(), genPhotons[0]->pt());
      hCutsDecayGammaGen_Pt -> Fill( genPhotons[0]->eta(), genPhotons[0]->pt());
    }

    //PHOTON MATHING //////////
    if (recoPhoton != nullptr && lowestRphoton < 0.03 ) { //reconstructed and matched photons

      //recoMatchedPhotons.push_back(static_cast<const reco::Photon*>(recoPhoton));
      hDecayGammaReco_Pt->Fill(recoPhoton->pt());
      hDecayGammaReco_Eta_Pt->Fill(recoPhoton->eta(), recoPhoton->pt());
      hDecayGamma_RecoPt_GenPt->Fill(recoPhoton->pt(), genPhotons[0]->pt());

      if( reconstructWithCuts == true && recoMatchedMu.size() == 2){

        //histograms with additional cuts
        hCutsDecayGamma_RecoGenDeltaR -> Fill(lowestRphoton);
        hCutsDecayGammaReco_Pt->Fill(recoPhoton->pt());
        hCutsDecayGammaReco_Eta_Pt->Fill(recoPhoton->eta(), recoPhoton->pt());
        hCutsDecayGamma_RecoPt_GenPt->Fill(recoPhoton->pt(), genPhotons[0]->pt());

        //generated photons
        hCutsDecayGammaGenMatched_Eta_Pt -> Fill( genPhotons[0]->eta(), genPhotons[0]->pt());
        hCutsDecayGammaGenMatched_Pt -> Fill( genPhotons[0]->eta(), genPhotons[0]->pt());

        cout << recoMatchedMu.size() << " reco muon size" << endl;
        if( recoMatchedMu[0] == nullptr) cout << "mu 1 nullprt ";
        if( recoMatchedMu[1] == nullptr) cout << "mu 2 nullprt " << endl;

        double bsMass = invariantMass(recoMatchedMu[0], recoMatchedMu[1], recoPhoton);
        double phiMass = invariantMass(recoMatchedMu[0], recoMatchedMu[1]);

        hBsMass -> Fill( bsMass );
        hPhiMass -> Fill( phiMass );
        cout << " invariant mass: " << bsMass << " " << phiMass << endl;
      }
    }
  }

  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(ControlChannelPhi);

