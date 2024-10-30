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


using namespace std;


//object definition
class Analysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit Analysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~Analysis();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool isSameChannel(const std::vector<int>&, const std::vector<int>&);
  //double delta_R(const reco::Candidate*)

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  TH1D *hBsPt;
  TH1D *hBPt;
  TH1D *hBPtFrag;
  TH1D *hBssPtFrag;
  TH1D *hBsPtFrag;
  TH1D *hAntiBsPtFrag;

  TH1D *hBsToMuMuG;
  TH1D *hBsToBd;
  TH1D *hBsToB;

  TH1D *hHad; 
  TH1D *hBsParents ;
  TH1D *hBsAncestor ; 
  TH1D *hBsProduct ; 

  TH1D *hGBsGDeltaR;
  TH2D *hGBsEta_GEta ;
  TH1D *hGammaBsPt ;

  TH1D *hMuPt;
  TH1D *hMuEta;
  TH2D *hMu1Eta_Mu2Eta;
  TH1D *hGammaPt;
  TH1D *hGammaEta;

  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;

  std::vector<int> MuMuG = {22, 13, -13};
  std::vector<int> BsStarG= {22, 533};
  std::vector<int> Bs= {531};

};


Analysis::Analysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles" ));

}

Analysis::~Analysis()
{
  cout <<" DTOR" << endl;
}

bool Analysis::isSameChannel(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}


void Analysis::beginJob()
{
  //create a histogram
  hBPt = new TH1D("hBPt","Transverse momentum of B mesons (multiplecounting); p_{t} [GeV]; Counts",1000, 0., 100.);
  hBsPt = new TH1D("hBsPt","Transverse momentum of B_{s}^{0}, #bar{B}_{s}^{0} (multiple counting) ; p_{T} [GeV]; Counts",1000, 0., 100.);
  hBPtFrag = new TH1D("hBPtFrag","Transverse momentum of B mesons ; p_{T} [GeV]; Counts",1000, 0., 100.);
  hBssPtFrag = new TH1D("hBssPtFrag","Transverse momentum of B_{s}^{0}, #bar{B}_{s}^{0} ; p_{T} [GeV]; Counts",1000, 0., 100.);
  hBsPtFrag = new TH1D("hBsPtFrag","Transverse momentum of B_{s}^{0} ; p_{T} [GeV]; Counts",1000, 0., 100.);
  hAntiBsPtFrag = new TH1D("hAntiBsPtFrag","Transverse momentum of #bar{B}_{s}^{0} ; p_{T} [GeV]; Counts",1000, 0., 100.);
  
  hBsToMuMuG = new TH1D("hBsToMuMuG ", "Number of B_{s}^{0}/#bar{B}_{s}^{0} #rightarrow #mu^{+}#mu^{-}#gamma decays; Decays per Event; Counts", 6, -0.5, 5.5);
  hBsToBd =  new TH1D("hBsToBd","Ratio of B_{s} to B_{d} mesons ; Ratio; Counts",56, -0.1, 5.5);
  hBsToB =  new TH1D("hBsToB","Ratio of B_{s} to all B mesons ; Ratio; Counts",120, -0.15, 1.05); 
  
  hHad =  new TH1D("hHad","Number of b quarks hadronizing into a B meson; B group; Events",45, 505.5, 550.5); 
  hBsParents =  new TH1D("hBsParents","Parents of B_{s}^{0} decaing into #mu^{+}#mu^{-}#gamma; B group; Events",45, 505.5, 550.5);
  hBsAncestor =  new TH1D("hBsAncestor","First produced B meson that evolved into B_{s}^{0} decaing into #mu^{+}#mu^{-}#gamma; ID; Events",45, 505.5, 550.5); 
  hBsProduct =  new TH1D("hBsProduct","Production of B_{s}^{0} ; Decay Channel; Events",7, -1.5, 5.5); 

  hGBsGDeltaR = new TH1D("hGBsGDeltaR","Transverse momentum of muons; #delta R; Events",200, 0., 1.);
  hGBsEta_GEta = new TH2D("hGBsEta_GEta","Pseudorapidity of #gamma; #eta_{B^{*}_{0}}; #eta_{B_{0}}; Events",100, -4., 4., 100, -4., 4.);
  hGammaBsPt = new TH1D("hGammaBsPt","Transverse momentum of #gamma produced with B_{s}^{0} ; p_{T} [GeV]; Events",1000, 0., 100.);


  hMuPt = new TH1D("hMuPt","Transverse momentum of muons; p_{T} [GeV]; Counts",1000, 0., 100.);
  hMuEta = new TH1D("hMuEta","Pseudorapidity of muons; #eta; Counts",1000, -5., 5.);
  hMu1Eta_Mu2Eta = new TH2D("hMu1Eta_Mu2Eta","Pseudorapidity of muons; #eta_{1}; #eta_{2}; Counts",100, -4., 4., 100, -4., 4.);
  hGammaPt = new TH1D("hGammaPt","Transverse momentum of photons; p_{T} [GeV]; Counts",1000, 0., 100.);
  hGammaEta = new TH1D("hGammaEta","Pseudorapidity of photons; #eta; Counts",1000, -5., 5.);

  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hBPt -> Write();
  hBsPt -> Write();
  hBPtFrag-> Write();
  hBssPtFrag-> Write();
  hBsPtFrag-> Write();
  hAntiBsPtFrag-> Write();
  
  hBsToMuMuG-> Write();
  hBsToBd-> Write();
  hBsToB-> Write();
  
  hHad-> Write();
  hBsParents -> Write();
  hBsAncestor -> Write();
  hBsProduct -> Write();

  hGBsGDeltaR-> Write();
  hGBsEta_GEta -> Write();
  hGammaBsPt -> Write();

  hMuPt-> Write();
  hMuEta-> Write();
  hMu1Eta_Mu2Eta-> Write();
  hGammaPt-> Write();
  hGammaEta-> Write();

  myRootFile.Close();

  delete hBPt;
  delete hBsPt;
  delete hBPtFrag;
  delete hBssPtFrag;
  delete hBsPtFrag;
  delete hAntiBsPtFrag;

  delete hBsToMuMuG;
  delete hBsToBd;
  delete hBsToB;
  
  delete hHad;
  delete hBsParents ;
  delete hBsAncestor ;
  delete hBsProduct ;

  delete hGBsGDeltaR;
  delete hGBsEta_GEta ;
  delete hGammaBsPt ;

  delete hMuPt;
  delete hMuEta;
  delete hMu1Eta_Mu2Eta;
  delete hGammaPt;
  delete hGammaEta;

  cout << "HERE Mgr::endJob()" << endl;
}


void Analysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE Analysis::analyze "<< std::endl;

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);

  int nBmesons = 0;
  int nBd = 0;
  int nBs = 0;
  int nMuMuGamma = 0;
  vector<int> initialBpdg;
  vector<const reco::Candidate*>  firstAncestor;
  vector<int>  firstProductionBs;

  const reco::Candidate* Gamma0;

  //std::cout << "Number of particles: " << genPar.size() << std::endl;
  
  for (const auto& part:genPar ){
    
    if ( abs(part.pdgId())/100 == 5){ //B mesons

      //double dr = reco::deltaR2(part, part);
      
      std::cout << std::endl;
      //std::cout << "B meseon:  " << part.pdgId() << std::endl;

      hBPt -> Fill(part.pt());

      int nMothers = part.numberOfMothers();
      int nDaughters = part.numberOfDaughters();
      bool isOldestB = true; //initial assumption -this B is the oldest one

      //Check the family tree

      std::vector<int> productionChannel( std::vector<int>(nMothers, 0)); //store information about the production channel

      std::cout << "Number of mothers: " << nMothers << std::endl;

      for (int i = 0; i < nMothers; ++i){ // loop over mothers
          std::cout << "Mother:  " ;
          std::cout << part.mother(i)->pdgId() ;
          //std::cout << "Grandmother:  " ; // 2 generations back
          productionChannel[i] = abs(part.mother(i)->pdgId());
          /*
          for(reco::Candidate::size_type k = 0; k < part.mother(i)->numberOfMothers();++k){ //loop over grandparents
            std::cout << part.mother(i)-> mother(k)->pdgId() << " ";
          }
*/
          std::cout << std::endl;

          if ( isOldestB == true && (part.mother(i)->pdgId())/100 == 5 ){
            isOldestB = false; //B among the parents - not the oldest B then
          }
      } //end of the loop over mothers

      if( isOldestB == true){
        
        initialBpdg.push_back( part.pdgId() ); //no abs

      } 
      //std::cout << "Is the oldest B? " << isOldestB << " particle: " <<  part.pdgId() << std::endl;


      //Check the descendants

      std::vector<std::vector<int>> decayChannel(2, std::vector<int>(nDaughters, 0)); //store information about the decay channel
                                                                                      // first row - IDs of the children; second row - if it's a B or not

      std::cout << "Number of daughters: " << nDaughters << std::endl;
      std::cout << "Daughters:  " ;

      for (int i = 0; i < nDaughters; ++i){ // loop over daughters
      
        const reco::Candidate* daughter = part.daughter(i);

        decayChannel[0][i] = daughter->pdgId();

        if ( abs(decayChannel[0][i])/100 != 5 ){ //not a B
           decayChannel[1][i] = 1;
        }
        std::cout << part.daughter(i)->pdgId() << "  ";
      } //end of the loop over daughters


      int sum = std::accumulate(decayChannel[1].begin(), decayChannel[1].end(), 0); //sum of the second row

      std::cout << std::endl;
      //std::cout <<"Number of not-B particles among the daughters: "<< sum << std::endl;

      if( sum == nDaughters) { 
        std::cout << "Last B in the chain" << std::endl;
        nBmesons++;

        hBPtFrag -> Fill(part.pt());

        if( abs(part.pdgId()) == 531){
          hBssPtFrag -> Fill(part.pt());
          nBs++;

          if(part.pdgId() == 531 ){
            hBsPtFrag -> Fill(part.pt());
          }else{
            hAntiBsPtFrag -> Fill(part.pt());
          }

        }else if(abs(part.pdgId()) < 530){
          //cout << "Bd Meson" << std::endl;
          nBd++;
        }
      }

      // fill muon and photon histograms from the Bs decay into Mu Mu gamma
      if(abs(part.pdgId()) == 531 && sum == nDaughters && isSameChannel(decayChannel[0], MuMuG)){ //this can be put in the line 301
        vector<float> muEta;
        
        for (int i = 0; i < nDaughters; ++i) {
          const reco::Candidate* daughter = part.daughter(i);
          if (abs(daughter->pdgId()) == 13){
            Gamma0 = daughter;
            hMuPt -> Fill(Gamma0->pt());
            hMuEta -> Fill(daughter->eta());
            muEta.push_back(daughter->eta());
          }
          if (abs(daughter->pdgId()) == 22){
            hGammaPt -> Fill(daughter->pt());
            hGammaEta -> Fill(daughter->eta());
          }
        }
        
        const reco::Candidate* ancestor = static_cast<const reco::Candidate*>(&part); //ancestor = currently analized B
        
        firstAncestor.push_back(ancestor);

        while (ancestor != nullptr) { 
          bool foundMotherB = false;
          cout << "Parent " << ancestor ->pdgId() <<std::endl;

          for (unsigned int i = 0; i < ancestor->numberOfMothers(); ++i) {
            //cout << ancestor->mother(i)->pdgId() << " " << std::endl;
/*
            if(firstAncestor.size() == 2 && isSameChannel(firstProductionBs, Bs)  ){
              firstProductionBs.push_back(abs(ancestor->mother(i) ->pdgId()));
              cout << "Rodzice" << abs(ancestor->mother(i) ->pdgId()) <<endl;
              cout << "First production of Bs!!!!!" << std::endl;
                  
            }
*/
            if (foundMotherB == false && (abs(ancestor->mother(i)->pdgId()) / 100) == 5) { //mother B found

                if(ancestor == &part) {
                  hBsParents -> Fill( abs(ancestor->mother(i) ->pdgId()) );
                  
                  if( isSameChannel(productionChannel, BsStarG)){
                    hBsProduct -> Fill(1) ;
                  }else if( isSameChannel(productionChannel, Bs)){
                    hBsProduct -> Fill(2) ;
                  }else{
                    hBsProduct -> Fill(0) ;
                  }
                }
                firstAncestor.push_back(ancestor->mother(i));
                ancestor = ancestor->mother(i); // new ancestor
                
                foundMotherB = true;
                
                break; // stop, and enter the while loop once again
            }
          }
          
          if (!foundMotherB) break;  // stop, if there is no B parent
        }


        cout << firstAncestor.size() << std::endl;

        cout << "The oldest mu mu gamma parent: " << firstAncestor.back() ->pdgId() << std::endl;
        hBsAncestor->Fill(abs(firstAncestor.back() ->pdgId()));
     
        hMu1Eta_Mu2Eta -> Fill(muEta[0], muEta[1]);
      }
      //end of the loop over the youngest Bs decaying into mu mu gamma

      if (abs(part.pdgId()) == 531){ //this can by put in some other if as well (?)

        hBsPt -> Fill(part.pt());
      }
      /*
      for(int i = 0; i < nDaughters; ++i){
        std::cout << decayChannel[1][i] << "  " ;

      }
*/
      if (isSameChannel(decayChannel[0], MuMuG)) nMuMuGamma++; //again (?) can be in some already existing if
      std::cout << std::endl;
    }//end of the "if" for B mesons
      
  } // end of the loop over gen particles
  


  if (initialBpdg.size() == 2 ){
    
    if(initialBpdg.at(0) ==firstAncestor.back()->pdgId() && initialBpdg.at(1) == firstAncestor.back()->pdgId()){ //if both are the same
      
      hHad ->Fill(firstAncestor.back()->pdgId());
    }else if (initialBpdg.at(0) ==firstAncestor.back()->pdgId()){ //if the first is the Bs ancestor, take the other one
      hHad ->Fill(abs(initialBpdg.at(1)));
    }else if (initialBpdg.at(1) ==firstAncestor.back()->pdgId()){
      hHad ->Fill(abs(initialBpdg.at(0)));
    }else{
      hHad ->Fill(506);
      std::cout<<"Something is wrong" << std::endl;
    }
  }

  if( firstAncestor.size() ==1 ){
    hBsParents ->Fill(506);
  }

  hBsToMuMuG -> Fill (nMuMuGamma);
  hBsToBd -> Fill(nBs*1.0/nBd);
  hBsToB -> Fill(nBs*1.0/nBmesons);
  
  cout << "Number of Bs:  " << nBs << " , Bd: " << nBd << " ,ratio: " << nBs*1.0/nBd << "  , B mesons " << nBmesons <<std::endl;
  cout << "Number of Bs -> Mu Mu Gamma: " << nMuMuGamma <<std::endl;
  cout <<"Number of Bs mesons: " << nBs << endl;
  for(const auto& m:firstAncestor  ){
    if(abs(m->pdgId()) == 533 ){
      for(reco::Candidate::size_type i=0; i<m->numberOfDaughters(); ++i){
        
        if(m->daughter(i)->pdgId() ==22){
          double deltaR = std::sqrt( reco::deltaR2(Gamma0->eta(), Gamma0->phi(), m->daughter(i)->eta(), m->daughter(i)->phi()) );
          hGBsGDeltaR ->Fill (deltaR);
          hGammaBsPt -> Fill(m->daughter(i)->pt());
          hGBsEta_GEta -> Fill(m->daughter(i)->eta(), Gamma0 -> eta());
          //cout << "deltaR" << deltaR <<endl;
        }
      }
    }
  }
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);

