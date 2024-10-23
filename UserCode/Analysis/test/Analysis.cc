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

  bool isSameDecay(const std::vector<int>&, const std::vector<int>&);
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

  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;

  std::vector<int> MuMuG = {22, 13, -13};

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

bool Analysis::isSameDecay(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
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
  hBsToBd =  new TH1D("hBsToBd","Ratio of B_{s} to B_{d} mesons ; Ratio; Counts",120, -0.15, 1.05);
  hBsToB =  new TH1D("hBsToB","Ratio of B_{s} to all B mesons ; Ratio; Counts",120, -0.15, 1.05); 


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

  std::cout << "Number of particles: " << genPar.size() << std::endl;
  //for (std::vector<reco::GenParticle>::const_iterator part = genPar.begin(); part < genPar.end(); part++) {}

  for (const auto& part:genPar ){
    
    if ( abs(part.pdgId())/100 == 5){

      std::cout << std::endl;
      std::cout << "B meseon:  " << part.pdgId() << std::endl;
      hBPt -> Fill(part.pt());
      int nDaughters = part.numberOfDaughters();

      std::vector<std::vector<int>> decayChannel(2, std::vector<int>(nDaughters, 0));

      std::cout << "Number of daughters: " << nDaughters << std::endl;
      std::cout << "Daughters:  " ;

      for (int i = 0; i < nDaughters; ++i) {
        const reco::Candidate* daughter = part.daughter(i);

        decayChannel[0][i] = daughter->pdgId();

        if ( abs(decayChannel[0][i])/100 != 5 ) decayChannel[1][i] = 1;

        std::cout << part.daughter(i)->pdgId() << "  ";
      }


      int sum = std::accumulate(decayChannel[1].begin(), decayChannel[1].end(), 0);

      std::cout << std::endl;
      std::cout <<"Number of not-B particles: "<< sum << std::endl;

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
          cout << "Bd Meson" << std::endl;
          nBd++;
        }
      }

      if (abs(part.pdgId()) == 531){

        hBsPt -> Fill(part.pt());
      }
      /*
      for(int i = 0; i < nDaughters; ++i){
        std::cout << decayChannel[1][i] << "  " ;

      }
*/
      if (isSameDecay(decayChannel[0], MuMuG)) nMuMuGamma++;
      std::cout << std::endl;
    }
      //std::cout << "Id: " << part.pdgId() << std::endl;
  }

  hBsToMuMuG -> Fill (nMuMuGamma);
  hBsToBd -> Fill(nBs*1.0/nBd);
  hBsToB -> Fill(nBs*1.0/nBmesons);
  
  cout << "Number of Bs:  " << nBs << " , Bd: " << nBd << " ,ratio: " << nBs*1.0/nBd << "  , B mesons " << nBmesons <<std::endl;
  cout << "Number of Bs -> Mu Mu Gamma: " << nMuMuGamma <<std::endl;
  cout <<"Number of Bs mesons: " << nBs << endl;
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);

