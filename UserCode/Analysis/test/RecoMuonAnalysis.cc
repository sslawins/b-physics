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
#include <DataFormats/MuonReco/interface/Muon.h>

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
class RecoMuonAnalysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit RecoMuonAnalysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~RecoMuonAnalysis();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool isSameDecay(const std::vector<int>&, const std::vector<int>&);
private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;
  edm::EDGetTokenT < vector<reco::Muon> > theMuonToken;

  // histograms
  TH1D* hMuPt;

  std::vector<int> MuMuG = {22, 13, -13};

};


RecoMuonAnalysis::RecoMuonAnalysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles"));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));

}

RecoMuonAnalysis::~RecoMuonAnalysis()
{
  cout <<" DTOR" << endl;
}

bool RecoMuonAnalysis::isSameDecay(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}


void RecoMuonAnalysis::beginJob()
{
  //create a histogram
  hMuPt = new TH1D("hRecoMuPt", "reco muon pT; pT; counts", 100, 0, 30);

  cout << "HERE RecoMuonAnalysis::beginJob()" << endl;
}

void RecoMuonAnalysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hMuPt -> Write();

  myRootFile.Close();

  delete hMuPt;

  cout << "HERE RecoMuonAnalysis::endJob()" << endl;
}


void RecoMuonAnalysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE RecoMuonAnalysis::analyze "<< std::endl;

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & muons = ev.get(theMuonToken);

  for (const auto& mu:muons ){
    hMuPt -> Fill(mu.pt());
  }


  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(RecoMuonAnalysis);

