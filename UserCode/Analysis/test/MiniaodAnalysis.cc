#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"

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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <sstream>
#include <iomanip> 
#include <utility>
#include <numeric>


using namespace std;


//object definition
class MiniaodAnalysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit MiniaodAnalysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~MiniaodAnalysis();

  //edm filter plugin specific functions
  virtual void beginJob();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
;
  bool isSameDecay(const std::vector<int>&, const std::vector<int>&);
private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  edm::EDGetTokenT < vector<pat::Photon> > thePhotonToken;
  edm::EDGetTokenT < vector<pat::PackedCandidate> > thePackedCandidateToken;
  edm::EDGetTokenT < vector<reco::Vertex> > theVertexToken;

  // histograms
  TH1D* hIso1;
  TH1D* hIso2;
  TH1D* hIso3;
  TH1D* hIso4;

  TProfile* hIso1VsNGamma;
  TProfile* hIso1VsNPV;

  std::vector<int> MuMuG = {22, 13, -13};

};


MiniaodAnalysis::MiniaodAnalysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  thePhotonToken = consumes< vector<pat::Photon>  >( edm::InputTag("slimmedPhotons"));
  thePackedCandidateToken = consumes< vector<pat::PackedCandidate>  >( edm::InputTag("packedPFCandidates"));
  theVertexToken = consumes< vector<reco::Vertex>  >( edm::InputTag("offlineSlimmedPrimaryVertices"));

}

MiniaodAnalysis::~MiniaodAnalysis()
{
  cout <<" DTOR" << endl;
}

bool MiniaodAnalysis::isSameDecay(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}


void MiniaodAnalysis::beginJob()
{
  //create a histogram
  hIso1 = new TH1D("hIso1", "hIso1; #frac{#sum pT_{PFCand}}{pT_{#gamma}}; counts}", 100, 0, 10);
  hIso2 = new TH1D("hIso2", "hIso2; #frac{#sum pT_{PFCand}}{pT_{#gamma}}; counts}", 100, 0, 10);
  hIso3 = new TH1D("hIso3", "hIso3; #frac{#sum pT_{PFCand}}{pT_{#gamma}}; counts}", 100, 0, 10);
  hIso4 = new TH1D("hIso4", "hIso4; #frac{#sum pT_{PFCand}}{pT_{#gamma}}; counts}", 100, 0, 10);

  hIso1VsNGamma = new TProfile("hIso1VsNGamma", "hIso1VsNGamma; #gamma multiplicity; #frac{#sum pT_{PFCand}}{pT_{#gamma}}", 10, 0, 10);
  hIso1VsNPV = new TProfile("hIso1VsNPV", "hIso1VsNPV; #PV multiplicity; #frac{#sum pT_{PFCand}}{pT_{#gamma}}", 10, 0, 10);


  cout << "HERE MiniaodAnalysis::beginJob()" << endl;
}

void MiniaodAnalysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hIso1->Write();
  hIso2->Write();
  hIso3->Write();
  hIso4->Write();

  hIso1VsNGamma->Write();
  hIso1VsNPV->Write();

  myRootFile.Close();

  delete hIso1;
  delete hIso2;
  delete hIso3;
  delete hIso4;

  delete hIso1VsNGamma;
  delete hIso1VsNPV;

  cout << "HERE MiniaodAnalysis::endJob()" << endl;
}


void MiniaodAnalysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE MiniaodAnalysis::analyze "<< std::endl;

  const std::vector<pat::Photon> & patPhotons = ev.get(thePhotonToken);
  const std::vector<pat::PackedCandidate> & patPackedCandidates = ev.get(thePackedCandidateToken);
  const std::vector<reco::Vertex> & patVertices = ev.get(theVertexToken);

  for (const auto& photon : patPhotons)
  {
    double sum = 0;
    for (const auto& packedCandidate : patPackedCandidates)
    {
      if(reco::deltaR(photon, packedCandidate) < 0.4)
      {
        sum += packedCandidate.pt();
      }
    }
    hIso1->Fill(sum/photon.pt());
    hIso1VsNGamma->Fill(patPhotons.size(), sum/photon.pt());
    hIso1VsNPV->Fill(patVertices.size(), sum/photon.pt());

    sum = 0;

    for (const auto& packedCandidate : patPackedCandidates)
    {
      if(reco::deltaR(photon, packedCandidate) < 0.1)
      {
        sum += packedCandidate.pt();
      }
    }
    hIso2->Fill(sum/photon.pt());

    sum = 0;

    for (const auto& packedCandidate : patPackedCandidates)
    {
      if(reco::deltaR(photon, packedCandidate) < 0.01)
      {
        sum += packedCandidate.pt();
      }
    }
    hIso3->Fill(sum/photon.pt());

    sum = 0;

    for (const auto& packedCandidate : patPackedCandidates)
    {
      if(reco::deltaR(photon, packedCandidate) < 0.001)
      {
        sum += packedCandidate.pt();
      }
    }
    hIso4->Fill(sum/photon.pt());

  }

  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(MiniaodAnalysis);

