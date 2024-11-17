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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/Math/interface/deltaR.h"

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
class DataAnalysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit DataAnalysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~DataAnalysis();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool isSameDecay(const std::vector<int>&, const std::vector<int>&);
private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  edm::EDGetTokenT < vector<reco::Muon> > theMuonToken;
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;

  // histograms
  TH1D* hMuPt;
  TH1D* hGammaPt;

  TH2D* hRecoMuPtVsEta;
  TH2D* hRecoGammaPtVsEta;
  std::vector<int> MuMuG = {22, 13, -13};

};


DataAnalysis::DataAnalysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));

}

DataAnalysis::~DataAnalysis()
{
  cout <<" DTOR" << endl;
}

bool DataAnalysis::isSameDecay(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}


void DataAnalysis::beginJob()
{
  //create a histogram
  hMuPt = new TH1D("hRecoMuPt", "reco muon pT; pT [GeV]; counts", 100, 0, 30);
  hGammaPt = new TH1D("hRecoGammaPt", "reco photon pT; pT [GeV]; counts", 100, 0, 30);

  hRecoMuPtVsEta = new TH2D("hRecoMuPtVsEta", "reco muon pT vs eta; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);
  hRecoGammaPtVsEta = new TH2D("hRecoGammaPtVsEta", "reco photon pT vs eta; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);

  cout << "HERE DataAnalysis::beginJob()" << endl;
}

void DataAnalysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hMuPt -> Write();
  hGammaPt -> Write();

  hRecoMuPtVsEta -> Write();
  hRecoGammaPtVsEta -> Write();

  myRootFile.Close();

  delete hMuPt;
  delete hGammaPt;

  delete hRecoMuPtVsEta;
  delete hRecoGammaPtVsEta;

  cout << "HERE DataAnalysis::endJob()" << endl;
}


void DataAnalysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE DataAnalysis::analyze "<< std::endl;

  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);

  for (const auto& recoMu : recoMuons)
  {
    hMuPt->Fill(recoMu.pt());
    hRecoMuPtVsEta->Fill(recoMu.pt(), recoMu.eta());
  }

  for (const auto& recoPh : recoPhotons)
  {
    hGammaPt->Fill(recoPh.pt());
    hRecoGammaPtVsEta->Fill(recoPh.pt(), recoPh.eta());
  }


  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(DataAnalysis);

