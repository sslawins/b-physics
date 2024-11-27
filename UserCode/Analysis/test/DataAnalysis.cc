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

#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

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
class DataAnalysis : public edm::stream::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit DataAnalysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~DataAnalysis();

  //edm filter plugin specific functions
  // virtual void beginJob();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  // virtual void endJob();
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  bool isSameDecay(const std::vector<int>&, const std::vector<int>&);
private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  edm::EDGetTokenT < vector<reco::Muon> > theMuonToken;
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;
  HLTConfigProvider hltConfig;
  std::string processName_;

  // histograms
  TH1D* hMuPt;
  TH1D* hGammaPt;

  TH2D* hRecoMuPtVsEta;
  TH2D* hRecoGammaPtVsEta;
  std::vector<int> MuMuG = {22, 13, -13};

};


DataAnalysis::DataAnalysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0), processName_(conf.getParameter<std::string>("processName"))
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


// void DataAnalysis::beginJob()
// {
//   //create a histogram
//   hMuPt = new TH1D("hRecoMuPt", "reco muon pT; pT [GeV]; counts", 100, 0, 30);
//   hGammaPt = new TH1D("hRecoGammaPt", "reco photon pT; pT [GeV]; counts", 100, 0, 30);

//   hRecoMuPtVsEta = new TH2D("hRecoMuPtVsEta", "reco muon pT vs eta; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);
//   hRecoGammaPtVsEta = new TH2D("hRecoGammaPtVsEta", "reco photon pT vs eta; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);

//   cout << "HERE DataAnalysis::beginJob()" << endl;
// }

// void DataAnalysis::endJob()
// {
//   //make a new Root file
//   TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

//   //write histogram data
//   hMuPt -> Write();
//   hGammaPt -> Write();

//   hRecoMuPtVsEta -> Write();
//   hRecoGammaPtVsEta -> Write();

//   myRootFile.Close();

//   delete hMuPt;
//   delete hGammaPt;

//   delete hRecoMuPtVsEta;
//   delete hRecoGammaPtVsEta;

//   cout << "HERE DataAnalysis::endJob()" << endl;
// }

void DataAnalysis::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig.init(iRun,iSetup,processName_,changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
     // The HLT config has actually changed wrt the previous Run, hence rebook your
     // histograms or do anything else dependent on the revised HLT config
     // hltConfig.dump("Streams");
     // hltConfig.dump("Datasets");
     // hltConfig.dump("Triggers");
     hltConfig.dump("PrescaleTable");
     // hltConfig.dump("ProcessPSet");
    }
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("MyAnalyzer") << " HLT config extraction failure with process name " << processName_;
    // In this case, all access methods will return empty values!
  }
}

void DataAnalysis::endRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  cout << "HERE DataAnalysis::endRun()" << endl;
}


void DataAnalysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
//   std::cout << " -------------------------------- HERE DataAnalysis::analyze "<< std::endl;

//   const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
//   const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);

//   // std::vector<std::string> tfiggerNames({"HLT_DoubleMu4_3_Bs_v19", "HLT_DoubleMu4_3_LowMass_v5", "HLT_DoubleMu4_LowMass_Displaced_v5", "HLT_DoubleMu4_3_Photon4_BsToMMG_v4", "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v4"});
//   // for (std::string triggerName : tfiggerNames)
//   // {
//   //   int prescaleValue = hltPrescale.prescaleValue(ev, es, triggerName);
//   //   cout << "Trigger: " << triggerName << " prescale: " << prescaleValue << endl;
//   // }
//   for (const auto& recoMu : recoMuons)
//   {
//     hMuPt->Fill(recoMu.pt());
//     hRecoMuPtVsEta->Fill(recoMu.pt(), recoMu.eta());
//   }

//   for (const auto& recoPh : recoPhotons)
//   {
//     hGammaPt->Fill(recoPh.pt());
//     hRecoGammaPtVsEta->Fill(recoPh.pt(), recoPh.eta());
//   }


//   cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(DataAnalysis);

