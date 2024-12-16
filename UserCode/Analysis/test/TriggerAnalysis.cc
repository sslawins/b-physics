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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

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
class TriggerAnalysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit TriggerAnalysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~TriggerAnalysis();

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
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;
  edm::EDGetTokenT < edm::TriggerResults > theTriggerResultsToken;

  // histograms


  TH1D* hHLT_DoubleMu4_3_Bs_v15;
  TH1D* hHLT_DoubleMu4_3_LowMass_v1;
  TH1D* hHLT_DoubleMu4_LowMass_Displaced_v1;
  TH1D* hHLT_DoubleMu4_3_Photon4_BsToMMG_v1;
  TH1D* hHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1;


  TH2D* hMuPtVsEtaWithTrigger;
  TH2D* hMuPtVsEta;

  TH1D* hMuPtWithRecoPhoton;

  TH2D* hMuPtVsEtaWithMuonConditionWithoutTrigger;

  std::vector<int> MuMuG = {22, 13, -13};

  int nPhotonsWithMuonCondition = 0;
  int nMuonsWithCondition = 0;
  int nTriggerMuMuG = 0;

  int nMuons34 = 0;
  int nMuons34_Photon = 0;

  int nHLT_DoubleMu4_3_Bs_v15 = 0;
  int nHLT_DoubleMu4_3_LowMass_v1 = 0;
  int nHLT_DoubleMu4_3_LowMass_v1_WithRecoMuons = 0;
  int nHLT_DoubleMu4_LowMass_Displaced_v1 = 0;
  int nHLT_DoubleMu4_3_Photon4_BsToMMG_v1 = 0;
  int nHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1 = 0;

  int nHLT_DoubleMu4_3_Bs_v15_RecoPhoton = 0;
  int nHLT_DoubleMu4_3_LowMass_v1_RecoPhoton = 0;
  int nHLT_DoubleMu4_LowMass_Displaced_v1_RecoPhoton = 0;
  int nHLT_DoubleMu4_3_Photon4_BsToMMG_v1_RecoPhoton = 0;
  int nHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1_RecoPhoton = 0;

  int nPhotonsWithMuonConditionWithoutTrigger = 0;

};


TriggerAnalysis::TriggerAnalysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles"));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));
  theTriggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));

}

TriggerAnalysis::~TriggerAnalysis()
{
  cout <<" DTOR" << endl;
}

bool TriggerAnalysis::isSameDecay(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}


void TriggerAnalysis::beginJob()
{
  //create a histogram
  hHLT_DoubleMu4_3_Bs_v15 = new TH1D("hHLT_DoubleMu4_3_Bs_v15", "HLT_DoubleMu4_3_Bs_v15; pT [GeV]; counts", 100, 0, 30);
  hHLT_DoubleMu4_3_LowMass_v1 = new TH1D("hHLT_DoubleMu4_3_LowMass_v1", "HLT_DoubleMu4_3_LowMass_v1; pT [GeV]; counts", 100, 0, 30);
  hHLT_DoubleMu4_LowMass_Displaced_v1 = new TH1D("hHLT_DoubleMu4_LowMass_Displaced_v1", "HLT_DoubleMu4_LowMass_Displaced_v1; pT [GeV]; counts", 100, 0, 30);
  hHLT_DoubleMu4_3_Photon4_BsToMMG_v1 = new TH1D("hHLT_DoubleMu4_3_Photon4_BsToMMG_v1", "HLT_DoubleMu4_3_Photon4_BsToMMG_v1; pT [GeV]; counts", 100, 0, 30);
  hHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1 = new TH1D("hHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1", "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1; pT [GeV]; counts", 100, 0, 30);

  hMuPtVsEtaWithTrigger = new TH2D("hMuPtVsEtaWithTrigger", "Muon pT vs Eta with trigger; pT [GeV]; Eta", 100, 0, 30, 100, -3, 3);
  hMuPtVsEta = new TH2D("hMuPtVsEta", "Muon pT vs Eta; pT [GeV]; Eta", 100, 0, 30, 100, -3, 3);

  hMuPtWithRecoPhoton = new TH1D("hMuPtWithRecoPhoton", "Muon pT with reco photon; pT [GeV]; counts", 100, 0, 30);

  hMuPtVsEtaWithMuonConditionWithoutTrigger = new TH2D("hMuPtVsEtaWithMuonConditionWithoutTrigger", "Muon pT vs Eta with muon condition without trigger; pT [GeV]; Eta", 100, 0, 30, 100, -3, 3);


  cout << "HERE TriggerAnalysis::beginJob()" << endl;
}

void TriggerAnalysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data

  hHLT_DoubleMu4_3_Bs_v15 -> Write();
  hHLT_DoubleMu4_3_LowMass_v1 -> Write();
  hHLT_DoubleMu4_LowMass_Displaced_v1 -> Write();
  hHLT_DoubleMu4_3_Photon4_BsToMMG_v1 -> Write();
  hHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1 -> Write();

  hMuPtVsEtaWithTrigger -> Write();
  hMuPtVsEta -> Write();

  hMuPtWithRecoPhoton -> Write();

  hMuPtVsEtaWithMuonConditionWithoutTrigger -> Write();

  myRootFile.Close();

  delete hHLT_DoubleMu4_3_Bs_v15;
  delete hHLT_DoubleMu4_3_LowMass_v1;
  delete hHLT_DoubleMu4_LowMass_Displaced_v1;
  delete hHLT_DoubleMu4_3_Photon4_BsToMMG_v1;
  delete hHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1;

  delete hMuPtVsEtaWithTrigger;
  delete hMuPtVsEta;

  delete hMuPtWithRecoPhoton;

  delete hMuPtVsEtaWithMuonConditionWithoutTrigger;

  cout << "nPhotonsWithMuonCondition: " << nPhotonsWithMuonCondition << endl;\
  cout << "nMuonsWithCondition: " << nMuonsWithCondition << endl;
  cout << "nTriggerMuMuG: " << nTriggerMuMuG << endl << endl;

  cout << "nMuons34: " << nMuons34 << endl;
  cout << "nMuons34_Photon: " << nMuons34_Photon << endl << endl;

  cout << "nHLT_DoubleMu4_3_Bs_v15: " << nHLT_DoubleMu4_3_Bs_v15 << endl;
  cout << "nHLT_DoubleMu4_3_LowMass_v1: " << nHLT_DoubleMu4_3_LowMass_v1 << endl;
  cout << "nHLT_DoubleMu4_LowMass_Displaced_v1: " << nHLT_DoubleMu4_LowMass_Displaced_v1 << endl;
  cout << "nHLT_DoubleMu4_3_Photon4_BsToMMG_v1: " << nHLT_DoubleMu4_3_Photon4_BsToMMG_v1 << endl;
  cout << "nHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1: " << nHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1 << endl << endl;

  cout << "nHLT_DoubleMu4_3_Bs_v15_RecoPhoton: " << nHLT_DoubleMu4_3_Bs_v15_RecoPhoton << endl;
  cout << "nHLT_DoubleMu4_3_LowMass_v1_RecoPhoton: " << nHLT_DoubleMu4_3_LowMass_v1_RecoPhoton << endl;
  cout << "nHLT_DoubleMu4_LowMass_Displaced_v1_RecoPhoton: " << nHLT_DoubleMu4_LowMass_Displaced_v1_RecoPhoton << endl;
  cout << "nHLT_DoubleMu4_3_Photon4_BsToMMG_v1_RecoPhoton: " << nHLT_DoubleMu4_3_Photon4_BsToMMG_v1_RecoPhoton << endl;
  cout << "nHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1_RecoPhoton: " << nHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1_RecoPhoton << endl << endl;

  cout << "nHLT_DoubleMu4_3_LowMass_v1_WithRecoMuons: " << nHLT_DoubleMu4_3_LowMass_v1_WithRecoMuons << endl;
  cout << "nPhotonsWithMuonConditionWithoutTrigger: " << nPhotonsWithMuonConditionWithoutTrigger << endl;

  cout << "HERE TriggerAnalysis::endJob()" << endl;
}


void TriggerAnalysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE TriggerAnalysis::analyze "<< std::endl;

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);

  vector<const reco::Candidate*> genMuons;
  vector<const reco::Muon*> recoMatchedMuons;
  vector<const reco::Candidate*> genMatchedMuons;

  vector<const reco::Candidate*> genPhotons;
  vector<const reco::Photon*> recoMatchedPhotons;
  vector<const reco::Candidate*> genMatchedPhotons;

  for(const auto& genP : genPar)
  {
    if (abs(genP.pdgId()) == 531)
    {
      vector<int> daughters;
      for(unsigned int i=0; i < genP.numberOfDaughters(); i++)
      {
        daughters.push_back(genP.daughter(i)->pdgId());
      }
      if(isSameDecay(daughters, MuMuG))
      {
        for(unsigned int i=0; i < genP.numberOfDaughters(); i++)
        {
          if(abs(genP.daughter(i)->pdgId()) == 13) genMuons.push_back(genP.daughter(i));
          if(abs(genP.daughter(i)->pdgId()) == 22) genPhotons.push_back(genP.daughter(i));
        }
      }
    }
  }


  // reco muon matching
  for (const reco::Candidate* genMu : genMuons)
  {
    float minDR = 10;
    const reco::Muon* bestMatchedMuon;
    bool matched = false;
    for (const auto& recoMu : recoMuons)
    {
      float dR = reco::deltaR(recoMu, *genMu);
      if (dR < minDR)
      {
        minDR = dR;
        bestMatchedMuon = &recoMu;
        matched = true;
      }
    }
    if (matched && minDR < 0.01)
    {
      recoMatchedMuons.push_back(bestMatchedMuon);
      genMatchedMuons.push_back(genMu);
    }
  }

  // reco photon matching
  for (const reco::Candidate* genPh : genPhotons)
  {
    float minDR = 10;
    const reco::Photon* bestMatchedPhoton;
    bool matched = false;
    for (const auto& recoPh : recoPhotons)
    {
      float dR = reco::deltaR(recoPh, *genPh);
      if (dR < minDR)
      {
        minDR = dR;
        bestMatchedPhoton = &recoPh;
        matched = true;
      }
    }
    if (matched && minDR < 0.02)
    {
      recoMatchedPhotons.push_back(bestMatchedPhoton);
      genMatchedPhotons.push_back(genPh);
    }
  }

  if(recoMatchedMuons.size() == 2)
  {
    if(recoMatchedMuons.at(0)->pt() > 4 && recoMatchedMuons.at(1)->pt() > 4 &&
    abs(recoMatchedMuons.at(0)->eta()) < 2.4 && abs(recoMatchedMuons.at(0)->eta()) < 2.4)
    {
      nMuonsWithCondition++;
      if (recoMatchedPhotons.size() > 0) nPhotonsWithMuonCondition++;
    }
  }

    if(recoMatchedMuons.size() == 2)
  {
    if(recoMatchedMuons.at(0)->pt() > 3 && recoMatchedMuons.at(1)->pt() > 3 && recoMatchedMuons.at(0)->pt() < 4 && recoMatchedMuons.at(1)->pt() < 4 &&
    abs(recoMatchedMuons.at(0)->eta()) < 2.4 && abs(recoMatchedMuons.at(0)->eta()) < 2.4)
    {
      nMuons34++;
      if (recoMatchedPhotons.size() > 0) nMuons34_Photon++;
    }
  }

  if(recoMatchedPhotons.size() > 0)
  {
    for(const auto recoMatchedMu : recoMatchedMuons)
    {
      hMuPtWithRecoPhoton->Fill(recoMatchedMu->pt());
    }
  }


  //trigger info
  const edm::TriggerResults & triggerResults = ev.get(theTriggerResultsToken);
  edm::TriggerNames triggerNames = ev.triggerNames(triggerResults);

  for (unsigned int i = 0; i < triggerResults.size(); i++) cout << triggerNames.triggerName(i) << endl;

  bool triggerFired = false;
  for (unsigned int i = 0; i < triggerResults.size(); i++)
  {
    TString name = triggerNames.triggerName(i);
    if(name == "HLT_DoubleMu4_3_Bs_v15" || name == "HLT_DoubleMu4_3_LowMass_v1" || name == "HLT_DoubleMu4_LowMass_Displaced_v1"  ||
    name == "HLT_DoubleMu4_3_Photon4_BsToMMG_v1" || name == "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1")
    {
      cout << "Trigger: " << name << "  " << triggerResults.accept(i) << endl;
      if(triggerResults.accept(i)) triggerFired = true;
    }

    if(name == "HLT_DoubleMu4_3_Bs_v15" && triggerResults.accept(i) == 1)
    {
      nHLT_DoubleMu4_3_Bs_v15++;
      if (recoMatchedPhotons.size() > 0) nHLT_DoubleMu4_3_Bs_v15_RecoPhoton++;
      for(const auto recoMatchedMu : recoMatchedMuons)
      {
        if(recoMatchedMu->charge() > 0) {hHLT_DoubleMu4_3_Bs_v15->Fill(recoMatchedMu->pt()); cout << "Filling" << endl;}
      }
    }
    if(name == "HLT_DoubleMu4_3_LowMass_v1" && triggerResults.accept(i) == 1)
    {
      if (recoMatchedMuons.size() == 2) nHLT_DoubleMu4_3_LowMass_v1_WithRecoMuons++;
      nHLT_DoubleMu4_3_LowMass_v1++;
      if (recoMatchedPhotons.size() > 0) nHLT_DoubleMu4_3_LowMass_v1_RecoPhoton++;
      for(const auto recoMatchedMu : recoMatchedMuons)
      {
        if(recoMatchedMu->charge() > 0) hHLT_DoubleMu4_3_LowMass_v1->Fill(recoMatchedMu->pt());
        hMuPtVsEtaWithTrigger->Fill(recoMatchedMu->pt(), recoMatchedMu->eta());
      }
    }

    if(recoMatchedMuons.size() == 2)
    {
      if(recoMatchedMuons.at(0)->pt() > 4 && recoMatchedMuons.at(1)->pt() > 4 &&
      abs(recoMatchedMuons.at(0)->eta()) < 2.4 && abs(recoMatchedMuons.at(0)->eta()) < 2.4)
      {
        if(name == "HLT_DoubleMu4_3_LowMass_v1" && triggerResults.accept(i) != 1)
        {
          hMuPtVsEtaWithMuonConditionWithoutTrigger->Fill(recoMatchedMuons.at(0)->pt(), recoMatchedMuons.at(0)->eta());
          if (recoMatchedPhotons.size() > 0) nPhotonsWithMuonConditionWithoutTrigger++;
        }
      }
    }


    if(name == "HLT_DoubleMu4_LowMass_Displaced_v1" && triggerResults.accept(i) == 1)
    {
      nHLT_DoubleMu4_LowMass_Displaced_v1++;
      if (recoMatchedPhotons.size() > 0) nHLT_DoubleMu4_LowMass_Displaced_v1_RecoPhoton++;
      for(const auto recoMatchedMu : recoMatchedMuons)
      {
        if(recoMatchedMu->charge() > 0) hHLT_DoubleMu4_LowMass_Displaced_v1->Fill(recoMatchedMu->pt());
      }
    }
    if(name == "HLT_DoubleMu4_3_Photon4_BsToMMG_v1" && triggerResults.accept(i) == 1)
    {
      nHLT_DoubleMu4_3_Photon4_BsToMMG_v1++;
      if (recoMatchedPhotons.size() > 0) nHLT_DoubleMu4_3_Photon4_BsToMMG_v1_RecoPhoton++;
      for(const auto recoMatchedMu : recoMatchedMuons)
      {
        if(recoMatchedMu->charge() > 0) hHLT_DoubleMu4_3_Photon4_BsToMMG_v1->Fill(recoMatchedMu->pt());
      }
    }
    if(name == "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1" && triggerResults.accept(i) == 1)
    {
      nHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1++;
      if (recoMatchedPhotons.size() > 0) nHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1_RecoPhoton++;
      for(const auto recoMatchedMu : recoMatchedMuons)
      {
        if(recoMatchedMu->charge() > 0) hHLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v1->Fill(recoMatchedMu->pt());
      }
    }
  }
  if(triggerFired) nTriggerMuMuG++;

  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(TriggerAnalysis);

