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
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;
  edm::EDGetTokenT < edm::TriggerResults > theTriggerResultsToken;

  // histograms
  TH1D* hBs0Pt;

  TH1D* hMuPt;
  TH1D* hGammaPt;

  TH1D* hGenMuPt;
  TH1D* hRecoMuWithMuonConditionPt;
  TH2D* hRecoVsGenMuPtWithMuonCondition;

  TH2D* hGenMuPtVsEta;
  TH2D* hGenGammaPtVsEta;

  TH2D* hRecoMuPtVsEta;
  TH2D* hRecoGammaPtVsEta;

  TH2D* hRecoVsGenMuPt;
  TH2D* hRecoVsGenGammaPt;

  TH1D* hMuPtError;
  TH1D* hGammaPtError;

  TH1D* hMuDeltaR;
  TH1D* hGammaDeltaR;
  TH1D* hMu1VsMu2DeltaR;

  TH1D* hGammaWithMuonConditionPt;
  TH1D* hMatchedGammaWithMuonConditionPt;
  TH1D* hRecoGammaWithMuonConditionPt;

  TH2D* hRecoVsGenGammaWithMuonConditionPt;

  TH2D* hGammaWithMuonConditionPtVsEta;
  TH2D* hMatchedGammaWithMuonConditionPtVsEta;

  TH1D* BsInvMass;

  std::vector<int> MuMuG = {22, 13, -13};

  int nPhotonsWithMuonCondition = 0;
  int nMuonsWithCondition = 0;
  int nTriggerMuMuG = 0;
};


RecoMuonAnalysis::RecoMuonAnalysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles"));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));
  //theTriggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));

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
  hBs0Pt = new TH1D("Bs0Pt", "Bs0 pT; pT [GeV]; counts", 1000, 0, 100);

  hMuPt = new TH1D("hRecoMuPt", "reco muon pT; pT [GeV]; counts", 100, 0, 30);
  hGammaPt = new TH1D("hRecoGammaPt", "reco photon pT; pT [GeV]; counts", 100, 0, 30);

  hGenMuPt = new TH1D("hGenMuPt", "gen muon pT; pT [GeV]; counts", 500, 0, 30);
  hRecoMuWithMuonConditionPt = new TH1D("hRecoMuWithMuonConditionPt", "reco muon with #mu pT > 4 |#eta| < 2.4; pT [GeV]; counts", 500, 0, 30);
  hRecoVsGenMuPtWithMuonCondition = new TH2D("hRecoVsGenMuPtWithMuonCondition", "reco vs gen muon with #mu pT > 4 |#eta| < 2.4; gen pT [GeV]; reco pT [GeV]", 500, 0, 30, 500, 0, 30);

  hGenMuPtVsEta = new TH2D("hGenMuPtVsEta", "gen muon pT vs eta; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);
  hGenGammaPtVsEta = new TH2D("hGenGammaPtVsEta", "gen photon pT vs eta; pT [GeV]; eta", 500, 0, 30, 1000, -10, 10);

  hRecoMuPtVsEta = new TH2D("hRecoMuPtVsEta", "reco muon pT vs eta; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);
  hRecoGammaPtVsEta = new TH2D("hRecoGammaPtVsEta", "reco photon pT vs eta; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);

  hRecoVsGenMuPt = new TH2D("hRecoVsGenMuPt", "reco vs gen muon pT; gen pT [GeV]; reco pT [GeV]", 100, 0, 30, 100, 0, 30);
  hRecoVsGenGammaPt = new TH2D("hRecoVsGenGammaPt", "reco vs gen photon pT; gen pT [GeV]; reco pT [GeV]", 100, 0, 30, 100, 0, 30);

  hMuPtError = new TH1D("hMuPtError", "reco muon pT error; pT error [GeV]; counts", 100, -1, 1);
  hGammaPtError = new TH1D("hGammaPtError", "reco photon pT error; pT error [GeV]; counts", 100, -1, 1);

  hMuDeltaR = new TH1D("hMuDeltaR", "reco muon deltaR; deltaR; counts", 1000, 0, 0.1);
  hGammaDeltaR = new TH1D("hGammaDeltaR", "reco photon deltaR; deltaR; counts", 1000, 0, 0.1);
  hMu1VsMu2DeltaR = new TH1D("hMu1VsMu2DeltaR", "reco muon1 vs muon2 deltaR; deltaR; counts", 1000, 0, 10);

  hGammaWithMuonConditionPt = new TH1D("hGammaWithMuonConditionPt", "gen #gamma with #mu pT > 4 |#eta| < 2.4; pT [GeV]; counts", 500, 0, 30);
  hMatchedGammaWithMuonConditionPt = new TH1D("hMatchedGammaWithMuonConditionPt", "matched gen #gamma with #mu pT > 4 |#eta| < 2.4; pT [GeV]; counts", 500, 0, 30);
  hRecoGammaWithMuonConditionPt = new TH1D("hRecoGammaWithMuonConditionPt", "reco #gamma with #mu pT > 4 |#eta| < 2.4; pT [GeV]; counts", 500, 0, 30);
  
  hRecoVsGenGammaWithMuonConditionPt = new TH2D("hRecoVsGenGammaWithMuonConditionPt", "reco vs gen #gamma with #mu pT > 4 |#eta| < 2.4; gen pT [GeV]; reco pT [GeV]", 500, 0, 30, 500, 0, 30);

  hGammaWithMuonConditionPtVsEta = new TH2D("hGammaWithMuonConditionPtVsEta", "gen #gamma with #mu pT > 4 |#eta| < 2.4; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);
  hMatchedGammaWithMuonConditionPtVsEta = new TH2D("hMatchedGammaWithMuonConditionPtVsEta", "matched gen #gamma with #mu pT > 4 |#eta| < 2.4; pT [GeV]; eta", 100, 0, 30, 100, -3, 3);

  BsInvMass = new TH1D("BsInvMass", "Bs invariant mass; mass [GeV]; counts", 48, 5, 6.2);

  cout << "HERE RecoMuonAnalysis::beginJob()" << endl;
}

void RecoMuonAnalysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hBs0Pt -> Write();

  hMuPt -> Write();
  hGammaPt -> Write();

  hGenMuPt -> Write();
  hRecoMuWithMuonConditionPt -> Write();
  hRecoVsGenMuPtWithMuonCondition -> Write();


  hGenMuPtVsEta -> Write();
  hGenGammaPtVsEta -> Write();

  hRecoMuPtVsEta -> Write();
  hRecoGammaPtVsEta -> Write();

  hRecoVsGenMuPt -> Write();
  hRecoVsGenGammaPt -> Write();

  hMuPtError -> Write();
  hGammaPtError -> Write();

  hMuDeltaR -> Write();
  hGammaDeltaR -> Write();
  hMu1VsMu2DeltaR -> Write();

  hGammaWithMuonConditionPt -> Write();
  hMatchedGammaWithMuonConditionPt -> Write();
  hRecoGammaWithMuonConditionPt -> Write();

  hRecoVsGenGammaWithMuonConditionPt -> Write();

  hGammaWithMuonConditionPtVsEta -> Write();
  hMatchedGammaWithMuonConditionPtVsEta -> Write();

  BsInvMass -> Write();

  myRootFile.Close();

  delete hBs0Pt;

  delete hMuPt;
  delete hGammaPt;

  delete hGenMuPt;
  delete hRecoMuWithMuonConditionPt;
  delete hRecoVsGenMuPtWithMuonCondition;

  delete hGenMuPtVsEta;
  delete hGenGammaPtVsEta;

  delete hRecoMuPtVsEta;
  delete hRecoGammaPtVsEta;

  delete hRecoVsGenMuPt;
  delete hRecoVsGenGammaPt;

  delete hMuPtError;
  delete hGammaPtError;

  delete hMuDeltaR;
  delete hGammaDeltaR;
  delete hMu1VsMu2DeltaR;

  delete hGammaWithMuonConditionPt;
  delete hMatchedGammaWithMuonConditionPt;
  delete hRecoGammaWithMuonConditionPt;

  delete hRecoVsGenGammaWithMuonConditionPt;

  delete hGammaWithMuonConditionPtVsEta;
  delete hMatchedGammaWithMuonConditionPtVsEta;

  delete BsInvMass;

  cout << "nPhotonsWithMuonCondition: " << nPhotonsWithMuonCondition << endl;\
  cout << "nMuonsWithCondition: " << nMuonsWithCondition << endl;
  cout << "nTriggerMuMuG: " << nTriggerMuMuG << endl;
  cout << "HERE RecoMuonAnalysis::endJob()" << endl;
}


void RecoMuonAnalysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE RecoMuonAnalysis::analyze "<< std::endl;

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);

  vector<const reco::Candidate*> genMuons;
  vector<const reco::Muon*> recoMatchedMuons;
  vector<const reco::Candidate*> genMatchedMuons;

  vector<const reco::Candidate*> genPhotons;
  vector<const reco::Photon*> recoMatchedPhotons;
  vector<const reco::Candidate*> genMatchedPhotons;

  // trigger info
  // const edm::TriggerResults & triggerResults = ev.get(theTriggerResultsToken);
  // edm::TriggerNames triggerNames = ev.triggerNames(triggerResults);

  // bool triggerFired = false;
  // for (unsigned int i = 0; i < triggerResults.size(); i++)
  // {
  //   TString name = triggerNames.triggerName(i);
  //   if(name == "HLT_DoubleMu4_3_Bs_v19" || name == "HLT_DoubleMu4_3_LowMass_v5" || name == "HLT_DoubleMu4_LowMass_Displaced_v5"  ||
  //   name == "HLT_DoubleMu4_3_Photon4_BsToMMG_v4" || name == "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG_v4")
  //   cout << "Trigger: " << name << "  " << triggerResults.accept(i) << endl;
  //   if(triggerResults.accept(i)) triggerFired = true;
  // }

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
        hBs0Pt->Fill(genP.pt());
        for(unsigned int i=0; i < genP.numberOfDaughters(); i++)
        {
          if(abs(genP.daughter(i)->pdgId()) == 13) genMuons.push_back(genP.daughter(i));
          if(abs(genP.daughter(i)->pdgId()) == 22) genPhotons.push_back(genP.daughter(i));
        }
      }
    }
  }

  hMu1VsMu2DeltaR->Fill(reco::deltaR(*(genMuons.at(0)), *(genMuons.at(1))));


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
    if (matched) hMuDeltaR->Fill(minDR);
    if (matched && minDR < 0.01)
    {
      recoMatchedMuons.push_back(bestMatchedMuon);
      genMatchedMuons.push_back(genMu);
      hRecoVsGenMuPt->Fill(genMu->pt(), bestMatchedMuon->pt());
      hMuPtError->Fill((bestMatchedMuon->pt() - genMu->pt())/genMu->pt());
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
    if (matched) hGammaDeltaR->Fill(minDR);
    if (matched && minDR < 0.03)
    {
      recoMatchedPhotons.push_back(bestMatchedPhoton);
      genMatchedPhotons.push_back(genPh);
      hRecoVsGenGammaPt->Fill(genPh->pt(), bestMatchedPhoton->pt());
      hGammaPtError->Fill((bestMatchedPhoton->pt() - genPh->pt())/genPh->pt());
    }
  }

  // if(recoMatchedMuons.size() == 2 && recoMatchedPhotons.size() == 1 && triggerFired)
  // {
  //   nTriggerMuMuG++;
  // }

  
  if(recoMatchedMuons.size() == 2 && recoMatchedPhotons.size() == 1)
  {
    if(recoMatchedMuons.at(0)->pt() > 4 && recoMatchedMuons.at(1)->pt() > 4 &&
    abs(recoMatchedMuons.at(0)->eta()) < 2.4 && abs(recoMatchedMuons.at(0)->eta()) < 2.4)
    {
      nPhotonsWithMuonCondition++;
      BsInvMass->Fill((recoMatchedMuons.at(0)->p4() + recoMatchedMuons.at(1)->p4() + recoMatchedPhotons.at(0)->p4()).M());
    }
  }


  // histograms for events with muon condition
  if(recoMatchedMuons.size() == 2)
  {
    if(recoMatchedMuons.at(0)->pt() > 4 && recoMatchedMuons.at(1)->pt() > 4 &&
    abs(recoMatchedMuons.at(0)->eta()) < 2.4 && abs(recoMatchedMuons.at(0)->eta()) < 2.4)
    {
      nMuonsWithCondition++;
      for (const auto genPh : genMatchedPhotons)
      {
        hMatchedGammaWithMuonConditionPt->Fill(genPh->pt());
        hMatchedGammaWithMuonConditionPtVsEta->Fill(genPh->pt(), genPh->eta());
      }

      for (const auto recoPh : recoMatchedPhotons)
      {
        hRecoGammaWithMuonConditionPt->Fill(recoPh->pt());
      }

      for (unsigned int i = 0; i < genMatchedPhotons.size(); i++)
      {
        hRecoVsGenGammaWithMuonConditionPt->Fill(genMatchedPhotons.at(i)->pt(), recoMatchedPhotons.at(i)->pt());
      }

      for (const auto genPh : genPhotons)
      {
        hGammaWithMuonConditionPt->Fill(genPh->pt());
        hGammaWithMuonConditionPtVsEta->Fill(genPh->pt(), genPh->eta());
      }

      for(unsigned int i = 0; i < recoMatchedMuons.size(); i++)
      {
        hRecoMuWithMuonConditionPt->Fill(recoMatchedMuons.at(i)->pt());
        hRecoVsGenMuPtWithMuonCondition->Fill(genMatchedMuons.at(i)->pt(), recoMatchedMuons.at(i)->pt());
      }
    }
  }

  for (const auto genMu : genMuons)
  {
    hGenMuPt->Fill(genMu->pt());
    hGenMuPtVsEta->Fill(genMu->pt(), genMu->eta());
  }

  for (const auto genPh : genPhotons)
  {
    hGenGammaPtVsEta->Fill(genPh->pt(), genPh->eta());
  }



  for (const reco::Muon* recoMatchedMu : recoMatchedMuons)
  {
    hMuPt->Fill(recoMatchedMu->pt());
    hRecoMuPtVsEta->Fill(recoMatchedMu->pt(), recoMatchedMu->eta());
  }

  for (const reco::Photon* recoMatchedPh : recoMatchedPhotons)
  {
    hGammaPt->Fill(recoMatchedPh->pt());
    hRecoGammaPtVsEta->Fill(recoMatchedPh->pt(), recoMatchedPh->eta());
  }


  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(RecoMuonAnalysis);

