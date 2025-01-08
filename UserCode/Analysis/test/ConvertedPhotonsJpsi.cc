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
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


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
class ConvertedPhotonsJpsi : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit ConvertedPhotonsJpsi(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~ConvertedPhotonsJpsi();

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
  edm::EDGetTokenT < vector<pat::CompositeCandidate> > theConversionsToken;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_fieldToken;
  edm::EDGetTokenT < edm::TriggerResults > theTriggerResultsToken;

  // histograms
  TH1D* nConvertedPhotons;

  TH1D* hBsMass;
  TH1D* hMuPt;
  TH1D* hGammaPt;
  TH1D* hGammaDeltaR;
  TH1D* hGammaPtWithTrigger;


  int nConvPhotons = 0;
  std::vector<int> MuMuG = {22, 13, -13};
};


ConvertedPhotonsJpsi::ConvertedPhotonsJpsi(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles"));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));
  theConversionsToken = consumes< vector<pat::CompositeCandidate> >( edm::InputTag("oniaPhotonCandidates","conversions"));
  m_fieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  theTriggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));

}

ConvertedPhotonsJpsi::~ConvertedPhotonsJpsi()
{
  cout <<" DTOR" << endl;
}

bool ConvertedPhotonsJpsi::isSameDecay(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}


void ConvertedPhotonsJpsi::beginJob()
{
  //create a histogram
  nConvertedPhotons = new TH1D("nConvertedPhotons", "nConvertedPhotons", 10, 0, 10);

  hBsMass = new TH1D("hBsMass", "hBsMass", 50, 3, 7);

  hMuPt = new TH1D("hMuPt", "hMuPt", 100, 0, 30);
  hGammaPt = new TH1D("hGammaPt", "hGammaPt", 100, 0, 30);
  hGammaDeltaR = new TH1D("hGammaDeltaR", "hGammaDeltaR", 100, 0, 0.05);
  hGammaPtWithTrigger = new TH1D("hGammaPtWithTrigger", "hGammaPtWithTrigger", 100, 0, 30);

  cout << "HERE ConvertedPhotonsJpsi::beginJob()" << endl;
}

void ConvertedPhotonsJpsi::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  nConvertedPhotons->Write();

  hBsMass->Write();

  hMuPt->Write();
  hGammaPt->Write();
  hGammaDeltaR->Write();
  hGammaPtWithTrigger->Write();

  myRootFile.Close();

  delete nConvertedPhotons;
  delete hBsMass;
  delete hMuPt;
  delete hGammaPt;
  delete hGammaDeltaR;
  delete hGammaPtWithTrigger;

  cout << "nConvertedPhotons: " << nConvPhotons << endl;


  cout << "HERE ConvertedPhotonsJpsi::endJob()" << endl;
}


void ConvertedPhotonsJpsi::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE ConvertedPhotonsJpsi::analyze "<< std::endl;

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);
  const pat::CompositeCandidateCollection * conversions = &(ev.get(theConversionsToken));
  auto const& field = es.getData(m_fieldToken);

  const edm::TriggerResults & triggerResults = ev.get(theTriggerResultsToken);
  edm::TriggerNames triggerNames = ev.triggerNames(triggerResults);


  vector<const reco::Candidate*> genMuons;
  vector<const reco::Muon*> recoMatchedMuons;
  vector<const reco::Candidate*> genMatchedMuons;

  vector<const reco::Candidate*> genPhotons;
  vector<RefCountedKinematicParticle> recoMatchedPhotons;
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

  vector<RefCountedKinematicParticle> convPhotons;
  for (pat::CompositeCandidateCollection::const_iterator conv = conversions->begin(); conv!= conversions->end(); ++conv) 
  {
    reco::TrackCollection convTracks;
    const reco::Track tk0=*conv->userData<reco::Track>("track0");
    const reco::Track tk1=*conv->userData<reco::Track>("track1");
    convTracks.push_back(tk0);
    convTracks.push_back(tk1);
    std::vector<reco::TransientTrack> EETT;
    reco::TransientTrack EETT1(convTracks[0], &field );
    reco::TransientTrack EETT2(convTracks[1], &field );
    EETT.push_back(EETT1);
    EETT.push_back(EETT2);
    const ParticleMass zero_mass(0);
    const ParticleMass PM_PDG_ELE_MASS(0.000511);
    // const ParticleMass zero_mass(2*PDG_ELE_MASS);
    float zero_sigma = 1E-6;
    float eleSigma = 1E-6;
    //
    std::vector<RefCountedKinematicParticle> PhotonParticles;
    KinematicParticleFactoryFromTransientTrack pFactory;
    PhotonParticles.push_back(pFactory.particle(EETT[0], PM_PDG_ELE_MASS, float(0), float(0), eleSigma));
    PhotonParticles.push_back(pFactory.particle(EETT[1], PM_PDG_ELE_MASS, float(0), float(0), eleSigma));
    KinematicParticleVertexFitter photonfitter;
    RefCountedKinematicTree photonVertexFitTree;
    photonVertexFitTree = photonfitter.fit(PhotonParticles);
    if (!photonVertexFitTree->isValid()) 
    {
      edm::ParameterSet pSet;
      pSet.addParameter<double>("maxDistance", 3);
      pSet.addParameter<int>("maxNbrOfIterations", 10000);
      KinematicParticleVertexFitter photonfitter2(pSet);
      photonVertexFitTree = photonfitter2.fit(PhotonParticles);
    }
    //
    if (!photonVertexFitTree->isValid()) continue;
    //
    KinematicParticleFitter csFitterPhoton;
    KinematicConstraint * pho_c = new MassKinematicConstraint(zero_mass, zero_sigma);
    // add mass constraint to the photon fit to do a constrained fit:
    photonVertexFitTree->movePointerToTheTop();
    photonVertexFitTree = csFitterPhoton.fit(pho_c, photonVertexFitTree);
    if (!photonVertexFitTree->isValid()) continue;

    RefCountedKinematicParticle fitPhoton = photonVertexFitTree->currentParticle();
    nConvPhotons++;
    convPhotons.push_back(fitPhoton);
  }

  nConvertedPhotons->Fill(convPhotons.size());

  // converted photon matching
  for (const reco::Candidate* genPh : genPhotons)
  {
    float minDR = 10;
    RefCountedKinematicParticle bestMatchedPhoton;
    bool matched = false;
    for (const auto& convPh : convPhotons)
    {
      float dR = reco::deltaR(convPh->currentState().globalMomentum().eta(), convPh->currentState().globalMomentum().phi(), genPh->eta(), genPh->phi());
      if (dR < minDR)
      {
        minDR = dR;
        bestMatchedPhoton = convPh;
        matched = true;
      }
    }
    if (matched) hGammaDeltaR->Fill(minDR);
    if (matched && minDR < 0.02)
    {
      recoMatchedPhotons.push_back(bestMatchedPhoton);
      genMatchedPhotons.push_back(genPh);
      // hRecoVsGenGammaPt->Fill(genPh->pt(), bestMatchedPhoton->pt());
      // hGammaPtError->Fill((bestMatchedPhoton->pt() - genPh->pt())/genPh->pt());
      hGammaPt->Fill(bestMatchedPhoton->currentState().globalMomentum().perp());
    }
  }

  for (unsigned int i = 0; i < triggerResults.size(); i++)
  {
    TString name = triggerNames.triggerName(i);
    if (name == "HLT_DoubleMu4_3_LowMass_v1" && triggerResults.accept(i) == 1 && recoMatchedPhotons.size() > 0)
    {
      hGammaPtWithTrigger->Fill(recoMatchedPhotons.at(0)->currentState().globalMomentum().perp());
    }
  }
  
  
  vector<RefCountedKinematicParticle> muonKinematicParticles;
  for(const auto& recoMu : recoMuons)
  {
    hMuPt->Fill(recoMu.pt());
    reco::TrackRef muTrack = recoMu.track();
    if(!muTrack) continue;
    reco::TransientTrack muonTT = reco::TransientTrack(muTrack, &field);

    const ParticleMass muon_mass(0.105658);
    float muon_sigma = 1E-6;

    KinematicParticleFactoryFromTransientTrack pFactory;
    muonKinematicParticles.push_back(pFactory.particle(muonTT, muon_mass, float(0), float(0), muon_sigma));
  }

  for (unsigned int i = 0; i < muonKinematicParticles.size(); i++)
  {
    for (unsigned int j = i+1; j < muonKinematicParticles.size(); j++)
    {
      for (unsigned int k = 0; k < convPhotons.size(); k++)
      {
        for (unsigned int l = k+1; l < convPhotons.size(); l++)
        {
          RefCountedKinematicParticle mu1 = muonKinematicParticles.at(i);
          RefCountedKinematicParticle mu2 = muonKinematicParticles.at(j);
          RefCountedKinematicParticle pho1 = convPhotons.at(k);
          RefCountedKinematicParticle pho2 = convPhotons.at(l);

          std::vector<RefCountedKinematicParticle> allParticlesMu;
          allParticlesMu.push_back(mu1);
          allParticlesMu.push_back(mu2);

          const ParticleMass jpsi_mass = 3.096916;
          float jpsiMsigma = 0.00004;


          KinematicParticleVertexFitter fitter;
          RefCountedKinematicTree jpsiFitTree = fitter.fit(allParticlesMu);

          if (!jpsiFitTree->isValid()) continue;

          KinematicConstraint* jpsi_mass_constraint = new MassKinematicConstraint(jpsi_mass, jpsiMsigma);
          KinematicParticleFitter csFitterJpsi;

          jpsiFitTree->movePointerToTheTop();
          jpsiFitTree = csFitterJpsi.fit(jpsi_mass_constraint, jpsiFitTree);

          if (!jpsiFitTree->isValid()) continue;

          jpsiFitTree->movePointerToTheTop();
          RefCountedKinematicParticle fitJpsi = jpsiFitTree->currentParticle();

          if (!fitJpsi->currentState().isValid()) continue;

          vector<RefCountedKinematicParticle> allParticles;
          allParticles.push_back(fitJpsi);
          allParticles.push_back(pho1);
          allParticles.push_back(pho2);

          KinematicParticleVertexFitter vertexFitter;
          RefCountedKinematicTree vertexFitTree = vertexFitter.fit(allParticles);

          if (!vertexFitTree->isValid()) continue;

          vertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle fitParticle = vertexFitTree->currentParticle();
          RefCountedKinematicVertex fitVertex = vertexFitTree->currentDecayVertex();

          if (!fitVertex->vertexIsValid()) continue;

          hBsMass->Fill(fitParticle->currentState().mass());
        }
      }
    }
  }

  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(ConvertedPhotonsJpsi);

