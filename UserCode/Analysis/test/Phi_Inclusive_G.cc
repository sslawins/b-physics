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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

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
#include <vector>
#include <stack>

#include "interface/BDecayAnalyzer.h"
#include "interface/HLTdecision.h"
#include "interface/MuonMatcher.h"
#include "interface/DecayTools.h"
#include "interface/ParticleMatcher.h"

using namespace std;


//object definition
class Phi_Inclusive_G : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit Phi_Inclusive_G(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~Phi_Inclusive_G();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  //Bs and phi decay - vectors for the products
  vector<const reco::Candidate*> genKaons; 
  vector<const reco::Candidate*> genPhotons;

  //tokens
  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;
  edm::EDGetTokenT < vector<reco::Muon> > theMuonToken;
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;
  edm::EDGetTokenT < vector<reco::PFCandidate> > theCandidateToken;
  //edm::EDGetTokenT < vector<pat::CompositeCandidate> > theConversionsToken;
  //edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_fieldToken;
  edm::EDGetTokenT < edm::TriggerResults > theTriggerResultsToken;
  //edm::EDGetTokenT < vector<reco::Vertex> > thePVToken;

  TH1D *hHLT;
  TH1D *hpTrecoMuons;
  TH1D *hPhiToKplusKminus;
  TH1D *hHLT_phiToKK;
  TH1D *hHLT_gamma;
  TH1D *hHLT_phiToKK_gamma;
  TH2D *hKKpT_recoVSgen;
  TH1D *hKK_deltaR;
  TH1D *hDelta_vz;

};

//constructor
Phi_Inclusive_G::Phi_Inclusive_G(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles" ));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  theCandidateToken = consumes< vector<reco::PFCandidate>  >( edm::InputTag("particleFlow"));
  thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));
  //theConversionsToken = consumes< vector<pat::CompositeCandidate> >( edm::InputTag("oniaPhotonCandidates","conversions"));
  //m_fieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  theTriggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
  //thePVToken = consumes< vector<reco::Vertex>  >( edm::InputTag("offlinePrimaryVertices"));
}

//destructor
Phi_Inclusive_G::~Phi_Inclusive_G()
{
  cout <<" DTOR" << endl;
}

void Phi_Inclusive_G::beginJob()
{
  hHLT = new TH1D("hHLT", "HLT decision",3,  -0.5, 2.5);
  hpTrecoMuons = new TH1D("hpTrecoMuons", "reco muon pt", 1000, 0, 50);
  hPhiToKplusKminus = new TH1D("hPhiToKplusKminus", "Phi to K+K- decay", 2, -0.5, 1.5);
  hHLT_phiToKK = new TH1D("hHLT_phiToKK", "HLT decision for Phi to K+K- decay", 3, -0.5, 2.5);
  hHLT_gamma = new TH1D("hHLT_gamma", "HLT decision for gamma", 3, -0.5, 2.5);
  hHLT_phiToKK_gamma = new TH1D("hHLT_phiToKK_gamma", "HLT decision for Phi to K+K- decay with gamma", 3, -0.5, 2.5);
  hKKpT_recoVSgen = new TH2D("hKKpT_recoVSgen", "K+K- reco-gen", 700, 0, 70, 100, 0, 70);
  hKK_deltaR = new TH1D("hKK_deltaR", "K+K- deltaR", 100, 0, 0.5);
  hDelta_vz = new TH1D("hDelta_vz", "delta vz", 500000, 0, 5);

  cout << "HERE Phi_Inclusive_G::beginJob()" << endl;
}

void Phi_Inclusive_G::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hHLT -> Write();
  hpTrecoMuons -> Write();
  hPhiToKplusKminus -> Write();
  hHLT_phiToKK -> Write();
  hHLT_gamma -> Write();
  hHLT_phiToKK_gamma -> Write();
  hKKpT_recoVSgen -> Write();
  hKK_deltaR -> Write();  
  hDelta_vz -> Write();

  myRootFile.Close();

  delete hHLT;
  delete hpTrecoMuons;
  delete hPhiToKplusKminus;
  delete hHLT_phiToKK;
  delete hHLT_gamma;
  delete hHLT_phiToKK_gamma;
  delete hKKpT_recoVSgen;
  delete hKK_deltaR;
  delete hDelta_vz;

  cout << "HERE Phi_Inclusive_G::endJob()" << endl;
}


void Phi_Inclusive_G::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{

  
  std::cout << " -------------------------------- HERE Phi_Inclusive_G::analyze "<< std::endl;

  genKaons.clear();
  genPhotons.clear();

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  //const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::PFCandidate> & recoCandidates = ev.get(theCandidateToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);
  
  std::cout << "PFCandidates: " << recoCandidates.size() << std::endl;
  
  BDecayAnalyzer bAnalyzer;
  std::vector<std::vector<const reco::Candidate*>> tree = bAnalyzer.analyzeBDecays(genPar, DecayTools::KK);
  
  genPhotons = bAnalyzer.getPhotons();
  genKaons = bAnalyzer.getPhiProducts();
  std::cout << "Phi products: " << genKaons.size() << std::endl;
  if( bAnalyzer.analyzeEvent(ev.id())){
    hPhiToKplusKminus->Fill(1);
    
    // print the family tree
    //bAnalyzer.printTheTree(tree);
    std::cout << "Number of kaons: " << genKaons.size() << std::endl;
    //Kaon matching

    ParticleMatcher<reco::PFCandidate> particleMatcher( recoCandidates, genKaons, 0.01, hKK_deltaR);
    particleMatcher.matchRecoToGen();
    particleMatcher.printMatchedParticles();
    std::cout << "Matching successful: " << particleMatcher.isSuccessful() << std::endl;
    std::vector< const reco::Candidate*> matchedParticles = particleMatcher.getMatched();
    if( particleMatcher.isSuccessful() ) {
      hKKpT_recoVSgen -> Fill(matchedParticles[0]->pt(), genKaons[0]->pt());
      hKKpT_recoVSgen -> Fill(matchedParticles[1]->pt(), genKaons[1]->pt());
      hDelta_vz -> Fill( abs(matchedParticles[0]->vz() - matchedParticles[1]->vz()));
    };

    //Photon matching
    ParticleMatcher<reco::Photon> photonMatcher( recoPhotons, genPhotons, 0.03);
    photonMatcher.matchRecoToGen();
    photonMatcher.printMatchedParticles();
    std::cout << "Matching successful: " << photonMatcher.isSuccessful() << std::endl;
    std::vector< const reco::Candidate*> matchedPhotons = photonMatcher.getMatched();
    
    HLTdecision HLT(theTriggerResultsToken, ev, theConfig);
    bool accepted = HLT.checkTriggers(ev, true); //print = true
    
    if(!accepted) hHLT->Fill(0);
    else{
      hHLT->Fill(1);
      if( particleMatcher.isSuccessful() && photonMatcher.isSuccessful()) hHLT_phiToKK_gamma->Fill(1);
      if (particleMatcher.isSuccessful()) hHLT_phiToKK->Fill(1);
      if (photonMatcher.isSuccessful()) hHLT_gamma->Fill(1);

    }
  }else hPhiToKplusKminus->Fill(0);

  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(Phi_Inclusive_G);

