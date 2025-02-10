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
#include <vector>
#include <stack>

#include "BDecayAnalyzer.h"


using namespace std;


//object definition
class PVdistance : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit PVdistance(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~PVdistance();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  //Bs and phi decay - vectors for the products
  vector<const reco::Candidate*> genMuons;
  vector<const reco::Muon*> recoMatchedMuons;
  vector<const reco::Muon*> recoMatchedMuonsOtherSide;
  vector<const reco::Candidate*> genPhotons;
  vector<const reco::Photon*> recoMatchedPhotons;

  //tokens
  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;
  edm::EDGetTokenT < vector<reco::Muon> > theMuonToken;
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;
  //edm::EDGetTokenT < vector<pat::CompositeCandidate> > theConversionsToken;
  //edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_fieldToken;
  edm::EDGetTokenT < edm::TriggerResults > theTriggerResultsToken;

};

//constructor
PVdistance::PVdistance(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles" ));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));
  //theConversionsToken = consumes< vector<pat::CompositeCandidate> >( edm::InputTag("oniaPhotonCandidates","conversions"));
  //m_fieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  theTriggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
}

//destructor
PVdistance::~PVdistance()
{
  cout <<" DTOR" << endl;
}


void PVdistance::beginJob()
{

  cout << "HERE PVdistance::beginJob()" << endl;
}

void PVdistance::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data

  myRootFile.Close();


  cout << "HERE PVdistance::endJob()" << endl;
}


void PVdistance::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE PVdistance::analyze "<< std::endl;

  genMuons.clear();
  recoMatchedMuons.clear();
  //genMuonsOtherSide.clear(); ----> this is exactly what cascadeMuons is
  recoMatchedMuonsOtherSide.clear();
  genPhotons.clear();
  recoMatchedPhotons.clear();

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);
  
  //const pat::CompositeCandidateCollection * conversions = &(ev.get(theConversionsToken));
  //auto const& field = es.getData(m_fieldToken);
  const edm::TriggerResults & triggerResults = ev.get(theTriggerResultsToken);
  edm::TriggerNames triggerNames = ev.triggerNames(triggerResults);

  cout <<"Number of triggers:   " <<triggerNames.size() << endl;
  for (unsigned int trgIter = 0; trgIter < triggerNames.size(); ++trgIter ){
    cout << triggerNames.triggerName(trgIter) << endl;
  }
  

  BDecayAnalyzer bAnalyzer;
  std::vector<std::vector<const reco::Candidate*>> tree = bAnalyzer.analyzeBDecays(genPar);
  //std::vector<const reco::Candidate*> cascadeMuons;
  genPhotons = bAnalyzer.getPhotons();
  genMuons = bAnalyzer.getMuons();
  //if(genMuons.size() == 2 && genPhotons.size() == 1) cout << genMuons[0]->pt() << "  " << genPhotons[0]->pt() << endl;
  
  if (genMuons.size() != 2 || genPhotons.size() != 1) {
    std::cout << "Skipping event: " << ev.id() << " (condition not met)" << std::endl;
    return;  
  }
  
  // print the family tree
  std::cout <<"Family tree: " << std::endl;
  for (const auto& lineage : tree) { 
    //if( abs(lineage.back()->pdgId()) == 531 ) nBs++;
    for (const auto* particle : lineage) { 
        std::cout << particle->pdgId() << " "; 
    }
    std::cout << std::endl; 
  }

  
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(PVdistance);

