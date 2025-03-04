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
#include "HLTdecision.h"
#include "MuonMatcher.h"


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
  vector<const reco::Candidate*> genPhotons;

  //tokens
  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;
  edm::EDGetTokenT < vector<reco::Muon> > theMuonToken;
  edm::EDGetTokenT < vector<reco::Photon> > thePhotonToken;
  //edm::EDGetTokenT < vector<pat::CompositeCandidate> > theConversionsToken;
  //edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_fieldToken;
  edm::EDGetTokenT < edm::TriggerResults > theTriggerResultsToken;
  edm::EDGetTokenT < vector<reco::Vertex> > thePVToken;

  TH1D *hPhiMass;
  TH1D *hPCAz;
  TH1D *hPvSvDistance;

};

//constructor
PVdistance::PVdistance(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles" ));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  //thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));
  //theConversionsToken = consumes< vector<pat::CompositeCandidate> >( edm::InputTag("oniaPhotonCandidates","conversions"));
  //m_fieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  theTriggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
  thePVToken = consumes< vector<reco::Vertex>  >( edm::InputTag("offlinePrimaryVertices"));
}

//destructor
PVdistance::~PVdistance()
{
  cout <<" DTOR" << endl;
}

void PVdistance::beginJob()
{

  hPhiMass = new TH1D("hPhiMass", "Reconstruction of #Phi ; M_{inv} [GeV]; Counts", 48 , 0.9, 1.2);
  hPCAz = new TH1D("hPCAz", "Distance between the point of closest approach and PV, z axis", 200, -10, 10);
  hPvSvDistance = new TH1D("hPvSvDistance", "Distance between PV and SV", 200, -10, 10);

  cout << "HERE PVdistance::beginJob()" << endl;
}

void PVdistance::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hPhiMass -> Write();
  hPCAz    -> Write();
  hPvSvDistance -> Write();

  myRootFile.Close();

  delete hPhiMass;
  delete hPCAz;
  delete hPvSvDistance;

  cout << "HERE PVdistance::endJob()" << endl;
}


void PVdistance::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{

  
  std::cout << " -------------------------------- HERE PVdistance::analyze "<< std::endl;

  genMuons.clear();
  genPhotons.clear();

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  //const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);
  //const pat::CompositeCandidateCollection * conversions = &(ev.get(theConversionsToken));
  //auto const& field = es.getData(m_fieldToken);
  
  HLTdecision HLTdecision(theTriggerResultsToken, ev, theConfig);
  bool accepted = HLTdecision.checkTriggers(ev, true); //print = true
  if(!accepted) return ;

  BDecayAnalyzer bAnalyzer;
  std::vector<std::vector<const reco::Candidate*>> tree = bAnalyzer.analyzeBDecays(genPar);
  genPhotons = bAnalyzer.getPhotons();
  genMuons = bAnalyzer.getMuons();
  if(! bAnalyzer.analyzeEvent(ev.id())) return;
  // print the family tree
  bAnalyzer.printTheTree(tree);

  MuonMatcher muonMatcher( recoMuons, genMuons, 0.1);
  muonMatcher.matchRecoToGen();
  muonMatcher.printMatchedMuons();
  std::cout << "Matching successful: " << muonMatcher.isSuccessful() << std::endl;
  std::vector< const reco::Candidate*> matchedMuons = muonMatcher.getMatched();

  for( const auto& muon : matchedMuons){
    std::cout << "Matched muons' SV: (" << muon->vertex().x() 
              << " , " << muon->vertex().Y()
              << " , " << muon->vertex().Z() << ")"
              << std::endl;
  }

  //ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>>  -->  math::XYZPoint
  math::XYZPoint sv = matchedMuons[0]->vertex();
  math::XYZPoint svG0 = genMuons[0]->vertex();
  math::XYZPoint svG1 = genMuons[1]->vertex();
  std::vector<math::XYZPoint> vertices = {sv,svG0, svG1};
  for(const auto& vertex : vertices){
      std::cout << "Matched muons' SV: (" << vertex.x() 
              << " , " << vertex.Y()
              << " , " << vertex.Z() << ")"
              << std::endl;
  }

  math::XYZVectorD pMuMu = matchedMuons[0]->momentum() + matchedMuons[1]->momentum();

  //ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > p_MuMu(0.,0.,0.);

  ////////////////// PV /////////////////
  const std::vector<reco::Vertex> & PVertices = ev.get(thePVToken);  

  std::cout << "Primary vertices: " << std::endl;
  for( const auto& vertex : PVertices){
    std::cout << vertex.position().Theta() << "  " << vertex.y() << "  " << vertex.z() << "  " << std::endl;
  }
  math::XYZPoint pv = PVertices[0].position();
  math::XYZPoint pca = sv;
  double s = ((pv - sv).Dot(pMuMu)) / pMuMu.Mag2();
  math::XYZVectorD y = s*pMuMu;
  pca += y;
  double minDistance = sqrt((pca - pv).Mag2());
  std::cout << "PV: (" << pv.X() << ", " << pv.Y() << ", " << pv.Z() << ")" << std::endl;
  std::cout << "PCA: (" << pca.X() << ", " << pca.Y() << ", " << pca.Z() << ")" << std::endl;
  std::cout << "Minimalna odległość: " << minDistance << std::endl;
  //cout <<"Number of triggers:   " <<triggerNames.size() << endl;
 
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(PVdistance);

