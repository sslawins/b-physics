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
class ConvertedPhotons : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit ConvertedPhotons(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~ConvertedPhotons();

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

  // histograms
  TH1D* nConvertedPhotons;



  std::vector<int> MuMuG = {22, 13, -13};
};


ConvertedPhotons::ConvertedPhotons(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles"));
  theMuonToken = consumes< vector<reco::Muon>  >( edm::InputTag("muons"));
  thePhotonToken = consumes< vector<reco::Photon>  >( edm::InputTag("photons"));
  theConversionsToken = consumes< vector<pat::CompositeCandidate> >( edm::InputTag("oniaPhotonCandidates","conversions"));
  m_fieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
}

ConvertedPhotons::~ConvertedPhotons()
{
  cout <<" DTOR" << endl;
}

bool ConvertedPhotons::isSameDecay(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}


void ConvertedPhotons::beginJob()
{
  //create a histogram
  nConvertedPhotons = new TH1D("nConvertedPhotons", "nConvertedPhotons", 10, 0, 10);

  cout << "HERE ConvertedPhotons::beginJob()" << endl;
}

void ConvertedPhotons::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  nConvertedPhotons->Write();

  myRootFile.Close();

  delete nConvertedPhotons;

  cout << "HERE ConvertedPhotons::endJob()" << endl;
}


void ConvertedPhotons::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE ConvertedPhotons::analyze "<< std::endl;

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);
  const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);
  const pat::CompositeCandidateCollection * conversions = &(ev.get(theConversionsToken));
  auto const& field = es.getData(m_fieldToken);

  nConvertedPhotons->Fill(conversions->size());

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
  }

  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(ConvertedPhotons);

