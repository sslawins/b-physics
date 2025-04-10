#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/Vector4D.h"

#include <vector>

#include "interface/HLTdecision.h"
#include "interface/DecayTools.h"

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

  template <typename T>
  double invariantMass(const std::vector<T>& particles, const std::vector<double>& masses) {
    if (particles.size() != masses.size()) {
        throw std::invalid_argument("Size mismatch: particles and masses vectors must have the same length");
    }

    ROOT::Math::PxPyPzEVector totalP4(0, 0, 0, 0);

    for (size_t i = 0; i < particles.size(); ++i) {
        ROOT::Math::PxPyPzEVector p4(particles[i].px(), particles[i].py(), particles[i].pz(),
                                     std::sqrt(masses[i] * masses[i] + particles[i].p() * particles[i].p()));
        totalP4 += p4;
    }

    return totalP4.M();
  }

  
private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;


  //tokens
  edm::EDGetTokenT < vector<pat::Muon> > theMuonToken;
  edm::EDGetTokenT < vector<pat::Photon> > thePhotonToken;
  edm::EDGetTokenT < vector<pat::PackedCandidate> > theCandidateToken;
  edm::EDGetTokenT<std::vector<reco::Vertex>> thePVToken;
  edm::EDGetTokenT<std::vector<reco::Vertex>> thePVWithBSToken;
  edm::EDGetTokenT < edm::TriggerResults > theTriggerResultsToken;
  
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTrackBuilderToken;

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
DataAnalysis::DataAnalysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theMuonToken = consumes< vector<pat::Muon>  >( edm::InputTag("slimmedDisplacedMuons"));
  theCandidateToken = consumes< vector<pat::PackedCandidate> >  ( edm::InputTag("packedPFCandidates"));
  thePhotonToken = consumes< vector<pat::Photon>  >( edm::InputTag("slimmedPhotons"));
  thePVToken = consumes<std::vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  thePVWithBSToken = consumes<std::vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVerticesWithBS"));
  theTriggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
  theTrackBuilderToken = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));
}
DataAnalysis::~DataAnalysis()
{
  cout <<" DTOR" << endl;
}

void DataAnalysis::beginJob()
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

  cout << "HERE DataAnalysis::beginJob()" << endl;
}

void DataAnalysis::endJob()
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

  cout << "HERE DataAnalysis::endJob()" << endl;
}

void DataAnalysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{

  
  std::cout << " -------------------------------- HERE DataAnalysis::analyze "<< std::endl;

  
  //const std::vector<reco::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<pat::PackedCandidate> & candidates = ev.get(theCandidateToken);
  const std::vector<pat::Photon> & recoPhotons = ev.get(thePhotonToken);
  const std::vector<pat::Muon> & recoMuons = ev.get(theMuonToken);
  const std::vector<reco::Vertex>& recoVertices = ev.get(thePVToken);
  const std::vector<reco::Vertex>& recoVerticesWithBS = ev.get(thePVWithBSToken);
  
  std::cout << "PFCandidates: " << candidates.size() << std::endl;
  std::cout << "Photons: " << recoPhotons.size() << std::endl;
  std::cout << "Muons: " << recoMuons.size() << std::endl;
  std::cout << "Primary Vertices: " << recoVertices.size() << std::endl;
  std::cout << "Primary Vertices with Beam Spot: " << recoVerticesWithBS.size() << std::endl;

  //if(recoPhotons == 0) return;

  for( unsigned int vtxIter = 0; vtxIter < recoVertices.size(); ++vtxIter){
    const reco::Vertex & vtx = recoVertices[vtxIter];
    cout << "Primary Vertex: " << vtx.position().x() << " " << vtx.position().y() << " " << vtx.position().z() << endl;
  }
  cout << "//////////////////////////////////////////////////" << endl;
  for( unsigned int vtxIter = 0; vtxIter < recoVerticesWithBS.size(); ++vtxIter){
    const reco::Vertex & vtx = recoVerticesWithBS[vtxIter];
    cout << "Primary Vertex: " << vtx.position().x() << " " << vtx.position().y() << " " << vtx.position().z() << endl;
  }
  double deltaRKaons = 0.03, deltaRPhotonVsKaons = 0.03;
  double deltaVzKaons = 0.5;
  double deltaMassPhi = 0.01;
  double deltaMassBs = 0.01;

  //double Kmass = 0.493677;
  //double PhiMass = 1.019455;
  //double BsMass = 5.36688;
  
  const auto & trackBuilder = es.getData(theTrackBuilderToken); 
  std::vector<reco::TransientTrack> kaonTTs;

  //Loop over the kaons
  for( std:: vector<pat::PackedCandidate>::const_iterator ic1 = candidates.begin(); ic1 < candidates.end(); ic1++){

    if( !ic1->hasTrackDetails() ) continue; // skip if a bestTrack can be extracted from this Candidate
    if( abs(ic1->pdgId()) ) continue; //K+/-
    //if( ic1->pt() < 2. ) continue; 
    

    for ( std:: vector<pat::PackedCandidate>::const_iterator ic2 = ic1+1; ic2 < candidates.end(); ic2++){
      if( !ic1->hasTrackDetails() ) continue; // skip if a bestTrack can be extracted from this Candidate
      if( abs(ic1->pdgId()) ) continue; //K+/-
      //if( ic1->pt() < 2. ) continue; 
      
      const reco::Track & rtrk1 = ic1 -> pseudoTrack();
      const reco::Track & rtrk2 = ic2 -> pseudoTrack();
      
      double vzBsPhi = (rtrk1.vz() + rtrk2.vz())/2;

      if( rtrk1.charge()*rtrk2.charge() != -1 ) continue;
      if( fabs(rtrk1.vz() - rtrk2.vz()) > deltaVzKaons) continue;
      if( reco::deltaR(rtrk1, rtrk2) > deltaRKaons) continue;
      //if( fabs(vzBsPhi))

    }

  }

  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(DataAnalysis);

