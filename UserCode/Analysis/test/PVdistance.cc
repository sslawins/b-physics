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

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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
#include "DecayTools.h"


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

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTrackBuilderToken;

  TH1D *hPhiMass;
  TH1D *hPCAz_reco0SV;
  TH1D *hPCAz;
  TH1D *hPCAz_genSV;
  TH1D *hPvSvDistance;
  TH1D *hMuMu_vz;

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

  theTrackBuilderToken = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));
}

//destructor
PVdistance::~PVdistance()
{
  cout <<" DTOR" << endl;
}

void PVdistance::beginJob()
{

  hPhiMass = new TH1D("hPhiMass", "Reconstruction of #Phi ; M_{inv} [GeV]; Counts", 48 , 0.9, 1.2);
  hPCAz_reco0SV = new TH1D("hPCAz_reco0SV", "Distance between the point of closest approach and PV, z axis", 200, 0., 0.2);
  hPCAz = new TH1D("hPCAz", "Distance between the point of closest approach and PV, z axis", 1000, 0., 1.);
  hPCAz_genSV = new TH1D("hPCAz_genSV", "Distance between the point of closest approach and PV, z axis", 200, 0., 0.2);
  hPvSvDistance = new TH1D("hPvSvDistance", "Distance between PV and SV", 4000, 0.,4.);
  hMuMu_vz = new TH1D("hMuMu_vz", "|delta vz|", 500000, 0.,5.);

  cout << "HERE PVdistance::beginJob()" << endl;
}

void PVdistance::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hPhiMass -> Write();
  hPCAz_genSV    -> Write();
  hPCAz    -> Write();
  hPCAz_reco0SV -> Write();
  hPvSvDistance -> Write();
  hMuMu_vz -> Write();

  myRootFile.Close();

  delete hPhiMass;
  delete hPCAz_genSV;
  delete hPCAz;
  delete hPCAz_reco0SV;
  delete hPvSvDistance;
  delete hMuMu_vz;

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
  
  const auto & trackBuilder = es.getData(theTrackBuilderToken);

  HLTdecision HLTdecision(theTriggerResultsToken, ev, theConfig);
  bool accepted = HLTdecision.checkTriggers(ev, true); //print = true
  if(!accepted) return ;

  BDecayAnalyzer bAnalyzer;
  std::vector<std::vector<const reco::Candidate*>> tree = bAnalyzer.analyzeBDecays(genPar, DecayTools::MuMu);
  genPhotons = bAnalyzer.getPhotons();
  genMuons = bAnalyzer.getPhiProducts();
  if(! bAnalyzer.analyzeEvent(ev.id())) return;
  // print the family tree
  bAnalyzer.printTheTree(tree);

  MuonMatcher muonMatcher( recoMuons, genMuons, 0.1);
  muonMatcher.matchRecoToGen();
  muonMatcher.printMatchedMuons();
  std::cout << "Matching successful: " << muonMatcher.isSuccessful() << std::endl;
  if(!muonMatcher.isSuccessful()) return;
  
  std::vector< const reco::Candidate*> matchedMuons = muonMatcher.getMatched();
  hMuMu_vz ->Fill(abs(matchedMuons[0]->vz() - matchedMuons[1]->vz()));
  for( const auto& muon : matchedMuons){
    std::cout << "Matched muons' SV: (" << muon->vertex().x() 
              << " , " << muon->vertex().Y()
              << " , " << muon->vertex().Z() << ")"
              << std::endl;
  }

  //ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>>  -->  math::XYZPoint
  math::XYZPoint sv = matchedMuons[0]->vertex();
  math::XYZPoint svG0 = genMuons[0]->vertex();
  //math::XYZPoint svG1 = genMuons[1]->vertex();
  //std::vector<math::XYZPoint> vertices = {sv,svG0, svG1};

  math::XYZVectorD pMuMu = matchedMuons[0]->momentum() + matchedMuons[1]->momentum();

  //ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > p_MuMu(0.,0.,0.);

  ////////////////// PV /////////////////
  const std::vector<reco::Vertex> & PVertices = ev.get(thePVToken);  

  std::cout << "Primary vertices: " << std::endl;
  for( const auto& vertex : PVertices){
    std::cout << vertex.position().x() << "  " << vertex.y() << "  " << vertex.z() << "  " << std::endl;
  }

  //PV - the first out of the list
  math::XYZPoint pv = PVertices[0].position();
  //s = (PV - SV) * (pMuMu) / |pMuMu|^2 
  
  double s = ((pv - sv).Dot(pMuMu)) / pMuMu.Mag2();
  //pca = sv + s*pMuMu
  math::XYZPoint pca = sv;
  pca += s*pMuMu;

  //double minDistance = sqrt((pca - pv).Mag2());
  std::cout << "PV: (" << pv.X() << ", " << pv.Y() << ", " << pv.Z() << ")" << std::endl;
  std::cout << "PCA: (" << pca.X() << ", " << pca.Y() << ", " << pca.Z() << ")" << std::endl;
  std::cout << "(PCA - PV)_z: " << abs(pca.z() - pv.z()) << std::endl;
  hPCAz_reco0SV->Fill(abs(pca.z() - pv.z()));
  hPvSvDistance->Fill(abs(pv.z() - sv.z()));
  
  double s_gen = ((pv - svG0).Dot(pMuMu)) / pMuMu.Mag2();
  //pca = sv + s*pMuMu
  pca = svG0;
  pca += s_gen*pMuMu;

  hPCAz_genSV -> Fill(abs(pca.z() - pv.z()));
  //cout <<"Number of triggers:   " <<triggerNames.size() << endl;

  vector<reco::TransientTrack> muonTTs;
  for (const auto& mu : matchedMuons)
  {
    const reco::Muon* muon = dynamic_cast<const reco::Muon*>(mu);
    if(muon)
    {
        reco::TrackRef muTrack = muon->track();
        if(!muTrack) continue;
        muonTTs.push_back(trackBuilder.build(muTrack));
    }
  }
  
  if(muonTTs.size() == 2)
  {
    KalmanVertexFitter kvf(true);

    reco::Vertex muonVertex = TransientVertex(kvf.vertex(muonTTs));

    math::XYZPoint fittedPoint = muonVertex.position();

    cout << "genBsDecayPoint: " << fittedPoint << " fittedPoint: " << fittedPoint << endl;

  }
 
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(PVdistance);

