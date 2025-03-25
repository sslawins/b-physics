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
#include "Math/SVector.h"

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
  TH1D *hPCAz_testSV;
  TH1D *hPCAz;
  TH1D *hPCAz_true;
  TH1D *hPvSv;
  TH1D *hMuMu_vz;
  TH1D *hMuMu_vz_allReco;
  TH1D *hPCA;
  TH1D *hPCA_T;
  TH1D *hPCA_true;
  TH1D *hPCA_T_true;
  TH1D *hPVequalToSV;
  TH1D *hNofPV;
  TH2D *h_pGT_pMMT;
  TH1D *h_pGT_pMMT_ratio;
  TH1D *h_alfa;
  TH1D *h_beta;

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

  hPhiMass = new TH1D("hPhiMass", "Reconstruction of #Phi ; M_{inv} [GeV]; Counts", 3000 , 0.9, 1.2);

  hPCAz_testSV = new TH1D("hPCAz_testSV", "|PCA-PV|_{z} calculated using Candidate::vertex() method; Distance [cm]; #Events", 200, 0., 0.2);
  
  hPCAz = new TH1D("hPCAz", "|PCA-PV|_{z}; Distance [cm]; #Events", 2000, 0., 1.);
  hPCA = new TH1D("hPCA", "|PCA-PV|; Distance [cm]; #Events", 2000, 0., 1.);
  hPCA_T = new TH1D("hPCA_T", "|PCA-PV|_{T}; Distance [cm]; #Events", 2000, 0., 1.);
  hPCAz_true = new TH1D("hPCAz_true", "|PCA-PV|_{z}; Distance [cm]; #Events", 2000, 0., 0.2);
  hPCA_true = new TH1D("hPCA_true", "|PCA-PV|; Distance [cm]; #Events", 2000, 0., 0.2);
  hPCA_T_true = new TH1D("hPCA_T_true", "|PCA-PV|_{T}; Distance [cm]; #Events", 2000, 0., 0.2);

  hPvSv = new TH1D("hPvSv", "|SV-PV|; Distance [cm]; #Events", 4000, 0.,4.);
  hMuMu_vz = new TH1D("hMuMu_vz", "|delta vz|; Distance [cm]; #Events", 500000, 0.,5.);
  hMuMu_vz_allReco = new TH1D("hMuMu_vz_allReco", "|v_{\\mu^+, z}-v_{\\mu^-, z}|; Distance [cm]; #Events", 500000, 0.,5.);

  hPVequalToSV = new TH1D("hPVequalToSV", "Number of events where PV^{gen}=SV^{gen}; N_{PV}; #Events", 3, -0.5, 2.5);
  hNofPV = new TH1D("hNofPV", "Number of PV^{reco}, fully reconstructed events; N_{PV}; #Events", 5, -0.5, 4.5);

  h_pGT_pMMT = new TH2D("h_pGT_pMMT", "Reconstructed transverse momentum of \\mu\\mu vs \\gamma; p_{\\mu\\mu} [GeV]; p_{\\gamma} [GeV]", 700, 0, 70, 700, 0, 70);
  h_pGT_pMMT_ratio = new TH1D("h_pGT_pMMT_ratio", "Ratio of (p_{\\mu\\mu}-p_{\\gamma})/p_{\\mu\\mu}; Ratio; #Events", 100, 0, 10);
  h_alfa = new TH1D("h_alfa", "Scaling factor; \\alfa; #Events", 2000, -10., 10.);
  h_beta = new TH1D("h_beta", "Scaling factor; \\beta; #Events", 2000, -10., 10.);

  cout << "HERE PVdistance::beginJob()" << endl;
}

void PVdistance::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hPhiMass      -> Write();
  hPCAz_true   -> Write();
  hPCAz         -> Write();
  hPCAz_testSV -> Write();
  hPvSv -> Write();
  hMuMu_vz      -> Write();
  hMuMu_vz_allReco -> Write();
  hPCA -> Write();
  hPCA_T -> Write();
  hPCA_true -> Write();
  hPCA_T_true -> Write();
  hPVequalToSV -> Write();
  hNofPV -> Write();
  h_pGT_pMMT -> Write();
  h_pGT_pMMT_ratio -> Write();
  h_alfa -> Write();
  h_beta -> Write();

  myRootFile.Close();

  delete hPhiMass;
  delete hPCAz_true;
  delete hPCAz;
  delete hPCAz_testSV;
  delete hPvSv;
  delete hMuMu_vz;
  delete hMuMu_vz_allReco;
  delete hPCA;
  delete hPCA_T;
  delete hPCA_true;
  delete hPCA_T_true;
  delete hPVequalToSV;
  delete hNofPV;
  delete h_pGT_pMMT;
  delete h_pGT_pMMT_ratio;
  delete h_alfa;
  delete h_beta;

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
  const std::vector<reco::Photon> & recoPhotons = ev.get(thePhotonToken);
  
  const auto & trackBuilder = es.getData(theTrackBuilderToken);

  HLTdecision HLTdecision(theTriggerResultsToken, ev, theConfig);
  bool accepted = HLTdecision.checkTriggers(ev, true); //print = true
  //if(!accepted) return ;

  ///////////////////////////////Trigger was fired
  BDecayAnalyzer bAnalyzer;
  std::vector<std::vector<const reco::Candidate*>> tree = bAnalyzer.analyzeBDecays(genPar, DecayTools::MuMu);
  genPhotons = bAnalyzer.getPhotons();
  genMuons = bAnalyzer.getPhiProducts();
  if(! bAnalyzer.analyzeEvent(ev.id())) return;
  // print the family tree
  bAnalyzer.printTheTree(tree);
  
  ///////////////////////////////Trigger was fired AND the desired decay occurred
  ParticleMatcher muonMatcher( recoMuons, genMuons, 0.01);
  muonMatcher.matchRecoToGen();
  muonMatcher.printMatchedParticles();
  std::cout << "Muons, Matching successful: " << muonMatcher.isSuccessful() << std::endl;
  if(!muonMatcher.isSuccessful()) return;
  

  //////////////////////////////Trigger was fired AND muons were matched
  ParticleMatcher photonMatcher( recoPhotons, genPhotons, 0.02);
  photonMatcher.matchRecoToGen();
  photonMatcher.printMatchedParticles();
  std::cout << "Photons, Matching successful: " << photonMatcher.isSuccessful() << std::endl;
  //if(!photonMatcher.isSuccessful()) return;

  std::vector< const reco::Candidate*> matchedMuons   = muonMatcher.getMatched();
  std::vector< const reco::Candidate*> matchedPhotons = photonMatcher.getMatched();

  //All muons, provided they are reconstructed
  hMuMu_vz_allReco->Fill (abs(matchedMuons[0]->vz() - matchedMuons[1]->vz()));

  /////////////////////////////Trigger was fired AND muons were matched AND photon was matched
  if( muonMatcher.isSuccessful() 
      && bAnalyzer.analyzeEvent(ev.id()) 
      && accepted 
      && photonMatcher.isSuccessful() ){
    
    //reconstructable events
    hMuMu_vz ->Fill(abs(matchedMuons[0]->vz() - matchedMuons[1]->vz()));

    ////////////////// PV /////////////////
    const std::vector<reco::Vertex> & PVertices = ev.get(thePVToken);  

    std::cout << "Primary vertices reco: " << std::endl;
    for( const auto& vertex : PVertices){
      std::cout << vertex.x() << "  " << vertex.y() << "  " << vertex.z() << "  " << std::endl;
    }
  
    math::XYZPoint pv     = PVertices[0].position();    // recoPV - the first out of the list
    math::XYZPoint pv_gen = tree[0].back()->vertex();   //vertex of the Bs0 meson
                                                        //the last one from the first row of the tree
    math::XYZPoint sv_test  = matchedMuons[0]->vertex();//Candidate::vertex() method (to check the difference)
    math::XYZPoint svG      = genMuons[0]->vertex();    //vertex of the muons from the phi decay
    
    math::XYZPoint svG_photon = genPhotons[0]->vertex(); //vertex of the photon from the Bs0 decay
    math::XYZPoint svG_phi1   = tree[0].back()->daughter(0)->vertex(); //vertex of the photon or phi from the Bs0 decay
    math::XYZPoint svG_phi2   = tree[0].back()->daughter(1)->vertex(); //vertex of the photon ot phi from the Bs0 decay

    math::XYZVectorD pMuMu      = matchedMuons[0]->momentum() + matchedMuons[1]->momentum();
    math::XYZVectorD pMuMu_gen  = genMuons[0]->momentum() + genMuons[1]->momentum();
    math::XYZVectorD pMuMuGamma = genMuons[0]->momentum() + genMuons[1]->momentum() + genPhotons[0]->momentum();

    std::cout<< "PV GEN:           " << "(" << pv_gen.X() << ", " << pv_gen.Y() << ", " << pv_gen.Z() << ")" << std::endl;

    math::XYZPoint pca_test       = DecayTools::pca(pv, sv_test, pMuMu);
    math::XYZPoint pca_gengen     = DecayTools::pca(pv_gen, svG, pMuMu_gen);
    math::XYZPoint pca_withGamma  = DecayTools::pca(pv_gen, svG, pMuMuGamma);

    vector<reco::TransientTrack> muonTTs;
    for (const auto& mu : matchedMuons){
      const reco::Muon* muon = dynamic_cast<const reco::Muon*>(mu);
      if(muon){
          reco::TrackRef muTrack = muon->track();
          if(!muTrack) continue;
          muonTTs.push_back(trackBuilder.build(muTrack));
      }
    }

    if(muonTTs.size() == 2){
      KalmanVertexFitter kvf(true);

      reco::Vertex muonVertex = TransientVertex(kvf.vertex(muonTTs));

      math::XYZPoint fittedSV = muonVertex.position();
      std::cout << "Gen muons' SV:    " << "(" << svG.x() 
                << " , " << svG.y()     << " , " << svG.z() << ")" << std::endl;
      
      std::cout << "Gen photon's SV:          " << svG_photon << std::endl;
      std::cout << "Gen Bs0's daughter[0] SV: " << svG_phi1 << std::endl;
      std::cout << "Gen Bs0's daughter[1] SV: " << svG_phi2 << std::endl;

      for( const auto& muon : matchedMuons){
        std::cout << "Matched muons' SV:" << "(" << muon->vertex().x() << ", " << muon->vertex().y()
                << ", " << muon->vertex().z() << ")"
                << std::endl;
      }

      math::XYZPoint pca_reco = DecayTools::pca(pv, fittedSV, pMuMu);

      hPCAz_testSV ->Fill(abs(pca_test.z() - pv.z()));
      
      hPCA -> Fill(sqrt((pca_reco - pv).Mag2()));
      hPCA_T -> Fill(sqrt((pca_reco - pv).Perp2()));
      hPCAz -> Fill(abs(pca_reco.z() - pv.z()));
      
      hPCA_true -> Fill(sqrt((pca_gengen - pv_gen).Mag2()));
      hPCA_T_true -> Fill(sqrt((pca_gengen - pv_gen).Perp2())); 
      hPCAz_true   -> Fill(abs(pca_gengen.z() - pv_gen.z()));
      
      hPvSv ->Fill(abs(pv.z() - fittedSV.z()));
      hNofPV -> Fill(PVertices.size());
      hPVequalToSV -> Fill(pv_gen == svG);

      std::cout<< "Fitted SV:" << "(" << fittedSV.X() << ", " << fittedSV.Y() << ", " << fittedSV.Z() << ")" << std::endl;
      std::cout<< "/////////////////////////////"  << std::endl;
      std::cout<< "True PCA: " << "(" << pca_gengen.X() << ", " << pca_gengen.Y() << ", " << pca_gengen.Z() << ")" << std::endl;
      std::cout<< "PCA:      " << "(" << pca_test.X() << ", " << pca_test.Y() << ", " << pca_test.Z() << ")" << std::endl;
      std::cout<< "PCA_reco: " << "(" << pca_reco.X() << ", " << pca_reco.Y() << ", " << pca_reco.Z() << ")" << std::endl;
      std::cout<< "/////////////////////////////"  << std::endl;
      std::cout<< "WithGamma:" << "(" << pca_withGamma.X() << ", " << pca_withGamma.Y() << ", " << pca_withGamma.Z() << ")" << std::endl;

      /////////////////////////////// p_{\mu\mu} vs p_{\gamma}

      math::XYZVectorD bsDirection = (fittedSV - pv).Unit(); //direction of the Bs0 meson
      math::XYZVectorD pMuMuT = pMuMu - bsDirection.Dot(pMuMu)*bsDirection;
      math::XYZVectorD pGammaT = matchedPhotons[0]->momentum() - bsDirection.Dot(matchedPhotons[0]->momentum())*bsDirection;

      h_pGT_pMMT -> Fill(pMuMuT.R(), pGammaT.R());
      h_pGT_pMMT_ratio -> Fill((pMuMuT.R() - pGammaT.R())/pMuMuT.R());

      // Calculate scaling factor alfa
      double alfa = -(pMuMuT.Dot(pGammaT)) / pGammaT.Mag2();
      double beta = - pMuMuT.Mag2() / pGammaT.Dot(pMuMuT);
      h_alfa->Fill(alfa);
      h_beta->Fill(beta);
      
    }

  }
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(PVdistance);

