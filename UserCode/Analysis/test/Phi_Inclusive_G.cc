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
  edm::EDGetTokenT < edm::TriggerResults > theTriggerResultsToken;
  edm::EDGetTokenT<std::vector<reco::Vertex>> thePVToken;

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTrackBuilderToken;

  TH1D *hHLT;
  TH1D *hHLT_phiToKK;
  TH1D *hHLT_gamma;
  TH1D *hHLT_phiToKK_gamma;

  TH1D *hPhiToKplusKminus;
  TH1D *hpTrecoMuons;
  TH2D *hKKGen_Eta_Pt;
  TH1D *hKKGen_pT;
  TH2D *hKKReco_Eta_Pt;
  TH2D *hKKpT_recoVSgen;
  TH1D *hKK_deltaR_allReco;
  TH1D *hKK_deltaR;
  TH1D *hKK_vz;
  TH1D *hKK_vz_allReco;

  TH1D *hPhiG_deltaR;
  TH1D *hBsKK_vz;

  TH1D *hPCAz;
  TH1D *hPCAz_true;
  TH1D *hPCA;
  TH1D *hPCA_T;
  TH1D *hPCA_true;
  TH1D *hPCA_T_true;

  //TH1D *h_pGT_pMMT_ratio;
  //TH1D *h_alfa;
  //TH1D *h_beta;
  TH1D *h_pGT_pMM_angle;

  TH1D *hPvSv;
  TH1D *hVtxProbKK;

  // Invariant masses
  TH1D *hPhiMass;
  TH1D *hRecoKaons_RecoPhotonDir_GenPhotonMag;
  TH1D *hRecoKaons_GenPhotonDir_RecoPhotonMag;
  TH1D *hRecoPhotonFourMomentaMass;
  TH1D *hGenPhotonFourMomentaMass;

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
  theTriggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
  thePVToken = consumes< vector<reco::Vertex>  >( edm::InputTag("offlinePrimaryVertices"));

  theTrackBuilderToken = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));
}

//destructor
Phi_Inclusive_G::~Phi_Inclusive_G()
{
  cout <<" DTOR" << endl;
}

void Phi_Inclusive_G::beginJob()
{
  hHLT = new TH1D("hHLT", "HLT decision",3,  -0.5, 2.5);
  hHLT_phiToKK = new TH1D("hHLT_phiToKK", "HLT decision for Phi to K+K- decay", 3, -0.5, 2.5);
  hHLT_gamma = new TH1D("hHLT_gamma", "HLT decision for gamma", 3, -0.5, 2.5);
  hHLT_phiToKK_gamma = new TH1D("hHLT_phiToKK_gamma", "HLT decision for Phi to K+K- decay with gamma", 3, -0.5, 2.5);
  
  hPhiToKplusKminus = new TH1D("hPhiToKplusKminus", "Phi to K+K- decay", 2, -0.5, 1.5);
  hpTrecoMuons = new TH1D("hpTrecoMuons", "reco muon pt", 1000, 0, 50); //NOT FILLED YET
  
  hKKGen_Eta_Pt = new TH2D("hKKGen_Eta_Pt", "K+K- gen; #eta; p_{T, gen}", 1000, -10., 10., 500, 0., 50.);
  hKKGen_pT = new TH1D("hKKGen_pT", "K+K- gen; p_{T, gen} [GeV]; Events", 500, 0, 50);
  hKKReco_Eta_Pt = new TH2D("hKKReco_Eta_Pt", "K+K- reco; #eta; p_{T, reco}", 1000, -10., 10., 500, 0., 50.);
  hKKpT_recoVSgen = new TH2D("hKKpT_recoVSgen", "K+K- reco-gen; p_{T, reco}; p_{T, gen}", 500, 0, 50, 500, 0, 50);
  hKK_deltaR_allReco = new TH1D("hKK_deltaR_allReco", "K+K- deltaR, all reconstructed particles; \\Delta R; Events", 100, 0, 0.5);
  hKK_deltaR = new TH1D("hKK_deltaR", "K+K- deltaR; \\Delta R; Events", 250, 0, 0.5);
  hKK_vz = new TH1D("hKK_vz", "|\\Delta v_z|; Distance [cm]; Events", 500000, 0, 5);
  hKK_vz_allReco = new TH1D("hKK_vz_allReco", "|v_{\\mu^+, z}-v_{\\mu^-, z}|; Distance [cm]; #Events", 3000, 0.,0.3);
  
  hPhiG_deltaR = new TH1D("hPhiG_deltaR", "Phi-Gamma deltaR; \\Delta R; Events", 250, 0, 5.0);
  hBsKK_vz = new TH1D("hBsKK_vz", "|v_{Bs, z} - v_{KK, z}|; Distance [cm]; Events", 500, 0, 5.0);

  /* FOR CASCADE MUONS - IF THERE IS TIME
  hMuMu_vz = new TH1D("hMuMu_vz", "|delta vz|; Distance [cm]; #Events", 500000, 0.,5.);
  hMuMu_vz_allReco = new TH1D("hMuMu_vz_allReco", "|v_{\\mu^+, z}-v_{\\mu^-, z}|; Distance [cm]; #Events", 500000, 0.,5.);
  */

  hPCAz = new TH1D("hPCAz", "|PCA-PV|_{z}; Distance [cm]; #Events", 2000, 0., 1.);
  hPCA = new TH1D("hPCA", "|PCA-PV|; Distance [cm]; #Events", 2000, 0., 1.);
  hPCA_T = new TH1D("hPCA_T", "|PCA-PV|_{T}; Distance [cm]; #Events", 2000, 0., 1.);
  hPCAz_true = new TH1D("hPCAz_true", "|PCA-PV|_{z}; Distance [cm]; #Events", 2000, 0., 0.2);
  hPCA_true = new TH1D("hPCA_true", "|PCA-PV|; Distance [cm]; #Events", 2000, 0., 0.2);
  hPCA_T_true = new TH1D("hPCA_T_true", "|PCA-PV|_{T}; Distance [cm]; #Events", 2000, 0., 0.2);

  h_pGT_pMM_angle = new TH1D("h_pGT_pMM_angle", "cos(p_{\\gamma,T}; p_{\\mu\\mu,T}); cos(\\alfa); Events", 1000, -1, 1);

  hPvSv = new TH1D("hPvSv", "|SV-PV|; Distance [cm]; #Events", 4000, 0.,4.);

  hVtxProbKK = new TH1D("hVtxProbKK", "Vtx prob of K^{+}K^{-}; Probability; Events", 1000, 0., 1.);
  
  hPhiMass = new TH1D("hPhiMass", "Reconstruction of #Phi ; M_{inv} [GeV]; Events", 3000 , 0.9, 1.2);
  hRecoKaons_RecoPhotonDir_GenPhotonMag = new TH1D("hRecoKaons_RecoPhotonDir_GenPhotonMag", 
    "Reco muons, reco photon direction, gen photon magnitude; M_{inv} [GeV]; Events", 400, 3, 7 );
  hRecoKaons_GenPhotonDir_RecoPhotonMag = new TH1D("hRecoKaons_GenPhotonDir_RecoPhotonMag", 
    "Reco muons, gen photon direction, reco photon magnitude; M_{inv} [GeV]; Events", 400, 3, 7 );
  hRecoPhotonFourMomentaMass = new TH1D("hRecoPhotonFourMomentaMass", 
    "Invariant mass of reco muons and reco photon; M_{inv} [GeV]; Events", 400, 3, 7 );
  hGenPhotonFourMomentaMass = new TH1D("hGenPhotonFourMomentaMass", 
    "Invariant mass of reco muons and gen photon; M_{inv} [GeV]; Events", 400, 3, 7 );

  cout << "HERE Phi_Inclusive_G::beginJob()" << endl;
}

void Phi_Inclusive_G::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hHLT -> Write();
  hHLT_phiToKK -> Write();
  hHLT_gamma -> Write();
  hHLT_phiToKK_gamma -> Write();
  
  hPhiToKplusKminus -> Write();
  hpTrecoMuons -> Write();
  hKKGen_Eta_Pt -> Write();
  hKKGen_pT -> Write();
  hKKReco_Eta_Pt -> Write();
  hKKpT_recoVSgen -> Write();
  hKK_deltaR_allReco -> Write();
  hKK_deltaR -> Write();  
  hKK_vz->Write();
  hKK_vz_allReco->Write();
  
  hPhiG_deltaR->Write();
  hBsKK_vz->Write();

  hPCAz->Write();
  hPCA->Write();
  hPCA_T->Write();
  hPCAz_true->Write();
  hPCA_true->Write();
  hPCA_T_true->Write();

  //h_pGT_pMMT_ratio->Write();
  //h_alfa->Write();
  //h_beta->Write();
  h_pGT_pMM_angle->Write();
  
  hPvSv->Write();
  hVtxProbKK->Write();

  //Invariant masses
  hPhiMass->Write();
  hRecoKaons_RecoPhotonDir_GenPhotonMag->Write();
  hRecoKaons_GenPhotonDir_RecoPhotonMag->Write();
  hRecoPhotonFourMomentaMass->Write();
  hGenPhotonFourMomentaMass->Write();
  
  myRootFile.Close();

  delete hHLT;
  delete hHLT_phiToKK;
  delete hHLT_gamma;
  delete hHLT_phiToKK_gamma;

  delete hPhiToKplusKminus;
  delete hpTrecoMuons;
  delete hKKGen_Eta_Pt;
  delete hKKGen_pT;
  delete hKKReco_Eta_Pt;
  delete hKKpT_recoVSgen;
  delete hKK_deltaR_allReco;
  delete hKK_deltaR;
  delete hKK_vz;
  delete hKK_vz_allReco;

  delete hPhiG_deltaR;
  delete hBsKK_vz;

  delete hPCAz;
  delete hPCA;
  delete hPCA_T;
  delete hPCAz_true;
  delete hPCA_true;
  delete hPCA_T_true;
  
  //delete h_pGT_pMMT_ratio;
  //delete h_alfa;
  //delete h_beta;
  delete h_pGT_pMM_angle;

  delete hPvSv;
  delete hVtxProbKK;

  delete hPhiMass;
  delete hRecoKaons_RecoPhotonDir_GenPhotonMag;
  delete hRecoKaons_GenPhotonDir_RecoPhotonMag;
  delete hRecoPhotonFourMomentaMass;
  delete hGenPhotonFourMomentaMass;
  
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

  const auto & trackBuilder = es.getData(theTrackBuilderToken);

  //////Trigger 
  HLTdecision HLT(theTriggerResultsToken, ev, theConfig);
  bool accepted = HLT.checkTriggers(ev, true); //print = true

  /////B quarks history
  BDecayAnalyzer bAnalyzer;
  std::vector<std::vector<const reco::Candidate*>> tree = bAnalyzer.analyzeBDecays(genPar, DecayTools::KK);
  genPhotons = bAnalyzer.getPhotons();
  genKaons = bAnalyzer.getPhiProducts();
  std::cout << "Phi daughters: " << genKaons.size() << std::endl;

  // there is no need to analyze further if there is no phi -> K+K- decay
  if( !bAnalyzer.analyzeEvent(ev.id())) hPhiToKplusKminus->Fill(0);
  else{
    hKKGen_Eta_Pt -> Fill(genKaons[0]->eta(), genKaons[0]->pt());
    hKKGen_Eta_Pt -> Fill(genKaons[1]->eta(), genKaons[1]->pt());
    hKKGen_pT -> Fill(genKaons[0]->pt());
    hKKGen_pT -> Fill(genKaons[1]->pt());
    hPhiToKplusKminus->Fill(1);
    // to find the PV
    
    /*const auto* particle = tree[0].back();

    while (particle->mother()) {
        particle = particle->mother();
        std::cout << "PDG ID mother: " << particle->pdgId() << std::endl;
    }*/

    // print the family tree
    bAnalyzer.printTheTree(tree);
    std::cout << "Number of kaons: " << genKaons.size() << std::endl;
    
    ///////////////////Kaon matching
    ParticleMatcher<reco::PFCandidate> particleMatcher( recoCandidates, genKaons, 0.01, hKK_deltaR_allReco);
    particleMatcher.matchRecoToGen();
    particleMatcher.printMatchedParticles();
    std::cout << "Matching successful: " << particleMatcher.isSuccessful() << std::endl;
    std::vector< const reco::Candidate*> matchedParticles = particleMatcher.getMatched();
    if( particleMatcher.isSuccessful() ) {
      // all the reco kaons provided they are reconstructable
      hKKpT_recoVSgen -> Fill(matchedParticles[0]->pt(), genKaons[0]->pt());
      hKKpT_recoVSgen -> Fill(matchedParticles[1]->pt(), genKaons[1]->pt());
      hKKReco_Eta_Pt -> Fill(matchedParticles[0]->eta(), matchedParticles[0]->pt());
      hKKReco_Eta_Pt -> Fill(matchedParticles[1]->eta(), matchedParticles[1]->pt());
      //hKK_deltaR -> Fill( reco::deltaR(matchedParticles[0]->momentum(), matchedParticles[1]->momentum()));
      hKK_vz_allReco -> Fill( abs(matchedParticles[0]->vz() - matchedParticles[1]->vz()));
    }

    //////////////////Photon matching
    ParticleMatcher<reco::Photon> photonMatcher( recoPhotons, genPhotons, 0.03);
    photonMatcher.matchRecoToGen();
    photonMatcher.printMatchedParticles();
    std::cout << "Matching successful: " << photonMatcher.isSuccessful() << std::endl;
    std::vector< const reco::Candidate*> matchedPhotons = photonMatcher.getMatched();
    
    
    //////////////////// HLT analysis
    if(!accepted){
      hHLT->Fill(0);
    }else if(accepted){  
      hHLT->Fill(1);
      if( particleMatcher.isSuccessful() && photonMatcher.isSuccessful()) hHLT_phiToKK_gamma->Fill(1);
      if (particleMatcher.isSuccessful()) hHLT_phiToKK->Fill(1);
      if (photonMatcher.isSuccessful()) hHLT_gamma->Fill(1);
    }

    std::cout <<"HLT/ Kaon matcher/ photon matcher: "<< accepted << particleMatcher.isSuccessful() << photonMatcher.isSuccessful() <<endl;
    
    if( particleMatcher.isSuccessful() && photonMatcher.isSuccessful()){
      std::cout << " muons were matched AND photon was matched" << std::endl;
      hKK_vz -> Fill( abs(matchedParticles[0]->vz() - matchedParticles[1]->vz()));
      hKK_deltaR -> Fill( reco::deltaR(matchedParticles[0]->momentum(), matchedParticles[1]->momentum()));
      const std::vector<reco::Vertex> & PVertices = ev.get(thePVToken); 
      
      /*for( const auto& vertex : PVertices){
        std::cout << vertex.x() << "  " << vertex.y() << "  " << vertex.z() << " ,  :" << vertex.tracksSize() <<std::endl;
      }*/

      math::XYZPoint pv = PVertices[0].position();
      math::XYZPoint pv_gen = tree[0][0]->vertex();
      math::XYZPoint svG = genKaons[0]->vertex();

      math::XYZVectorD pKK = matchedParticles[0]->momentum() + matchedParticles[1]->momentum();
      math::XYZVectorD pKK_gen = genKaons[0]->momentum() + genKaons[1]->momentum();

      math::XYZPoint pca_true = DecayTools::pca(pv_gen, svG, pKK_gen);

      hPCA_true->Fill(sqrt((pca_true - pv_gen).Mag2()));
      hPCA_T_true->Fill(sqrt((pca_true - pv_gen).Perp2()));
      hPCAz_true->Fill(abs(pca_true.z() - pv_gen.z()));

      vector<reco::TransientTrack> kaonTTs;

      for (const auto& kk : matchedParticles){
        const reco::Track * rtrk = kk -> bestTrack();
        if(!rtrk) continue;
        kaonTTs.push_back(trackBuilder.build(rtrk));
      
      }

      if(kaonTTs.size() == 2){
        KalmanVertexFitter kvf(true);

        reco::Vertex kaonVertex = TransientVertex(kvf.vertex(kaonTTs));
        double probability = TMath::Prob(kaonVertex.chi2(),kaonVertex.ndof());
        hVtxProbKK ->Fill( probability );
        std::cout << "Vtx prob of kaons: " << probability << std::endl;
        math::XYZPoint fittedSV = kaonVertex.position();
        std::cout << "Fitted SV: " << fittedSV << std::endl;
        math::XYZPoint pca_reco = DecayTools::pca(pv, fittedSV, pKK);
        std::cout << "PCA: " << pca_reco << std::endl;

        hPCA->Fill(sqrt((pca_reco - pv).Mag2()));
        hPCA_T->Fill(sqrt((pca_reco - pv).Perp2()));
        hPCAz->Fill(abs(pca_reco.z() - pv.z()));
        hPvSv->Fill(abs(pv.z() - fittedSV.z()));
      }
      /////////////////////////////////////
      math::XYZVectorD bsDirection = (svG - pv).Unit();
      math::XYZVectorD pKK_T = pKK - bsDirection.Dot(pKK) * bsDirection;

      double alfa = -(pKK_T.Dot(pKK)) / pKK.Mag2();
      //h_alfa->Fill(alfa);

      
      // Add reco kaons
      std::vector<const reco::Candidate*> recoKaonsWithRecoPhoton;
      std::vector<const reco::Candidate*> recoKaonsWithGenPhoton;


      recoKaonsWithRecoPhoton.insert(recoKaonsWithRecoPhoton.end(), matchedParticles.begin(), matchedParticles.end());
      recoKaonsWithGenPhoton.insert(recoKaonsWithGenPhoton.end(), matchedParticles.begin(), matchedParticles.end());

      // Reco kaons reco photon
      const reco::Candidate* recoPhotonAsKaon = dynamic_cast<const reco::Candidate*>(matchedPhotons[0]);
      if (recoPhotonAsKaon) {
        recoKaonsWithRecoPhoton.push_back(recoPhotonAsKaon);
      }else{
        std::cerr << "Error: dynamic_cast failed for reco photon." << std::endl;
      }
      std::cout << "RecoPhotonFourMomenta: " << recoKaonsWithRecoPhoton.size()<< std::endl;
      auto recoPhotonFourMomenta = DecayTools::fourMomenta(recoKaonsWithRecoPhoton, DecayTools::KKGmasses);

      // Reco kaons gen photon
      recoKaonsWithGenPhoton.push_back(genPhotons[0]);
      std::cout << "GenPhotonFourMomenta: " << recoKaonsWithGenPhoton.size()<< std::endl;
      auto genPhotonFourMomenta = DecayTools::fourMomenta(recoKaonsWithGenPhoton, DecayTools::KKGmasses);

      // Reco kaons, reco photon direction, gen photon magnitude
      double P4recoDirGenMag = DecayTools::scaledInvariant(matchedParticles, genPhotons[0], matchedPhotons[0], DecayTools::KKGmasses);
      // Reco kaons, gen photon direction, reco photon magnitude
      double P4recoMagGenDir = DecayTools::scaledInvariant(matchedParticles, matchedPhotons[0], genPhotons[0], DecayTools::KKGmasses);

      std::cout << "Invariant mass of the Bs meson: " << recoPhotonFourMomenta.M() << std::endl;
      std::cout << "Invariant mass of the Bs meson (gen): " << genPhotonFourMomenta.M() << std::endl;
      std::cout << "Invariant mass of the Bs meson (reco kaons, reco photon direction, gen photon magnitude): " << P4recoDirGenMag << std::endl;
      std::cout << "Invariant mass of the Bs meson (reco kaons, gen photon direction, reco photon magnitude): " << P4recoMagGenDir << std::endl;
      std::cout << "Invariant mass of the Phi meson (reco muons): " << DecayTools::invariantMass(matchedParticles, DecayTools::KKMasses) << std::endl;

      hPhiMass ->Fill(DecayTools::invariantMass(matchedParticles, DecayTools::KKMasses));
      hRecoPhotonFourMomentaMass->Fill(recoPhotonFourMomenta.M());
      hGenPhotonFourMomentaMass->Fill(genPhotonFourMomenta.M());
      hRecoKaons_GenPhotonDir_RecoPhotonMag->Fill(P4recoMagGenDir);
      hRecoKaons_RecoPhotonDir_GenPhotonMag->Fill(P4recoDirGenMag);

      hPhiG_deltaR->Fill(reco::deltaR(matchedPhotons[0]->momentum(), pKK));
      hBsKK_vz->Fill(abs(pv.z() - (matchedParticles[0]->vz() + matchedParticles[1]->vz())/2.));
    }
  }
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(Phi_Inclusive_G);

