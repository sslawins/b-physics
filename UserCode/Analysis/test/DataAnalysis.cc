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
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/Vector4D.h"

#include <vector>
#include <numeric>
#include <algorithm>

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

  //function to calculate the PCA
  math::XYZPoint pca(math::XYZPoint pv, math::XYZPoint sv, math::XYZVectorD pMuMu);

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

  //histograms

  TH1D *hpTMuons;
  TH1D *hMM_deltaR;
  TH1D *hMMDelta_vz;
  TH1D *hMuMu_mass;
  TH1D *hMuMu_mass_disp;
  TH1D *hMM_prob;

  TH1D *hpTKaons;
  TH1D *hKK_deltaR;
  TH1D *hKKDelta_vz;
  TH1D *hKK_mass;
  TH1D *hKKmass_disp;
  TH1D *hKK_prob;

  TH1D *hKKGmass;
  TH1D *hMuMuGmass;

  TH1D *hNoPV;
  TH1D *hPV_distance;
  TH1D *hpTG;
  TH2D *hpT_etaG;

  TH1D *hPCA1smallest; // the smallest distance along the z axis between PV and the PCA of the two kaons (-1 for cases when it is not 1st PV)
  TH1D *hPCA1vs2; //the difference between the smallest and the second smallest distance along the z axis between PV and the PCA of the two kaons
  TH1D *hPCAvsVz; 
  TH1D *h3Ddisp;
  TH1D *h3DdispErr;
  TH1D *h3DdispSignificance;

  TH1D *hPCA1smallest_mu; // the smallest distance along the z axis between PV and the PCA of the two muons (-1 for cases when it is not 1st PV)
  TH1D *hPCA1vs2_mu; //the difference between the smallest and the second smallest distance along the z axis between PV and the PCA of the two muons
  TH1D *hPCAvsVz_mu;
  TH1D *h3Ddisp_mu;
  TH1D *h3DdispErr_mu;
  TH1D *h3DdispSignificance_mu;
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
  hpTMuons = new TH1D("hpTMuons", "reco muons pt; p_T [GeV]; Events ", 3000, 0, 30);
  hMM_deltaR = new TH1D("hMM_deltaR", "Mu+Mu- deltaR; \\Delta R; Events", 1100, 0, 11);
  hMMDelta_vz = new TH1D("hMMDelta_vz", "delta vz; \\Delta v_{z} [cm]; Events", 50000, 0, 5);
  hMuMu_mass = new TH1D("hMuMu_mass", "MuMu mass; M_{inv} [GeV]; Events", 80000, 0.8, 1.6);
  hMuMu_mass_disp = new TH1D("hMuMu_mass_disp", "Displaced Mu+Mu-Photon mass; M_{inv} [GeV]; Events", 80000, 0.8, 1.6);
  hMM_prob = new TH1D("hMM_prob", "MuMu vertex prob; Probability; Events", 1000, 0, 1);

  hpTKaons = new TH1D("hpTKaons", "reco kaons pt; p_T [GeV]; Events ", 3000, 0, 30);
  hKK_deltaR = new TH1D("hKK_deltaR", "K+K- deltaR; \\Delta R; Events", 1100, 0, 11);
  hKKDelta_vz = new TH1D("hDelta_vz", "delta vz; \\Delta v_{z} [cm]; Events", 50000, 0, 5);
  hKK_mass = new TH1D("hKK_mass", "K+K- mass; M_{inv} [GeV]; Events", 80000, 0.8, 1.6);
  hKKmass_disp = new TH1D("hKKmass_disp", "Displaced K+K-Photon mass; M_{inv} [GeV]; Events", 80000, 0.8, 1.6);
  hKK_prob = new TH1D("hKK_prob", "K+K- vertex prob; Probability; Events", 1000, 0, 1);

  hKKGmass = new TH1D("hKKGmass", "K+K-Photon mass; M_{inv} [GeV]; Events", 8000, 0.0, 8.0);
  hMuMuGmass = new TH1D("hMuMuGmass", "MuMuPhoton mass; M_{inv} [GeV]; Events", 8000, 0.0, 8.0);

  hNoPV = new TH1D("hNoPV", "Number of PVs;", 101, -0.5, 100.5);
  //hPV_distance = new TH1D("hPV_distance", "Distance beetween PVs_{z}", 4000, 0, 40);
  hpTG = new TH1D("hpTG", "gamma pT; p_T [GeV]; Events", 10000, 0, 100);
  hpT_etaG = new TH2D("hpT_etaG", "gamma pT vs eta; \\eta; p_{T} [GeV]", 1000, 0, 100, 1000, -5, 5);

  hPCA1smallest = new TH1D("hPCA1smallest", "min(|PCA-PV|)_{z}; Distance [cm]; Events", 1000, 0 , 0.5);
  hPCA1vs2 = new TH1D("hPCA1vs2", "(min(|PCA-PV|)-second_min(|PCA-PV|)_{z} Distance [cm]; Events", 5000, 0, 0.5);
  hPCAvsVz = new TH1D("hPCAvsVz", "min(|PCA-PV|)_{z} vs min(|PV-PV|)_{z}; Distance [cm]; Events", 3, -0.5, 2.5);
  h3Ddisp = new TH1D("h3Ddisp", "3D displacement; Distance [cm]; Events", 2000, 0, 5);
  h3DdispErr = new TH1D("h3DdispErr", "3D displacement error; Distance [cm]; Events", 2000, 0, 5);
  h3DdispSignificance = new TH1D("h3DdispSignificance", "3D displacement significance; Significance; Events", 2000, 0, 20);

  hPCA1smallest_mu = new TH1D("hPCA1smallest_mu", "min(|PCA-PV|)_{z}; Distance [cm]; Events", 1000, 0 , 0.5);
  hPCA1vs2_mu = new TH1D("hPCA1vs2_mu", "(min(|PCA-PV|)-second_min(|PCA-PV|)_{z} Distance [cm]; Events", 5000, 0, 0.5);
  hPCAvsVz_mu = new TH1D("hPCAvsVz_mu", "min(|PCA-PV|)_{z} vs min(|PV-PV|)_{z}; Distance [cm]; Events", 3, -0.5, 2.5);
  h3Ddisp_mu = new TH1D("h3Ddisp_mu", "3D displacement; Distance [cm]; Events", 2000, 0, 5);
  h3DdispErr_mu = new TH1D("h3DdispErr_mu", "3D displacement error; Distance [cm]; Events", 2000, 0, 5);
  h3DdispSignificance_mu = new TH1D("h3DdispSignificance_mu", "3D displacement significance; Significance; Events", 2000, 0, 20);
  
  cout << "HERE DataAnalysis::beginJob()" << endl;
}

void DataAnalysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  hpTMuons ->Write();
  hMM_deltaR -> Write();
  hMMDelta_vz -> Write();
  hMuMu_mass -> Write();
  hMuMu_mass_disp -> Write();
  hMM_prob -> Write();

  hpTKaons -> Write();
  hKK_deltaR -> Write();
  hKKDelta_vz -> Write();
  hKK_mass -> Write();
  hKKmass_disp -> Write();
  hKK_prob -> Write();

  hKKGmass -> Write();
  hMuMuGmass -> Write();

  hNoPV -> Write();
  //hPV_distance -> Write();
  hpTG -> Write();
  hpT_etaG -> Write();

  hPCA1smallest -> Write(); // the smallest distance along the z axis between PV and the PCA of the two kaons (-1 for cases when it is not 1st PV)
  hPCA1vs2 -> Write(); //the difference between the smallest and the second smallest distance along the z axis between PV and the PCA of the two kaons
  hPCAvsVz -> Write();
  h3Ddisp -> Write();
  h3DdispErr -> Write();
  h3DdispSignificance -> Write();

  hPCA1smallest_mu -> Write(); // the smallest distance along the z axis between PV and the PCA of the two muons (-1 for cases when it is not 1st PV)
  hPCA1vs2_mu -> Write(); //the difference between the smallest and the second smallest distance along the z axis between PV and the PCA of the two muons
  hPCAvsVz_mu -> Write();
  h3Ddisp_mu -> Write();
  h3DdispErr_mu -> Write();
  h3DdispSignificance_mu -> Write();

  myRootFile.Close();

  delete hpTMuons;
  delete hMM_deltaR;
  delete hMMDelta_vz;
  delete hMuMu_mass;
  delete hMuMu_mass_disp;
  delete hMM_prob;

  delete hpTKaons;
  delete hKK_deltaR;
  delete hKKDelta_vz;
  delete hKK_mass;
  delete hKKmass_disp;
  delete hKK_prob;

  delete hKKGmass;
  delete hMuMuGmass;

  delete hNoPV;
  //delete hPV_distance;
  delete hpTG;
  delete hpT_etaG;

  delete hPCA1smallest; // the smallest distance along the z axis between PV and the PCA of the two kaons (-1 for cases when it is not 1st PV)
  delete hPCA1vs2; //the difference between the smallest and the second smallest distance along the z axis between PV and the PCA of the two kaons
  delete hPCAvsVz;
  delete h3Ddisp;
  delete h3DdispErr;
  delete h3DdispSignificance;

  delete hPCA1smallest_mu; // the smallest distance along the z axis between PV and the PCA of the two muons (-1 for cases when it is not 1st PV)
  delete hPCA1vs2_mu; //the difference between the smallest and the second smallest distance along the z axis between PV and the PCA of the two muons
  delete hPCAvsVz_mu;
  delete h3Ddisp_mu;
  delete h3DdispErr_mu;
  delete h3DdispSignificance_mu;


  cout << "HERE DataAnalysis::endJob()" << endl;
}

math::XYZPoint DataAnalysis::pca(math::XYZPoint pv, math::XYZPoint sv, math::XYZVectorD pMuMu){
  //s = (PV - SV) * (pMuMu) / |pMuMu|^2 
  double s = ((pv - sv).Dot(pMuMu)) / pMuMu.Mag2();
  //pca = sv + s*pMuMu
  math::XYZPoint PCA = sv + s*pMuMu;

  return PCA;
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
  //const std::vector<reco::Vertex>& recoVertices = ev.get(thePVWithBSToken);
  
  std::cout << "PFCandidates: " << candidates.size() << std::endl;
  std::cout << "Photons: " << recoPhotons.size() << std::endl;
  std::cout << "Muons: " << recoMuons.size() << std::endl;
  //std::cout << "Primary Vertices: " << recoVertices.size() << std::endl;
  std::cout << "Primary Vertices with Beam Spot: " << recoVertices.size() << std::endl;

  //if(recoPhotons == 0) return;
  hNoPV->Fill(recoVertices.size());
/*
  std::vector<double> eventPt2;

  for( unsigned int vtxIter = 0; vtxIter < recoVertices.size(); ++vtxIter){
    const reco::Vertex & vtx = recoVertices[vtxIter];
    //cout << "Primary Vertex: " << vtx.position().x() << " " << vtx.position().y() << " " << vtx.position().z() << endl;
    for( unsigned int vtxIter2 = vtxIter+1; vtxIter2 < recoVertices.size(); ++vtxIter2){
      const reco::Vertex & vtx2 = recoVertices[vtxIter2];
      double deltaVZ = fabs(vtx.position().z() - vtx2.position().z());
      hPV_distance->Fill(deltaVZ);
    }
    //std::cout << vtx.tracksSize() << " tracks in the vertex" << std::endl;
    double vtxPt2 = 0.0;
    for (auto iTrack = vtx.tracks_begin(); iTrack != vtx.tracks_end(); ++iTrack) {
      double pt = (*iTrack)->pt();
      std::cout << "Track pT: " << pt << std::endl;
      vtxPt2 += pt * pt;
      
    }
    eventPt2.push_back(vtxPt2);
  }
  hPVpT2->Fill(eventPt2[0]/std::accumulate(eventPt2.begin(), eventPt2.end(), 0.0));
  //std::cout<< "1st PV pT^2: " << eventPt2[0] << std::endl;
  //std::cout << "Event pT^2: " << std::accumulate(eventPt2.begin(), eventPt2.end(), 0.0) << std::endl;
  */

  cout << "//////////////////////////////////////////////////" << endl;

  double deltaRKaons = 0.25, deltaRPhotonVsKaons = 0.03, 
          deltaRMuons = 0.3, deltaRMuonsVsPhoton = 1.5;
  double deltaVzKaons = 0.05;
  double deltaVzMuons = 0.05;
  double deltaMassPhi = 0.010;
  //double deltaMassBs = 0.01;
  
  std::vector<double> massesKK = {0.493677, 0.493677};
  std::vector<double> massesMuMu = {0.105658, 0.105658};
  //std::vector<double> massesKKG = {0.493677, 0.493677, 0.0};
  //std::vector<double> massesMuMuG = {0.105658, 0.105658, 0.0};

  double PhiMass = 1.019455;
  //double BsMass = 5.36688;
  
  const auto & trackBuilder = es.getData(theTrackBuilderToken); 

  //////////////////Muons
  for(std::vector<pat::Muon>::const_iterator im1 = recoMuons.begin(); im1 < recoMuons.end(); im1++){
    //primary vertex ....
    
    if( !im1->isGlobalMuon() ) continue;
    
    hpTMuons -> Fill(im1->pt());
    
    if( im1->pt() < 1. ) continue;

    for (std::vector<pat::Muon>::const_iterator im2 = im1+1; im2 < recoMuons.end(); im2++) {
      if (!im2->isGlobalMuon() ) continue;
      if ( im2-> charge() * im1->charge() >= 0) continue;
      if( im2->pt() < 1. ) continue;
      //reco::TrackRef mu2Ref = im1->get<reco::TrackRef>();

      hMM_deltaR -> Fill( reco::deltaR(*im1, *im2));
      hMMDelta_vz -> Fill( abs(im1->vz() - im2->vz()));
      if(abs(im1->vz()-im2->vz()) > deltaVzMuons) continue;
      if( reco::deltaR(*im1, *im2) > deltaRMuons) continue;
      
      reco::TrackRef mu1Ref = im1->track();
      if (!mu1Ref) continue;
      reco::TrackRef mu2Ref = im2->track();
      if (!mu2Ref)continue;
      
      std::vector<reco::TransientTrack> muonTTs;
      
      muonTTs.push_back(trackBuilder.build(mu1Ref));
      muonTTs.push_back(trackBuilder.build(mu2Ref));
      KalmanVertexFitter kvf(true); //true?

      TransientVertex tvertex = kvf.vertex(muonTTs);
      reco::Vertex tvMM(TransientVertex(kvf.vertex(muonTTs)));

      double probability = TMath::Prob(tvMM.chi2(),tvMM.ndof());
      std::cout << "Muon vertex probability: " << probability << std::endl;
      hMM_prob -> Fill(probability);
      if(probability < 0.1) continue;

      std::vector<pat::Muon> mm = {*im1, *im2};
      double mMuMu = invariantMass(mm, massesMuMu);
      hMuMu_mass -> Fill(mMuMu); 
      // di-muon mass
      std::cout << "Di-muon mass: " << mMuMu << std::endl;
      

      if( fabs(mMuMu - PhiMass) > deltaMassPhi ) continue;
      double phiVz = (im1->vz() + im2->vz()) / 2.;

      //pv
      std::vector<double> distPCAs_mu, distVzs_mu;
      for( std::vector<reco::Vertex>::const_iterator ivtx = recoVertices.begin(); ivtx < recoVertices.end(); ivtx++){
        const reco::Vertex & vtx = *ivtx;
        
        math::XYZPoint cloasestApproach = pca(vtx.position(), tvMM.position(), (im1->momentum() + im2->momentum()));
        double deltaPCA = fabs(cloasestApproach.z() - vtx.position().z());
        distPCAs_mu.push_back(deltaPCA);

        double deltaVz = fabs(vtx.position().z() - phiVz);
        distVzs_mu.push_back(deltaVz);

      }
      // the smallest element in the vector
      auto minPCA_it = std::min_element(distPCAs_mu.begin(), distPCAs_mu.end());
      int minPCA_index = std::distance(distPCAs_mu.begin(), minPCA_it);
      std::cout << "minPCA_index: " << minPCA_index << "; minPCA distance: " << *minPCA_it 
                << "PCA: " << recoVertices[minPCA_index].position() << std::endl;
      auto minVz_it = std::min_element(distVzs_mu.begin(), distVzs_mu.end());
      int minVz_index = std::distance(distVzs_mu.begin(), minVz_it);
      std::cout << "minVz_index: " << minVz_index << "; minVz distance: " << *minVz_it 
                << "Vz: " << recoVertices[minVz_index].position() << std::endl;
      // do both methods lead to the same vertex?
      hPCAvsVz_mu -> Fill(minPCA_index == minVz_index);
      std::cout << "same PCA and Vz: " << (minPCA_index == minVz_index) << std::endl;
      std::vector<double> sortedPCAdist_mu = distPCAs_mu;
      std::sort(sortedPCAdist_mu.begin(), sortedPCAdist_mu.end());
      std::cout << "sorted PCA distances, first: " << sortedPCAdist_mu[0] 
                << " second: " << sortedPCAdist_mu[1] << std::endl;
      // difference between the smallest and second smallest distance between PCA and PV
      hPCA1vs2_mu -> Fill(abs(sortedPCAdist_mu[0] - sortedPCAdist_mu[1]));
      // the smallest distance between PCA and PV
      hPCA1smallest_mu -> Fill(sortedPCAdist_mu[0]);
      // if the chosen PV (PCA method) is not the first from the list
      
      /////////// DISPLACEMENT
      // PCA method used
      math::XYZPoint pvMinPCA_mu = recoVertices[minPCA_index].position();  // pozycja PV
      // Uzyskaj pozycję wierzchołka SV
      GlobalPoint svPos_mu = tvertex.position();
      // Oblicz wektor przemieszczenia (SV - PV)
      math::XYZPoint svPosXYZ_mu(svPos_mu.x(), svPos_mu.y(), svPos_mu.z());
      math::XYZVector displacement_mu = svPosXYZ_mu - pvMinPCA_mu;
      // Macierz kowariancji pozycji SV
      AlgebraicSymMatrix33 svCov_mu = tvertex.positionError().matrix();
      AlgebraicSymMatrix33 pvCov_mu = recoVertices[minPCA_index].covariance();
      AlgebraicSymMatrix33 totalCov_mu = svCov_mu + pvCov_mu;
      // Obliczanie niepewności przemieszczenia
      AlgebraicVector3 dispVec_mu(displacement_mu.x(), displacement_mu.y(), displacement_mu.z());
      double displacementError_mu = std::sqrt(ROOT::Math::Similarity(dispVec_mu, totalCov_mu));
      std::cout << "Displacement: (" 
                << displacement_mu.x() << ", "
                << displacement_mu.y() << ", "
                << displacement_mu.z() << "), lenght: " << displacement_mu.R()<< std::endl;
      std::cout << "Displacement error: " << displacementError_mu << std::endl;
      double displacementSignificance_mu = displacement_mu.r() / displacementError_mu;
      
      h3Ddisp_mu -> Fill(displacement_mu.R());
      h3DdispErr_mu -> Fill(displacementError_mu);
      h3DdispSignificance_mu -> Fill(displacementSignificance_mu);
      if(displacementSignificance_mu < 4) continue;
      hMuMu_mass_disp -> Fill(mMuMu);

      /////////// PHOTON
      for( std::vector<pat::Photon>::const_iterator ipho = recoPhotons.begin(); ipho < recoPhotons.end(); ipho++){
        auto phi_p4 = im1->p4() + im2->p4();
        double deltaRPhoton = reco::deltaR(phi_p4, ipho->p4());
        //if( deltaRPhoton > deltaRPhotonVsMuons) continue;
        //if(fabs(ipho->vz() - phiVz) > deltaVzMuons) continue;

        auto total_p4 = phi_p4 + ipho->p4();
        double mMMG = total_p4.M();
        hMuMuGmass -> Fill(mMMG);
      }
    }
  }

  
  /////////////////Loop over the kaons
  for( std:: vector<pat::PackedCandidate>::const_iterator ic1 = candidates.begin(); ic1 < candidates.end(); ic1++){
    
    if( !ic1->hasTrackDetails() ) continue; // skip if a bestTrack cannot be extracted from this Candidate
    if( abs(ic1->pdgId()) != 211 ) continue; // hadron (includes K+/-)
    hpTKaons -> Fill(ic1->pt());
    if( ic1->pt() < 1. ) continue; 
    
    for ( std:: vector<pat::PackedCandidate>::const_iterator ic2 = ic1+1; ic2 < candidates.end(); ic2++){
      if(ic1->charge() * ic2->charge() >= 0) continue; // opposite charges
      if( !ic2->hasTrackDetails() ) continue; // skip if a bestTrack can be extracted from this Candidate
      if( abs(ic2->pdgId()) != 211) continue; //K+/-
      if( ic2->pt() < 1. ) continue;

      //std::cout << "DeltaR: " << reco::deltaR(*ic1, *ic2) << std::endl;
      //std::cout << "DeltaVz: " << abs(ic1->vz() - ic2->vz()) << std::endl;
      hKK_deltaR -> Fill( reco::deltaR(*ic1, *ic2));
      hKKDelta_vz -> Fill( abs(ic1->vz() - ic2->vz()));

      if(abs(ic1->vz()-ic2->vz()) > deltaVzKaons) continue;
      if( reco::deltaR(*ic1, *ic2) > deltaRKaons) continue;
      
      const reco::Track & rtrk1 = ic1 -> pseudoTrack();
      const reco::Track & rtrk2 = ic2 -> pseudoTrack();
      if( rtrk1.charge()*rtrk2.charge() != -1 ) continue; //opposite charges
      std::vector<reco::TransientTrack> kaonTTs;
      kaonTTs.push_back(trackBuilder.build(rtrk1));
      kaonTTs.push_back(trackBuilder.build(rtrk2));

      KalmanVertexFitter kvf(true);
      TransientVertex tvertex = kvf.vertex(kaonTTs);
      reco::Vertex tvKK(TransientVertex(kvf.vertex(kaonTTs)));

      double probability = TMath::Prob(tvKK.chi2(),tvKK.ndof());
      std::cout << "Kaon vertex probability: " << probability << std::endl;
      hKK_prob -> Fill(probability);
      if(probability < 0.1) continue;

      std::vector<pat::PackedCandidate> kk = { *ic1, *ic2 };
      double mKK = invariantMass(kk, massesKK);
      hKK_mass -> Fill(mKK); // di-kaon mass
      
      std::cout << "Di-kaon mass: " << mKK << std::endl;
      std::cout << std::endl;
      if( fabs(mKK - PhiMass) > deltaMassPhi ) continue;
      //commonVtxKaons.push_back(kk);
      
      double phiVz = (ic1->vz() + ic2->vz()) / 2.;

      //pv
      
      std::vector<double> distPCAs, distVzs;

      for( std::vector<reco::Vertex>::const_iterator ivtx = recoVertices.begin(); ivtx < recoVertices.end(); ivtx++){
        const reco::Vertex & vtx = *ivtx;
        
        math::XYZPoint cloasestApproach = pca(vtx.position(), tvKK.position(), (ic1->momentum() + ic2->momentum()));
        double deltaPCA = fabs(cloasestApproach.z() - vtx.position().z());
        distPCAs.push_back(deltaPCA);

        double deltaVz = fabs(vtx.position().z() - phiVz);
        distVzs.push_back(deltaVz);

      }
      // the smallest element in the vector
      auto minPCA_it = std::min_element(distPCAs.begin(), distPCAs.end());
      int minPCA_index = std::distance(distPCAs.begin(), minPCA_it);
      std::cout << "minPCA_index: " << minPCA_index << "; minPCA distance: " << *minPCA_it 
                << "PCA: " << recoVertices[minPCA_index].position() << std::endl;

      auto minVz_it = std::min_element(distVzs.begin(), distVzs.end());
      int minVz_index = std::distance(distVzs.begin(), minVz_it);
      std::cout << "minVz_index: " << minVz_index << "; minVz distance: " << *minVz_it 
                << "Vz: " << recoVertices[minVz_index].position() << std::endl;

      // do both methods lead to the same vertex?
      hPCAvsVz -> Fill(minPCA_index == minVz_index);

      std::cout << "same PCA and Vz: " << (minPCA_index == minVz_index) << std::endl;
      std::vector<double> sortedPCAdist = distPCAs;
      std::sort(sortedPCAdist.begin(), sortedPCAdist.end());

      std::cout << "sorted PCA distances, first: " << sortedPCAdist[0] 
                << " second: " << sortedPCAdist[1] << std::endl;
      // difference between the smallest and second smallest distance between PCA and PV
      hPCA1vs2 -> Fill(abs(sortedPCAdist[0] - sortedPCAdist[1]));

      // the smallest distance between PCA and PV
      hPCA1smallest -> Fill(sortedPCAdist[0]);

      /////////// DISPLACEMENT

      // PCA method used
      math::XYZPoint pvMinPCA = recoVertices[minPCA_index].position();  // pozycja PV

      // Uzyskaj pozycję wierzchołka SV
      GlobalPoint svPos = tvertex.position();

      // Oblicz wektor przemieszczenia (SV - PV)
      math::XYZPoint svPosXYZ(svPos.x(), svPos.y(), svPos.z());
      math::XYZVector displacement = svPosXYZ - pvMinPCA;

      // Macierz kowariancji pozycji SV
      AlgebraicSymMatrix33 svCov = tvertex.positionError().matrix();
      AlgebraicSymMatrix33 pvCov = recoVertices[minPCA_index].covariance();
      AlgebraicSymMatrix33 totalCov = svCov + pvCov;
      
      // Obliczanie niepewności przemieszczenia
      AlgebraicVector3 dispVec(displacement.x(), displacement.y(), displacement.z());
      double displacementError = std::sqrt(ROOT::Math::Similarity(dispVec, totalCov));

      std::cout << "Displacement: (" 
                << displacement.x() << ", "
                << displacement.y() << ", "
                << displacement.z() << "), lenght: " << displacement.R()<< std::endl;

      std::cout << "Displacement error: " << displacementError << std::endl;

      double displacementSignificance = displacement.r() / displacementError;
      
      h3Ddisp -> Fill(displacement.R());
      h3DdispErr -> Fill(displacementError);
      h3DdispSignificance -> Fill(displacementSignificance);
      if(displacementSignificance < 4) continue;

      hKKmass_disp -> Fill(mKK);

      

      ///// PHOTON 
      for( std::vector<pat::Photon>::const_iterator ipho = recoPhotons.begin(); ipho < recoPhotons.end(); ipho++){

        auto phi_p4 = ic1->p4() + ic2->p4();
        double deltaR_photon_phi = reco::deltaR(phi_p4, ipho->p4());
        //if( deltaR_photon_phi > deltaRPhotonVsKaons ) continue;
        
        //if( fabs(ipho->vz() - phiVz) > deltaVzKaons ) continue;
        
        std::cout << "Photon pT: " << ipho->pt() << std::endl;
        hpTG -> Fill(ipho->pt());
        hpT_etaG -> Fill(ipho->eta(), ipho->pt() );

        auto p4_total = ic1->p4() + ic2->p4() + ipho->p4();
        hKKGmass -> Fill(p4_total.M()); // di-kaon-photon mass

        //std::cout << "Photon mass: " << invariantMass(kk, massesKKG) << std::endl;
      }

    }

  }
  
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(DataAnalysis);

