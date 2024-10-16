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

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <sstream>
#include <iomanip> 
#include <utility>
#include <typeinfo>


using namespace std;


//object definition
class Analysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit Analysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~Analysis();

  //edm filter plugin specific functions
  virtual void beginJob();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  TH1D *BspT;
  TH1D *BpT;
  TH1D *BToBs;

  TH1D *mupT;
  TH1D *gammapT;

  TH1D *muEta;
  TH1D *gammaEta;

  TH2D *mu1Eta_mu2Eta;

  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;


};


Analysis::Analysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles" ));

}

Analysis::~Analysis()
{
  cout <<" DTOR" << endl;
}

void Analysis::beginJob()
{
  //create a histogram
  BspT = new TH1D("BspT","Bs pT; pT [GeV]; #events",100, 0., 50.);
  BpT = new TH1D("BpT","B (not Bs); pT [GeV]; #events",100, 0., 50.);
  BToBs = new TH1D("BsToB","Bs to B; pT [GeV]; #events",100, -1, 4);

  mupT = new TH1D("mupT","mu pT; pT [GeV]; #events",100, 0., 50.);
  gammapT = new TH1D("gammapT","gamma pT; pT [GeV]; #events",100, 0., 50.);

  muEta = new TH1D("muEta","mu eta; eta; #events",100, -4., 4.);
  gammaEta = new TH1D("gammaEta","gamma eta; eta; #events",100, -4., 4.);

  mu1Eta_mu2Eta = new TH2D("mu1Eta_mu2Eta","mu1 eta vs mu2 eta; mu1 eta; mu2 eta; #events",100, -4., 4., 100, -4., 4.);

  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  BspT -> Write();
  BpT -> Write();
  BToBs -> Write();

  mupT -> Write();
  gammapT -> Write();

  muEta -> Write();
  gammaEta -> Write();

  mu1Eta_mu2Eta -> Write();

  myRootFile.Close();


  cout << "HERE Analysis::endJob()" << endl;
}


void Analysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE Analysis::analyze "<< std::endl;
  //bool debug = true;

  const vector<reco::GenParticle> & genParv = ev.get(theGenParticleToken);
  int nBs  = 0;
  int nB   = 0;

  for(auto genPar : genParv)
  {
    if(abs(genPar.pdgId())/100 == 5)
    {
      cout << "Particle PDGid: " << genPar.pdgId() << endl;
      cout << "\t daughter PDGid:";

      bool isNotFinal = false;

      vector<float> muEtaV;
      muEtaV.clear();
      for(unsigned int i=0; i < genPar.numberOfDaughters(); i++)
      {
        int d_pdgid = genPar.daughterRef(i)->pdgId();
        cout << " " << d_pdgid;
        if (abs(d_pdgid)/100 == 5) isNotFinal = true;

        if (abs(genPar.pdgId() == 531) && abs(d_pdgid) == 13)
        {
          mupT->Fill(genPar.daughterRef(i)->pt());
          muEta->Fill(genPar.daughterRef(i)->eta());
          muEtaV.push_back(genPar.daughterRef(i)->eta());
        }

        if (abs(genPar.pdgId() == 531) && abs(d_pdgid) == 22)
        {
          gammapT->Fill(genPar.daughterRef(i)->pt());
          gammaEta->Fill(genPar.daughterRef(i)->eta());
        }

      }
      cout << endl << endl;

      if (muEtaV.size() == 2)
      {
        mu1Eta_mu2Eta->Fill(muEtaV[0], muEtaV[1]);
      }

      if (abs(genPar.pdgId()) == 531 && !isNotFinal)
      {
        BspT->Fill(genPar.pt());
        nBs++;
      }
      if (abs(genPar.pdgId()) == 511 && !isNotFinal)
      {
        BpT->Fill(genPar.pt());
        nB++;
      }
    }


  } // particle loop

  BToBs->Fill((float)nB/nBs);
  cout <<"*** Analyzed event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl << endl;
}

DEFINE_FWK_MODULE(Analysis);

