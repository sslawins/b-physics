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

  TH1D *pT;

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
  pT = new TH1D("pT","pT; pT [GeV]; #events",100, 0., 50.);


  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  pT -> Write();

  myRootFile.Close();

  delete pT;


  cout << "HERE Analysis::endJob()" << endl;
}


void Analysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE Analysis::analyze "<< std::endl;
  //bool debug = true;

  const vector<reco::GenParticle> & genParv = ev.get(theGenParticleToken);

  for(auto genPar : genParv)
  {
    if(abs(genPar.pdgId())/100 == 5)
    {
      pT->Fill(genPar.pt());
      cout << "Particle id: " << genPar.pdgId() << endl;
      cout << "\t daughter id:";
      for(unsigned int i=0; i < genPar.numberOfDaughters(); i++)
      {
        cout << " " << genPar.daughterRef(i)->pdgId();
      }
      cout << endl << endl;
    }

  }
  cout <<"*** Analyzed event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl << endl;
}

DEFINE_FWK_MODULE(Analysis);

