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
#include <numeric>


using namespace std;


//object definition
class HadronizationAnalysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit HadronizationAnalysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~HadronizationAnalysis();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool isSameDecay(const std::vector<int>&, const std::vector<int>&);
private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  TH1D *hPdgidPositive;
  TH1D *hPdgidNegative;

  edm::EDGetTokenT < vector<reco::GenParticle> > theGenParticleToken;

  std::vector<int> MuMuG = {22, 13, -13};

};


HadronizationAnalysis::HadronizationAnalysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;

  theGenParticleToken = consumes< vector<reco::GenParticle>  >( edm::InputTag("genParticles" ));

}

HadronizationAnalysis::~HadronizationAnalysis()
{
  cout <<" DTOR" << endl;
}

bool HadronizationAnalysis::isSameDecay(const std::vector<int>& dec1, const std::vector<int>& dec2) {
    
    if (dec1.size() != dec2.size()) {
        return false; 
    }

    std::set<int> dec1Set(dec1.begin(), dec1.end());
    std::set<int> dec2Set(dec2.begin(), dec2.end());

    return dec1Set == dec2Set;
}


void HadronizationAnalysis::beginJob()
{
  //create a histogram

  hPdgidPositive = new TH1D("hPdgidPos","Pdgid of particles; pdgid; Counts", 100, 499.5, 599.5);
  hPdgidNegative = new TH1D("hPdgidNeg","Pdgid of particles; pdgid; Counts", 100, -599.5, -499.5);


  cout << "HERE HadronizationAnalysis::beginJob()" << endl;
}

void HadronizationAnalysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data

  hPdgidPositive-> Write();
  hPdgidNegative-> Write();

  myRootFile.Close();

  delete hPdgidPositive;
  delete hPdgidNegative;

  cout << "HERE HadronizationAnalysis::endJob()" << endl;
}


void HadronizationAnalysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE HadronizationAnalysis::analyze "<< std::endl;

  const std::vector<reco::GenParticle> & genPar = ev.get(theGenParticleToken);

  int nBmesons = 0;
  int nBd = 0;
  int nBs = 0;
  int nMuMuGamma = 0;

  int BmesonPdgId = 0;
  vector<vector<const reco::Candidate*>> ancestors;

  std::cout << "Number of particles: " << genPar.size() << std::endl;
  //for (std::vector<reco::GenParticle>::const_iterator part = genPar.begin(); part < genPar.end(); part++) {}

  for (const auto& part:genPar ){
    
    if ( abs(part.pdgId())/100 == 5){

      //std::cout << std::endl;
      //std::cout << "B meseon:  " << part.pdgId() << std::endl;
      int nDaughters = part.numberOfDaughters();

      std::vector<std::vector<int>> decayChannel(2, std::vector<int>(nDaughters, 0));

      //std::cout << "\tDaughters:  " ;

      for (int i = 0; i < nDaughters; ++i) {
        const reco::Candidate* daughter = part.daughter(i);

        decayChannel[0][i] = daughter->pdgId();

        if ( abs(decayChannel[0][i])/100 != 5 ) decayChannel[1][i] = 1;

        //std::cout << part.daughter(i)->pdgId() << "  ";
      }
      bool isFinal=false;
      int sum = std::accumulate(decayChannel[1].begin(), decayChannel[1].end(), 0);

      if (sum == nDaughters) isFinal = true;

      if(isFinal) nBmesons++;


      // find the B mesons that are not Bs forced to decay to mu mu gamma
      if(isFinal && !(isSameDecay(decayChannel[0], MuMuG) && abs(part.pdgId()) == 531))
      {
        const reco::Candidate* mother = &part;
        while (abs(mother->mother(0)->pdgId())/100 == 5){
          mother = mother->mother(0);
        }
        BmesonPdgId = mother->pdgId();

        // get the ancestors until the proton
        vector<const reco::Candidate*> ancestorLevel;
        ancestorLevel.push_back(mother);
        ancestors.push_back(ancestorLevel);
        while (true)
        {
          ancestorLevel.clear();
          bool endLoop = true;
          for(auto ancestor : ancestors.back())
          {
            if (ancestor->numberOfMothers() > 0) endLoop = false;
            for(unsigned int i = 0; i < ancestor->numberOfMothers(); i++)
            {
              ancestorLevel.push_back(ancestor->mother(i));
            }
          }
          if (endLoop) break;
          ancestors.push_back(ancestorLevel);
        }
      }
    } // end of loop over B mesons
  } // end of loop over gen particles

  if(nBmesons == 2 && BmesonPdgId > 0) hPdgidPositive->Fill(BmesonPdgId);
  if(nBmesons == 2 && BmesonPdgId < 0) hPdgidNegative->Fill(BmesonPdgId);

  if(nBmesons ==2)
  {
    cout << "B meson pdgId: " << BmesonPdgId << endl;
    cout << "Ancestors: " << endl;
    for (const auto& ancestorLevel:ancestors)
    {
      for (const auto& ancestor:ancestorLevel)
      {
        cout << "  " << ancestor->pdgId();
      }
      cout << endl;
    }
  }

  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(HadronizationAnalysis);

