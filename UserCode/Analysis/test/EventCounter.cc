#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include <sstream>
#include <iomanip> 
#include <utility>
#include <numeric>
#include <iostream>


using namespace std;


//object definition
class EventCounter : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit EventCounter(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~EventCounter();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

};


EventCounter::EventCounter(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;
}

EventCounter::~EventCounter()
{
  cout <<" DTOR" << endl;
}



void EventCounter::beginJob()
{
  cout << "HERE EventCounter::beginJob()" << endl;
}

void EventCounter::endJob()
{

  cout << "HERE EventCounter::endJob()" << endl;
}


void EventCounter::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE EventCounter::analyze "<< std::endl;
  cout <<"*** Analyze event: " << ev.id() <<" analysed event count:" << ++theEventCount << endl;
}

DEFINE_FWK_MODULE(EventCounter);

