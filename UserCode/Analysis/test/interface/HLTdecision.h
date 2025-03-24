#ifndef HLTDECISION_H
#define HLTDECISION_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include <vector>
#include <string>

class HLTdecision {
    const edm::TriggerResults& triggerResults;
    const edm::TriggerNames& triggerNames;
    const std::vector<std::string> trgOfInterest; 

public:
    HLTdecision(edm::EDGetTokenT<edm::TriggerResults>, const edm::Event&, const edm::ParameterSet&);
    ~HLTdecision();

    bool Passed(const std::string& );

    bool checkTriggers( const edm::Event& , bool print = false );
    void printAllPaths( std::ofstream&);

};

#endif // HLTDECISION_H
