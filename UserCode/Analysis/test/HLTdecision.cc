#include <iostream>
#include <fstream>
#include "HLTdecision.h"


HLTdecision::HLTdecision(edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken,
                         const edm::Event& ev_,
                         const edm::ParameterSet& conf)
    : triggerResults(ev_.get(triggerResultsToken)),
      triggerNames(ev_.triggerNames(triggerResults)),
      trgOfInterest(std::move(conf.getParameter<std::vector<std::string>>("trg"))) {}

HLTdecision::~HLTdecision(){
    // Destructor implementation
}

bool HLTdecision::Passed(const std::string& triggerName) {
    unsigned int trgIndex = triggerNames.triggerIndex(triggerName);
    return triggerResults.accept(trgIndex);
}

bool HLTdecision::checkTriggers(const edm::Event& ev,
                                bool print) {
    bool accepted = false;

    if (print) std::cout << "Trigger results: ";

    for (const auto& trgIter : trgOfInterest){
        bool decision = Passed(trgIter);
        if ( decision ) accepted = true;
        if (print) std::cout << trgIter << ": " << decision << std::endl;
    }

    if (!accepted) {
        std::cout << "Skipping event: " << ev.id() << " No trigger fired" << std::endl;
        return false;
    }
    return accepted;
}
/*
void HLTdecision::printAllPaths(){
    std::cout << "Fired trigger paths:" << std::endl;
    for (const auto& trgIter : triggerNames) {
        
        if (Passed(trgIter)){
            std::cout << "      " << trgIter << std::endl;
        }
        
    }
}*/

void HLTdecision::printAllPaths(std::ofstream& file) {
    if (!file) { 
        std::cerr << "Error: File stream is not open!" << std::endl;
        return;
    }

    file << "Fired trigger paths:\n";
    std::cout << "Fired trigger paths:" << std::endl;

    for (unsigned int i = 0; i < triggerNames.size(); ++i) {
        std::string trgName = triggerNames.triggerName(i);
        if (Passed(trgName)) {
            file << "      " << trgName << "\n";
            std::cout << "      " << trgName << std::endl;
        }
    }
    file.close();
}
