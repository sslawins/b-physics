#include <iostream>
#include "HLTdecision.h"


HLTdecision::HLTdecision(edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken,
                         const edm::Event& ev_,
                         const edm::ParameterSet& conf)
    : triggerResults(ev_.get(triggerResultsToken)),
      triggerNames(ev_.triggerNames(triggerResults)),
      trgOfInterest(std::move(conf.getParameter<std::vector<std::string>>("trg"))) {}

HLTdecision::~HLTdecision() {
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

    for (const auto& trgIter : trgOfInterest) {
        bool decision = Passed(trgIter);
        if (decision) accepted = true;
        if (print) std::cout << trgIter << ": " << decision << std::endl;
    }
    if (!accepted) {
        std::cout << "Skipping event: " << ev.id() << " No trigger fired" << std::endl;
        return false;
    }
    return accepted;
}

