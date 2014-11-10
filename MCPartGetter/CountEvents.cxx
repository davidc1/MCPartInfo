#ifndef COUNTEVENTS_CXX
#define COUNTEVENTS_CXX

#include "CountEvents.h"

namespace larlite {

  bool CountEvents::initialize() {

    _hEvents = new TH1I("hEvents","Number of Events in File",30,0,30);

    _evtN = 0;
    
    return true;
  }

  bool CountEvents::analyze(storage_manager* storage) {

    _evtN += 1;

    return true;
  }

  bool CountEvents::finalize() {

    _hEvents->Fill(_evtN);
    _hEvents->Write();

    return true;
  }

}
#endif
