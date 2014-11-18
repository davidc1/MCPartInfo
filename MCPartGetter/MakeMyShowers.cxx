#ifndef MAKEMYSHOWERS_CXX
#define MAKEMYSHOWERS_CXX

#include "MakeMyShowers.h"

namespace larlite {

  bool MakeMyShowers::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool MakeMyShowers::analyze(storage_manager* storage) {


    // get MCShowers
    auto evt_mcshower = storage->get_data<event_mcshower>("mcreco");  
    
    // make my MCShowers
    auto my_mcshower = storage->get_data<event_mcshower>("davidc1");  

    *(my_mcshower) = *(evt_mcshower);

    my_mcshower->set_subrun(evt_mcshower->subrun());
    my_mcshower->set_run(evt_mcshower->run());
    my_mcshower->set_event_id(evt_mcshower->event_id());

  
    return true;
  }

  bool MakeMyShowers::finalize() {

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }

}
#endif
