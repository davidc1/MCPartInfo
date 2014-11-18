#ifndef MAKEMYTRACKS_CXX
#define MAKEMYTRACKS_CXX

#include "MakeMyTracks.h"

namespace larlite {

  bool MakeMyTracks::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool MakeMyTracks::analyze(storage_manager* storage) {
  
    // get MCTracks
    auto evt_mctracks = storage->get_data<event_mctrack>("mcreco");
    
    // make my MCTracks
    auto my_mctracks = storage->get_data<event_mctrack>("davidc1");

    *(my_mctracks) = *(evt_mctracks);

    my_mctracks->set_subrun(evt_mctracks->subrun());
    my_mctracks->set_run(evt_mctracks->run());
    my_mctracks->set_event_id(evt_mctracks->event_id());


    return true;
  }

  bool MakeMyTracks::finalize() {

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
