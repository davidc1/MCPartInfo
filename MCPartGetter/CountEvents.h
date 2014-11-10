/**
 * \file CountEvents.h
 *
 * \ingroup MCInfo
 * 
 * \brief Class def header for a class CountEvents
 *
 * @author David Caratelli
 */

/** \addtogroup MCInfo

    @{*/

#ifndef COUNTEVENTS_H
#define COUNTEVENTS_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class CountEvents
     User custom analysis class made by david
   */
  class CountEvents : public ana_base{
  
  public:

    /// Default constructor
    CountEvents(){ _name="CountEvents"; _fout=0; };

    /// Default destructor
    virtual ~CountEvents(){};

    /** IMPLEMENT in CountEvents.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CountEvents.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CountEvents.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  private:
      
    TH1I* _hEvents;
    int _evtN;

  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
