/**
 * \file MakeMyTracks.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class MakeMyTracks
 *
 * @author david
 */

/** \addtogroup MCPartGetter

    @{*/

#ifndef LARLITE_MAKEMYTRACKS_H
#define LARLITE_MAKEMYTRACKS_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MakeMyTracks
     User custom analysis class made by david
   */
  class MakeMyTracks : public ana_base{
  
  public:

    /// Default constructor
    MakeMyTracks(){ _name="MakeMyTracks"; _fout=0;};

    /// Default destructor
    virtual ~MakeMyTracks(){};

    /** IMPLEMENT in MakeMyTracks.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MakeMyTracks.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MakeMyTracks.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    protected:

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
