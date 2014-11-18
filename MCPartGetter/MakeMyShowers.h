/**
 * \file MakeMyShowers.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class MakeMyShowers
 *
 * @author david
 */

/** \addtogroup MCPartGetter

    @{*/

#ifndef LARLITE_MAKEMYSHOWERS_H
#define LARLITE_MAKEMYSHOWERS_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MakeMyShowers
     User custom analysis class made by david
   */
  class MakeMyShowers : public ana_base{
  
  public:

    /// Default constructor
    MakeMyShowers(){ _name="MakeMyShowers"; _fout=0;};

    /// Default destructor
    virtual ~MakeMyShowers(){};

    /** IMPLEMENT in MakeMyShowers.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MakeMyShowers.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MakeMyShowers.cc! 
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
