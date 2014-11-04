/**
 * \file sample_ana.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class sample_ana
 *
 * @author david
 */

/** \addtogroup MCPartGetter

    @{*/

#ifndef LARLITE_MCPARTGETTER_SAMPLE_ANA_H
#define LARLITE_MCPARTGETTER_SAMPLE_ANA_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class sample_ana
     User custom analysis class made by david
   */
  class sample_ana : public ana_base{
  
  public:

    /// Default constructor
    sample_ana() : ana_base() { _name="sample_ana"; }

    /// Default destructor
    virtual ~sample_ana(){};

    /** IMPLEMENT in sample_ana.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in sample_ana.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in sample_ana.cc! 
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
