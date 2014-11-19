/**
 * \file GetPi0Distributions.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class GetPi0Distributions
 *
 * @author david
 */

/** \addtogroup MCPartGetter

    @{*/

#ifndef LARLITE_GETPI0DISTRIBUTIONS_H
#define LARLITE_GETPI0DISTRIBUTIONS_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class GetPi0Distributions
     User custom analysis class made by david
   */
  class GetPi0Distributions : public ana_base{
  
  public:

    /// Default constructor
    GetPi0Distributions(){ _name="GetPi0Distributions"; _fout=0;};

    /// Default destructor
    virtual ~GetPi0Distributions(){};

    /** IMPLEMENT in GetPi0Distributions.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GetPi0Distributions.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GetPi0Distributions.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clearTree();

    protected:

    TTree *_tree;

    int _Nshowers;
    double _totE;
    double _ShowerE1;
    double _ShowerE2;
    double _Angle;
    double _VtxDist1;
    double _VtxDist2;

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
