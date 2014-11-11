/**
 * \file ComptonBackground.h
 *
 * \ingroup MCInfo
 * 
 * \brief Class def header for a class ComptonBackground
 *
 * @author David Caratelli
 */

/** \addtogroup MCInfo
    
    @{*/

#ifndef MAKESHOWERS_H
#define MAKESHOWERS_H

#include "Analysis/ana_base.h"
#include "MCgetter.h"
#include "BasicTool/GeoAlgo/TrajectoryInVolume.h"
#include "BasicTool/GeoAlgo/PointToLineDist.h"
#include "BasicTool/GeoAlgo/TwoLineIntersection.h"
#include "BasicTool/GeoAlgo/SegmentPoCA.h"
#include "LArUtil/Geometry.h"
#include <string>

namespace larlite {
  /**
     \class MakeShowers
     User custom analysis class made by david
  */
  class MakeShowers : public ana_base{
    
  public:

/// Default constructor
    MakeShowers(){ _name="MakeShowers"; _fout=0; _verbose=false; SetProperties(); };

    /// Default destructor
    virtual ~MakeShowers(){};

    /** IMPLEMENT in MakeShowers.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MakeShowers.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MakeShowers.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetVerbose(bool on) { _verbose = on; }

    void SetProperties();

    /// Set MCgetter object
    void SetMCgetter(MCgetter mcgetter) { _MCgetter = mcgetter; }

    /// prepare TTree
    void prepareTree();

    void ResetTree();

    protected:

    /// Process to be searched
    std::vector<std::pair<int,std::string> > _process;

    /// verbose
    bool _verbose;

    /// double Energy cut
    double _Ecut;

    /// Event number
    int _evtN;

    /// MCgetter to make mc particle map
    MCgetter _MCgetter;

    /// GeoAlg for TPC containment
    geoalgo::TrajectoryInVolume _inTPCAlgo;
    /// GeoAlg for point to line dist
    geoalgo::PointToLineDist _pointDist;
    /// GeoAlg for poka cut
    geoalgo::TwoLineIntersection _lineIntersection;
    /// GeoAlg for PoCA cut
    geoalgo::SegmentPoCA _PoCA;


    // shower tree
    TTree *_showertree;
    // variables
    int _showerDaughters;
    double _showerE;
    int _showerPDG;
    double _showerStartX;
    double _showerStartY;
    double _showerStartZ;
    std::vector<std::vector<std::vector<double> > > ShowerTraj;
    int _inTPC;

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
