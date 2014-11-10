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

#ifndef COSMICTRACKS_H
#define COSMICTRACKS_H

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
     \class CosmicTracks
     User custom analysis class made by david
  */
  class CosmicTracks : public ana_base{
    
  public:

/// Default constructor
    CosmicTracks(){ _name="CosmicTracks"; _fout=0; _verbose=false; SetProperties(); };

    /// Default destructor
    virtual ~CosmicTracks(){};

    /** IMPLEMENT in CosmicTracks.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CosmicTracks.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CosmicTracks.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void getMaxAngle(std::vector<std::vector<double> > *track, std::vector<double> *energies);

    double getAngle(std::vector<double> dir1, std::vector<double> dir2);

    double getAngle(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2);

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


    // muon tree
    TTree *_muontree;
    // variables
    double _muonE;
    int _muonPDG;
    double _muonStartX;
    double _muonStartY;
    double _muonStartZ;
    double _muonEndX;
    double _muonEndY;
    double _muonEndZ;
    std::vector<std::vector<double> > MuonTraj;
    std::vector<double> MuonEnergies;
    double _maxDeflection;
    double _maxDeflectionX;
    double _maxDeflectionY;
    double _maxDeflectionZ;
    double _deflectionE;
    double _deflectionDistToEnd;
    double _dirChange;
    int _inTPC;

    //histogram for muon track length
    TH1D *_hMuonTotLen;
     

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
