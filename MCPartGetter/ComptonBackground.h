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

#ifndef COMPTONBACKGROUND_H
#define COMPTONBACKGROUND_H

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
     \class ComptonBackground
     User custom analysis class made by david
  */
  class ComptonBackground : public ana_base{
    
  public:

/// Default constructor
    ComptonBackground(){ _name="ComptonBackground"; _fout=0; _verbose=false; SetProperties(); };

    /// Default destructor
    virtual ~ComptonBackground(){};

    /** IMPLEMENT in ComptonBackground.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ComptonBackground.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ComptonBackground.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetVerbose(bool on) { _verbose = on; }

    void fillPoCAParams(std::vector<double> ePoCA, std::vector<double> eStart, std::vector<double> eDir);

    /// Function to set physics process to be analyzed
    void SetProcess(std::vector<int> PDGs, std::vector<std::string> procs);

    /// Set MCgetter object
    void SetMCgetter(MCgetter mcgetter) { _MCgetter = mcgetter; }

    /// Set Distance for Cut for muons
    void SetCutDist(double d) { _cutDist = d; }

    /// prepare TTree
    void prepareTree();

    /// Set All Properties
    void SetProperties();

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



    //Cut Distance
    double _cutDist;

    /// Tree for analysis
    TTree* _tree;
    //variables for analysis
    int    _isPrimary;
    double _Energy;
    int    _inTPC;
    double _AncestorE;
    int    _AncestorPDG;
    double _AncestorDist;
    double _MotherE;
    double _MotherEndE;
    int    _MotherPDG;
    double _MotherDist;
    double _StartX;
    double _StartY;
    double _StartZ;
    double _PX;
    double _PY;
    double _PZ;
    double _StartT;
    double _EndX;
    double _EndY;
    double _EndZ;
    double _PoCAtoAncestor;
    double _PoCAtoAncestorDist;
    double _minMuonDist;
    double _minMuonPoka;
    double _PoCADist;
    int _PoCADistAfterStart;
    std::vector<std::vector<double> > PartTraj;
    std::vector<std::vector<double> > MotherTraj;
    std::vector<std::vector<double> > AncestorTraj;
    std::string Process;
    std::string ProcHist;
    int _PDG;

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
