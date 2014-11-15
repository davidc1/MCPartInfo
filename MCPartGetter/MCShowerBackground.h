/**
 * \file MCShowerBackground.h
 *
 * \ingroup MCInfo
 * 
 * \brief Class def header for a class MCShowerBackground
 *
 * @author David Caratelli
 */

/** \addtogroup MCInfo
    
    @{*/

#ifndef MCSHOWERBACKGROUND_H
#define MCSHOWERBACKGROUND_H

#include "Analysis/ana_base.h"
#include "LArUtil/Geometry.h"
#include "ShowerCutCalculator.h"
#include <vector>
#include <string>
#include <time.h>

namespace larlite {
  /**
     \class MCShowerBackground
     User custom analysis class made by david
  */
  class MCShowerBackground : public ana_base{
    
  public:

/// Default constructor
    MCShowerBackground(){ _name="MCShowerBackground"; _fout=0; _verbose=false; };

    /// Default destructor
    virtual ~MCShowerBackground(){};

    /** IMPLEMENT in MCShowerBackground.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MCShowerBackground.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MCShowerBackground.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetVerbose(bool on) { _verbose = on; }

    /// prepare TTree
    void prepareTree();

    void resetTree();

    double addTrack(mctrack track);

    protected:

    /// Process to be searched
    std::vector<std::pair<int,std::string> > _process;

    /// verbose
    bool _verbose;

    /// double Energy cut
    double _Ecut;

    /// Event number
    int _evtN;

    /// Instance of ShowerCutCalculator
    ShowerCutCalculator _cutParamCalculator;

    /// All muons tracks
    std::vector<std::vector<std::vector<double> > > _allTracks;
    


    //Cut Distance
    double _cutDist;

    /// Tree for analysis
    TTree* _tree;
    //variables for analysis
    int    _isPrimary;       ///< 1 = is Primary. 0 = not.
    double _E;               ///< shower energy
    int    _inTPC;           ///< shower is in TPC
    double _aE;              ///< Energy of ancestor
    int    _aPdg;            ///< PDG of ancestor
    double _mE;              ///< Energy of mother
    int    _mPdg;            ///< PDG of mother
    double _xStart;          ///< shower start point x
    double _yStart;          ///< shower start point y
    double _zStart;          ///< shower start point z
    double _px;              ///< shower initial Px
    double _py;              ///< shower initial Py
    double _pz;              ///< shower initial Pz
    double _tStart;          ///< shower start time
    std::string proc;        ///< shower process
    int _pdg;                ///< PDG code of shower
    double _muDist;          ///< Distance to closest muon
    double _muIP;            ///< Impact Param to closest muon
    double _distToIP;        ///< Start to IP on track distance
    double _forwardToWall;   ///< shower forward distance to TPC wall
    double _backToWall;      ///< shower backwards distance to TPC wall

    std::vector<std::vector<double> > MuonTraj;

    //histogram for muon track length
    TH1D *_hMuonTotLen;

    // evaluate time-performance
    clock_t t;
     

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
