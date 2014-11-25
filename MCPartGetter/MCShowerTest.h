/**
 * \file MCShowerTest.h
 *
 * \ingroup MCInfo
 * 
 * \brief Class def header for a class MCShowerTest
 *
 * @author David Caratelli
 */

/** \addtogroup MCInfo
    
    @{*/

#ifndef MCSHOWERTEST_H
#define MCSHOWERTEST_H

#include "Analysis/ana_base.h"
#include "LArUtil/Geometry.h"
#include "ShowerCutCalculator.h"
#include <vector>
#include <string>
#include <time.h>
#include <stdlib.h>

namespace larlite {
  /**
     \class MCShowerTest
     User custom analysis class made by david
  */
  class MCShowerTest : public ana_base{
    
  public:

/// Default constructor
    MCShowerTest(){ _name="MCShowerTest"; _fout=0; _verbose=false; };

    /// Default destructor
    virtual ~MCShowerTest(){};

    /** IMPLEMENT in MCShowerTest.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MCShowerTest.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MCShowerTest.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetVerbose(bool on) { _verbose = on; }

    /// prepare TTree
    void prepareTree();

    void resetTree();

    double addTrack(mctrack track);

    protected:

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
    std::vector<int> _allTrackIDs;


    double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;


    //Cut Distance
    double _cutDist;


	TTree * _ana_tree ;
		
	int _run ;
	int _subrun ;
	int _event ;

	std::string _process ;
	int _PDG ;
	int _trackID ;

	double _X ;
	double _Y ;
	double _Z ;
	double _Xmc ;
	double _Ymc ;
	double _Zmc ;

	double _T ;

	double _Px ;
	double _Py ;
	double _Pz ;
	double _Pxmc ;
	double _Pymc ;
	double _Pzmc ;

	double _E ;

	int _inActiveVolume ;

	double _distAlongTraj ;
	double _distBackAlongTraj ;

	double _minMuDist;
	double _minMuIP;
	double _distToIP;

	double _minMuDistExceptAncestor;
	double _minMuIPExceptAncestor;
	double _distToIPExceptAncestor;

	double _ancDist;
	double _ancIP;
	double _ancToIP;

	double _ancDistmc;
	double _ancIPmc;
	double _ancToIPmc;

	//Save info about parent as well
	int _parentPDG ;
	double _parentX ;
	double _parentY ;
	double _parentZ ;
	double _parentT ;

    double _parentPx;
	double _parentPy;
	double _parentPz ; 
	double _parentE ; 

	int _parentInActiveVolume ;


	//Save info about ancestor as well
	int _ancestorPDG ;
	double _ancestorX ;
	double _ancestorY ;
	double _ancestorZ ;
	double _ancestorT ;

	double _ancestorPx;
	double _ancestorPy;
	double _ancestorPz ; 
	double _ancestorE ; 

	int _ancTrackFound;

	int _ancestorInActiveVolume ;

	std::vector<std::vector<double> > ShowerTraj;
	std::vector<std::vector<double> > MotherTraj;
	std::vector<std::vector<double> > AncestorTraj;

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
