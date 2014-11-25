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
#include <string>
#include <algorithm>

namespace larlite {
  /**
     \class MakeShowers
     User custom analysis class made by david
  */
  class MakeShowers : public ana_base{
    
  public:

/// Default constructor
    MakeShowers(){ _name="MakeShowers"; _fout=0; _verbose=false; _Ecut = 0; };

    /// Default destructor
    virtual ~MakeShowers(){};

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetVerbose(bool on) { _verbose = on; }

    void SetEcut(double E) { _Ecut = E; }

    void findMCShowers(treenode tree,
		       event_mcpart *evt_part,
		       event_mctree *evt_tree,
		       event_mcshower *evt_mcshower);

    void makeMCShower(treenode tree,
		      event_mcpart *evt_part,
		      event_mctree *evt_tree,
		      event_mcshower *evt_mcshower);

    void getAllShowerParticles(treenode tree,
			       std::vector<unsigned int> &trackIDs,
			       event_mcpart *evt_part,
			       event_mctree *evt_tree,
			       treenode toptree);

    void PrepareTree();
    void ResetTree();
    void resetPartTree();
    

    protected:

    /// verbose
    bool _verbose;

    /// Event number
    int _evtN;

    /// Energy cut for shower creation [ GeV ]
    double _Ecut;


    // Shower Tree information
    TTree *_showertree;
    // variables for tree
    int    _showerTrackID;
    double _showerE;
    int    _showerPDG;
    double _showerStartX;
    double _showerStartY;
    double _showerStartZ;
    std::vector<std::vector<std::vector<double> > > ShowerTraj;
    std::vector<std::vector<double> > AncestorTraj;
    int    _inTPC;
    int    _eventN;
    std::string _showerProcess;
    int    _numEl;
    int    _numCompt;
    int    _numConv;
    int    _numComptE;
    int    _numConvE;

    // Part Tree information
    TTree *_parttree;
    // variables for tree
    int    _partTrackID;
    double _partE;
    int    _partPDG;
    std::string _partProcess;
    std::vector<std::vector<double> > PartTraj;
    int    _shrID;
    int    _shrPDG;
    double _shrE;
    std::string _shrProc;

    //mcgetter
    MCgetter _MCgetter;

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
