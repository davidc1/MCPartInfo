/**
 * \file MCPartExample.h
 *
 * \ingroup MCInfo
 * 
 * \brief Class def header for a class MCPartExample
 *
 * @author David Caratelli
 */

/** \addtogroup MCInfo

    @{*/

#ifndef MCPARTEXAMPLE_H
#define MCPARTEXAMPLE_H

#include "Analysis/ana_base.h"
#include "MCgetter.h"
#include "LArUtil/Geometry.h"
#include <string>

namespace larlite {
  /**
     \class MCPartExample
     User custom analysis class made by david
   */
  class MCPartExample : public ana_base{
  
  public:

    /// Default constructor
    MCPartExample(){ _name="MCPartExample"; _fout=0; _verbose=false; SetProperties(); };

    /// Default destructor
    virtual ~MCPartExample(){};

    /** IMPLEMENT in MCPartExample.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MCPartExample.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MCPartExample.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetVerbose(bool on) { _verbose = on; }

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
    int    _MotherPDG;
    double _MotherDist;
    double _StartX;
    double _StartY;
    double _StartZ;
    double _StartT;
    double _EndX;
    double _EndY;
    double _EndZ;
    double _minMuonDist;
    std::vector<std::vector<double> > PartTraj;
    std::vector<std::vector<double> > MotherTraj;
    std::vector<std::vector<double> > AncestorTraj;
    std::string Process;
    std::string ProcHist;

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
