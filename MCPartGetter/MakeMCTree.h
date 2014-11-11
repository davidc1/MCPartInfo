/**
 * \file MakeMCTree.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class MakeMCTree
 *
 * @author david
 */

/** \addtogroup MCPartGetter

    @{*/

#ifndef LARLITE_MAKEMCTREE_H
#define LARLITE_MAKEMCTREE_H

#include "Analysis/ana_base.h"
#include "DataFormat/treenode.h"

namespace larlite {
  /**
     \class MakeMCTree
     User custom analysis class made by david
   */
  class MakeMCTree : public ana_base{
  
  public:

    /// Default constructor
    MakeMCTree(){ _name="MakeMCTree"; _fout=0;};

    /// Default destructor
    virtual ~MakeMCTree(){};

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setParticleMap(event_mctree *event_tree, event_mcpart *event_part);

    int searchParticleMap(std::set<int>::iterator it){ return _ParticleMap.find(*it)->second; }

    int searchParticleMap(const int trackID);

    void setTrees(event_mctree *event_tree, event_mcpart *event_part);

    void AddNodes(mcpart part, treenode& parentnode, int ancestor, event_mcpart *event_part);

    protected:

    std::map<int, int> _ParticleMap;

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
