#ifndef MAKESHOWERS_CXX
#define MAKESHOWERS_CXX

#include "MakeShowers.h"
#include <time.h>

namespace larlite {

  bool MakeShowers::initialize() {

    if (_verbose) { _MCgetter.SetVerbose(true); }

    /// Prepare Tree
    prepareTree();

    _evtN = 0;

    return true;
  }


  void MakeShowers::prepareTree(){

    if(_showertree) delete _showertree;
    _showertree = new TTree("shower_tree","");
    _showertree->Branch("_showerDaughters",&_showerDaughters,"showerDaughters/I");
    _showertree->Branch("_showerE",&_showerE,"showerE/D");
    _showertree->Branch("_showerPDG",&_showerPDG,"showerPDG/I");
    _showertree->Branch("_showerStartX",&_showerStartX,"showerStartX/D");
    _showertree->Branch("_showerStartY",&_showerStartY,"showerStartY/D");
    _showertree->Branch("_showerStartZ",&_showerStartZ,"showerStartZ/D");
    _showertree->Branch("ShowerTraj",&ShowerTraj);
    _showertree->Branch("_inTPC",&_inTPC,"inTPC/I");

  }

  void MakeShowers::SetProperties(){

    /// set volume for TrajectoryInVolume algorithm
    _inTPCAlgo.SetVolume( 0, 
			  2*(::larutil::Geometry::GetME()->DetHalfWidth()),
			  -(::larutil::Geometry::GetME()->DetHalfHeight()),
			  ::larutil::Geometry::GetME()->DetHalfHeight(),
			  0,
			  ::larutil::Geometry::GetME()->DetLength());
  }

  
  bool MakeShowers::analyze(storage_manager* storage) {

    auto *event_part = (event_mcpart*)(storage->get_data(data::kMCParticle,"largeant"));

    if (!event_part) {
      std::cout << "Noooo! " << std::endl;
      return false;
    }

    // make the particle map & Tree-structure
    _MCgetter.Reset(event_part);

    // Now Call Shower Maker from MCgetter
    // Needs to be called from each primary (i.e. TreeTop)
    // because it digs into the tree structure and finds showers
    std::cout << "Number of TreeTops: " << _MCgetter.getTreeTops().size() << std::endl;
    for (size_t i=0; i < _MCgetter.getTreeTops().size(); i++)
      _MCgetter.findMCShowers(_MCgetter.getTreeTops().at(i));

    // now get all showers
    std::vector<std::vector<int> > showers = _MCgetter.getAllShowers();
    std::cout << "Number of showers in event: " << showers.size() << std::endl;

    // loop over showers and fill tree
    for (size_t j=0; j < showers.size(); j++){
      ResetTree();
      std::vector<int> thisshower = showers.at(j);
      if ( _MCgetter.searchParticleMap(thisshower.at(0)) > 0 ){
	mcpart showerTop = event_part->at(_MCgetter.searchParticleMap(thisshower.at(0)));
	_showerDaughters = thisshower.size();
	_inTPC = 0;
	_showerE      = showerTop.Trajectory().at(0).E();
	_showerStartX = showerTop.Trajectory().at(0).X();
	_showerStartY = showerTop.Trajectory().at(0).Y();
	_showerStartZ = showerTop.Trajectory().at(0).Z();
	_showerPDG    = showerTop.PdgCode();
	// now get trajectories (in TPC only)
	for (size_t u=0; u < thisshower.size(); u++){
	  mcpart tmppart = event_part->at(_MCgetter.searchParticleMap(thisshower.at(u)));
	  std::vector<std::vector<double> > thistrack = _MCgetter.getTrajectoryPoints(&tmppart);
	  if (thistrack.size() > 1){
	    _inTPC = 1;
	    ShowerTraj.push_back(thistrack);
	  }
	}//for all particles in the shower
	std::cout << "Energy : " << _showerE << "\tdaughters: " << _showerDaughters 
		  << "\tinTPC: " << _inTPC << "\tMCDaughters: " << showerTop.Daughters().size()
		  << std::endl;
	_showertree->Fill();
      }
    }
      
    
    _evtN += 1;
    
    return true;
  }
  
  bool MakeShowers::finalize() {

    _showertree->Write();
    return true;
  }

  void MakeShowers::ResetTree(){

    _showerDaughters = -1;
    _showerE = -1;
    _showerPDG = -1;
    _showerStartX = 0;
    _showerStartY = 0;
    _showerStartZ = 0;
    ShowerTraj.clear();
    _inTPC = -1;

    return;
  }

}
#endif
