#ifndef MAKESHOWERS_CXX
#define MAKESHOWERS_CXX

#include "MakeShowers.h"
#include <time.h>

namespace larlite {

  bool MakeShowers::initialize() {

    _evtN = 0;

    PrepareTree();

    return true;
  }


  void MakeShowers::PrepareTree(){

    if(_showertree) delete _showertree;
    _showertree = new TTree("shower_tree","");
    _showertree->Branch("_showerE",&_showerE,"showerE/D");
    _showertree->Branch("_showerPDG",&_showerPDG,"showerPDG/I");
    _showertree->Branch("_showerStartX",&_showerStartX,"showerStartX/D");
    _showertree->Branch("_showerStartY",&_showerStartY,"showerStartY/D");
    _showertree->Branch("_showerStartZ",&_showerStartZ,"showerStartZ/D");
    _showertree->Branch("ShowerTraj",&ShowerTraj);
    _showertree->Branch("_inTPC",&_inTPC,"inTPC/I");
    _showertree->Branch("Process",&Process);
    _showertree->Branch("_eventN",&_eventN,"eventN/I");
    
  }
  


  bool MakeShowers::analyze(storage_manager* storage) {

    // get mcpart data. This is the input
    auto evt_part = storage->get_data<event_mcpart>("largeant");
    // get mctree data
    auto evt_tree = storage->get_data<event_mctree>("davidc");
    // create new mcshowers
    auto evt_mcshower = storage->get_data<event_mcshower>("davidc");


    evt_mcshower->set_subrun(evt_part->subrun());
    evt_mcshower->set_run(evt_part->run());
    evt_mcshower->set_event_id(evt_part->event_id());


    if (!evt_part) {
      std::cout << "Noooo! " << std::endl;
      return false;
    }

    // Now Call Shower Maker from MCgetter
    // Needs to be called from each primary (i.e. TreeTop)
    // because it digs into the tree structure and finds showers
    for (size_t i=0; i < evt_tree->size(); i++)
      findMCShowers(evt_tree->at(i), evt_part, evt_tree, evt_mcshower);

    _evtN += 1;

    evt_part->clear();

    return true;
  }
  
  bool MakeShowers::finalize() {

    _showertree->Write();

    return true;
  }

  void MakeShowers::findMCShowers(treenode tree, event_mcpart *evt_part, event_mctree *evt_tree, event_mcshower *evt_mcshower){


    // if this particle has PDG == 22 AND
    // if any of the daughters of a particle are PDG == 11 or -11
    // then this particle should make a shower
    if ( tree.getNodeIndex() == 129035 ){
      std::cout << "found it! " << std::endl;
      std::cout << "TrackID: " << evt_part->at( evt_tree->searchParticleMap(tree.getNodeIndex()) ).TrackId() << std::endl;
      std::cout << "PDG: " << evt_part->at( evt_tree->searchParticleMap(tree.getNodeIndex()) ).PdgCode() << std::endl;
      std::cout << "Energy: " << evt_part->at( evt_tree->searchParticleMap(tree.getNodeIndex()) ).Trajectory().at(0).E() << std::endl;
    }
    if ( (abs(evt_part->at( evt_tree->searchParticleMap(tree.getNodeIndex()) ).PdgCode()) == 11) or
	 (abs(evt_part->at( evt_tree->searchParticleMap(tree.getNodeIndex()) ).PdgCode()) == 22) )
      {
	
	if ((evt_part->at( evt_tree->searchParticleMap(tree.getNodeIndex()) ).Trajectory().at(0).E() > _Ecut) ){
	  makeMCShower(tree, evt_part, evt_tree, evt_mcshower);
	}
      } // if particle's PDG == 22 and above energy cut
    
    // otherwise...keep on looking
    else{
      std::vector<::treenode> children = tree.getChildren();
      for (size_t i=0; i < children.size(); i++)
	findMCShowers(children.at(i), evt_part, evt_tree, evt_mcshower);
    }
    
    return;
  }

  void MakeShowers::makeMCShower(treenode tree, event_mcpart *evt_part, event_mctree *evt_tree, event_mcshower *evt_mcshower){

    //first clear tree
    ResetTree();

    // object to be filled
    mcshower shr;

    //get all daughter charged particles to make up the shower
    std::vector<unsigned int> showerPartTrackIDs;
    getAllShowerParticles(tree, showerPartTrackIDs, evt_part, evt_tree);
    // add al daughter tracks
    shr.DaughterTrackID(showerPartTrackIDs);

    //get particle that generates shower
    mcpart thispart = evt_part->at(evt_tree->searchParticleMap(tree.getNodeIndex()));
    
    if ( abs(thispart.PdgCode()) == 11){
      std::vector<std::vector<double> > thistrack = _MCgetter.getTrajectoryPointsInTPC(&thispart,0);
      if (thistrack.size() > 0)
	ShowerTraj.push_back(thistrack);
    }


    shr.PdgCode(thispart.PdgCode());
    shr.Process(thispart.Process());
    shr.TrackID(thispart.TrackId());
    shr.Start(thispart.Trajectory().at(0));
    shr.End(thispart.Trajectory().back());

    _showerE = thispart.Trajectory().at(0).E();
    _showerStartX = thispart.Trajectory().at(0).X();
    _showerStartY = thispart.Trajectory().at(0).Y();
    _showerStartZ = thispart.Trajectory().at(0).Z();
    _eventN = _evtN;
    if (ShowerTraj.size() == 0)
      _inTPC = 0;
    else{
      _inTPC = 1;
      int totpoints = 0;
      for (size_t u=0; u < ShowerTraj.size(); u++)
	totpoints += ShowerTraj.at(u).size();

      _showertree->Fill();
    }
    
    _showerPDG = thispart.PdgCode();
    Process = thispart.Process();


    // mother of showering particle
    if (evt_tree->searchParticleMap(tree.getParentId()) > 0){
      mcpart mother = evt_part->at(evt_tree->searchParticleMap(tree.getParentId()));
      shr.MotherPdgCode(mother.PdgCode());
      shr.MotherProcess(mother.Process());
      shr.MotherTrackID(mother.TrackId());
      shr.MotherStart(mother.Trajectory().at(0));
      shr.MotherEnd(mother.Trajectory().back());
    }

    // mother of showering particle
    if (evt_tree->searchParticleMap(tree.getAncestorId())){
      mcpart ancestor = evt_part->at(evt_tree->searchParticleMap(tree.getAncestorId()));
      shr.AncestorPdgCode(ancestor.PdgCode());
      shr.AncestorProcess(ancestor.Process());
      shr.AncestorTrackID(ancestor.TrackId());
      shr.AncestorStart(ancestor.Trajectory().at(0));
      shr.AncestorEnd(ancestor.Trajectory().back());
    }
    evt_mcshower->push_back(shr);
    
    return;
  }


  void MakeShowers::getAllShowerParticles(treenode tree, std::vector<unsigned int> &trackIDs, event_mcpart *evt_part, event_mctree *evt_tree){
    
    //loop over all children of node and add charged particles to shower vector
    std::vector<treenode> daughters = tree.getChildren();
    if (daughters.size() == 0) { return; }
    for (size_t i=0; i < daughters.size(); i++){
      trackIDs.push_back(daughters.at(i).getNodeIndex());
      //if e-/e+ add track info to ShowerTraj
      if ( abs(evt_part->at(evt_tree->searchParticleMap(tree.getNodeIndex())).PdgCode()) == 11 ){
	mcpart thispart = evt_part->at(evt_tree->searchParticleMap(tree.getNodeIndex()));
	std::vector<std::vector<double> > thistrack = _MCgetter.getTrajectoryPointsInTPC(&thispart,0);
	if (thistrack.size() > 0)
	  ShowerTraj.push_back(thistrack);
      }//if charged particle
      getAllShowerParticles(daughters.at(i), trackIDs, evt_part, evt_tree);
    }
    
    return;
  }

  // reset tere
  void MakeShowers::ResetTree(){

    if (_showertree){

      _inTPC = -1;
      _showerStartX = -1;
      _showerStartY = -1;
      _showerStartZ = -1;
      _showerE = -1;
      ShowerTraj.clear();
      _showerPDG = -1;
      Process = "";
      _eventN = -1;

    }

    return;
  }

}
#endif
