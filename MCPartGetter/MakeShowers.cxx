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
    _showertree->Branch("_showerTrackID",&_showerTrackID,"showerTrackID/I");
    _showertree->Branch("_showerE",&_showerE,"showerE/D");
    _showertree->Branch("_showerPDG",&_showerPDG,"showerPDG/I");
    _showertree->Branch("_showerStartX",&_showerStartX,"showerStartX/D");
    _showertree->Branch("_showerStartY",&_showerStartY,"showerStartY/D");
    _showertree->Branch("_showerStartZ",&_showerStartZ,"showerStartZ/D");
    _showertree->Branch("ShowerTraj",&ShowerTraj);
    _showertree->Branch("AncestorTraj",&AncestorTraj);
    _showertree->Branch("_inTPC",&_inTPC,"inTPC/I");
    _showertree->Branch("_showerProcess",&_showerProcess);
    _showertree->Branch("_eventN",&_eventN,"eventN/I");
    _showertree->Branch("_numEl",&_numEl,"numEl/i");
    _showertree->Branch("_numCompt",&_numCompt,"numCompt/I");
    _showertree->Branch("_numConv",&_numConv,"numConv/I");
    _showertree->Branch("_numComptE",&_numComptE,"numComptE/I"); //above some energy threshold (200 MeV)
    _showertree->Branch("_numConvE",&_numConvE,"numConvE/I");

    if (_parttree) delete _parttree;
    _parttree = new TTree("part_tree","");
    _parttree->Branch("_partTrackID",&_partTrackID,"partTrackID/I");
    _parttree->Branch("_partE",&_partE,"partE/D");
    _parttree->Branch("_partPDG",&_partPDG,"partPDG/I");
    _parttree->Branch("_partProcess",&_partProcess);
    _parttree->Branch("PartTraj",&PartTraj);
    _parttree->Branch("_shrID",&_shrID,"shrID/I");
    _parttree->Branch("_shrPDG",&_shrPDG,"shrPDG/I");
    _parttree->Branch("_shrE",&_shrE,"shrE/D");
    _parttree->Branch("_shrProc",&_shrProc);
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
    _parttree->Write();

    return true;
  }

  void MakeShowers::findMCShowers(treenode tree, event_mcpart *evt_part, event_mctree *evt_tree, event_mcshower *evt_mcshower){


    // if this particle has PDG == 22 AND
    // if any of the daughters of a particle are PDG == 11 or -11
    // then this particle should make a shower
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
    getAllShowerParticles(tree, showerPartTrackIDs, evt_part, evt_tree, tree);
    // add al daughter tracks
    shr.DaughterTrackID(showerPartTrackIDs);

    //get particle that generates shower
    mcpart thispart = evt_part->at(evt_tree->searchParticleMap(tree.getNodeIndex()));
    
    if ( abs(thispart.PdgCode()) == 11){
      std::vector<std::vector<double> > thistrack = _MCgetter.getTrajectoryPointsInTPC(&thispart,0);
      if (thistrack.size() > 0)
	ShowerTraj.push_back(thistrack);
    }

    //find ancestor track
    int ancID = tree.getAncestorId();
    if (evt_tree->searchParticleMap(ancID) > 0){
      mcpart anc = evt_part->at(evt_tree->searchParticleMap(ancID));
      AncestorTraj = _MCgetter.getTrajectoryPointsInTPC(&anc,0);
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
    _showerProcess = thispart.Process();


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


  void MakeShowers::getAllShowerParticles(treenode tree, std::vector<unsigned int> &trackIDs, event_mcpart *evt_part, event_mctree *evt_tree, treenode toptree){
    
    //loop over all children of node and add charged particles to shower vector
    std::vector<treenode> daughters = tree.getChildren();
    if (daughters.size() == 0) { return; }
    for (size_t i=0; i < daughters.size(); i++){
      
      //find this trackID in vector: if not found go on
      std::vector<unsigned int>::iterator it;
      it = find(trackIDs.begin(), trackIDs.end(), daughters.at(i).getNodeIndex());
      // if not found
      if ( it == trackIDs.end() ){
      
	trackIDs.push_back(daughters.at(i).getNodeIndex());
	//if e-/e+ add track info to ShowerTraj
	if ( abs(evt_part->at(evt_tree->searchParticleMap(daughters.at(i).getNodeIndex())).PdgCode()) == 11 ){
	  mcpart thispart = evt_part->at(evt_tree->searchParticleMap(daughters.at(i).getNodeIndex()));
	  std::vector<std::vector<double> > thistrack = _MCgetter.getTrajectoryPointsInTPC(&thispart,0);
	  if (thistrack.size() > 0)
	    ShowerTraj.push_back(thistrack);
	  //keep track of stats for shower (num of particles)
	  _numEl += 1;
	  if (thispart.Process() == "compt"){
	    _numCompt += 1;
	  if (thispart.Trajectory().at(0).E() > 0.2)
	    _numComptE += 1;
	  }
	  if (thispart.Process() == "conv"){
	    _numConv += 1;
	    if (thispart.Trajectory().at(0).E() > 0.2)
	      _numConvE += 1;
	  }
	  //Fill stats for this particle
	  resetPartTree();
	  _partTrackID = thispart.TrackId();
	  _partE = thispart.Trajectory().at(0).E();
	  _partPDG = thispart.PdgCode();
	  _partProcess = thispart.Process();
	  if (thistrack.size() > 0)
	    PartTraj = thistrack;
	  mcpart top = evt_part->at(evt_tree->searchParticleMap(toptree.getNodeIndex()));
	  _shrID = top.TrackId();
	  _shrPDG = top.PdgCode();
	  _shrE = top.Trajectory().at(0).E();
	  _shrProc = top.Process();
	  _parttree->Fill();
	  
	}//if charged particle
      }//if this particle was not already added
      getAllShowerParticles(daughters.at(i), trackIDs, evt_part, evt_tree, toptree);
    }
    
    return;
  }


  void MakeShowers::resetPartTree(){

    if (_parttree){
      _partTrackID = -1;
      _partE = -1;
      _partPDG = -1;
      _partProcess = "";
      PartTraj.clear();
      _shrID = -1;
      _shrPDG = -1;
      _shrE = -1;
      _shrProc = "";
    }

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
      AncestorTraj.clear();
      _showerPDG = -1;
      _showerProcess = "";
      _eventN = -1;
      _numEl = 0;
      _numCompt = 0;
      _numConv = 0;
      _numComptE = 0;
      _numConvE = 0;

    }

    return;
  }

}
#endif
