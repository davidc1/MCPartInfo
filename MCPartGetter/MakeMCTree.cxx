#ifndef MAKEMCTREE_CXX
#define MAKEMCTREE_CXX

#include "MakeMCTree.h"

namespace larlite {

  bool MakeMCTree::initialize() {


    return true;
  }
  
  bool MakeMCTree::analyze(storage_manager* storage) {

    // get mcpart data. This is the input
    auto *event_part = (event_mcpart*)(storage->get_data(data::kMCParticle,"largeant"));

    // create new event_tree object...will fill it later
    auto *event_tree = (event_mctree*)(storage->get_data(data::kMCTree,""));

    // first need to make map
    setParticleMap(event_tree, event_part);

    // now can go through the map and add trees to event_tree
    setTrees(event_tree, event_part);
  
    event_part->clear();

    std::cout << event_tree->size() << std::endl;

    return true;
  }

  bool MakeMCTree::finalize() {

    return true;
  }

  void MakeMCTree::setParticleMap(event_mctree *event_tree, event_mcpart *event_part){

    // if event_part empty return
    if (event_part->size() == 0){
      std::cout << "Event Part not filled...exit" << std::endl;
      return;
    }

    //clear module's map
    _ParticleMap.clear();

    //clear map before filling it
    event_tree->clearMap();
     //Fill the Map for particles and trackIDs
    for (size_t i=0; i < event_part->size(); i++){
      _ParticleMap[event_part->at(i).TrackId()] = i;
    }//for all particles in event_part
    
    // set this map for the event_tree object
    event_tree->setMap(_ParticleMap);

    return;
  } 


  int MakeMCTree::searchParticleMap(const int trackID){
    
    try{
      return _ParticleMap.at(trackID);
    }
    catch (const std::out_of_range& oor){
      //       std::cerr << "Out of Range Error: Mother TrackID not found. Likely a Primary particle." << std::endl;
    }
    return -1;
  }


  void MakeMCTree::setTrees(event_mctree *event_tree, event_mcpart *event_part){

     // if event_part empty return
     if (event_part->size() == 0){
       std::cout << "Event Part not filled...exit" << std::endl;
       return;
     }

     for (size_t i=0; i < event_part->size(); i++){
       mcpart part = event_part->at(i);
       if (part.Process() == "primary"){
	 mctree node(part.TrackId());
	 node.setPrimary(true);
	 node.setAncestorId(part.TrackId());
	 // add all sub-nodes for this root node
	 //	 AddNodes(part, node, part.TrackId(), event_part);
	 event_tree->push_back(node);
       }// if the process is primary -> new root node
     }// for all particles

    return;
  }


  void MakeMCTree::AddNodes(mcpart part, treenode& parentnode, int ancestor, event_mcpart *event_part){

    std::set<Int_t> daughters = part.Daughters();
    //loop over daugher iterator -> add daughters to node
    for (std::set<int>::iterator it=daughters.begin(); it!=daughters.end(); ++it){
      mcpart tmpPart = event_part->at(this->searchParticleMap(it));
      ::treenode tmpnode(tmpPart.TrackId());
      tmpnode.setParentId(part.TrackId());
      tmpnode.setAncestorId(ancestor);
      tmpnode.setPrimary(false);
      //use recursion to fill this node's children as well
      parentnode.AddChild(tmpnode);
      AddNodes(tmpPart, tmpnode, ancestor, event_part);
    }

    return;
  }
  
  
}
#endif
