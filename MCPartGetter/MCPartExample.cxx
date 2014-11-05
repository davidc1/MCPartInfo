#ifndef CALCULATEBACKGROUND_CXX
#define CALCULATEBACKGROUND_CXX

#include "MCPartExample.h"
#include <time.h>

namespace larlite {

  bool MCPartExample::initialize() {

    if (_verbose) { _MCgetter.SetVerbose(true); }

    /// Prepare Tree
    prepareTree();

    _evtN = 0;

    return true;
  }


  void MCPartExample::prepareTree(){

    if (_tree) delete _tree;
    _tree = new TTree("ana_tree","");
    _tree->Branch("_isPrimary",&_isPrimary,"_isPrimary/I");
    _tree->Branch("_Energy",&_Energy,"Energy/D");
    _tree->Branch("_inTPC",&_inTPC,"inTPC/I");
    _tree->Branch("_StartX",&_StartX,"StartX/D"); 
    _tree->Branch("_StartY",&_StartY,"StartY/D"); 
    _tree->Branch("_StartZ",&_StartZ,"StartZ/D"); 
    _tree->Branch("_StartT",&_StartT,"StartT/D");
    _tree->Branch("_EndX",&_EndX,"EndX/D"); 
    _tree->Branch("_EndY",&_EndY,"EndY/D"); 
    _tree->Branch("_EndZ",&_EndZ,"EndZ/D"); 
    _tree->Branch("PartTraj",&PartTraj);
    _tree->Branch("ProcHist",&ProcHist);
    _tree->Branch("Process",&Process);
    _tree->Branch("_MotherPDG",&_MotherPDG,"MotherPDG/I");
    _tree->Branch("_MotherE",&_MotherE,"MotherE/D");
    _tree->Branch("_MotherDist",&_MotherDist,"MotherDist/D");
    _tree->Branch("MotherTraj",&MotherTraj);
    _tree->Branch("_AncestorPDG",&_AncestorPDG,"AncestorPDG/I");
    _tree->Branch("_AncestorE",&_AncestorE,"AncestorE/D");
    _tree->Branch("_AncestorDist",&_AncestorDist,"AncestorDist/D");
    _tree->Branch("AncestorTraj",&AncestorTraj);
    _tree->Branch("_minMuonDist",&_minMuonDist,"MinMuonDist/D");


  }

  void MCPartExample::SetProperties(){

    /// set volume for TrajectoryInVolume algorithm
    _inTPCAlgo.SetVolume( 0, 
			  2*(::larutil::Geometry::GetME()->DetHalfWidth()),
			  -(::larutil::Geometry::GetME()->DetHalfHeight()),
			  ::larutil::Geometry::GetME()->DetHalfHeight(),
			  0,
			  ::larutil::Geometry::GetME()->DetLength());
  }

  
  bool MCPartExample::analyze(storage_manager* storage) {

    auto *event_part = (event_mcpart*)(storage->get_data(data::kMCParticle,"largeant"));

    if (!event_part) {
      std::cout << "Noooo! " << std::endl;
      return false;
    }

    // make the particle map & Tree-structure
    _MCgetter.Reset(event_part);
    
    // If the _MCgetter.getAllPDGs() function was called somewhere
    // before the Reset function, then while looping through particles
    // _MCgetter will also find ones that have the right PDG and are
    // above the enrergy cut.
    // The particles that succesfully make the cut are placed in a
    // vector (_TreeNodes in _MCgetter) which holds vectors of TrackId
    // values. Each element in _TreeNodes is a list of different PDG
    // particles, the order defined one the list of PDGs that we are
    // searching for is set. See mac/mcparexample.py to see where
    // this list is set.
    // This list of trackId values can be retrieved with the
    // _MCgetter.getTreeNodelist() function as shown below
    std::vector<TreeNode> result = _MCgetter.getTreeNodelist().at(0);
    std::vector<TreeNode> muon = _MCgetter.getTreeNodelist().at(1);
    
    std::cout << "MCParticles:   " << event_part->size() << std::endl;
    std::cout << "Matches found: " << result.size() << std::endl;
    std::cout << "Muons found:   " << muon.size() << std::endl;

    //prepare list of muon tracks (each track a list of 3D points)
    std::vector< std::vector< std::vector<double> > > muonTracks;
    muonTracks.clear();
    for (size_t h=0; h < muon.size(); h++){
      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsInTPC(&(event_part->at(_MCgetter.searchParticleMap(muon.at(h).getNodeIndex()))),0);
      if (muonTraj.size() > 1)
	muonTracks.push_back(muonTraj);
    }//for all muons
    
    for (size_t j=0; j < result.size(); j++){

      if (_MCgetter.searchParticleMap(result.at(j).getNodeIndex()) >= 0){
	mcpart part = event_part->at(_MCgetter.searchParticleMap( result.at(j).getNodeIndex() ));
	// Now get particle track
	// Trajectory consisting only of start & end points
	std::vector<std::vector<double> > partTrajStartEnd = _MCgetter.getStartEndTrajectory(&part);
	std::vector<std::vector<double> > partTrajAll = _MCgetter.getTrajectoryPoints(&part);
	std::vector<std::vector<double> > partTrajInTPC = _MCgetter.getTrajectoryPointsInTPC(&part,0); //0 cm buffer...up to TPC boundaries
	std::vector<std::vector<double> > partTrajInTPCExtra = _MCgetter.getTrajectoryPointsInTPC(&part,100);//100 cm buffer...track up to 1 meter away from all directions
	std::string process = _MCgetter.getFullProcess(part);
	
	if ( !(result.at(j).isPrimary()) ){
	  if (_MCgetter.searchParticleMap(result.at(j).getParentId()) >= 0)
	    mcpart mother = event_part->at(_MCgetter.searchParticleMap( result.at(j).getParentId() ));
	  if (_MCgetter.searchParticleMap(result.at(j).getAncestorId()) >= 0)
	    mcpart ancestor = event_part->at(_MCgetter.searchParticleMap( result.at(j).getAncestorId() ));
	}//if particle not primary
	_tree->Fill();
	
      }//if searchparticlemap returns ok value
    }//for all particles
    _evtN += 1;
    
    return true;
  }

  bool MCPartExample::finalize() {

    _tree->Write();

    return true;
  }


  void MCPartExample::SetProcess(std::vector<int> PDGs, std::vector<std::string> procs){

    _process.clear();

    if ( PDGs.size() != procs.size() )
      std::cout << "Error: PDG and Process list must have same length. No process added." << std::endl;
    else{
      std::pair<int,std::string> pairTmp;
      for (size_t i=0; i < PDGs.size(); i++){
	pairTmp = std::make_pair(PDGs.at(i),procs.at(i));
	_process.push_back(pairTmp);
      }
    }
      
    return;
  }

}
#endif
