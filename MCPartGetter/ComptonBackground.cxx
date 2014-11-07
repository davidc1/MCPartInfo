#ifndef COMPTONBACKGROUND_CXX
#define COMPTONBACKGROUND_CXX

#include "ComptonBackground.h"
#include <time.h>

namespace larlite {

  bool ComptonBackground::initialize() {

    if (_verbose) { _MCgetter.SetVerbose(true); }

    /// Prepare Tree
    prepareTree();

    _evtN = 0;

    return true;
  }


  void ComptonBackground::prepareTree(){

    if (_tree) delete _tree;
    _tree = new TTree("ana_tree","");
    _tree->Branch("_isPrimary",&_isPrimary,"isPrimary/I");
    _tree->Branch("_Energy",&_Energy,"Energy/D");
    _tree->Branch("_inTPC",&_inTPC,"inTPC/I");
    _tree->Branch("_StartX",&_StartX,"StartX/D"); 
    _tree->Branch("_StartY",&_StartY,"StartY/D"); 
    _tree->Branch("_StartZ",&_StartZ,"StartZ/D"); 
    _tree->Branch("_PX",&_PX,"PX/D"); 
    _tree->Branch("_PY",&_PY,"PY/D"); 
    _tree->Branch("_PZ",&_PZ,"PZ/D"); 
    _tree->Branch("_StartT",&_StartT,"StartT/D");
    _tree->Branch("_EndX",&_EndX,"EndX/D"); 
    _tree->Branch("_EndY",&_EndY,"EndY/D"); 
    _tree->Branch("_EndZ",&_EndZ,"EndZ/D"); 
    _tree->Branch("_PDG",&_PDG,"PDG/I");
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
    _tree->Branch("_minMuonDist",&_minMuonDist,"minMuonDist/D");
    _tree->Branch("_minMuonPoka",&_minMuonPoka,"minMuonPoka/D");
    _tree->Branch("_PoCADist",&_PoCADist,"PoCADist/D");


  }

  void ComptonBackground::SetProperties(){

    /// set volume for TrajectoryInVolume algorithm
    _inTPCAlgo.SetVolume( 0, 
			  2*(::larutil::Geometry::GetME()->DetHalfWidth()),
			  -(::larutil::Geometry::GetME()->DetHalfHeight()),
			  ::larutil::Geometry::GetME()->DetHalfHeight(),
			  0,
			  ::larutil::Geometry::GetME()->DetLength());
  }

  
  bool ComptonBackground::analyze(storage_manager* storage) {

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
    std::vector<TreeNode> antimuon = _MCgetter.getTreeNodelist().at(2);
    
    std::cout << "MCParticles:   " << event_part->size() << std::endl;
    std::cout << "Matches found: " << result.size() << std::endl;
    std::cout << "Muons found:   " << muon.size()+antimuon.size() << std::endl;

    // prepare list of muon tracks (each track a list of 3D points)
    // This will be used for the background cuts
    std::vector< std::vector< std::vector<double> > > muonTracks;
    muonTracks.clear();
    for (size_t h=0; h < muon.size(); h++){
      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsInTPC(&(event_part->at(_MCgetter.searchParticleMap(muon.at(h).getNodeIndex()))),0);
      if (muonTraj.size() > 1)
	muonTracks.push_back(muonTraj);
    }//for all muons
    for (size_t h=0; h < antimuon.size(); h++){
      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsInTPC(&(event_part->at(_MCgetter.searchParticleMap(antimuon.at(h).getNodeIndex()))),0);
      if (muonTraj.size() > 1)
	muonTracks.push_back(muonTraj);
    }//for all anti-muons

    for (size_t j=0; j < result.size(); j++){

      //Reset Branch variables
      ResetTree();

      if (_MCgetter.searchParticleMap(result.at(j).getNodeIndex()) >= 0){
	mcpart part = event_part->at(_MCgetter.searchParticleMap( result.at(j).getNodeIndex() ));
	// Now get particle track
	// Trajectory consisting only of start & end points
	PartTraj = _MCgetter.getTrajectoryPointsInTPC(&part,0); //0 cm buffer...up to TPC boundaries
	if (PartTraj.size() > 0){
	  _inTPC = 1;
	  _StartX = PartTraj.at(0).at(0);
	  _StartY = PartTraj.at(0).at(1);
	  _StartZ = PartTraj.at(0).at(2);
	  std::vector<double> partStart = {_StartX,_StartY,_StartZ};
	  _StartT = part.Trajectory().at(0).T();
	  _PX = part.Trajectory().at(0).Px();
	  _PY = part.Trajectory().at(0).Py();
	  _PZ = part.Trajectory().at(0).Pz();
	  _Energy = part.Trajectory().at(0).E();
	  Process = part.Process();
	  ProcHist = _MCgetter.getFullProcess(part);
	  _PDG = part.PdgCode();
	  
	  if ( !(result.at(j).isPrimary()) ){
	    _isPrimary = 0;
	    if (_MCgetter.searchParticleMap(result.at(j).getParentId()) >= 0){
	      mcpart mother = event_part->at(_MCgetter.searchParticleMap( result.at(j).getParentId() ));
	      MotherTraj = _MCgetter.getTrajectoryPointsInTPC(&mother,0);
	      if (MotherTraj.size() > 0)
		_MotherDist = _pointDist.DistanceToTrack(partStart,MotherTraj);
	      _MotherPDG = mother.PdgCode();
	      _MotherE   = mother.Trajectory().at(0).E();
	    }//Mother
	    if (_MCgetter.searchParticleMap(result.at(j).getAncestorId()) >= 0){
	      mcpart ancestor = event_part->at(_MCgetter.searchParticleMap( result.at(j).getAncestorId() ));
	      AncestorTraj = _MCgetter.getTrajectoryPointsInTPC(&ancestor,0);
	      if (AncestorTraj.size() > 0)
	      _AncestorDist = _pointDist.DistanceToTrack(partStart,AncestorTraj);
	      _AncestorPDG = ancestor.PdgCode();
	      _AncestorE   = ancestor.Trajectory().at(0).E();
	    }//Ancestor
	  }//if particle not primary
	  else { _isPrimary = 1; _MotherPDG = part.Mother(); _AncestorPDG = -1; }

	  
	  // Figure out distance to nearest muon
	  double minDist = 9999999.;
	  double minPoka = 9999999.;
	  std::vector<double> c1, c2;
	  double t1,t2;
	  double partMom = sqrt(_PX*_PX+_PY*_PY+_PZ*_PZ);
	  std::vector<double> partDir = {_PX/partMom,_PY/partMom,_PZ/partMom};
	  std::vector<double> partOrigin = { partStart.at(0)-partDir.at(0)*300,
					     partStart.at(1)-partDir.at(1)*300,
					     partStart.at(2)-partDir.at(2)*300 };

	  for (size_t y=0; y < muonTracks.size(); y++){

	    double tmpPoka = _PoCA.ClosestApproachToTrajectory(muonTracks.at(y),partOrigin,partStart,c1,c2,t1,t2);
	    //calculate distance from PoCA point to e- start point
	    if (c2.size()==3)
	      _PoCADist = sqrt( (c2.at(0)-partStart.at(0))*(c2.at(0)-partStart.at(0)) +
				(c2.at(1)-partStart.at(1))*(c2.at(1)-partStart.at(1)) +
				(c2.at(2)-partStart.at(2))*(c2.at(2)-partStart.at(2)) );
	    double tmpDist = _pointDist.DistanceToTrack(partStart,muonTracks.at(y));
	    if (tmpDist < minDist) { minDist = tmpDist; }
	    if (tmpPoka < minPoka) { minPoka = tmpPoka; }
	  }
	  if (minDist == 9999999.) { minDist = -1; }
	  if (minPoka == 9999999.) { minPoka = -1; }
	  _minMuonDist = minDist;
	  _minMuonPoka = minPoka;
	  
	}//if in TPC
	else { _inTPC = 0; }
	
	
	// Now Fill Tree!
	_tree->Fill();
	
      }//if searchparticlemap returns ok value
    }//for all particles
    _evtN += 1;
    
    return true;
  }
  
  bool ComptonBackground::finalize() {

    _tree->Write();

    return true;
  }


  void ComptonBackground::SetProcess(std::vector<int> PDGs, std::vector<std::string> procs){

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

  void ComptonBackground::ResetTree(){
    
    if (_tree){
      
      _isPrimary = -1;
      _Energy    = -1;
      _inTPC     = -1;
      _StartX    = -1;
      _StartY    = -1;
      _StartZ    = -1;
      _StartT    = -1;
      _EndX      = -1;
      _EndY      = -1;
      _EndZ      = -1;
      _PX        = -1;
      _PY        = -1;
      _PZ        = -1;
      PartTraj.clear();
      ProcHist   = "";
      Process    = "";
      _MotherPDG = -1;
      _MotherE   = -1;
      _MotherDist = -1;
      MotherTraj.clear();
      _AncestorPDG = -1;
      _AncestorE   = -1;
      _AncestorDist = -1;
      AncestorTraj.clear();
      _minMuonDist = -1;
      _minMuonPoka = -1;
      _PoCADist    = -1;
    }      
  }

}
#endif
