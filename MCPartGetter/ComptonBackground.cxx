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

    _hMuonTotLen = new TH1D("hMuonTotLen","Summed Length of All Muons in one Event; Sum length [meters]", 100, 0, 100);

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
    _tree->Branch("_MotherEndE",&_MotherEndE,"MotherEndE/D");
    _tree->Branch("_MotherDist",&_MotherDist,"MotherDist/D");
    _tree->Branch("MotherTraj",&MotherTraj);
    _tree->Branch("_AncestorPDG",&_AncestorPDG,"AncestorPDG/I");
    _tree->Branch("_AncestorE",&_AncestorE,"AncestorE/D");
    _tree->Branch("_AncestorDist",&_AncestorDist,"AncestorDist/D");
    _tree->Branch("AncestorTraj",&AncestorTraj);
    _tree->Branch("_PoCAtoAncestor",&_PoCAtoAncestor,"PoCAtoAncestor/D");
    _tree->Branch("_PoCAtoAncestorDist",&_PoCAtoAncestorDist,"PoCAtoAncestorDist/D");
    _tree->Branch("_minMuonDist",&_minMuonDist,"minMuonDist/D");
    _tree->Branch("_minMuonPoka",&_minMuonPoka,"minMuonPoka/D");
    _tree->Branch("_PoCADist",&_PoCADist,"PoCADist/D");
    _tree->Branch("_PoCADistAfterStart",&_PoCADistAfterStart,"POCADistAfterStart/I");


    if(_muontree) delete _muontree;
    _muontree = new TTree("muon_tree","");
    _muontree->Branch("_muonE",&_muonE,"muonE/D");
    _muontree->Branch("_muonPDG",&_muonPDG,"muonPDG/I");
    _muontree->Branch("_muonStartX",&_muonStartX,"muonStartX/D");
    _muontree->Branch("_muonStartY",&_muonStartY,"muonStartY/D");
    _muontree->Branch("_muonStartZ",&_muonStartZ,"muonStartZ/D");
    _muontree->Branch("_muonEndX",&_muonEndX,"muonEndX/D");
    _muontree->Branch("_muonEndY",&_muonEndY,"muonEndY/D");
    _muontree->Branch("_muonEndZ",&_muonEndZ,"muonEndZ/D");
    _muontree->Branch("MuonTraj",&MuonTraj);

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
    /*    
    std::cout << "MCParticles:   " << event_part->size() << std::endl;
    std::cout << "Matches found: " << result.size() << std::endl;
    std::cout << "Muons found:   " << muon.size()+antimuon.size() << std::endl;
    */
    //keep track of total lenght of all muon tracks in event
    double totMuonLen = 0;

    // prepare list of muon tracks (each track a list of 3D points)
    // This will be used for the background cuts
    std::vector< std::vector< std::vector<double> > > muonTracks;
    muonTracks.clear();
    for (size_t h=0; h < muon.size(); h++){
      mcpart mu = event_part->at(_MCgetter.searchParticleMap(muon.at(h).getNodeIndex()));
      _muonE = mu.Trajectory().at(0).E();
      _muonPDG = mu.PdgCode();
      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsInTPC(&mu,0);
      MuonTraj = muonTraj;
      if (muonTraj.size() > 1){
	muonTracks.push_back(muonTraj);
	_muonStartX = muonTraj.at(0).at(0);
	_muonStartY = muonTraj.at(0).at(1);
	_muonStartZ = muonTraj.at(0).at(2);
	_muonEndX = muonTraj.back().at(0);
	_muonEndY = muonTraj.back().at(1);
	_muonEndZ = muonTraj.back().at(2);
	totMuonLen += sqrt( (_muonEndX - _muonStartX)*(_muonEndX - _muonStartX) +
			    (_muonEndY - _muonStartY)*(_muonEndY - _muonStartY) +
			    (_muonEndZ - _muonStartZ)*(_muonEndZ - _muonStartZ) );
      }//if trajectory size > 1
      _muontree->Fill();
    }//for all muons
    for (size_t h=0; h < antimuon.size(); h++){
      mcpart mu = event_part->at(_MCgetter.searchParticleMap(antimuon.at(h).getNodeIndex()));
      _muonE = mu.Trajectory().at(0).E();
      _muonPDG = mu.PdgCode();
      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsInTPC(&mu,0);
      MuonTraj = muonTraj;
      if (muonTraj.size() > 1){
	muonTracks.push_back(muonTraj);
	_muonStartX = muonTraj.at(0).at(0);
	_muonStartY = muonTraj.at(0).at(1);
	_muonStartZ = muonTraj.at(0).at(2);
	_muonEndX = muonTraj.back().at(0);
	_muonEndY = muonTraj.back().at(1);
	_muonEndZ = muonTraj.back().at(2);
	totMuonLen += sqrt( (_muonEndX - _muonStartX)*(_muonEndX - _muonStartX) +
			    (_muonEndY - _muonStartY)*(_muonEndY - _muonStartY) +
			    (_muonEndZ - _muonStartZ)*(_muonEndZ - _muonStartZ) );

      }//if trajectory size > 1
      _muontree->Fill();
    }//for all anti-muons
    _hMuonTotLen->Fill(totMuonLen/100.);

    for (size_t j=0; j < result.size(); j++){

      //Reset Branch variables
      ResetTree();

      if (_MCgetter.searchParticleMap(result.at(j).getNodeIndex()) >= 0){
	mcpart part = event_part->at(_MCgetter.searchParticleMap( result.at(j).getNodeIndex() ));

	//used for PoCA
	std::vector<double> c1 = {-1000,-1000,-1000};
	std::vector<double> c2 = {-1000,-1000,-1000};
	std::vector<double> PoCAPointMU = {-1000,-1000,-1000};
	std::vector<double> PoCAPointE = {-1000,-1000,-1000};

	// Now get particle track
	// Trajectory consisting only of start & end points
	PartTraj = _MCgetter.getTrajectoryPointsInTPC(&part,0); //0 cm buffer...up to TPC boundaries
	if (PartTraj.size() > 0){
	  _inTPC = 1;
	  _StartX = PartTraj.at(0).at(0);
	  _StartY = PartTraj.at(0).at(1);
	  _StartZ = PartTraj.at(0).at(2);
	  std::vector<double> partStart = {_StartX,_StartY,_StartZ};
	  _PX = part.Trajectory().at(0).Px();
	  _PY = part.Trajectory().at(0).Py();
	  _PZ = part.Trajectory().at(0).Pz();
	  double partMom = sqrt(_PX*_PX+_PY*_PY+_PZ*_PZ);
	  std::vector<double> partDir = {_PX/partMom,_PY/partMom,_PZ/partMom};
	  std::vector<double> partOrigin = { partStart.at(0)-partDir.at(0)*300,
					     partStart.at(1)-partDir.at(1)*300,
					     partStart.at(2)-partDir.at(2)*300 };
	  std::vector<double> partEnd = { partStart.at(0)+partDir.at(0)*10,
					  partStart.at(1)+partDir.at(1)*10,
					  partStart.at(2)+partDir.at(2)*10 };
	  
	  _StartT = part.Trajectory().at(0).T();
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
		_MotherDist = _pointDist.DistanceToTrack(&partStart,&MotherTraj);
	      _MotherPDG = mother.PdgCode();
	      _MotherE   = mother.Trajectory().at(0).E();
	      //given time of interaction that produced electron, find step of mother right before and get energy
	      double tmin = 0;
	      for (size_t m=0; m < mother.Trajectory().size(); m++){
		if ( mother.Trajectory().at(m).T() < _StartT )
		  tmin = m;
	      }
	      _MotherEndE   = mother.Trajectory().at(tmin).E();
	    }//Mother
	    if (_MCgetter.searchParticleMap(result.at(j).getAncestorId()) >= 0){
	      mcpart ancestor = event_part->at(_MCgetter.searchParticleMap( result.at(j).getAncestorId() ));
	      AncestorTraj = _MCgetter.getTrajectoryPointsInTPC(&ancestor,0);
	      if (AncestorTraj.size() > 0){
		_AncestorDist = _pointDist.DistanceToTrack(&partStart,&AncestorTraj);
		_PoCAtoAncestor = _PoCA.ClosestApproachToTrajectory(&AncestorTraj,&partOrigin,&partStart,c1,c2);
		_PoCAtoAncestorDist = sqrt ( (c2.at(0)-partStart.at(0))*(c2.at(0)-partStart.at(0)) +
					     (c2.at(1)-partStart.at(1))*(c2.at(1)-partStart.at(1)) +
					     (c2.at(2)-partStart.at(2))*(c2.at(2)-partStart.at(2)) );
	      }
	      _AncestorPDG = ancestor.PdgCode();
	      _AncestorE   = ancestor.Trajectory().at(0).E();
	    }//Ancestor
	  }//if particle not primary
	  else { _isPrimary = 1; _MotherPDG = part.Mother(); _AncestorPDG = -1; }

	  
	  // Figure out distance to nearest muon
	  double minDist = 9999999.;
	  double minPoka = 9999999.;
	  c1 = {-1000,-1000,-1000};
	  c2 = {-1000,-1000,-1000};
	  PoCAPointMU = {-1000,-1000,-1000};
	  PoCAPointE = {-1000,-1000,-1000};

	  for (size_t y=0; y < muonTracks.size(); y++){

	    double tmpPoka = _PoCA.ClosestApproachToTrajectory(&muonTracks.at(y),&partOrigin,&partEnd,c1,c2);
	    //calculate distance from PoCA point to e- start point
	    double tmpDist = _pointDist.DistanceToTrack(&partStart,&muonTracks.at(y));
	    if (tmpDist < minDist) { minDist = tmpDist; }
	    if (tmpPoka < minPoka) {
	      minPoka = tmpPoka;
	      PoCAPointE = c2;
	      PoCAPointMU  = c1;
	    }
	  }
	  if (minDist == 9999999.) { minDist = -1; }
	  if (minPoka == 9999999.) { minPoka = -1; }
	  _minMuonDist = minDist;
	  _minMuonPoka = minPoka;
	  fillPoCAParams(PoCAPointE, partStart, partDir);

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
    _muontree->Write();
    _hMuonTotLen->Write();

    return true;
  }


  void ComptonBackground::fillPoCAParams(std::vector<double> ePoCA, std::vector<double> eStart, std::vector<double> eDir){

    _PoCADist = sqrt( (ePoCA.at(0)-eStart.at(0))*(ePoCA.at(0)-eStart.at(0)) +
		      (ePoCA.at(1)-eStart.at(1))*(ePoCA.at(1)-eStart.at(1)) +
		      (ePoCA.at(2)-eStart.at(2))*(ePoCA.at(2)-eStart.at(2)) );
    std::vector<double> vec = { ePoCA.at(0)-eStart.at(0),
				ePoCA.at(1)-eStart.at(1),
				ePoCA.at(2)-eStart.at(2) };
    double vecmag = sqrt( (vec.at(0)*vec.at(0)) +
			  (vec.at(1)*vec.at(1)) +
			  (vec.at(2)*vec.at(2)) );
    double vec_dir = (vec.at(0)*eDir.at(0) + vec.at(1)*eDir.at(1) + vec.at(2)*eDir.at(2))/vecmag;
    if (vec_dir > 0 ) { _PoCADistAfterStart = 1; _PoCADist *= -1; }
    else { _PoCADistAfterStart = 0; }


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
      _MotherEndE= -1;
      _MotherDist = -1;
      MotherTraj.clear();
      _AncestorPDG = -1;
      _AncestorE   = -1;
      _AncestorDist = -1;
      AncestorTraj.clear();
      _PoCAtoAncestor = -1;
      _PoCAtoAncestorDist = -1;
      _minMuonDist = -1;
      _minMuonPoka = -1;
      _PoCADist    = -1;
      _PoCADistAfterStart = -1;
    }      
  }

}
#endif
