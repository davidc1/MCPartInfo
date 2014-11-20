#ifndef MCSHOWERBACKGROUND_CXX
#define MCSHOWERBACKGROUND_CXX

#include "MCShowerBackground.h"
#include <time.h>

namespace larlite {

  bool MCShowerBackground::initialize() {

    /// Prepare Tree
    prepareTree();

    /// count events
    _evtN = 0;

    /// setup _cutParamCalculator (set up TPC boundaries)
    _cutParamCalculator.SetAlgoProperties();

    /// prepare total muon length histogram
    _hMuonTotLen = new TH1D("hMuonTotLen","Summed Length of All Muons in one Event; Sum length [meters]", 100, 0, 100);

    // detector boundaries
    _xmin = 0;
    _xmax = ::larutil::Geometry::GetME()->DetHalfWidth()*2;
    _ymin = -::larutil::Geometry::GetME()->DetHalfHeight();
    _ymax = ::larutil::Geometry::GetME()->DetHalfHeight();
    _zmin = 0;
    _zmax = ::larutil::Geometry::GetME()->DetLength();

    return true;
  }


  
  bool MCShowerBackground::analyze(storage_manager* storage) {

    srand (time(NULL));
    
    // get MCShowers
    auto evt_mcshower = storage->get_data<event_mcshower>("davidc1");
    // get MCTracks
    auto evt_mctracks = storage->get_data<event_mctrack>("davidc1");

    //keep track of total lenght of all muon tracks in event
    double totMuonLen = 0;
    // make a vector of all tracks. Do this only once
    _allTracks.clear();
    _allTrackIDs.clear();
    double frac=0;
    for (size_t m=0; m < evt_mctracks->size(); m++){
      if ( (evt_mctracks->at(m).Start().T() > -0.8E6) and (evt_mctracks->at(m).Start().T() < 0.8E6) ){
      frac += 1;
      totMuonLen += addTrack(evt_mctracks->at(m));
      }
    }
    _hMuonTotLen->Fill(totMuonLen/100.);
    // now loop over all showers

    for (size_t s=0; s < evt_mcshower->size(); s++){

      //get current shower
      mcshower shr = evt_mcshower->at(s);


      _run = evt_mctracks->run() ;
      _subrun = evt_mctracks->subrun();
      _event = evt_mctracks->event_id(); 


      // Now get particle track
      // Trajectory consisting only of start & end points
      /*
      if ( (shr.DetProfile().X() > _xmin) and (shr.DetProfile().X() < _xmax) and
	   (shr.DetProfile().Y() > _ymin) and (shr.DetProfile().Y() < _ymax) and
	   (shr.DetProfile().Z() > _zmin) and (shr.DetProfile().Z() < _zmax) )
      */
      if (shr.DetProfile().X() == 0)
	_inActiveVolume = 0;
      else{

	_trackID = shr.TrackID();
	_inActiveVolume = 1;
	/*
	_X = rand() % 256;
	_Y = rand() % 116;
	_Z = rand() % 1036;
	*/
	_X = shr.DetProfile().X();
	_Y = shr.DetProfile().Y();
	_Z = shr.DetProfile().Z();

	std::vector<double> shrStart = {_X, _Y, _Z};
	//	_Px = shr.DetProfile().Px();
	//	_Py = shr.DetProfile().Py();
	//	_Pz = shr.DetProfile().Pz();
	_Px = shr.Start().Px();
	_Py = shr.Start().Py();
	_Pz = shr.Start().Pz();

	double shrMom = sqrt(_Px*_Px+_Py*_Py+_Pz*_Pz);
	std::vector<double> shrDir = {_Px/shrMom,_Py/shrMom,_Pz/shrMom};
	std::vector<double> shrOrigin = { shrStart.at(0)-shrDir.at(0)*300,
					  shrStart.at(1)-shrDir.at(1)*300,
					  shrStart.at(2)-shrDir.at(2)*300 };
	std::vector<double> shrEnd = { shrStart.at(0)+shrDir.at(0)*10,
				       shrStart.at(1)+shrDir.at(1)*10,
				       shrStart.at(2)+shrDir.at(2)*10 };

	_T          = shr.DetProfile().T();
	_E          = shr.Start().E();
	_process    = shr.Process();
	_PDG        = shr.PdgCode();

	// get mother information
	// if shower starts with photon -> mother is the photon
	if (_PDG == 22){
	  _parentPDG = shr.PdgCode();
	  _parentX   = shr.Start().X();
	  _parentY   = shr.Start().Y();
	  _parentZ   = shr.Start().Z();
	  _parentT   = shr.Start().T();
	  _parentPx  = shr.Start().Px();
	  _parentPy  = shr.Start().Py();
	  _parentPz  = shr.Start().Pz();
	  _parentE   = shr.Start().E();
	}	  
	// otherwise -> electron/positron's mother is the "mother"
	else{
	  _parentPDG = shr.MotherPdgCode();
	  _parentX = shr.MotherStart().X();
	  _parentY = shr.MotherStart().Y();
	  _parentZ = shr.MotherStart().Z();
	  _parentT = shr.MotherStart().T();
	  _parentPx = shr.MotherStart().Px();
	  _parentPy = shr.MotherStart().Py();
	  _parentPz = shr.MotherStart().Pz();
	  _parentE = shr.MotherStart().E();
	}

	std::vector<double> pVtx = { _parentX, _parentY, _parentZ } ; 
	if(_cutParamCalculator.isInVolume(pVtx))
	  _parentInActiveVolume = 1; 
	else
	  _parentInActiveVolume = 0;
	
	

	// get ancestor information
	_ancestorPDG = shr.AncestorPdgCode();
	_ancestorX   = shr.AncestorStart().X();
	_ancestorY   = shr.AncestorStart().Y();
	_ancestorZ   = shr.AncestorStart().Z();
	_ancestorT   = shr.AncestorStart().T();
	_ancestorPx  = shr.AncestorStart().Px();
	_ancestorPy  = shr.AncestorStart().Py();
	_ancestorPz  = shr.AncestorStart().Pz();
	_ancestorE   = shr.AncestorStart().E();

	_ancestorInActiveVolume = 0;
	
	//if ancestor is pi+/pi-/mu+/mu-/proton/e+/e- then find distance to track
	if ( (abs(_ancestorPDG) == 11) or (abs(_ancestorPDG) == 13) or
	     (abs(_ancestorPDG == 211)) or (_ancestorPDG == 2212) ){
	  
	  // try and find ancestor in mctracks and find IP/Dist to it
	  for (size_t m=0; m < evt_mctracks->size(); m++)
	    if (evt_mctracks->at(m).TrackID() == shr.AncestorTrackID()){

	      //get trajectory points
	      std::vector<std::vector<std::vector<double> > > tracks;
	      tracks.clear();
	      std::vector<std::vector<double> > thistrack;
	      std::vector<int> ID;
	      ID.clear();
	      AncestorTraj.clear();
	      for (size_t i=0; i < evt_mctracks->at(m).size(); i++){
		thistrack.push_back( {evt_mctracks->at(m).at(i).X(), evt_mctracks->at(m).at(i).Y(), evt_mctracks->at(m).at(i).Z()} );
		AncestorTraj.push_back( {evt_mctracks->at(m).at(i).X(), evt_mctracks->at(m).at(i).Y(), evt_mctracks->at(m).at(i).Z()} );
	      }
	      if (thistrack.size() > 1){
		// then the ancestor is a track that leaves charge in the detector
		_ancestorInActiveVolume = 1;
		tracks.push_back(thistrack);
		ID.push_back(evt_mctracks->at(m).TrackID());
		_cutParamCalculator.getNearestMuonParams(&shrStart, &shrDir, &tracks, &ID, 0, _ancDist, _ancIP, _ancToIP);
	      }//if track > 1 in length
	    }//if ancestor mctrack is found
	}//for correct PDGs

	
	// get results from algorithms
	_cutParamCalculator.getNearestMuonParams(&shrStart, &shrDir, &_allTracks, &_allTrackIDs, 0,  _minMuDist, _minMuIP, _distToIP);
	_cutParamCalculator.getNearestMuonParams(&shrStart, &shrDir, &_allTracks, &_allTrackIDs, shr.AncestorTrackID(), 
						 _minMuDistExceptAncestor, _minMuIPExceptAncestor, _distToIPExceptAncestor);
	_cutParamCalculator.getDistanceToWall(shrStart, shrDir, _distAlongTraj, _distBackAlongTraj);
	

	// Now Fill Tree!
	// Fill only if inActiveVolume
	_ana_tree->Fill();
      }
    
    }//for all particles
    
    _evtN += 1;
    
    return true;
  }

  bool MCShowerBackground::finalize() {
    
    _ana_tree->Write();
    _hMuonTotLen->Write();

    return true;
  }


  double MCShowerBackground::addTrack(mctrack track){

    double totLen = 0;
    if ( (abs(track.PdgCode()) == 13) and (track.size() > 1) ){
      std::vector<std::vector<double> > thisTrack;    

      for (size_t i=0; i < track.size(); i++){
	thisTrack.push_back( {track.at(i).X(), track.at(i).Y(), track.at(i).Z()} );
	if (i > 0)
	  totLen += pow ( (track.at(i-1).X()-track.at(i).X())*(track.at(i-1).X()-track.at(i).X()) +
			  (track.at(i-1).Y()-track.at(i).Y())*(track.at(i-1).Y()-track.at(i).Y()) +
			  (track.at(i-1).Z()-track.at(i).Z())*(track.at(i-1).Z()-track.at(i).Z()), 0.5);
      }// for all track steps
      /*
      thisTrack.push_back( {track.at(0).X(), track.at(0).Y(), track.at(0).Z()} ); 
      thisTrack.push_back( {track.back().X(), track.back().Y(), track.back().Z()} ); 
      */
      _allTrackIDs.push_back(track.TrackID());
      _allTracks.push_back(thisTrack);
    }// if a muon
    return totLen;
  }


  void MCShowerBackground::prepareTree(){


    if(!_ana_tree) {
      _ana_tree = new TTree("ana_tree","");
      
      _ana_tree->Branch("_run",&_run,"run/I");//---------------Event Run Number
      _ana_tree->Branch("_subrun",&_subrun,"subrun/I");//------Event SubRun Number
      _ana_tree->Branch("_event",&_event,"event/I");//------------------Event ID
      
      _ana_tree->Branch("_process","std::string",&_process);//-G4 Process (deltas are 'muIoni')
      _ana_tree->Branch("_PDG",&_PDG,"PDG/I");//---------------PDG code of particle/shower
      _ana_tree->Branch("_trackID",&_trackID,"trackID/I");//---G4 Track ID
      
      _ana_tree->Branch("_X",&_X,"X/D");//---------------------G4 start X point of mcparticle / first energy dep. point of shower
      _ana_tree->Branch("_Y",&_Y,"Y/D");//---------------------G4 start Y point of mcparticle / first energy dep. point of shower
      _ana_tree->Branch("_Z",&_Z,"Z/D");//---------------------G4 start Z point of mcparticle / first energy dep. point of shower
      _ana_tree->Branch("_T",&_T,"T/D");//---------------------G4 start T point of mcparticle / first energy dep. point of shower
      
      _ana_tree->Branch("_Px",&_Px,"Px/D");//------------------G4 start Px point of mcparticle / first energy dep. point of shower
      _ana_tree->Branch("_Py",&_Py,"Py/D");//------------------G4 start Py point of mcparticle / first energy dep. point of shower
      _ana_tree->Branch("_Pz",&_Pz,"Pz/D");//------------------G4 start Pz point of mcparticle / first energy dep. point of shower
      _ana_tree->Branch("_E",&_E,"E/D");//---------------------G4 start Pz point of mcparticle / G4 Energy of shower @ start point
      
      _ana_tree->Branch("_inActiveVolume",&_inActiveVolume,"inActiveVolume/I");//----------{_X,_Y,_Z} in TPC
      
      _ana_tree->Branch("_distAlongTraj",&_distAlongTraj,"distAlongTraj/D") ;//--------------Forward distance from {_X,_Y,Z} to Wall
      _ana_tree->Branch("_distBackAlongTraj",&_distBackAlongTraj,"distBackAlongTraj/D") ;//--Backwards distance from {_X,_Y,Z} to Wall
      
      _ana_tree->Branch("_minMuDist",&_minMuDist,"minMuDist/D");//-----Distance from {_X,_Y,_Z} to nearest muon track
      _ana_tree->Branch("_minMuIP",&_minMuIP,"minMuIP/D");//-----------Impact parameter to nearest muon
      _ana_tree->Branch("_distToIP",&_distToIP,"distToIP/D");//--------dist between IP point on e-/shr track and {_X,_Y,_Z}
      
      _ana_tree->Branch("_minMuDistExceptAncestor",&_minMuDistExceptAncestor,"minMuDistExceptAncestor/D");//--not looking at ancestor track
      _ana_tree->Branch("_minMuIPExceptAncestor",&_minMuIPExceptAncestor,"minMuIPExceptAncestor/D");//--------not looking at ancestor track
      _ana_tree->Branch("_distToIPExceptAncestor",&_distToIPExceptAncestor,"distToIPExceptAncestor/D");//-----not looking at ancestor track
      
      ////PARENT INFO
      _ana_tree->Branch("_parentPDG",&_parentPDG,"parentPDG/I");//---Mother PDG code
      _ana_tree->Branch("_parentX",&_parentX,"parentX/D");//---------G4 start X of mother
      _ana_tree->Branch("_parentY",&_parentY,"parentY/D");//---------G4 start Y of mother
      _ana_tree->Branch("_parentZ",&_parentZ,"parentZ/D");//---------G4 start Z of mother
      _ana_tree->Branch("_parentT",&_parentT,"parentT/D");//---------G4 start T of mother
      
      _ana_tree->Branch("_parentPx",&_parentPx,"parentPx/D");//------G4 start Px of mother
      _ana_tree->Branch("_parentPy",&_parentPy,"parentPy/D");//------G4 start Py of mother
      _ana_tree->Branch("_parentPz",&_parentPz,"parentPz/D");//------G4 start Pz of mother
      _ana_tree->Branch("_parentE",&_parentE,"parentE/D");//---------G4 start E of mother
      
      _ana_tree->Branch("_parentInActiveVolume",&_parentInActiveVolume,"parentInActiveVolume/I");//---{_parentX,_parentY,_parentZ} in TPC
      
      _ana_tree->Branch("_ancDist",&_ancDist,"ancDist/D");//----Distance from {_X,_Y,_Z} to ancestor (if ancestor is e+/e-/mu+/mu-/p/pi+/pi-)
      _ana_tree->Branch("_ancIP",&_ancIP,"ancIP/D");//----------Impact Param to ancestor (if ancestor is e+/e-/mu+/mu-/p/pi+/pi-)
      _ana_tree->Branch("_ancToIP",&_ancToIP,"ancToIP/D");//----dist between IP and start (if ancestor is e+/e-/mu+/mu-/p/pi+/pi-)
      
      ////ANCESTOR INFO
      _ana_tree->Branch("_ancestorPDG",&_ancestorPDG,"ancestorPDG/I");//---Ancestor PDG code
      _ana_tree->Branch("_ancestorX",&_ancestorX,"ancestorX/D");//---------G4 start X of ancestor
      _ana_tree->Branch("_ancestorY",&_ancestorY,"ancestorY/D");//---------G4 start Y of ancestor
      _ana_tree->Branch("_ancestorZ",&_ancestorZ,"ancestorZ/D");//---------G4 start Z of ancestor
      _ana_tree->Branch("_ancestorT",&_ancestorT,"ancestorT/D");//---------G4 start T of ancestor
      
      _ana_tree->Branch("_ancestorPx",&_ancestorPx,"ancestorPx/D");//------G4 start Px of ancestor
      _ana_tree->Branch("_ancestorPy",&_ancestorPy,"ancestorPy/D");//------G4 start Py of ancestor
      _ana_tree->Branch("_ancestorPz",&_ancestorPz,"ancestorPz/D");//------G4 start Pz of ancestor
      _ana_tree->Branch("_ancestorE",&_ancestorE,"ancestorE/D");//---------G4 start E of ancestor
      
      _ana_tree->Branch("_ancestorInActiveVolume",&_ancestorInActiveVolume,"ancestorInActiveVolume/I");//---{_ancX,_ancY,_ancZ} in TPC

      _ana_tree->Branch("ShowerTraj",&ShowerTraj);
      _ana_tree->Branch("MotherTraj",&MotherTraj);
      _ana_tree->Branch("AncestorTraj",&AncestorTraj);
      
    }
  }
  

    
  
  void MCShowerBackground::resetTree(){

   _run     = -1;
   _subrun  = -1;
   _event 	= -1;

   _process = "NONE";
   _PDG  	= -1;
   _trackID = -1;

   _X  		= -9999999;
   _Y 		= -9999999;
   _Z 		= -9999999;
   _T 		= -9999999;

   _Px 		= -9999999;
   _Py 		= -9999999;
   _Pz 		= -9999999;
   _E 		= -9999999;

   _distAlongTraj     = -9999999;
   _distBackAlongTraj = -9999999;

   _minMuDist = -1;
   _minMuIP = -1;
   _distToIP = -1;

   _minMuDistExceptAncestor = -1;
   _minMuIPExceptAncestor = -1;
   _distToIPExceptAncestor = -1;

   _ancDist = -1;
   _ancIP = -1;
   _ancToIP = -1;

   //parent
   _parentPDG = -1;
   _parentX   = -9999999;
   _parentY   = -9999999;
   _parentZ   = -9999999;
   _parentT   = -9999999;

   _parentPx  = -9999999;
   _parentPy  = -9999999;
   _parentPz  = -9999999;
   _parentE   = -9999999;


   //ancestor
   _ancestorPDG = -1;
   _ancestorX   = -9999999;
   _ancestorY   = -9999999;
   _ancestorZ   = -9999999;
   _ancestorT   = -9999999;

   _ancestorPx  = -9999999;
   _ancestorPy  = -9999999;
   _ancestorPz  = -9999999;
   _ancestorE   = -9999999;


   _inActiveVolume = -99 ;
   _parentInActiveVolume = -99 ;
   _ancestorInActiveVolume = -99 ;

   ShowerTraj.clear();
   MotherTraj.clear();
   AncestorTraj.clear();

  }

}
#endif
