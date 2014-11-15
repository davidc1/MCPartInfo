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

    return true;
  }


  
  bool MCShowerBackground::analyze(storage_manager* storage) {
    
    // get MCShowers
    auto evt_mcshower = storage->get_data<event_mcshower>("mcreco");
    // get MCTracks
    auto evt_mctracks = storage->get_data<event_mctrack>("mcreco");

    //keep track of total lenght of all muon tracks in event
    double totMuonLen = 0;
    // make a vector of all tracks. Do this only once
    _allTracks.clear();
    for (size_t m=0; m < evt_mctracks->size(); m++)
      totMuonLen += addTrack(evt_mctracks->at(m));
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
      _inActiveVolume = 1;
      if ( shr.DetProfile().X() == 0 )
	_inActiveVolume = 0;
      else{

	_trackID = shr.TrackID();
	_inActiveVolume = 1;
	_X = shr.DetProfile().X();
	_Y = shr.DetProfile().Y();
	_Z = shr.DetProfile().Z();
	std::vector<double> shrStart = {_X, _Y, _Z};
	_Px = shr.DetProfile().Px();
	_Py = shr.DetProfile().Py();
	_Pz = shr.DetProfile().Pz();
	double shrMom = sqrt(_Px*_Px+_Py*_Py+_Pz*_Pz);
	std::vector<double> shrDir = {_Px/shrMom,_Py/shrMom,_Pz/shrMom};
	std::vector<double> partOrigin = { shrStart.at(0)-shrDir.at(0)*300,
					   shrStart.at(1)-shrDir.at(1)*300,
					   shrStart.at(2)-shrDir.at(2)*300 };
	std::vector<double> partEnd = { shrStart.at(0)+shrDir.at(0)*10,
					shrStart.at(1)+shrDir.at(1)*10,
					shrStart.at(2)+shrDir.at(2)*10 };

	_T      = shr.DetProfile().T();
	_E      = shr.Start().E();
	_process    = shr.Process();
	_PDG    = shr.PdgCode();

	// get mother information
	_parentPDG = shr.MotherPdgCode();
	_parentX = shr.MotherStart().X();
	_parentY = shr.MotherStart().Y();
	_parentZ = shr.MotherStart().Z();
	_parentT = shr.MotherStart().T();
	_parentPx = shr.MotherStart().Px();
	_parentPy = shr.MotherStart().Py();
	_parentPz = shr.MotherStart().Pz();
	_parentE = shr.MotherStart().E();
	std::vector<double> pVtx = { _parentX, _parentY, _parentZ } ; 
	if(_cutParamCalculator.isInVolume(pVtx))
	  _parentInActiveVolume = 1; 
	else
	  _parentInActiveVolume = 0;
	
	

	// get ancestory information
	_ancestorPDG = shr.AncestorPdgCode();
	_ancestorX = shr.AncestorStart().X();
	_ancestorY = shr.AncestorStart().Y();
	_ancestorZ = shr.AncestorStart().Z();
	_ancestorT = shr.AncestorStart().T();
	_ancestorPx = shr.AncestorStart().Px();
	_ancestorPy = shr.AncestorStart().Py();
	_ancestorPz = shr.AncestorStart().Pz();
	_ancestorE = shr.AncestorStart().E();
	std::vector<double> aVtx = { _ancestorX, _ancestorY, _ancestorZ } ; 
	if(_cutParamCalculator.isInVolume(aVtx))
	  _ancestorInActiveVolume = 1;  
	else
	  _ancestorInActiveVolume = 0;


	
	// get results from algorithms
	_cutParamCalculator.getNearestMuonParams(&shrStart, &shrDir, &_allTracks, _minMuDist, _minMuIP, _distToIP);
	_cutParamCalculator.getDistanceToWall(shrStart, shrDir, _distAlongTraj, _distBackAlongTraj);

      }
      // Now Fill Tree!
      _ana_tree->Fill();
      
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
      _allTracks.push_back(thisTrack);
    }// if a muon
    return totLen;
  }
    

  void MCShowerBackground::prepareTree(){


    if(!_ana_tree) {
      _ana_tree = new TTree("ana_tree","");
      
	  _ana_tree->Branch("_run",&_run,"run/I");
	  _ana_tree->Branch("_subrun",&_subrun,"subrun/I");
	  _ana_tree->Branch("_event",&_event,"event/I");

	  _ana_tree->Branch("_process","std::string",&_process);
	  _ana_tree->Branch("_PDG",&_PDG,"PDG/I");
	  _ana_tree->Branch("_trackID",&_trackID,"trackID/I");

	  _ana_tree->Branch("_X",&_X,"X/D");
	  _ana_tree->Branch("_Y",&_Y,"Y/D");
	  _ana_tree->Branch("_Z",&_Z,"Z/D");
	  _ana_tree->Branch("_T",&_T,"T/D");

	  _ana_tree->Branch("_Px",&_Px,"Px/D");
	  _ana_tree->Branch("_Py",&_Py,"Py/D");
	  _ana_tree->Branch("_Pz",&_Pz,"Pz/D");
	  _ana_tree->Branch("_E",&_E,"E/D");

	  _ana_tree->Branch("_inActiveVolume",&_inActiveVolume,"inActiveVolume/I");
	  
	  _ana_tree->Branch("_distAlongTraj",&_distAlongTraj,"distAlongTraj/D") ;
	  _ana_tree->Branch("_distBackAlongTraj",&_distBackAlongTraj,"distBackAlongTraj/D") ;

	  _ana_tree->Branch("_minMuDist",&_minMuDist,"minMuDist/D");
	  _ana_tree->Branch("_minMuIP",&_minMuIP,"minMuIP/D");
	  _ana_tree->Branch("_distToIP",&_distToIP,"distToIP/D");

	  _ana_tree->Branch("_minMuDistExceptAncestor",&_minMuDistExceptAncestor,"minMuDistExceptAncestor/D");
	  _ana_tree->Branch("_minMuIPExceptAncestor",&_minMuIPExceptAncestor,"minMuIPExceptAncestor/D");
	  _ana_tree->Branch("_distToIPExceptAncestor",&_distToIPExceptAncestor,"distToIPExceptAncestor/D");


	  ////ANCESTOR INFO
	  //
	  _ana_tree->Branch("_parentPDG",&_parentPDG,"parentPDG/I");
	  _ana_tree->Branch("_parentX",&_parentX,"parentX/D");
	  _ana_tree->Branch("_parentY",&_parentY,"parentY/D");
	  _ana_tree->Branch("_parentZ",&_parentZ,"parentZ/D");
	  _ana_tree->Branch("_parentT",&_parentT,"parentT/D");

	  _ana_tree->Branch("_parentPx",&_parentPx,"parentPx/D");
	  _ana_tree->Branch("_parentPy",&_parentPy,"parentPy/D");
	  _ana_tree->Branch("_parentPz",&_parentPz,"parentPz/D");
	  _ana_tree->Branch("_parentE",&_parentE,"parentE/D");

	  _ana_tree->Branch("_parentInActiveVolume",&_parentInActiveVolume,"parentInActiveVolume/I");

	  _ana_tree->Branch("_ancDist",&_ancDist,"ancDist/D");
	  _ana_tree->Branch("_ancIP",&_ancIP,"ancIP/D");
	  _ana_tree->Branch("_ancToIP",&_ancToIP,"ancToIP/D");



	  ////PARENT INFO
	  //
	  _ana_tree->Branch("_ancestorPDG",&_ancestorPDG,"ancestorPDG/I");
	  _ana_tree->Branch("_ancestorX",&_ancestorX,"ancestorX/D");
	  _ana_tree->Branch("_ancestorY",&_ancestorY,"ancestorY/D");
	  _ana_tree->Branch("_ancestorZ",&_ancestorZ,"ancestorZ/D");
	  _ana_tree->Branch("_ancestorT",&_ancestorT,"ancestorT/D");

	  _ana_tree->Branch("_ancestorPx",&_ancestorPx,"ancestorPx/D");
	  _ana_tree->Branch("_ancestorPy",&_ancestorPy,"ancestorPy/D");
	  _ana_tree->Branch("_ancestorPz",&_ancestorPz,"ancestorPz/D");
	  _ana_tree->Branch("_ancestorE",&_ancestorE,"ancestorE/D");

	  _ana_tree->Branch("_ancestorInActiveVolume",&_ancestorInActiveVolume,"ancestorInActiveVolume/I");


//	  _ana_tree->Branch("MuonTraj",&MuonTraj) ;
	
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
  }

}
#endif
