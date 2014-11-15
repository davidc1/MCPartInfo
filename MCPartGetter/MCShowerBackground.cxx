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

    //    std::cout << "number of showers: " << evt_mcshower->size() << std::endl;
    //    std::cout << "number of tracks: " << evt_mctracks->size() << std::endl;

    //keep track of total lenght of all muon tracks in event
    double totMuonLen = 0;
    // make a vector of all tracks. Do this only once
    _allTracks.clear();
    //clock_t j;
    //j = clock();
    for (size_t m=0; m < evt_mctracks->size(); m++)
      totMuonLen += addTrack(evt_mctracks->at(m));
    //j = clock() - j;
    //std::cout << "get muon tracks: " << 1000*((float)j)/CLOCKS_PER_SEC << " ms" << std::endl;
    _hMuonTotLen->Fill(totMuonLen/100.);
    std::cout << "Total length for muons in event: " << totMuonLen/100. << std::endl;
    std::cout << "number of tracks: " << _allTracks.size() << std::endl;
    // now loop over all showers

    //clock_t u;
    //u = clock();

    for (size_t s=0; s < evt_mcshower->size(); s++){

      //get current shower
      mcshower shr = evt_mcshower->at(s);

      // Now get particle track
      // Trajectory consisting only of start & end points
      _inTPC = 1;
      if ( shr.DetProfile().X() == 0 )
	_inTPC = 0;
      else{


	_inTPC = 1;
	_xStart = shr.DetProfile().X();
	_yStart = shr.DetProfile().Y();
	_zStart = shr.DetProfile().Z();
	std::vector<double> shrStart = {_xStart, _yStart, _zStart};
	_px = shr.DetProfile().Px();
	_py = shr.DetProfile().Py();
	_pz = shr.DetProfile().Pz();
	double shrMom = sqrt(_px*_px+_py*_py+_pz*_pz);
	std::vector<double> shrDir = {_px/shrMom,_py/shrMom,_pz/shrMom};
	std::vector<double> partOrigin = { shrStart.at(0)-shrDir.at(0)*300,
					   shrStart.at(1)-shrDir.at(1)*300,
					   shrStart.at(2)-shrDir.at(2)*300 };
	std::vector<double> partEnd = { shrStart.at(0)+shrDir.at(0)*10,
					shrStart.at(1)+shrDir.at(1)*10,
					shrStart.at(2)+shrDir.at(2)*10 };

	_tStart = shr.DetProfile().T();
	_E      = shr.Start().E();
	//proc    = shr.Process();
	_pdg    = shr.PdgCode();

	
	// get results from algorithms
	//	t = clock();
	_cutParamCalculator.getNearestMuonParams(&shrStart, &shrDir, &_allTracks, _muDist, _muIP, _distToIP);
	//	t = clock() - t;
	//	std::cout << "Nearest Muon: " << 1000*((float)t)/CLOCKS_PER_SEC << " ms" << std::endl;
	//	clock_t k;
	//	k = clock();
	_cutParamCalculator.getDistanceToWall(shrStart, shrDir, _forwardToWall, _backToWall);
	//	k = clock() - k;
	//	std::cout << "Wall Dist: " << 1000*((float)k)/CLOCKS_PER_SEC << " ms" << std::endl;

      }
      // Now Fill Tree!
      _tree->Fill();
      
    }//for all particles
    //u = clock() - u;
    //std::cout << "Tot to analyze this shower (inTPC): " << 1000*((float)u)/CLOCKS_PER_SEC << " ms" << std::endl;
    _evtN += 1;
    
    return true;
  }
  
  bool MCShowerBackground::finalize() {

    _tree->Write();
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

    if (_tree) delete _tree;
    _tree = new TTree("ana_tree","");
    _tree->Branch("_isPrimary",&_isPrimary,"isPrimary/I");
    _tree->Branch("_E",&_E,"Energy/D");
    _tree->Branch("_inTPC",&_inTPC,"inTPC/I");
    _tree->Branch("_xStart",&_xStart,"xStart/D"); 
    _tree->Branch("_yStart",&_yStart,"yStart/D"); 
    _tree->Branch("_zStart",&_zStart,"zStart/D"); 
    _tree->Branch("_px",&_px,"Showerpx/D"); 
    _tree->Branch("_py",&_py,"Showerpy/D"); 
    _tree->Branch("_pz",&_pz,"Showerpz/D"); 
    _tree->Branch("_tStart",&_tStart,"tStart/D");
    _tree->Branch("_pdg",&_pdg,"PDG/I");
    //_tree->Branch("proc",&proc);
    _tree->Branch("_mPdg",&_mPdg,"MotherPDG/I");
    _tree->Branch("_mE",&_mE,"mE/D");
    _tree->Branch("_aPdg",&_aPdg,"AncestorPDG/I");
    _tree->Branch("_aE",&_aE,"AncestorE/D");
    _tree->Branch("_muDist",&_muDist,"DistanceToClosestMuon/D");
    _tree->Branch("_muIP",&_muIP,"ImpactParamToClosestMuon/D");
    _tree->Branch("_forwardToWall",&_forwardToWall,"FowardDistToWall/D");
    _tree->Branch("_backToWall",&_backToWall,"BackwardsDistToWall/D");
    
  }


    
  
  void MCShowerBackground::resetTree(){
    
    if (_tree){
      
      _isPrimary = 0;
      _E         = -1;
      _inTPC     = -1;
      _xStart    = -1;
      _yStart    = -1;
      _zStart    = -1;
      _tStart    = -1;
      _px        = -1;
      _py        = -1;
      _pz        = -1;
      proc       = "";
      _mPdg      = 0;
      _mE        = -1;
      _aPdg      = -1;
      _aE        = -1;
      _muDist    = -1;
      _muIP      = -1;
      _distToIP  = -1;
      _forwardToWall = -1;
      _backToWall    = -1;
    }      
  }

}
#endif
