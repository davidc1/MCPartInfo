#ifndef COSMICTRACKS_CXX
#define COSMICTRACKS_CXX

#include "CosmicTracks.h"
#include <time.h>

namespace larlite {

  bool CosmicTracks::initialize() {

    if (_verbose) { _MCgetter.SetVerbose(true); }

    /// Prepare Tree
    prepareTree();

    _evtN = 0;

    _hMuonTotLen = new TH1D("hMuonTotLen","Summed Length of All Muons in one Event; Sum length [meters]", 100, 0, 100);

    return true;
  }


  void CosmicTracks::prepareTree(){

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
    _muontree->Branch("MuonEnergies",&MuonEnergies);
    _muontree->Branch("_maxDeflection",&_maxDeflection,"maxDeflection/D");
    _muontree->Branch("_maxDeflectionX",&_maxDeflectionX,"maxDeflectionX/D");
    _muontree->Branch("_maxDeflectionY",&_maxDeflectionY,"maxDeflectionY/D");
    _muontree->Branch("_maxDeflectionZ",&_maxDeflectionZ,"maxDeflectionZ/D");
    _muontree->Branch("_deflectionE",&_deflectionE,"deflectionE/D");
    _muontree->Branch("_deflectionDistToEnd",&_deflectionDistToEnd,"deflectionDistToEnd/D");
    _muontree->Branch("_dirChange",&_dirChange,"dirChange/D");
    _muontree->Branch("_inTPC",&_inTPC,"inTPC/I");

  }

  void CosmicTracks::SetProperties(){}

  
  bool CosmicTracks::analyze(storage_manager* storage) {

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
    std::vector<TreeNode> muon = _MCgetter.getTreeNodelist().at(0);
    std::vector<TreeNode> antimuon = _MCgetter.getTreeNodelist().at(1);
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
      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsInTPC(&mu,2000);
      //      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsBeforeTPC(&mu,2000,_inTPC);
      //      if (_inTPC==0 )
      //	continue;
      _muonE = mu.Trajectory().at(0).E();
      _muonPDG = mu.PdgCode();
      std::vector<double> muonEnergies = _MCgetter.getEnergyPointsInTPC(&mu,2000);
      MuonTraj = muonTraj;
      MuonEnergies = muonEnergies;
      if (muonTraj.size() > 2){
	_inTPC = 1;
	getMaxAngle(&muonTraj, &muonEnergies);
	std::vector<double> beginDir = {muonTraj.at(1).at(0)-muonTraj.at(0).at(0),
					muonTraj.at(1).at(0)-muonTraj.at(0).at(0),
					muonTraj.at(1).at(0)-muonTraj.at(0).at(0) };
	std::vector<double> endDir = {muonTraj.at(muonTraj.size()-1).at(0)-muonTraj.at(muonTraj.size()-2).at(0),
				      muonTraj.at(muonTraj.size()-1).at(0)-muonTraj.at(muonTraj.size()-2).at(0),
				      muonTraj.at(muonTraj.size()-1).at(0)-muonTraj.at(muonTraj.size()-2).at(0) };
	_dirChange = (180./3.14)*acos(getAngle(beginDir,endDir));

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
      else { _inTPC = 0; }
      _muontree->Fill();
    }//for all muons
    for (size_t h=0; h < antimuon.size(); h++){
      mcpart mu = event_part->at(_MCgetter.searchParticleMap(antimuon.at(h).getNodeIndex()));
      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsInTPC(&mu,2000);
      //      std::vector<std::vector<double> > muonTraj = _MCgetter.getTrajectoryPointsBeforeTPC(&mu,2000,_inTPC);
      //      if (_inTPC == 0)
      //	continue;
      _muonE = mu.Trajectory().at(0).E();
      _muonPDG = mu.PdgCode();
      std::vector<double> muonEnergies = _MCgetter.getEnergyPointsInTPC(&mu,2000);
      MuonTraj = muonTraj;
      MuonEnergies = muonEnergies;
      if (muonTraj.size() > 2){
	_inTPC = 1;
	getMaxAngle(&muonTraj, &muonEnergies);
	std::vector<double> beginDir = {muonTraj.at(1).at(0)-muonTraj.at(0).at(0),
					muonTraj.at(1).at(0)-muonTraj.at(0).at(0),
					muonTraj.at(1).at(0)-muonTraj.at(0).at(0) };
	std::vector<double> endDir = {muonTraj.at(muonTraj.size()-1).at(0)-muonTraj.at(muonTraj.size()-2).at(0),
				      muonTraj.at(muonTraj.size()-1).at(0)-muonTraj.at(muonTraj.size()-2).at(0),
				      muonTraj.at(muonTraj.size()-1).at(0)-muonTraj.at(muonTraj.size()-2).at(0) };
	_dirChange = (180./3.14)*acos(getAngle(beginDir,endDir));

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
      else { _inTPC = 0; }
      _muontree->Fill();
    }//for all anti-muons
    _hMuonTotLen->Fill(totMuonLen/100.);

    _evtN += 1;
    
    return true;
  }
  
  bool CosmicTracks::finalize() {

    _muontree->Write();
    _hMuonTotLen->Write();

    return true;
  }

  void CosmicTracks::getMaxAngle(std::vector<std::vector<double> > *track, std::vector<double> *energies){

    double maxangle = 1;
    int maxtick = 0;
    double energy = 0;

    for (size_t i=0; i < (track->size()-2); i++){
      double thisangle = getAngle(track->at(i),track->at(i+1),track->at(i+2));
      if ( thisangle < maxangle ){
	maxangle = thisangle;
	maxtick = i+1;
	energy = energies->at(i+1);
      }
    }

    _maxDeflection = (180./3.14)*acos(maxangle);
    _maxDeflectionX = track->at(maxtick).at(0);
    _maxDeflectionY = track->at(maxtick).at(1);
    _maxDeflectionZ = track->at(maxtick).at(2);
    _deflectionE = energy;
    _deflectionDistToEnd = pow( (track->at(maxtick).at(0)-track->back().at(0))*(track->at(maxtick).at(0)-track->back().at(0)) +
				(track->at(maxtick).at(1)-track->back().at(1))*(track->at(maxtick).at(1)-track->back().at(1)) +

				(track->at(maxtick).at(2)-track->back().at(2))*(track->at(maxtick).at(2)-track->back().at(2)) , 0.5 );

    return;
  }
      

  double CosmicTracks::getAngle(std::vector<double> dir1, std::vector<double> dir2){

    double mag1 = pow( dir1.at(0)*dir1.at(0) + dir1.at(1)*dir1.at(1) + dir1.at(2)*dir1.at(2), 0.5 );
    double mag2 = pow( dir2.at(0)*dir2.at(0) + dir2.at(1)*dir2.at(1) + dir2.at(2)*dir2.at(2), 0.5 );

    double angle = 1;

    if ( (mag1 > 0) && (mag2 > 0) )
      angle = (dir1.at(0)*dir2.at(0) + dir1.at(1)*dir2.at(1) + dir1.at(2)*dir2.at(2)) / (mag1*mag2);

    return angle;
  }
  
  double CosmicTracks::getAngle(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2){
    
    std::vector<double> dir1 = { p1.at(0)-p0.at(0), p1.at(1)-p0.at(1), p1.at(2)-p0.at(2) };
    std::vector<double> dir2 = { p2.at(0)-p1.at(0), p2.at(1)-p1.at(1), p2.at(2)-p1.at(2) };

    double mag1 = pow( dir1.at(0)*dir1.at(0) + dir1.at(1)*dir1.at(1) + dir1.at(2)*dir1.at(2), 0.5 );
    double mag2 = pow( dir2.at(0)*dir2.at(0) + dir2.at(1)*dir2.at(1) + dir2.at(2)*dir2.at(2), 0.5 );

    double angle = 1;

    if ( (mag1 > 0) && (mag2 > 0) )
      angle = (dir1.at(0)*dir2.at(0) + dir1.at(1)*dir2.at(1) + dir1.at(2)*dir2.at(2)) / (mag1*mag2);

    return angle;
  }


  void CosmicTracks::ResetTree(){

    _muonE = -1;
    _muonPDG = -1;
    _muonStartX = 0;
    _muonStartY = 0;
    _muonStartZ = 0;
    _muonEndX = 0;
    _muonEndY = 0;
    _muonEndZ = 0;
    MuonTraj.clear();
    MuonEnergies.clear();
    _maxDeflection = -1;
    _maxDeflectionX = -1;
    _maxDeflectionY = -1;
    _maxDeflectionZ = -1;
    _deflectionE = -1;
    _deflectionDistToEnd = -1;
    _dirChange = -1;
    _inTPC = -1;

    return;
  }

}
#endif
