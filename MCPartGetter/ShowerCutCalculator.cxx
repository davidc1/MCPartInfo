#ifndef SHOWERCUTCALCULATOR_CXX
#define SHOWERCUTCALCULATOR_CXX

#include "ShowerCutCalculator.h"


void ShowerCutCalculator::SetAlgoProperties(){
  
  /// set volume for TrajectoryInVolume algorithm
  _inTPCAlgo.SetVolume( 0, 
			2*(::larutil::Geometry::GetME()->DetHalfWidth()),
			-(::larutil::Geometry::GetME()->DetHalfHeight()),
			::larutil::Geometry::GetME()->DetHalfHeight(),
			0,
			::larutil::Geometry::GetME()->DetLength());

  _DistToBoxWall.SetXYZMin( 0,
			    -(::larutil::Geometry::GetME()->DetHalfHeight()),
			    0);

  _DistToBoxWall.SetXYZMax( 2*(::larutil::Geometry::GetME()->DetHalfWidth()),
			    ::larutil::Geometry::GetME()->DetHalfHeight(),
			    ::larutil::Geometry::GetME()->DetLength());
  
}


void ShowerCutCalculator::getNearestMuonParams(std::vector<double> *shrStart,
					       std::vector<double> *shrDir,
					       std::vector<std::vector<std::vector<double> > > *muonTracks,
					       double &Dist,
					       double &IP,
					       double &DistToIP){
  
  // use shower's start point + direction and list of muon trajectory to find:
  // 1) PoCA and PoCA distance to shower start point w/ PoCA GeoAlgo class
  // 2) MuonDist for cylinder cut
  
  // initialize values for PoCA cut
  double minIP = 10000;
  std::vector<double> IP_onShr = {-1000, -1000, -1000};
  std::vector<double> c1 = {-1000,-1000,-1000};
  std::vector<double> c2 = {-1000,-1000,-1000};
  std::vector<double> PoCAPointMU = {-1000,-1000,-1000};
  std::vector<double> PoCAPointE = {-1000,-1000,-1000};
  // initialize values for muon proximity
  double minDist = 10000;

  // use shower start & momentum to define a segment
  // which starts 3 meters before start point
  // and ends 10 cm after start point
  // aligned with shr momentum. Use for Impact Parameter
  std::vector<double> shrOrigin = { shrStart->at(0)-300*shrDir->at(0),
				    shrStart->at(1)-300*shrDir->at(1),
				    shrStart->at(2)-300*shrDir->at(2) };

  std::vector<double> shrEnd = { shrStart->at(0)+10*shrDir->at(0),
				 shrStart->at(1)+10*shrDir->at(1),
				 shrStart->at(2)+10*shrDir->at(2) };

  // loop over all muon tracks and calculate value per-muon
  for (size_t u=0; u < muonTracks->size(); u++){

    if (muonTracks->at(u).size() > 1){
      //std::cout << "Points in Track: " << muonTracks->at(u).size() << std::endl;
      // distance to muon track
      //clock_t j;
      //j = clock();
      double tmpDist = _pointDist.DistanceToTrack(shrStart, &(muonTracks->at(u)));
      //std::cout << "tmpDist = " << tmpDist << std::endl;
      //j = clock() - j;
      //std::cout << "Muon Cylidner Time: " << 1000*((float)j)/CLOCKS_PER_SEC << " ms" << std::endl;
      // Impact parameter
      //clock_t k;
      //k = clock();
      double tmpIP = _PoCA.ClosestApproachToTrajectory(&muonTracks->at(u), &shrOrigin, &shrEnd, c1, c2);
      //k = clock() - k;
      //std::cout << "Impact Parameter Time: " << 1000*((float)k)/CLOCKS_PER_SEC << " ms" << std::endl;
      
      if (tmpDist < minDist) { minDist = tmpDist; }
      if (tmpIP < minIP) { 
	minIP = tmpIP; 
	IP_onShr = c2;
      }// if best muon so far

    }// if the muon's track is > 1

  }// for all muon tracks

  DistToIP = pow( (IP_onShr.at(0)-shrStart->at(0))*(IP_onShr.at(0)-shrStart->at(0)) +
		  (IP_onShr.at(1)-shrStart->at(1))*(IP_onShr.at(1)-shrStart->at(1)) +
		  (IP_onShr.at(2)-shrStart->at(2))*(IP_onShr.at(2)-shrStart->at(2)), 0.5 );

  
  Dist = sqrt(minDist);
  IP = sqrt(minIP);
  
  return;
}// function    



void ShowerCutCalculator::getDistanceToWall(std::vector<double> shrStart,
					    std::vector<double> shrDir,
					    double &distToWallForwards,
					    double &distToWallBackwards){

  distToWallForwards = _DistToBoxWall.DistanceToWall(shrStart,shrDir,1);
  distToWallBackwards = _DistToBoxWall.DistanceToWall(shrStart,shrDir,0);

  return;
}
  
#endif
