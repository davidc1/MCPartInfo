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

  _inVol.SetVolume(0,256.35,-116.5,116.5,0,1036.8) ;
  
}

bool ShowerCutCalculator::isInVolume(std::vector<double> point){
  
  return _inVol.PointInVolume(point);
  
}


void ShowerCutCalculator::getNearestMuonParams(std::vector<double> *shrStart,
					       std::vector<double> *shrDir,
					       std::vector<std::vector<std::vector<double> > > *muonTracks,
					       std::vector<int> *muonIDs,
					       int ancestorID,
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

    if ( (muonTracks->at(u).size() > 1) and (muonIDs->at(u) != ancestorID) ){
      double tmpDist = _pointDist.DistanceToTrack(*shrStart, muonTracks->at(u));
      double tmpIP = _PoCA.ClosestApproachToTrajectory(muonTracks->at(u), shrOrigin, shrEnd, c1, c2);
      
      if (tmpDist < minDist) { minDist = tmpDist; }
      if (tmpIP < minIP) { 
	minIP = tmpIP; 
	IP_onShr = c2;
      }// if best muon so far

    }// if the muon's track is > 1
    
  }// for all muon tracks

  if (minDist != 10000){
    DistToIP = pow( (IP_onShr.at(0)-shrStart->at(0))*(IP_onShr.at(0)-shrStart->at(0)) +
		    (IP_onShr.at(1)-shrStart->at(1))*(IP_onShr.at(1)-shrStart->at(1)) +
		    (IP_onShr.at(2)-shrStart->at(2))*(IP_onShr.at(2)-shrStart->at(2)), 0.5 );
    
    
    Dist = sqrt(minDist);
    IP = sqrt(minIP);
    
    //need to figure out if IP point is "before" or "after" start point w.r.t. momentum direction
    if (DistToIP > 0.001){
      std::vector<double> vec = { IP_onShr.at(0)-shrStart->at(0),
				  IP_onShr.at(1)-shrStart->at(1),
				  IP_onShr.at(2)-shrStart->at(2) };
      double vecmag = sqrt( (vec.at(0)*vec.at(0)) + (vec.at(1)*vec.at(1)) + (vec.at(2)*vec.at(2)) );
      double vec_dir = (vec.at(0)*shrDir->at(0) + vec.at(1)*shrDir->at(1) + vec.at(2)*shrDir->at(2))/vecmag;
      if (vec_dir == 1 ) { DistToIP *= -1; }
    }
  }

  return;
}// function    


void ShowerCutCalculator::getAncestorMuonParams(std::vector<double> *shrStart,
						std::vector<double> *shrDir,
						std::vector<std::vector<double> >  *muonTrack,
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
  
  if ( muonTrack->size() > 1 ){

    double tmpDist = _pointDist.DistanceToTrack(*shrStart, *muonTrack);
    double tmpIP = _PoCA.ClosestApproachToTrajectory(*muonTrack, shrOrigin, shrEnd, c1, c2);
    
    if (tmpDist < minDist) { minDist = tmpDist; }
    if (tmpIP < minIP) { 
      minIP = tmpIP; 
      IP_onShr = c2;
    }// if best muon so far

    DistToIP = pow( (IP_onShr.at(0)-shrStart->at(0))*(IP_onShr.at(0)-shrStart->at(0)) +
		    (IP_onShr.at(1)-shrStart->at(1))*(IP_onShr.at(1)-shrStart->at(1)) +
		    (IP_onShr.at(2)-shrStart->at(2))*(IP_onShr.at(2)-shrStart->at(2)), 0.5 );
    
    
    Dist = sqrt(minDist);
    IP = sqrt(minIP);
    
    //need to figure out if IP point is "before" or "after" start point w.r.t. momentum direction
    if (DistToIP > 0.001){
      std::vector<double> vec = { IP_onShr.at(0)-shrStart->at(0),
				  IP_onShr.at(1)-shrStart->at(1),
				  IP_onShr.at(2)-shrStart->at(2) };
      double vecmag = sqrt( (vec.at(0)*vec.at(0)) + (vec.at(1)*vec.at(1)) + (vec.at(2)*vec.at(2)) );
      double vec_dir = (vec.at(0)*shrDir->at(0) + vec.at(1)*shrDir->at(1) + vec.at(2)*shrDir->at(2))/vecmag;
      if (vec_dir == 1 ) { DistToIP *= -1; }
    }
  }// if the muon's track is > 1
  
  return;
}// function    



void ShowerCutCalculator::getDistanceToWall(std::vector<double> shrStart,
					    std::vector<double> shrDir,
					    double &distToWallForwards,
					    double &distToWallBackwards){
  //  std::cout << "Shower Position: [" << shrStart.at(0) << ", "
  //	    << shrStart.at(1) << ", " << shrStart.at(2) << "]" << std::endl;
  //  std::cout << "Shower Direction: [" << shrDir.at(0) << ", "
  //	    << shrDir.at(1) << ", " << shrDir.at(2) << "]" << std::endl;
  distToWallForwards = _DistToBoxWall.DistanceToWall(shrStart,shrDir,1);
  distToWallBackwards = _DistToBoxWall.DistanceToWall(shrStart,shrDir,0);
  //  std::cout << "Dist Back to Wall: " << distToWallBackwards << std::endl;

  return;
}
  
#endif
