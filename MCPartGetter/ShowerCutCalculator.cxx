#ifndef SHOWERCUTCALCULATOR_CXX
#define SHOWERCUTCALCULATOR_CXX

#include "ShowerCutCalculator.h"


void ShowerCutCalculator::SetAlgoProperties(){

  /// Set TPC
  _TpcBox = geoalgo::AABox(0, -(::larutil::Geometry::GetME()->DetHalfHeight()), 0,
			   2*(::larutil::Geometry::GetME()->DetHalfWidth()),
			   ::larutil::Geometry::GetME()->DetHalfHeight(),
			   ::larutil::Geometry::GetME()->DetLength());
  
}


bool ShowerCutCalculator::isInVolume(const std::vector<double>& point){

  geoalgo::Point_t p(point);

  return _TpcBox.Contain(p);
  
}


void ShowerCutCalculator::getNearestMuonParams(const geoalgo::Point_t& shrStart,
					       const geoalgo::Vector_t& shrDir,
					       const std::vector<geoalgo::Trajectory_t >& muonTracks,
					       const std::vector<int>& muonIDs,
					       const int ancestorID,
					       double &Dist,
					       double &IP,
					       double &DistToIP) const
{
  
  // use shower's start point + direction and list of muon trajectory to find:
  // 1) PoCA and PoCA distance to shower start point w/ PoCA GeoAlgo class
  // 2) MuonDist for cylinder cut
  
  // initialize values for PoCA cut
  double minIP = std::numeric_limits<double>::max();
  geoalgo::Point_t IP_onShr(3);
  geoalgo::Vector c1(3);
  geoalgo::Vector c2(3);
  //  std::vector<double> c1 = {-1000,-1000,-1000};
  //  std::vector<double> c2 = {-1000,-1000,-1000};
  // initialize values for muon proximity
  double minDist = std::numeric_limits<double>::max();

  // use shower start & momentum to define a segment
  // which starts 3 meters before start point
  // and ends 10 cm after start point
  // aligned with shr momentum. Use for Impact Parameter
  geoalgo::LineSegment shrSeg(shrStart-shrDir*300,shrStart+shrDir*10);

  // loop over all muon tracks and calculate value per-muon
  for (size_t u=0; u < muonTracks.size(); u++){

    if ( (muonTracks.at(u).size() > 1) and (muonIDs.at(u) != ancestorID) ){

      double tmpDist = _geoAlgo.SqDist(shrStart, muonTracks.at(u));
      double tmpIP = _geoAlgo.SqDist(muonTracks.at(u), shrSeg, c1, c2);
      if (tmpDist < minDist) { minDist = tmpDist; }
      if (tmpIP < minIP) { 
	minIP = tmpIP; 
	IP_onShr = c2;
      }// if best muon so far

    }// if the muon's track is > 1
    
  }// for all muon tracks

  if (minDist != std::numeric_limits<double>::max()){
    DistToIP = IP_onShr.Dist(shrStart);
    Dist = sqrt(minDist);
    IP = sqrt(minIP);
    
    //need to figure out if IP point is "before" or "after" start point w.r.t. momentum direction
    if (DistToIP > 0.001){
      double dot = ( (IP_onShr-shrStart)/((IP_onShr-shrStart).Length()) )*shrDir;
      if (dot == 1 ) { DistToIP *= -1; }
    }
  }

  return;
}// function    


void ShowerCutCalculator::getAncestorMuonParams(const geoalgo::Point_t& shrStart,
						const geoalgo::Vector_t& shrDir,
						const geoalgo::Trajectory_t& muonTrack,
						double &Dist,
						double &IP,
						double &DistToIP) const
{
  
  // use shower's start point + direction and list of muon trajectory to find:
  // 1) PoCA and PoCA distance to shower start point w/ PoCA GeoAlgo class
  // 2) MuonDist for cylinder cut
  
  // initialize values for PoCA cut
  double minIP = std::numeric_limits<double>::max();
  geoalgo::Point_t IP_onShr(3);
  geoalgo::Point_t c1(3);
  geoalgo::Point_t c2(3);
  // initialize values for muon proximity
  double minDist = std::numeric_limits<double>::max();

  // use shower start & momentum to define a segment
  // which starts 3 meters before start point
  // and ends 10 cm after start point
  // aligned with shr momentum. Use for Impact Parameter
  geoalgo::LineSegment_t shrSeg(shrStart-shrDir*300,shrStart+shrDir*10);

  if ( muonTrack.size() > 1 ){

    double tmpDist = _geoAlgo.SqDist(shrStart,muonTrack);//DistanceToTrack(*shrStart, *muonTrack);
    double tmpIP = _geoAlgo.SqDist(muonTrack,shrSeg,c1,c2);//_PoCA.ClosestApproachToTrajectory(*muonTrack, shrOrigin, shrEnd, c1, c2);
    
    if (tmpDist < minDist) { minDist = tmpDist; }
    if (tmpIP < minIP) { 
      minIP = tmpIP; 
      IP_onShr = c2;
    }// if best muon so far

    DistToIP = IP_onShr.Dist(shrStart);
    Dist = sqrt(minDist);
    IP = sqrt(minIP);
    
    //need to figure out if IP point is "before" or "after" start point w.r.t. momentum direction
    if (DistToIP > 0.001){
      double dot = ( (IP_onShr-shrStart)/((IP_onShr-shrStart).Length()) )*shrDir;
      if (dot == 1 ) { DistToIP *= -1; }
    }
  }// if the muon's track is > 1
  
  return;
}// function    



void ShowerCutCalculator::getDistanceToWall(const geoalgo::Point_t& shrStart,
					    const geoalgo::Vector_t& shrDir,
					    double &distToWallForwards,
					    double &distToWallBackwards) const
{
  geoalgo::HalfLine_t sDir(shrStart,shrDir);

  distToWallForwards = sqrt(shrStart.SqDist(_geoAlgo.Intersection(_TpcBox,sDir)));//_DistToBoxWall.DistanceToWall(shrStart,shrDir,1);
  distToWallBackwards = sqrt(shrStart.SqDist(_geoAlgo.Intersection(_TpcBox,sDir,true)));//_DistToBoxWall.DistanceToWall(shrStart,shrDir,0);
  return;
}
  
#endif
