/**
 * \file ShowerCutCalculator.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class ShowerCutCalculator
 *
 * @author David Caratelli
 */

/** \addtogroup MCPartGetter

    @{*/
#ifndef LARLITE_SHOWERCUTCALCULATOR_H
#define LARLITE_SHOWERCUTCALCULATOR_H

#include <iostream>
#include "LArUtil/Geometry.h"
#include <time.h>
// include GeoAlgo functions
#include "BasicTool/GeoAlgo/TrajectoryInVolume.h"
#include "BasicTool/GeoAlgo/PointToLineDist.h"
#include "BasicTool/GeoAlgo/TwoLineIntersection.h"
#include "BasicTool/GeoAlgo/SegmentPoCA.h"
#include "BasicTool/GeoAlgo/DistToBoxWall.h"
#include "BasicTool/GeoAlgo/TrajectoryInVolume.h"

/**
   \class ShowerCutCalculator
   User defined class ShowerCutCalculator ... these comments are used to generate
   doxygen documentation!
 */
class ShowerCutCalculator{

 public:
  
  /// Default constructor
  ShowerCutCalculator(){};
  
  /// Default destructor
  virtual ~ShowerCutCalculator(){};

  void SetAlgoProperties();

  bool isInVolume(std::vector<double> point);

  void getNearestMuonParams(std::vector<double> *shrStart,
			    std::vector<double> *shrDir,
			    std::vector<std::vector<std::vector<double> > > *muonTracks,
			    std::vector<int> *muonIDs,
			    int ancestorID,
			    double &Dist,
			    double &IP,
			    double &DistToIP);

  void getAncestorMuonParams(std::vector<double> *shrStart,
			     std::vector<double> *shrDir,
			     std::vector<std::vector<double> > *muonTrack,
			     double &Dist,
			     double &IP,
			     double &DistToIP);


  void getDistanceToWall(std::vector<double> shrStart,
			 std::vector<double> shrDir,
			 double &distToWallForwards,
			 double &distToWallBackwards);


 private:

  // GeoAlgo algorithm used
  geoalgo::DistToBoxWall _DistToBoxWall;
  /// GeoAlg for TPC containment
  geoalgo::TrajectoryInVolume _inTPCAlgo;
  /// GeoAlg for point to line dist
  geoalgo::PointToLineDist _pointDist;
  /// GeoAlg for poka cut
  geoalgo::TwoLineIntersection _lineIntersection;
  /// GeoAlg for PoCA cut
  geoalgo::SegmentPoCA _PoCA;
  /// GeoAlg for InVolume
  geoalgo::TrajectoryInVolume _inVol ;  
  

};

#endif
/** @} */ // end of doxygen group 

