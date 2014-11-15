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

  void getNearestMuonParams(std::vector<double> *shrStart,
			    std::vector<double> *shrDir,
			    std::vector<std::vector<std::vector<double> > > *muonTracks,
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
  

};

#endif
/** @} */ // end of doxygen group 

