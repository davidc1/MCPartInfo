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
#include "BasicTool/GeoAlgo/GeoAlgoConstants.h"
#include "BasicTool/GeoAlgo/DistanceAlgo.h"
#include "BasicTool/GeoAlgo/IntersectAlgo.h"

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

  bool isInVolume(const std::vector<double>& point);

  void getNearestMuonParams(const geoalgo::Point_t& shrStart,
			    const geoalgo::Vector_t& shrDir,
			    const std::vector<geoalgo::Trajectory_t >& muonTracks,
			    const std::vector<int>& muonIDs,
			    const int ancestorID,
			    double &Dist,
			    double &IP,
			    double &DistToIP) const;

  void getAncestorMuonParams(const geoalgo::Point_t& shrStart,
			     const geoalgo::Vector_t& shrDir,
			     const geoalgo::Trajectory_t& muonTrack,
			     double &Dist,
			     double &IP,
			     double &DistToIP) const;


  void getDistanceToWall(const geoalgo::Point_t& shrStart,
			 const geoalgo::Vector_t& shrDir,
			 double &distToWallForwards,
			 double &distToWallBackwards) const;


 private:

  // GeoAlgo Distance Algo
  geoalgo::DistanceAlgo _dAlgo;
  // geoalgo Intersection Algo
  geoalgo::IntersectAlgo _iAlgo;

  // TPC AABos
  geoalgo::AABox _TpcBox;


};

#endif
/** @} */ // end of doxygen group 

