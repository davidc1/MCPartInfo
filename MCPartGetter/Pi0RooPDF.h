/**
 * \file Pi0RooPDF.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class Pi0RooPDF
 *
 * @author David Caratelli
 */

/** \addtogroup MCPartGetter

    @{*/
#ifndef LARLITE_PI0ROOPDF_H
#define LARLITE_PI0ROOPDF_H

#include <iostream>
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooGenericPdf.h"
#include <iostream>

/**
   \class Pi0RooPDF
   User defined class Pi0RooPDF ... these comments are used to generate
   doxygen documentation!
 */
class Pi0RooPDF: public RooGenericPdf {
  
 public:
  
  /// Default constructor
  Pi0RooPDF(){};
  
 Pi0RooPDF(std::string n) : RooGenericPdf(n.c_str(),"exp(x/lambda)",
					  RooArgSet(RooRealVar("x","Distance [cm]",0,100),
						    RooRealVar("lambda","Radiation Length [cm]",14.) ) )
   
    {};
  
  /// Default destructor
  virtual ~Pi0RooPDF(){};
  
  /// Initialize variables
  void InitializeVariables();

  /// generate data sample from PDG
  void GenerateData(int nSamples);

 private:
  
  RooDataSet *_data;
  
  RooRealVar _x;
  
  RooRealVar _lambda;

};

#endif
/** @} */ // end of doxygen group 

