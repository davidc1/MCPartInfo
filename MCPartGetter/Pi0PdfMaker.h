/**
 * \file Pi0PdfMaker.h
 *
 * \ingroup MCPartGetter
 * 
 * \brief Class def header for a class Pi0PdfMaker
 *
 * @author David Caratelli
 */

/** \addtogroup MCPartGetter
    David Caratelli
    @{*/
#ifndef LARLITE_PI0PDFMAKER_H
#define LARLITE_PI0PDFMAKER_H

#include <iostream>
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include <iostream>

/**
   \class Pi0PdfMaker
   User defined class Pi0PdfMaker ... these comments are used to generate
   doxygen documentation!
 */
class Pi0PdfMaker {
  
 public:
  
  /// Default constructor
  Pi0PdfMaker(){};
  
  /// Default destructor
  virtual ~Pi0PdfMaker(){};

  RooAbsPdf* RadLenPdf(double xmin, double xmax, double lambda);

  RooAbsPdf* dEdxPdf(int pdg);

  /// Set "data" TTree
  void setTree(TTree *t) { tree = t; }
  
 private:

  // TTree with data sets in branches
  TTree *tree;

  
};

#endif
/** @} */ // end of doxygen group 

