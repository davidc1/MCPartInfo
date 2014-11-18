#ifndef PI0ROOPDF_CXX
#define PI0ROOPDF_CXX

#include "Pi0RooPDF.h"

void Pi0RooPDF::InitializeVariables(){

  _x = RooRealVar("x","Distance [cm]",0,100);

  _lambda = RooRealVar("lambda","Radiation Length [1/cm]",1./14.0);

  return;
}


void Pi0RooPDF::GenerateData(int nSamples){

  _data = this->generate(_x,nSamples);

  return;
}

#endif
