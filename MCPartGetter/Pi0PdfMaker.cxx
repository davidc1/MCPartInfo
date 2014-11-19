#ifndef PI0PDFMAKER_CXX
#define PI0PDFMAKER_CXX

#include "Pi0PdfMaker.h"

RooAbsPdf* Pi0PdfMaker::RadLenPdf(double xmin, double xmax, double lambda){

  RooRealVar x("x","Distance [cm]",xmin,xmax);
  RooRealVar l("l","Radiation Length [cm]",lambda);
  RooExponential *expo = new RooExponential("expo","Radiation Length", x, l);

  return expo;
}


RooAbsPdf* Pi0PdfMaker::dEdxPdf(int pdg){
  
  RooRealVar dedx("dEdx","dEdx [MeV/cm]",0,10);
  RooRealVar part("part","Particle PDG",-1000,1000);
  RooRealVar sigma("sigma","Spread in dEdx",1);
  if (pdg==22)
    part.setVal(4);
  else if ( (pdg==11) or (pdg==-11) )
    part.setVal(2);
  else
    part.setVal(0);
  RooGaussian *gauss = new RooGaussian("gauss","dEdx Gauss Pdf",dedx,part,sigma);

  return gauss;
}

#endif
