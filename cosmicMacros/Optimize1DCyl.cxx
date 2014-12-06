#include <iostream>

void Optimize1DCyl(){

  TCanvas *canv = new TCanvas("cCuts","cCuts",800,600);

  // specify number of bins
  const Int_t nBins = 18;
  // specify min E to plot
  const Double_t Emin = 100;
  // specify max E to plot
  const Double_t Emax = 1000;

  std::string SiglFile = "mcshowerSigFix.root";
  std::string BkgdFile = "mcshowerBkgd.root";

  //std::string fidVol = " && (X > 25) && (X < 231) && (Y > -91) && (Y < 91) && (Z > 5) && (Z < 936)";
  std::string fidVol = " && distAlongTraj > 70 && (X > 5) && (X < 251) && (Y > -111) && (Y < 111) && (Z > 5) && (Z < 1031)";

  std::string OptOver = "distBackAlongTraj";

  TFile sigl = TFile(SiglFile.c_str());
  TFile bkgd = TFile(BkgdFile.c_str());

  TTree *siglTree = sigl.Get("ana_tree");
  TTree *bkgdTree = bkgd.Get("ana_tree");

  double scanMinCut = 0;
  double scanMaxCut = 75;
  int Nscans = 500;

  double min = 1000;
  double max = -1000;



  TH1D *hEAllSigl;
  TH1D *hECylBackSigl;
  TH1D *hEAllBkgd;
  TH1D *hECylBackBkgd;

  Int_t Nevents = siglTree->GetEntries("inActiveVolume==1 && (Z>17) && (Z<1019) && (X>17) && (X<239) && (Y>-99) && (Y<99)");
  std::cout << "Number of Signal Showers: " << Nevents << std::endl;
  Double_t CCInclusive = 820;

  // Energy spectrum of all showers
  TH1D *hEAllSigl = new TH1D("hEAllSigl","Energy Spectrum",nBins,Emin,Emax);
  siglTree->Draw("E>>hEAllSigl",Form("inActiveVolume==1 %s",fidVol.c_str()));
  hEAllSigl->Scale(CCInclusive/Nevents);
  //fill array
  Double_t SignalNoCuts[nBins];
  for (int s=0; s < nBins; s++){
    SignalNoCuts[s] = hEAllSigl->GetBinContent(s+1);
  }
  //significance of signal w/ only fiducial volume cuts
  double sensOptimal = 0;
  for (int u=0; u < nBins; u++){
    if (SignalNoCuts[u] != 0)
      sensOptimal += sqrt(SignalNoCuts[u]);
  }
  sensOptimal /= nBins;
  std::cout << "Best chi^2: " << sensOptimal << std::endl;

  TH1D *hOpt = new TH1D("hOpt",Form("Significance of Signal - Diff. w.r.t. Signal Only Significance [%.2lf]; Cut on Cylinder Dist. [cm]",sensOptimal),
			Nscans, scanMinCut, scanMaxCut);

  //loop for optimization
  for (int dCut = 0; dCut < Nscans; dCut++){

    double Cut = scanMinCut+(dCut/(float)Nscans)*(scanMaxCut-scanMinCut);

    std::cout << "Cut value: " << Cut << std::endl;
    
    //----------------------------//
    //-----------SIGNAL-----------//
    //----------------------------//
    
    
    if (hECylBackSigl)
      delete hECylBackSigl;
    // All showers with Cyl > 10 cm & Back > 20 cm
    TH1D *hECylBackSigl = new TH1D("hECylBackSigl","Energy Spectrum - Dist to Wall Back and Cyl",nBins,Emin,Emax);
    siglTree->Draw("E>>hECylBackSigl",Form("inActiveVolume==1 && minMuDist>=%g %s",Cut,fidVol.c_str()));

    // scale histograms to 6.6E20 POT
    hECylBackSigl->Scale(CCInclusive/Nevents);

    Double_t SignalEvts[nBins];
    for (int s=0; s < nBins; s++)
      SignalEvts[s] = hECylBackSigl->GetBinContent(s+1);
    
    //----------------------------//
    //---------BACKGROUND---------//
    //----------------------------//
    
    // get number of events from histogram
    TH1D *muons = bkgd.Get("hTrackTotLen");
    double evts = muons->GetEntries();
    Double_t Time = evts*(6.4E-3); //in seconds
    Double_t ExperimentTime = (1.6E-6)*(6.6E20)/(5E12);
    
    if (hECylBackBkgd)
      delete hECylBackBkgd;
    
    // All showers with Cyl > 10 cm & Back > 20 cm
    TH1D *hECylBackBkgd = new TH1D("hECylBackBkgd","Energy Spectrum - Dist to Wall Back and Cyl",nBins,Emin,Emax);
    
    bkgdTree->Draw("E>>hECylBackBkgd",Form("inActiveVolume==1 && minMuDist>=%g %s",Cut,fidVol.c_str()));
    //bkgdTree->Draw("E>>hECylBackBkgd",Form("inActiveVolume==1  && (X > 25) && (X < 231) && (Y > -91) && (Y < 91) && (Z > 5) && (Z < 936)",Cut));

    //scale histogram assuming that surviving showers are all photon
    // *(0.04+0.96*0.06)
    //    hECylBackBkgd->Scale(0.04+0.96*0.06);

    
    // scale histograms to 6.6E20 POT
    hECylBackBkgd->Scale(ExperimentTime/Time);

    Double_t BackgroundEvts[nBins];
    for (int s=0; s < nBins; s++)
      BackgroundEvts[s] = hECylBackBkgd->GetBinContent(s+1);
    /*
    std::cout << "Cut Values. " << OptOver.c_str() << ": " << Cut << std::endl;
    std::cout << "Surviving Events. Signal: " << SignalEvts[1]
	      << ", Background: " << BackgroundEvts[1] << std::endl;
    std::cout << std::endl;
    */
    
    //loop over histogram bins
    double opt = 0;
    for (int b=0; b < nBins; b++){
      if (SignalEvts[b] != 0){
	double err = sqrt( BackgroundEvts[b] + SignalEvts[b] );
	opt += SignalEvts[b]/err;
      }
    }
    opt /= nBins;
    if (opt-sensOptimal < min) { min = opt-sensOptimal; }
    if (opt-sensOptimal > max) { max = opt-sensOptimal; }
    std::cout << "OPT: " << opt << std::endl;
    hOpt->SetBinContent(dCut+1,opt-sensOptimal);
  }//loop over cut parameters


  gStyle->SetOptStat(0);
  hOpt->Draw();
  hOpt->GetYaxis()->SetRangeUser(max-1,max+0.1);
  hOpt->SetLineWidth(2);
  hOpt->SetLineColor(4);
  hOpt->SetFillColor(4);
  hOpt->SetFillStyle(3003);
  canv->Modified();
  canv->Update();
  canv->SaveAs("optCylinderBigFV.pdf");

  return;
}
