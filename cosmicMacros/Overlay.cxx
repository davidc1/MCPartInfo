#include <iostream>

void Overlay(){

  TCanvas *canv = new TCanvas("cCuts","cCuts",800,600);

  // specify number of bins
  const Int_t nBins = 18;
  // specify min E to plot
  const Double_t Emin = 100;
  // specify max E to plot
  const Double_t Emax = 1000;

  //cut value
  double cutvalCyl = 13;
  double cutvalBack = 34;

  std::string SiglFile = "mcshowerSigFix.root";
  std::string BkgdFile = "mcshowerBkgd.root";

  std::string OptOver = "distBackAlongTraj";

  TFile sigl = TFile(SiglFile.c_str());
  TFile bkgd = TFile(BkgdFile.c_str());

  TTree *siglTree = sigl.Get("ana_tree");
  TTree *bkgdTree = bkgd.Get("ana_tree");

  //std::string fidVol = " && (X > 25) && (X < 231) && (Y > -91) && (Y < 91) && (Z > 5) && (Z < 936)";
  //std::string fidVol = " && (X > 25) && (X < 231) && (Y > -91) && (Y < 91) && (Z > 30) && (Z < 936)";
  std::string fidVol = " && distAlongTraj > 70 && (X > 5) && (X < 251) && (Y > -111) && (Y < 111) && (Z > 5) && (Z < 1031)";


  double min = 1000;
  double max = 0;

  TH1D *hEAllSigl;
  TH1D *hECylBackSigl;
  TH1D *hECylBackBkgd;

  //----------------------------//
  //-----------SIGNAL-----------//
  //----------------------------//
  
  Int_t Nevents = siglTree->GetEntries("inActiveVolume==1 && (Z>17) && (Z<1019) && (X>17) && (X<239) && (Y>-99) && (Y<99)");
  Double_t CCInclusive = 820;
  
  if (hEAllSigl)
    delete hEAllSigl;
  if (hECylBackSigl)
    delete hECylBackSigl;
  // Energy spectrum of all showers
  TH1D *hEAllSigl = new TH1D("hEAllSigl","Energy Spectrum",nBins,Emin,Emax);
  // All showers with Cyl > 10 cm & Back > 20 cm
  TH1D *hECylBackSigl = new TH1D("hECylBackSigl","Energy Spectrum - Dist to Wall Back and Cyl",nBins,Emin,Emax);

  siglTree->Draw("E>>hEAllSigl",Form("inActiveVolume==1 %s",fidVol.c_str()));
  siglTree->Draw("E>>hECylBackSigl",Form("inActiveVolume==1 && minMuDist>%g && distBackAlongTraj>%g %s",cutvalCyl,cutvalBack,fidVol.c_str()));

  // scale histograms to 6.6E20 POT
  hEAllSigl->Scale(CCInclusive/Nevents);
  hECylBackSigl->Scale(CCInclusive/Nevents);


  Double_t SignalNoCuts[nBins];
  for (int s=0; s < nBins; s++){
    SignalNoCuts[s] = hEAllSigl->GetBinContent(s+1);
  }

  double tot_sig = 0;
  Double_t SignalEvts[nBins];
  for (int s=0; s < nBins; s++){
    SignalEvts[s] = hECylBackSigl->GetBinContent(s+1);
    tot_sig += SignalEvts[s];
  }
  std::cout << "Signal Events: " << tot_sig << std::endl;


  //significance of signal w/ only fiducial volume cuts
  double sensOptimal = 0;
  for (int u=0; u < nBins; u++){
    if (SignalNoCuts[u] != 0)
      sensOptimal += sqrt(SignalNoCuts[u]);
  }
  sensOptimal /= nBins;
    
  //----------------------------//
  //---------BACKGROUND---------//
  //----------------------------//
  
  // get number of events from histogram
  TH1D *muons = bkgd.Get("hTrackTotLen");
  double evts = muons->GetEntries();
  Double_t Time = evts*(6.4E-3); //in seconds
  Double_t ExperimentTime = (1.6E-6)*(6.6E20)/(5E12);
  
  // All showers with Cyl > 10 cm & Back > 20 cm
  TH1D *hECylBackBkgd = new TH1D("hECylBackBkgd","Energy Spectrum - Dist to Wall Back and Cyl",nBins,Emin,Emax);
  //  TH1D *h22ECylBackBkgd = new TH1D("h22ECylBackBkgd","Energy Spectrum - Dist to Wall Back and Cyl",nBins,Emin,Emax);
  
  //  bkgdTree->Draw("E>>hECylBackBkgd","inActiveVolume==1 && (X > 25) && (X < 231) && (Y > -91) && (Y < 91) && (Z > 5) && (Z < 936)");
  bkgdTree->Draw("E>>hECylBackBkgd",Form("inActiveVolume==1 && minMuDist>%g && distBackAlongTraj>%g %s",cutvalCyl,cutvalBack,fidVol.c_str()));
  //bkgdTree->Draw("E>>hECylBackBkgd",Form("inActiveVolume==1 && minMuDist>%g && distBackAlongTraj>%g && (X > 25) && (X < 231) && (Y > -91) && (Y < 91) && (Z > 5) && (Z < 936)",0.,0.));
  //  bkgdTree->Draw("E>>h22ECylBackBkgd",Form("inActiveVolume==1 && _process==\"conv\" && minMuDist>%g && distBackAlongTraj>%g && (X > 25) && (X < 231) && (Y > -91) && (Y < 91) && (Z > 5) && (Z < 936)",cutvalCyl,cutvalBack));

  //  for (int x=0; x < nBins; x++)
  //    hECylBackBkgd->SetBinContent(x+1,hECylBackBkgd->GetBinContent(x+1)-0.94*h22ECylBackBkgd->GetBinContent(x+1));

  // scale histograms to 6.6E20 POT
  hECylBackBkgd->Scale(ExperimentTime/Time);

  Double_t BackgroundEvts[nBins];
  double bkgd_evts = 0;
  for (int s=0; s < nBins; s++){
    BackgroundEvts[s] = hECylBackBkgd->GetBinContent(s+1);
    bkgd_evts += BackgroundEvts[s];
  }
  std::cout << "Background events: " << bkgd_evts << std::endl;

  //Calculate Statistica Significance:
  double sig = 0;
  for (int h=0; h < nBins; h++){
    if (SignalEvts[h] != 0){
      double err = sqrt( SignalEvts[h] + BackgroundEvts[h] );
      sig += SignalEvts[h] / err;
    }
  }
  sig /= nBins;
      
  std::cout << "Loss in Significance [chi^2/Bin]: " << sig-sensOptimal << std::endl;

  //Stacked histogram
  THStack *plot = new THStack("plot",Form("BNB Ve Signal w/ Cosmic Bkgd - Loss in Significance: %.2lf; Electron Shower E [ GeV ]; Events - for 6.6E20 POT",sig-sensOptimal));

  //  THStack *plot = new THStack("plot","BNB Ve Signal - No Cosmic Bkgd; Electron Shower E [ GeV ]; Events / 50 MeV - for 6.6E20 POT");

  // Histogram to be plotted
  //TH1D *hSig = new TH1D("hSig",Form("BNB Ve Signal w/ Cosmic Bkgd - Loss in Significance: %.2lf; Electron Shower E [ GeV ]; Events - for 6.6E20 POT",sig-sensOptimal)
  TH1D *hSig = new TH1D("hSig",Form("BNB Ve Signal w/o Cosmic Bkgd - Significance [Chi^2/Bin]: %.2lf; Electron Shower E [ GeV ]; Events - for 6.6E20 POT",sensOptimal)
			,nBins,Emin/1000.,Emax/1000.);
  for (int s=0; s < nBins; s++){
    //hSig->SetBinContent(s+1,SignalEvts[s]+BackgroundEvts[s]);
    //hSig->SetBinError(s+1,sqrt(SignalEvts[s]+BackgroundEvts[s]));
    hSig->SetBinContent(s+1,SignalEvts[s]);
    hSig->SetBinError(s+1,sqrt(SignalEvts[s]));
  }
  hSig->SetLineColor(4);
  hSig->SetLineWidth(2);
  hSig->SetFillColor(4);
  hSig->SetFillStyle(3002);
  // Histogram to be plotted
  TH1D *hBkgd = new TH1D("hBkgd","Background Events",nBins,Emin/1000.,Emax/1000.);
  for (int s=0; s < nBins; s++){
    hBkgd->SetBinContent(s+1,BackgroundEvts[s]);
  }
  hBkgd->SetLineColor(2);
  hBkgd->SetLineWidth(2);
  hBkgd->SetFillColor(2);
  hBkgd->SetFillStyle(3002);

  leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry(hSig,"Signal Events");
  //leg->AddEntry(hBkgd,"Background Events");
  leg->SetFillColor(kWhite);

  hSig->Draw("E");
  hSig->GetYaxis()->SetRangeUser(0,60);
  gStyle->SetOptStat(0);
  //hBkgd->Draw("same");
  //  plot->Add(hBkgd);
  //  plot->Add(hSig);
  //  plot->Draw("E");
  //  plot->GetYaxis()->SetRangeUser(0,75);
  leg->Draw();
  /*
  gStyle->SetOptStat(0);
  hOpt->Draw();
  hOpt->GetYaxis()->SetRangeUser(max-1,max+0.1);
  hOpt->SetLineWidth(2);
  hOpt->SetFillStyle(3003);
  */
  canv->Modified();
  canv->Update();
  canv->SaveAs("SignalBigFV.pdf");

  return;
}
