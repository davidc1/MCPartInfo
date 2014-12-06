#include <iostream>

void Optimize(){

  TCanvas *canv = new TCanvas("cCuts","cCuts",900,600);

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

  TFile sigl = TFile(SiglFile.c_str());
  TFile bkgd = TFile(BkgdFile.c_str());

  TTree *siglTree = sigl.Get("ana_tree");
  TTree *bkgdTree = bkgd.Get("ana_tree");

  Int_t Nevents = siglTree->GetEntries("inActiveVolume==1 && (Z>17) && (Z<1019) && (X>17) && (X<239) && (Y>-99) && (Y<99)");
  Double_t CCInclusive = 820;

  double scanMinCyl = 0;
  double scanMinBk  = 0;
  double scanMaxCyl = 50;
  double scanMaxBk  = 100;
  int Nscans = 50;

  double min = 1000;
  double max = -1000;
  Double_t bestX[1];
  Double_t bestY[1];

  //First get chi^2 for signal only
  TH1D *hSignal = new TH1D("hSignal","Energy Spectrum",nBins,Emin,Emax);
  double chi = 0;
  siglTree->Draw("E>>hSignal",Form("inActiveVolume==1 %s",fidVol.c_str()));
  hSignal->Scale(CCInclusive/Nevents);
  for (int h=0; h < nBins; h++)
    chi += sqrt(hSignal->GetBinContent(h+1));
  chi /= nBins;
 std:cout << "Optimal chi^2: " << chi << std::endl;


  TH2D *hOpt = new TH2D("hOpt",Form("Significance of Signal - Diff. w.r.t. Signal Only [%.2lf]; Cut on Muon-Cylinder Dist [cm]; Cut on Back-To-Wall [cm]",chi),
			Nscans-1, scanMinCyl, scanMaxCyl, Nscans-1, scanMinBk, scanMaxBk);

  TH1D *hECylBackSigl;
  TH1D *hECylBackBkgd;

  //loop for optimization
  for (int dCyl = 0; dCyl < Nscans; dCyl++){

    double cylCut = scanMinCyl+(dCyl/(float)Nscans)*(scanMaxCyl-scanMinCyl);

    std::cout << "Iteration: " << dCyl << " / " << Nscans << std::endl; 
    
    for (int dBack = 0; dBack < Nscans; dBack++){

      double backCut = scanMinBk+(dBack/(float)Nscans)*(scanMaxBk-scanMinBk);


      //----------------------------//
      //-----------SIGNAL-----------//
      //----------------------------//
      

      
      if (hECylBackSigl)
	delete hECylBackSigl;
      // All showers with Cyl > 10 cm & Back > 20 cm
      TH1D *hECylBackSigl = new TH1D("hECylBackSigl","Energy Spectrum - Dist to Wall Back and Cyl",nBins,Emin,Emax);

      siglTree->Draw("E>>hECylBackSigl",Form("inActiveVolume==1 && distBackAlongTraj>%g && minMuDist>=%g %s",backCut,cylCut,fidVol.c_str()));

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
      
      bkgdTree->Draw("E>>hECylBackBkgd",Form("inActiveVolume==1 && distBackAlongTraj>%g && minMuDist>=%g %s",backCut,cylCut,fidVol.c_str()));

      // scale histograms to 6.6E20 POT
      hECylBackBkgd->Scale(ExperimentTime/Time);
      
      Double_t BackgroundEvts[nBins];
      for (int s=0; s < nBins; s++)
	BackgroundEvts[s] = hECylBackBkgd->GetBinContent(s+1);

      /*
      std::cout << "Cut Values. BDtW: " << backCut << ", Cylinder: " << cylCut << std::endl;
      std::cout << "Surviving Events. Signal: " << SignalEvts[1]
		<< ", Background: " << BackgroundEvts[1] << std::endl;
      std::cout << std::endl;
      */

      //loop over histogram bins
      double opt = 0;
      for (int b=0; b < nBins; b++){
	if ( (BackgroundEvts[b]==0) && (SignalEvts[b]== 0))
	  std::cout << "Divide by 0!" << std::endl;
	double err = sqrt( BackgroundEvts[b] + SignalEvts[b] );
	opt += SignalEvts[b]/err;
      }
      opt /= nBins;
      if ((opt-chi) < min) { min = opt-chi; }
      if ((opt-chi) > max) { max = opt-chi; bestX[0] = cylCut; bestY[0] = backCut; }
      hOpt->SetBinContent(dCyl,dBack,opt-chi);
    }//for back-to-wall loop
  }//for cylinder loop


  TGraph *optg = new TGraph(1,bestX,bestY);
  optg->SetMarkerStyle(20);
  optg->SetMarkerSize(3);
  optg->SetMarkerColor(3);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBorderMode(0);
  hOpt->Draw("COLZ");
  optg->Draw("P");
  leg = new TLegend(0.6,0.7,0.9,0.8);
  leg->AddEntry(optg,Form("Optimal Value: %g, %g",bestX[0],bestY[0]));
  leg->Draw();
  leg->SetFillColor(kWhite);
  hOpt->GetZaxis()->SetTitleOffset(0.1);
  //  hOpt->GetZaxis()->SetTitle("Chi^2 / Bin");
  hOpt->GetZaxis()->SetRangeUser(max-1,max+0.1);
  //canv->SetLogz();
  canv->Modified();
  canv->Update();
  canv->SaveAs("optBigFV.pdf");

  return;
}
