#include <iostream>

void MakeCosmogenicPlotsSignal(){

  TCanvas *cCuts = new TCanvas("cCuts","cCuts",800,600);

  // specify number of bins
  const Int_t nBins = 18;
  // specify min E to plot
  const Double_t Emin = 100;
  // specify max E to plot
  const Double_t Emax = 1000;
  // specify cylinder cut distance [ cm ]
  Int_t CylDist = 13;
  // specify IP cut distance [ cm ]
  Int_t IPDist = 13;
  // specify back to wall cut distance [ cm ]
  Int_t BackDist = 34;

  std::string inFile = "mcshowerSigFix.root";

  //std::string fidVol = " && (X > 25) && (X < 231) && (Y > -91) && (Y < 91) && (Z > 5) && (Z < 936)";
  std::string fidVol = " && distAlongTraj > 70 && (X > 5) && (X < 251) && (Y > -111) && (Y < 111) && (Z > 5) && (Z < 1031)";

  std::cout << "Opening File: " <<  inFile << std::endl;
  TFile f = TFile(inFile.c_str());

  TTree *tree = f.Get("ana_tree");

  // get number of events from histogram
  Int_t Nevents = tree->GetEntries("inActiveVolume==1 && (Z>17) && (Z<1019) && (X>17) && (X<239) && (Y>-99) && (Y<99)");
  TH1D *hMuonLength = f.Get("hTrackTotLen");
  double avgMuonLen = hMuonLength->GetMean();
  Double_t CCInclusive = 820;

  // Energy spectrum of all showers
  TH1D *h = new TH1D("h","Energy Spectrum",nBins,Emin,Emax);
  // Energy spectrum of all showers
  TH1D *hEAll = new TH1D("hEAll","Energy Spectrum",nBins,Emin,Emax);
  // All showers > 10 cm Away from a track
  TH1D *hECyl = new TH1D("hECyl","Energy Spectrum - Cylinder Cut",nBins,Emin,Emax);
  // All showers with IP > 10 cm
  TH1D *hEIP = new TH1D("hEIP","Energy Spectrum - IP Cut",nBins,Emin,Emax);
  // All showers with distBackToWall > 20 cm
  TH1D *hEBack = new TH1D("hEBack","Energy Spectrum - Dist to Wall Back",nBins,Emin,Emax);
  // All showers with Cyl > 10 cm & Back > 20 cm
  TH1D *hECylBack = new TH1D("hECylBack","Energy Spectrum - Dist to Wall Back and Cyl",nBins,Emin,Emax);
  // All showers with IP > 10 cm & Back > 20 cm
  TH1D *hEIPBack = new TH1D("hEIPBack","Energy Spectrum - Dist to Wall Back and IP",nBins,Emin,Emax);

  tree->Draw("E>>h",Form("inActiveVolume==1 %s",fidVol.c_str()));
  tree->Draw("E>>hEAll",Form("inActiveVolume==1 %s",fidVol.c_str()));
  tree->Draw("E>>hECyl",Form("inActiveVolume==1 && minMuDist>%d %s",CylDist,fidVol.c_str()));
  tree->Draw("E>>hEIP",Form("inActiveVolume==1 && minMuIP>%d %s ",IPDist,fidVol.c_str()));
  tree->Draw("E>>hEBack",Form("inActiveVolume==1 && distBackAlongTraj>%d %s",BackDist,fidVol.c_str()));
  tree->Draw("E>>hECylBack",Form("inActiveVolume==1 && distBackAlongTraj>%d && minMuDist>%d %s",BackDist,CylDist,fidVol.c_str()));
  tree->Draw("E>>hEIPBack",Form("inActiveVolume==1 && distBackAlongTraj>%d && minMuIP>%d %s",BackDist,IPDist,fidVol.c_str()));


  //get distance-based performance of cylinder cut
  TH1D *hCyl = new TH1D("hCyl","hCyl",10,0,1);
  TGraph CylPerformance[5];
  TLine CylBackOfEnvelope[5];
  for (Int_t i=0; i < 5; i++){
    if (hCyl)
      delete hCyl;
    TH1D *hCyl = new TH1D("hCyl","Energy Spectrum - Cylinder Cut",nBins,Emin,Emax);
    Int_t CutVal = i*6+1;
    tree->Draw("E>>hCyl",Form("inActiveVolume==1 && minMuDist>%d %s",CutVal,fidVol.c_str()));
    Double_t Espectrum[nBins];
    Double_t Eff[nBins];
    for (int n=0; n < hEAll->GetNbinsX(); n++){
      Espectrum[n] = hEAll->GetBinCenter(n+1)/1000.;
      Eff[n] = 1-hCyl->GetBinContent(n+1)/hEAll->GetBinContent(n+1);
    }
    CylPerformance[i] = new TGraph(nBins, Espectrum, Eff);
  }
  cCuts->Clear();
  legCyl = new TLegend(0.5,0.6,0.9,0.9);
  legCyl->SetFillColor(kWhite);
  for (Int_t u=0; u < 5; u++){
    if (u==0){
      CylPerformance[u].Draw();
      CylPerformance[u].SetMarkerColor(2+u);
      CylPerformance[u].SetMarkerStyle(20);
      CylPerformance[u].GetXaxis()->SetTitle("Energy [ GeV ]");
      CylPerformance[u].GetYaxis()->SetTitle("Cut Efficiency [ Frac. ]");
      CylPerformance[u].GetYaxis()->SetRangeUser(0,0.5);
      CylPerformance[u].SetTitle("Eff. of Cylinder Cut on Signal");
      double th = ((u*6+1)*(u*6+1)*3.14*avgMuonLen*100)/(256*232*1036.);
      CylBackOfEnvelope[u] = new TLine(Emin/1000.,th,Emax/1000.,th);
      std::cout << "theory: " << th << std::endl;
      CylBackOfEnvelope[u].SetLineColor(2+u);
      CylBackOfEnvelope[u].SetLineStyle(2);
      //CylBackOfEnvelope[u].Draw();
      legCyl->AddEntry(&CylPerformance[u],Form("Cylinder Cut Eff. on Signal - %d cm",u*6+1));
    }
    else{
      CylPerformance[u].Draw("PL");
      CylPerformance[u].SetMarkerColor(2+u);
      CylPerformance[u].SetMarkerStyle(20);
      double th = ((u*6+1)*(u*6+1)*3.14*avgMuonLen*100)/(256*232*1036.);
      CylBackOfEnvelope[u] = new TLine(Emin/1000.,th,Emax/1000.,th);
      std::cout << "theory: " << th << std::endl;
      CylBackOfEnvelope[u].SetLineColor(2+u);
      //      CylBackOfEnvelope[u].SetLineStyle(2);
      //CylBackOfEnvelope[u].Draw();
      legCyl->AddEntry(&CylPerformance[u],Form("Cylinder Cut Eff. on Signal. - %d cm", u*6+1));
    }
  }
  legCyl->Draw();
  cCuts->SaveAs("CyleffSignalOneFrameBigFV.pdf");
  cCuts->Clear();


  //get distance-based performance of cylinder cut
  TH1D *hBack = new TH1D("hBack","hBack",10,0,1);
  TGraph BackPerformance[5];
  for (Int_t i=0; i < 5; i++){
    if (hBack)
      delete hBack;
    TH1D *hBack = new TH1D("hBack","Energy Spectrum - Backinder Cut",nBins,Emin,Emax);
    Int_t CutVal = i*20+1;
    tree->Draw("E>>hBack",Form("inActiveVolume==1 && distBackAlongTraj>%d %s",CutVal,fidVol.c_str()));
    Double_t Espectrum[nBins];
    Double_t Eff[nBins];
    std::cout << "No fid vol: " << h->GetBinContent(1)
	      << "\tFid Vol: " << hEAll->GetBinContent(1)
	      << "\tCut: " << hBack->GetBinContent(1) << std::endl;
    for (int n=0; n < hEAll->GetNbinsX(); n++){
      Espectrum[n] = hEAll->GetBinCenter(n+1)/1000.;
      Eff[n] = 1-hBack->GetBinContent(n+1)/hEAll->GetBinContent(n+1);
    }
    BackPerformance[i] = new TGraph(nBins, Espectrum, Eff);
  }
  cCuts->Clear();
  legBack = new TLegend(0.5,0.6,0.9,0.9);
  legBack->SetFillColor(kWhite);
  for (Int_t u=0; u < 5; u++){
    if (u==0){
      BackPerformance[u].Draw();
      BackPerformance[u].SetMarkerColor(2+u);
      BackPerformance[u].SetMarkerStyle(20);
      BackPerformance[u].GetXaxis()->SetTitle("Energy [ GeV ]");
      BackPerformance[u].GetYaxis()->SetTitle("Cut Efficiency [ Frac. ]");
      BackPerformance[u].GetYaxis()->SetRangeUser(0,0.5);
      BackPerformance[u].SetTitle("Eff. of Back-To-Wall Cut on Signal");
      legBack->AddEntry(&BackPerformance[u],Form("Back Dist. Cut Eff. on Signal - %d cm", u*20+1));
    }
    else{
      BackPerformance[u].Draw("PL");
      BackPerformance[u].SetMarkerColor(2+u);
      BackPerformance[u].SetMarkerStyle(20);
      legBack->AddEntry(&BackPerformance[u],Form("Back Dist. Cut Eff. on Signal - %d cm", u*20+1));
    }
  }
  legBack->Draw();
  cCuts->SaveAs("BackeffSignalOneFrameBigFV.pdf");
  cCuts->Clear();



  //get distance-based performance of cylinder cut
  TH1D *hIP = new TH1D("hIP","hIP",10,0,1);
  TGraph IPPerformance[5];
  for (Int_t i=0; i < 5; i++){
    if (hIP)
      delete hIP;
    TH1D *hIP = new TH1D("hIP","Energy Spectrum - IP Cut",nBins,Emin,Emax);
    Int_t CutVal = i*6+1;
    tree->Draw("E>>hIP",Form("inActiveVolume==1 && minMuIP>%d %s",CutVal,fidVol.c_str()));
    Double_t Espectrum[nBins];
    Double_t Eff[nBins];
    for (int n=0; n < hEAll->GetNbinsX(); n++){
      Espectrum[n] = hEAll->GetBinCenter(n+1)/1000.;
      Eff[n] = 1-hIP->GetBinContent(n+1)/hEAll->GetBinContent(n+1);
    }
    IPPerformance[i] = new TGraph(nBins, Espectrum, Eff);
  }
  cCuts->Clear();
  legIP = new TLegend(0.5,0.6,0.9,0.9);
  legIP->SetFillColor(kWhite);
  for (Int_t u=0; u < 5; u++){
    if (u==0){
      IPPerformance[u].SetMarkerColor(2+u);
      IPPerformance[u].SetMarkerStyle(20);
      IPPerformance[u].GetXaxis()->SetTitle("Energy [ GeV ]");
      IPPerformance[u].GetYaxis()->SetTitle("Cut Efficiency [ Frac. ]");
      IPPerformance[u].GetYaxis()->SetRangeUser(0,1);
      IPPerformance[u].SetTitle("Eff. of Impact Param. Cut on Signal");
      IPPerformance[u].Draw();
      legIP->AddEntry(&IPPerformance[u],Form("Impact Param Cut Eff. on Signal - %d cm", u*6+1));
    }
    else{
      IPPerformance[u].SetMarkerColor(2+u);
      IPPerformance[u].SetMarkerStyle(20);
      IPPerformance[u].Draw("PL");
      legIP->AddEntry(&IPPerformance[u],Form("Impact Param Cut Eff. on Signal - %d cm", u*6+1));
    }
  }
  legIP->Draw();
  cCuts->SaveAs("IPeffSignalOneFrameBigFV.pdf");
  cCuts->Clear();
    


  // get fractional efficiency of various cuts
  Double_t Espectrum[nBins];
  Double_t EffCyl[nBins];
  Double_t EffIP[nBins];
  Double_t EffBack[nBins];
  Double_t EffCylBack[nBins];
  Double_t EffIPBack[nBins];

  for (int n=0; n < hEAll->GetNbinsX(); n++){
    Espectrum[n] = hEAll->GetBinCenter(n+1)/1000.;
    EffCyl[n] = 1-hECyl->GetBinContent(n+1)/hEAll->GetBinContent(n+1);
    EffIP[n] = 1-hEIP->GetBinContent(n+1)/hEAll->GetBinContent(n+1);
    EffBack[n] = 1-hEBack->GetBinContent(n+1)/hEAll->GetBinContent(n+1);
    EffCylBack[n] = 1-hECylBack->GetBinContent(n+1)/hEAll->GetBinContent(n+1);
    EffIPBack[n] = 1-hEIPBack->GetBinContent(n+1)/hEAll->GetBinContent(n+1);
  }



  effCyl = new TGraph(nBins, Espectrum, EffCyl);
  effIP = new TGraph(nBins, Espectrum, EffIP);
  effBack = new TGraph(nBins, Espectrum, EffBack);
  effCylBack = new TGraph(nBins, Espectrum, EffCylBack);
  effIPBack = new TGraph(nBins, Espectrum, EffIPBack);

  effBack->SetMarkerColor(2);
  effBack->SetMarkerStyle(20);

  effCyl->SetMarkerColor(3);
  effCyl->SetMarkerStyle(20);

  effIP->SetMarkerColor(3);
  effIP->SetMarkerStyle(22);

  effCylBack->SetMarkerColor(4);
  effCylBack->SetMarkerStyle(20);

  effIPBack->SetMarkerColor(4);
  effIPBack->SetMarkerStyle(22);

  effCyl->GetXaxis()->SetTitle("Energy [ GeV ]");
  effCyl->GetYaxis()->SetTitle("Cut Efficiency [ Frac. ]");
  effCyl->GetYaxis()->SetRangeUser(0,1);
  effCyl->SetTitle("Cut Efficiency on Signal Showers");
  effCyl->Draw();

  effIP->Draw("PL");
  effCylBack->Draw("PL");
  effIPBack->Draw("PL");
  effBack->Draw("PL");
  
  leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetFillColor(kWhite);
  leg->AddEntry(effCyl,Form("Cylinder - %d cm",CylDist));
  leg->AddEntry(effIP,Form("IP - %d cm",IPDist));
  leg->AddEntry(effBack,Form("Back Dist - %d cm",BackDist));
  leg->AddEntry(effCylBack,Form("Cyl - %d cm & Back - %d cm",CylDist,BackDist));
  leg->AddEntry(effIPBack,Form("IP - %d cm & Back - %d cm",IPDist,BackDist));
  leg->Draw();

  cCuts->SaveAs("effSignalOneFrameBigFV.pdf");

  cCuts->Clear();

  // scale histograms to 6.6E20 POT
  hEAll->Scale(CCInclusive/Nevents);
  hECyl->Scale(CCInclusive/Nevents);
  hEIP->Scale(CCInclusive/Nevents);
  hEBack->Scale(CCInclusive/Nevents);
  hECylBack->Scale(CCInclusive/Nevents);
  hEIPBack->Scale(CCInclusive/Nevents);

  gStyle->SetOptStat(0);
  hEAll->SetTitle("Signal Events Energy Spectrum [6.6E20 POT] - Cut Performance");
  hEAll->GetXaxis()->SetTitle("Energy [MeV]");
  hEAll->GetYaxis()->SetTitle("Events");
  //hEAll->GetYaxis()->SetRangeUser(10,3E6);
  hEAll->SetLineWidth(2);
  //  hEAll->SetLineStyle(2);
  hEAll->SetLineColor(4);

  hECylBack->SetLineWidth(2);
  //  hECylBack->SetLineStyle(2);
  hECylBack->SetLineColor(3);

  hECyl->SetLineWidth(2);
  //  hECyl->SetLineStyle(2);
  hECyl->SetLineColor(1);

  hEBack->SetLineWidth(2);
  //  hEBack->SetLineStyle(2);
  hEBack->SetLineColor(2);

  hEAll->Draw();
  hECylBack->Draw("same");
  hECyl->Draw("same");
  hEBack->Draw("same");

  leg2 = new TLegend(0.5,0.7,0.9,0.9);
  leg2->SetFillColor(kWhite);
  leg2->AddEntry(hEAll,"All EM Showers","LF");
  leg2->AddEntry(hECyl,Form("After %d cm Cylinder Cut",CylDist));
  leg2->AddEntry(hEBack,Form("After %d cm Back-To-Wall Cut",BackDist));
  leg2->AddEntry(hECylBack,Form("After %d cm Cyl & %d cm Back Cut",CylDist,BackDist));
  leg2->Draw();

  //  cCuts->SetLogy();
  cCuts->Modified();
  cCuts->Update();
  cCuts->SaveAs("histoPerformanceSignalOneFrameBigFV.pdf");
  return;
}
