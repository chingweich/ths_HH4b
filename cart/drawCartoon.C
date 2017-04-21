#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <TLatex.h>
#include <string>
#include <sstream>
#include "TF1.h"
void setTH(TH1D * h_data){
	h_data->SetTitle("");
	h_data->SetLineWidth(2);
	//h_data->SetMaximum(3000);
 h_data->SetMarkerSize(1.6);

 h_data->GetYaxis()->SetTitleOffset(1.7);
 h_data->GetZaxis()->SetTitleOffset(0.75);
 // size of axis labels
 h_data->GetXaxis()->SetTitleSize(0.04);
 h_data->GetYaxis()->SetTitleSize(0.04);
 h_data->GetZaxis()->SetTitleSize(0.035);
 h_data->GetXaxis()->SetLabelSize(0.04);
 h_data->GetYaxis()->SetLabelSize(0.04); 
 h_data->GetZaxis()->SetLabelSize(0.025);
 
//h_data->SetPadRightMargin(0.029);
//h_data->SetPadLeftMargin(0.1509);
h_data->SetNdivisions(000, "XYZ");
}


void drawCartoon(){
	
	
	 //TStyle* ts=new TStyle();
	 /*
	  ts
 ts->SetPaintTextFormat("2.1f");
 //ts->SetPalette(57);
 ts->SetFrameLineWidth(3);
ts->SetPadRightMargin(0.029);
ts->SetPadLeftMargin(0.1509);
ts->SetNdivisions(605, "XYZ");
ts->cd();

*/

gStyle->SetPadRightMargin(0.039);
gStyle->SetPadLeftMargin(0.1509);
gStyle->SetNdivisions(605, "XYZ");
gStyle->SetOptStat(0);
 TCanvas* c=new TCanvas("c","",0,0,600,600);
 
gStyle->SetPadRightMargin(0.029);
gStyle->SetPadLeftMargin(0.1509);
//cout<<gStyle->GetPadLeftMargin();
// gStyle->SetPadLeftMargin(1.559);
 //c->Update();
 
 TFile* tf2[10];
   tf2[0]=TFile::Open("../N1/1.root");
   
    string hadflv[4]={"bb","b","cc","udcsg"};
	
	
 //h_name.push_back("totalMass");
 //h_name.push_back("totalMassRed");
 
 TH1D * th1[2];
 th1[0]=(TH1D*) tf2[0]->Get("totalMass_bb");
 for(int i=0;i<3;i++){
	 TH1D * th2=(TH1D*) tf2[0]->Get(Form("totalMass_%s",hadflv[i+1].data()));
	 th2->Sumw2();
	 th1[0]->Add(th2);
 }
 th1[1]=(TH1D*) tf2[0]->Get("totalMassRed_bb");
 for(int i=0;i<3;i++){
	 TH1D * th2=(TH1D*) tf2[0]->Get(Form("totalMassRed_%s",hadflv[i+1].data()));
	 th2->Sumw2();
	 th1[1]->Add(th2);
 }
 th1[0]->GetXaxis()->SetRangeUser(1000,1600);
 th1[0]->SetMaximum(th1[0]->GetMaximum()*1.3);
 th1[0]->SetYTitle("Events");
 th1[0]->SetXTitle("M_{jj} [GeV]");
 setTH(th1[0]);
 setTH(th1[1]);
 th1[0]->Draw("hist");
 th1[1]->SetLineColor(2);
 th1[1]->Draw("samehist");
 
 TLegend *leg = new TLegend(0.2, 0.58, 0.3, 0.77);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(th1[0],"M_{jj}");
  TF1* tf1[2];
  tf1[0]=new TF1("fa1","gaus()",th1[0]->GetMaximumBin()*0.7,th1[0]->GetMaximumBin()*1.3);
  tf1[1]=new TF1("fa2","gaus()",th1[0]->GetMaximumBin()*0.7,th1[0]->GetMaximumBin()*1.3);
  
 leg->Draw("same");
 TLatex *lar = new TLatex();
 
 TF1* tf4;
 tf4=new TF1("fa3","[0]+x*[1]",1000,1800);
 tf4->SetParameter(0,1000);
 tf4->SetParameter(1,-0.2);
 
 TH1D* bkg= (TH1D*)th1[0]->Clone();
 bkg->Clear();
 for(int i=0;i<bkg->GetNbinsX();i++){
	 TRandom *eventGenerator = new TRandom();
	 bkg->SetBinContent(i,eventGenerator->Poisson(tf4->Eval(bkg->GetBinCenter(i))));
	 bkg->SetBinError(i,0.001);
	 th1[0]->SetBinError(i,0);
 }
 bkg->Draw("same");
 
// tf4->SetParameter(1,1000+0.3*1400);
tf4->Draw("same");
  lar->SetTextSize(0.044);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
 // lar->DrawLatex(0.7, 0.94, " 35.9 fb^{-1} (13 TeV)");
  
  lar->SetTextSize(0.070);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(62);

		
		//lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
		
		
		TLatex *lar2 = new TLatex();

		lar2->SetNDC(kTRUE);
		lar2->SetTextSize(0.04);
		lar2->SetLineWidth(5);
		lar2->SetTextAlign(12);
		lar->DrawLatex(0.2, 0.86, "CMS");
		lar2->DrawLatex(0.34, 0.83, "#it{#bf{Simulation}} ");
		th1[0]->Scale(0.03);
		bkg->Add(th1[0]);
		bkg->Draw("e");
		
  bkg->SetLineColor(kBlack);
  bkg->SetMarkerStyle(8);
  bkg->SetMarkerSize(1);
  bkg->Fit(tf4,"","",1000,1200);
  tf4->Draw("same");
	c->Print("cart.pdf");
}