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
#include <TMatrixDSym.h>
#include <TSystem.h>
#include <TFitResult.h>
#include <TLatex.h>
#include <string>
#include <sstream>
#include "TF1.h"
void setTH(TGraph * h_data){
	h_data->SetTitle("");
	h_data->SetLineWidth(2);
	//h_data->SetMaximum(3000);
 h_data->SetMarkerSize(1.6);

 h_data->GetYaxis()->SetTitleOffset(1.5);
 //h_data->GetZaxis()->SetTitleOffset(0.75);
 // size of axis labels
 h_data->GetXaxis()->SetTitleSize(0.04);
 h_data->GetYaxis()->SetTitleSize(0.04);
 //h_data->GetZaxis()->SetTitleSize(0.035);
 h_data->GetXaxis()->SetLabelSize(0.04);
 h_data->GetYaxis()->SetLabelSize(0.04); 
 //h_data->GetZaxis()->SetLabelSize(0.025);
 
//h_data->SetPadRightMargin(0.029);
//h_data->SetPadLeftMargin(0.1509);
h_data->GetXaxis()->SetNdivisions(605, "XYZ");
h_data->GetYaxis()->SetNdivisions(605, "XYZ");
}


void drawUncertBase(string st){


gStyle->SetPadRightMargin(0.039);
gStyle->SetPadLeftMargin(0.129);
gStyle->SetNdivisions(605, "XYZ");
gStyle->SetOptStat(0);
TCanvas* c=new TCanvas("c","",0,0,600,600);
gStyle->SetPadRightMargin(0.029);
gStyle->SetPadLeftMargin(0.1509);

TFile * tf1[3];
tf1[0]= TFile::Open("central_noUncert/limit.root");
tf1[1]= TFile::Open(Form("%sUp_noUncert/limit.root",st.data()));
tf1[2]= TFile::Open(Form("%sDn_noUncert/limit.root",st.data()));
TGraphAsymmErrors* tg1[3];
tg1[0]=(TGraphAsymmErrors*)tf1[0]->Get("LimitExpectedCLs");
tg1[1]=(TGraphAsymmErrors*)tf1[1]->Get("LimitExpectedCLs");
tg1[2]=(TGraphAsymmErrors*)tf1[2]->Get("LimitExpectedCLs");


for(int i=0;i<tg1[0]->GetN();i++){
	double x1,y1,x2,y2,x3,y3;
	tg1[0]->GetPoint(i,x1,y1);
	tg1[1]->GetPoint(i,x2,y2);
	tg1[2]->GetPoint(i,x3,y3);
	
	tg1[1]->SetPoint(i,x2,y2/y1);
	tg1[2]->SetPoint(i,x3,y3/y1);
}

//tg1[0]->Draw("APL");
setTH(tg1[1]);
tg1[1]->SetMaximum(1.4);
tg1[1]->SetMinimum(0.6);
tg1[1]->SetLineStyle(1);
tg1[2]->SetLineStyle(1);
tg1[1]->SetMarkerStyle(20);
tg1[2]->SetMarkerStyle(21);
tg1[1]->SetLineColor(1);
tg1[2]->SetLineColor(2);
tg1[1]->SetMarkerColor(1);
tg1[2]->SetMarkerColor(2);
tg1[1]->GetXaxis()->SetTitle("M_{X} [GeV]");
tg1[1]->GetYaxis()->SetTitle("95% upper limit on #sigma: ratio _{scaled / central}");
tg1[1]->Draw("PLA");
tg1[2]->Draw("PLsame");
TLegend *leg = new TLegend(0.7, 0.68, 0.88, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  //leg->AddEntry(tg1[0]," cerntral","pl");
  leg->AddEntry(tg1[1]," up","pl");
  leg->AddEntry(tg1[2]," down","pl");
  
leg->Draw("same");

TLatex *lar = new TLatex();

  lar->SetTextSize(0.044);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.66, 0.94, " 35.9 fb^{-1} (13 TeV)");
  
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
		lar->DrawLatex(0.17, 0.85, "CMS");
		//lar2->DrawLatex(0.31, 0.77, "#it{#bf{Simulation}} ");
		lar2->DrawLatex(0.31, 0.83, "#it{#bf{Preliminary}} ");

	c->Print(Form("plots_uncert/%s.pdf",st.data()));
}
void drawUncertBaseLL(string st){


gStyle->SetPadRightMargin(0.039);
gStyle->SetPadLeftMargin(0.129);
gStyle->SetNdivisions(605, "XYZ");
gStyle->SetOptStat(0);
TCanvas* c=new TCanvas("c","",0,0,600,600);
gStyle->SetPadRightMargin(0.029);
gStyle->SetPadLeftMargin(0.1509);

TFile * tf1[3];
tf1[0]= TFile::Open("central_noUncert/LL/limit.root");
tf1[1]= TFile::Open(Form("%sUp_noUncert/LL/limit.root",st.data()));
tf1[2]= TFile::Open(Form("%sDn_noUncert/LL/limit.root",st.data()));
TGraphAsymmErrors* tg1[3];
tg1[0]=(TGraphAsymmErrors*)tf1[0]->Get("LimitExpectedCLs");
tg1[1]=(TGraphAsymmErrors*)tf1[1]->Get("LimitExpectedCLs");
tg1[2]=(TGraphAsymmErrors*)tf1[2]->Get("LimitExpectedCLs");


for(int i=0;i<tg1[0]->GetN();i++){
	double x1,y1,x2,y2,x3,y3;
	tg1[0]->GetPoint(i,x1,y1);
	tg1[1]->GetPoint(i,x2,y2);
	tg1[2]->GetPoint(i,x3,y3);
	
	tg1[1]->SetPoint(i,x2,y2/y1);
	tg1[2]->SetPoint(i,x3,y3/y1);
}

//tg1[0]->Draw("APL");
setTH(tg1[1]);
tg1[1]->SetMaximum(1.4);
tg1[1]->SetMinimum(0.6);
tg1[1]->SetLineStyle(1);
tg1[2]->SetLineStyle(1);
tg1[1]->SetMarkerStyle(20);
tg1[2]->SetMarkerStyle(21);
tg1[1]->SetLineColor(1);
tg1[2]->SetLineColor(2);
tg1[1]->SetMarkerColor(1);
tg1[2]->SetMarkerColor(2);
tg1[1]->GetXaxis()->SetTitle("M_{X} [GeV]");
tg1[1]->GetYaxis()->SetTitle("95% upper limit on #sigma: ratio _{scaled / central}");
tg1[1]->Draw("PLA");
tg1[2]->Draw("PLsame");
TLegend *leg = new TLegend(0.7, 0.68, 0.88, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  //leg->AddEntry(tg1[0]," cerntral","pl");
  leg->AddEntry(tg1[1]," up","pl");
  leg->AddEntry(tg1[2]," down","pl");
  
leg->Draw("same");

TLatex *lar = new TLatex();

  lar->SetTextSize(0.044);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.66, 0.94, " 35.9 fb^{-1} (13 TeV)");
  
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
		lar->DrawLatex(0.17, 0.85, "CMS");
		//lar2->DrawLatex(0.31, 0.77, "#it{#bf{Simulation}} ");
		lar2->DrawLatex(0.31, 0.83, "#it{#bf{Preliminary}} ");

	c->Print(Form("plots_uncert/%s_LL.pdf",st.data()));
}
void drawUncertBaseTT(string st){


gStyle->SetPadRightMargin(0.039);
gStyle->SetPadLeftMargin(0.129);
gStyle->SetNdivisions(605, "XYZ");
gStyle->SetOptStat(0);
TCanvas* c=new TCanvas("c","",0,0,600,600);
gStyle->SetPadRightMargin(0.029);
gStyle->SetPadLeftMargin(0.1509);

TFile * tf1[3];
tf1[0]= TFile::Open("central_noUncert/TT/limit.root");
tf1[1]= TFile::Open(Form("%sUp_noUncert/TT/limit.root",st.data()));
tf1[2]= TFile::Open(Form("%sDn_noUncert/TT/limit.root",st.data()));
TGraphAsymmErrors* tg1[3];
tg1[0]=(TGraphAsymmErrors*)tf1[0]->Get("LimitExpectedCLs");
tg1[1]=(TGraphAsymmErrors*)tf1[1]->Get("LimitExpectedCLs");
tg1[2]=(TGraphAsymmErrors*)tf1[2]->Get("LimitExpectedCLs");


for(int i=0;i<tg1[0]->GetN();i++){
	double x1,y1,x2,y2,x3,y3;
	tg1[0]->GetPoint(i,x1,y1);
	tg1[1]->GetPoint(i,x2,y2);
	tg1[2]->GetPoint(i,x3,y3);
	
	tg1[1]->SetPoint(i,x2,y2/y1);
	tg1[2]->SetPoint(i,x3,y3/y1);
}

//tg1[0]->Draw("APL");
setTH(tg1[1]);
tg1[1]->SetMaximum(1.4);
tg1[1]->SetMinimum(0.6);
tg1[1]->SetLineStyle(1);
tg1[2]->SetLineStyle(1);
tg1[1]->SetMarkerStyle(20);
tg1[2]->SetMarkerStyle(21);
tg1[1]->SetLineColor(1);
tg1[2]->SetLineColor(2);
tg1[1]->SetMarkerColor(1);
tg1[2]->SetMarkerColor(2);
tg1[1]->GetXaxis()->SetTitle("M_{X} [GeV]");
tg1[1]->GetYaxis()->SetTitle("95% upper limit on #sigma: ratio _{scaled / central}");
tg1[1]->Draw("PLA");
tg1[2]->Draw("PLsame");
TLegend *leg = new TLegend(0.7, 0.68, 0.88, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  //leg->AddEntry(tg1[0]," cerntral","pl");
  leg->AddEntry(tg1[1]," up","pl");
  leg->AddEntry(tg1[2]," down","pl");
  
leg->Draw("same");

TLatex *lar = new TLatex();

  lar->SetTextSize(0.044);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.66, 0.94, " 35.9 fb^{-1} (13 TeV)");
  
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
		lar->DrawLatex(0.17, 0.85, "CMS");
		//lar2->DrawLatex(0.31, 0.77, "#it{#bf{Simulation}} ");
		lar2->DrawLatex(0.31, 0.83, "#it{#bf{Preliminary}} ");

	c->Print(Form("plots_uncert/%s_TT.pdf",st.data()));
}
void drawUncert_pu(){
	drawUncertBase("pu");
	drawUncertBaseLL("pu");
	drawUncertBaseTT("pu");
}