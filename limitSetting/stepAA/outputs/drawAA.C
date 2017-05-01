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
void setTH(TH1D * h_data){
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

void setRatio(TH1D* h_ratio){
	h_ratio->GetYaxis()->SetTitle("#frac{Data-Fit}{#sigma}");
	h_ratio->GetXaxis()->SetTitle("M^{reduced}_{jj} [GeV]");
  h_ratio->GetYaxis()->SetTitleOffset(0.45);
  h_ratio->GetXaxis()->SetLabelSize(0.125);
  h_ratio->GetXaxis()->SetLabelOffset(0.005);
  h_ratio->GetXaxis()->SetTitleSize(0.08);
  h_ratio->GetXaxis()->SetTitleOffset(1.2);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.12);
  h_ratio->GetYaxis()->SetNdivisions(505);
  //h_ratio->GetYaxis()->SetRangeUser(0,2);
}


void drawAlphaBase(string st){


gStyle->SetPadRightMargin(0.039);
gStyle->SetPadLeftMargin(0.129);
gStyle->SetNdivisions(605, "XYZ");
gStyle->SetOptStat(0);
TCanvas* c=new TCanvas("c","",0,0,600,600);
 c->Divide(1,2);
Float_t up_height     = 0.8;
  Float_t dw_correction = 1.375;
  Float_t dw_height= (1-up_height)*dw_correction;

  //TPad* c_up = (TPad*) c.GetListOfPrimitives()->FindObject("c_1");
  TPad* c_up = (TPad*) c->cd(1);
  //c_up->SetLogy(1);
  TPad* c_dw = (TPad*) c->GetListOfPrimitives()->FindObject("c_2"); 

  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);
  
  c_up->SetLogy(1);
gStyle->SetPadRightMargin(0.029);
gStyle->SetPadLeftMargin(0.1509);
c->SetLogy(1);
TFile * tf1= TFile::Open(Form("HHSR_%s_Data.root",st.data()));
TH1D* th2[2];
th2[0]=(TH1D*)tf1->Get("data");
th2[1]=(TH1D*)tf1->Get("antitag");
th2[0]->Rebin(50);
th2[1]->Rebin(50);
TF1* tf4[2];
//(1 + @0*@3)*exp(-@0*@2/(1+(@0*@1*@2)))
//tf4[0]=new TF1("L1", "(1+x*[0])",-1000,5000);
tf4[0]=new TF1("L1", "(1+x*[2])*[3]*exp(-x*[1]/(1+(x*[0]*[1])))",1100,3000);
tf4[1]=new TF1("L2", "(1+x*[2])*[3]*exp(-x*[1]/(1+(x*[0]*[1])))",1100,3000);
tf4[0]->SetLineColor(2);
tf4[1]->SetLineColor(3);

if(st.find("TT")!= std::string::npos){
	tf4[0]->SetParameter(0,0.0250418);
tf4[0]->SetParameter(1,0.0445427);
tf4[0]->SetParameter(2,0.0052978);
tf4[1]->SetParameter(0,0.0306015);
tf4[1]->SetParameter(1,0.0182569);
tf4[1]->SetParameter(2,0.0052978);
	
}
else {
	
	tf4[0]->SetParameter(0,0.0284413);
tf4[0]->SetParameter(1,0.0256102);
tf4[0]->SetParameter(2,0.0052978);
tf4[1]->SetParameter(0,0.0298375);
tf4[1]->SetParameter(1,0.0199731);
tf4[1]->SetParameter(2,0.00312665);
}

tf4[0]->SetParameter(3,1);
tf4[0]->SetParameter(3,th2[0]->Integral()*50/tf4[0]->Integral(1100,2800));

tf4[1]->SetParameter(3,1);
tf4[1]->SetParameter(3,th2[0]->Integral()*50/tf4[1]->Integral(1100,2800));
setTH(th2[0]);


th2[0]->GetXaxis()->SetRangeUser(1100,2800);
th2[0]->SetMinimum(0.5);

if(st.find("LL")!= std::string::npos)th2[1]->SetMaximum(1e4);
else th2[0]->SetMaximum(100);
th2[0]->SetYTitle("Events / 50 GeV");
th2[1]->SetYTitle("Events / 50 GeV");
th2[0]->Draw("e");
th2[0]->SetMarkerStyle(20);
th2[0]->SetLineColor(1);
th2[1]->SetLineColor(1);
th2[1]->SetMarkerStyle(20);
tf4[0]->Draw("same");
tf4[1]->Draw("same");

TLegend *leg = new TLegend(0.7, 0.68, 0.9, 0.87);
leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(th2[0]," data");
  leg->AddEntry(tf4[0]," pre-fit");
  leg->AddEntry(tf4[1]," post-fit");
  tf4[0]->SetFillColor(0);
  tf4[1]->SetFillColor(0);
leg->Draw("same");
TLatex *lar = new TLatex();

  lar->SetTextSize(0.044);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.66, 0.94, " 35.9 fb^{-1} (13 TeV)");
  
  lar->SetTextSize(0.09);
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

		for(int i=1;i<th2[0]->GetNbinsX();i++){
	if(th2[0]->GetBinContent(i)<1){
		th2[0]->SetBinContent(i,0);
		th2[0]->SetBinContent(i,0);
	}
}
	
		
c_dw->cd();
TH1D* th3[2];
th3[0]=(TH1D*)th2[0]->Clone("a");
for(int i=1;i<th2[0]->GetNbinsX();i++){
	if(th2[0]->GetBinContent(i))th3[0]->SetBinContent(i,(th2[0]->GetBinContent(i)-tf4[1]->Eval(th2[0]->GetBinCenter(i)))/th2[0]->GetBinError(i));
	else th3[0]->SetBinContent(i,0);
	th3[0]->SetBinError(i,0.000001);
	//else th3[0]->SetBinContent(i,0);
}
th3[0]->SetMaximum(2);
th3[0]->SetMinimum(-2);
setRatio(th3[0]);
th3[0]->Draw("e");



c->Print(Form("aa/data_%s.pdf",st.data()));
c_up->cd();
TF1* tf5[2];
//exp(-@0*@2/(1+(@0*@1*@2)))
tf5[0]=new TF1("LinearFit_","[2]*exp(-x*[1]/(1+(x*[0]*[1])))",1100,2800);
tf5[1]=new TF1("LinearFit_","[2]*exp(-x*[1]/(1+(x*[0]*[1])))",1100,2800);
tf5[0]->SetLineColor(2);
tf5[1]->SetLineColor(3);
//tf4[0]->SetParameter(0,0.0250418);
setTH(th2[1]);
if(st.find("TT")!= std::string::npos){
	tf5[0]->SetParameter(0,0.0250418);
tf5[0]->SetParameter(1,0.0445427);
//tf4[0]->SetParameter(2,0.0052978);
tf5[1]->SetParameter(0,0.0306015);
tf5[1]->SetParameter(1,0.0182569);
//tf4[1]->SetParameter(2,0.0052978);
}
else {
	tf5[0]->SetParameter(0,0.0284413);
tf5[0]->SetParameter(1,0.0256102);
tf5[1]->SetParameter(0,0.0298375);
tf5[1]->SetParameter(1,0.0199731);
}

tf5[0]->SetParameter(2,1);
tf5[0]->SetParameter(2,th2[1]->Integral()*50/tf5[0]->Integral(1100,2800));

tf5[1]->SetParameter(2,1);
tf5[1]->SetParameter(2,th2[1]->Integral()*50/tf5[1]->Integral(1100,2800));

//th2[0]->Rebin(50);
th2[1]->GetXaxis()->SetRangeUser(1100,2800);
th2[1]->SetMinimum(0.5);
if(st.find("LL")!= std::string::npos)th2[1]->SetMaximum(1e5);
else th2[1]->SetMaximum(1e3);

th2[1]->Draw();
tf5[0]->Draw("same");
tf5[1]->Draw("same");
leg->Draw("same");
for(int i=1;i<th2[1]->GetNbinsX();i++){
	if(th2[1]->GetBinContent(i)<1){
		th2[1]->SetBinContent(i,0);
		th2[1]->SetBinContent(i,0);
	}
}
 lar->SetTextSize(0.044);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.66, 0.94, " 35.9 fb^{-1} (13 TeV)");
  
  lar->SetTextSize(0.09);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(62);
		lar2->SetNDC(kTRUE);
		lar2->SetTextSize(0.04);
		lar2->SetLineWidth(5);
		lar2->SetTextAlign(12);
		lar->DrawLatex(0.17, 0.85, "CMS");
		//lar2->DrawLatex(0.31, 0.77, "#it{#bf{Simulation}} ");
		lar2->DrawLatex(0.31, 0.83, "#it{#bf{Preliminary}} ");
c_dw->cd();
th3[1]=(TH1D*)th2[1]->Clone("a");
for(int i=1;i<th2[1]->GetNbinsX();i++){
	if(th2[1]->GetBinContent(i))th3[1]->SetBinContent(i,(th2[1]->GetBinContent(i)-tf5[1]->Eval(th2[1]->GetBinCenter(i)))/th2[1]->GetBinError(i));
	else th3[1]->SetBinContent(i,0);
	th3[1]->SetBinError(i,0.000001);
}
th3[1]->SetMaximum(2);
th3[1]->SetMinimum(-2);
setRatio(th3[1]);
th3[1]->Draw("e");

c->Print(Form("aa/anti_%s.pdf",st.data()));
}
void drawAA(){
	drawAlphaBase("LL");
	drawAlphaBase("TT");
}