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


void drawAlphaBase(string st){


gStyle->SetPadRightMargin(0.039);
gStyle->SetPadLeftMargin(0.129);
gStyle->SetNdivisions(605, "XYZ");
gStyle->SetOptStat(0);
TCanvas* c=new TCanvas("c","",0,0,600,600);
gStyle->SetPadRightMargin(0.029);
gStyle->SetPadLeftMargin(0.1509);

TFile * tf1= TFile::Open(Form("th_%s.root",st.data()));
TH1D* th2=(TH1D*)tf1->Get("ths_N");
TH1D* th4=(TH1D*)tf1->Get("ths_A");
TH1D* th3=(TH1D*)tf1->Get("ths_SIG0");
setTH(th2);

th2->GetXaxis()->SetRangeUser(600,1800);
th2->SetXTitle("M^{reduced}_{jj} [GeV]");
th2->SetYTitle(Form("Events / %2.0f GeV",th2->GetBinWidth(1)));
th2->Draw();
th2->SetMaximum(th2->GetMaximum()*0.6);
th3->Scale(10);
//th3->Draw("same");
//th4->Draw("same");
TLegend *leg = new TLegend(0.5, 0.68, 0.8, 0.87);
leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(th2," Predicted background");
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

c->Print(Form("al/mjj_%s.pdf",st.data()));
}
void drawAlpha3(){
	drawAlphaBase("LL");
	drawAlphaBase("TT");
}