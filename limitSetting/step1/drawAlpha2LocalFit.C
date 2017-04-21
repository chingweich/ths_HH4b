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


void drawAlphaBase(string st){


gStyle->SetPadRightMargin(0.039);
gStyle->SetPadLeftMargin(0.129);
gStyle->SetNdivisions(605, "XYZ");
gStyle->SetOptStat(0);
TCanvas* c=new TCanvas("c","",0,0,600,600);
gStyle->SetPadRightMargin(0.029);
gStyle->SetPadLeftMargin(0.1509);

TFile * tf1= TFile::Open(Form("th_%s.root",st.data()));
TGraph* th2=(TGraph*)tf1->Get("ths_1D");

th2->GetXaxis()->SetTitle("M_{leading AK8} -125 [GeV]");
th2->GetXaxis()->SetRangeUser(-81,81);
if(st.find("LL")!= std::string::npos)th2->GetYaxis()->SetRangeUser(0.18,0.3);
else th2->GetYaxis()->SetRangeUser(0.04,0.12);
th2->GetYaxis()->SetTitle("R_{p/f}");
th2->Draw("AP");
TF1* tf4=new TF1("LinearFit_", "[0]+ [1]*x + [2]*x^2 ",-81,81);
TF1* err[2];
err[0]=new TF1("errUp", "[0]+ [1]*x + [2]*x*x + sqrt(([3]*[3]) + (2*x*[6]) + (x*x*[4]*[4]) + (2*x*x*[7]) + (2*x*x*x*[8]) + (x*x*x*x*[5]*[5])) ",-81,81);
err[1]=new TF1("errUp", "[0]+ [1]*x + [2]*x*x - sqrt(([3]*[3]) + (2*x*[6]) + (x*x*[4]*[4]) + (2*x*x*[7]) + (2*x*x*x*[8]) + (x*x*x*x*[5]*[5]))",-81,81);
tf4->SetParameter(0,0.2);
//tf4->SetParameter(1,-2.84377e-04);
tf4->SetParameter(2,7.83062e-06);
TFitResultPtr r =th2->Fit(tf4,"S","",-81,81);
TMatrixDSym cov = r->GetCovarianceMatrix();
cov.Print();
setTH(th2);
for(int i=0;i<2;i++){
	err[i]->SetParameter(0,tf4->GetParameter(0));
	err[i]->SetParameter(1,tf4->GetParameter(1));
	err[i]->SetParameter(2,tf4->GetParameter(2));
	err[i]->SetParameter(3,tf4->GetParErrors()[0]);
	err[i]->SetParameter(4,tf4->GetParErrors()[1]);
	err[i]->SetParameter(5,tf4->GetParErrors()[2]);
	err[i]->SetParameter(6,cov(0,1));
	err[i]->SetParameter(7,cov(0,2));
	err[i]->SetParameter(8,cov(1,2));
	err[i]->SetLineStyle(7);
	err[i]->SetLineColor(2);
	err[i]->Draw("same");
}
tf4->Draw("same");
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

c->Print(Form("al/1d_%s.pdf",st.data()));
}
void drawAlpha2(){
	drawAlphaBase("LL");
	drawAlphaBase("TT");
}