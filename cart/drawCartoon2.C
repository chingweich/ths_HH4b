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
#include <TLine.h>
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
void setTH(TH2D * h_data){
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


void drawCartoon2(){
	
	
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

TH2D * th2=new TH2D("","",2,0,2,2,0,2);
TLine* tl1[3];
tl1[0]=new TLine(0,1,2,1);
tl1[1]=new TLine(0.6,0,0.6,2);
tl1[2]=new TLine(1.3,0,1.3,2);
setTH(th2);

th2->Draw();
th2->SetXTitle("M_{leading AK8}");
th2->SetYTitle("double-b tagger_{leading AK8}");
tl1[0]->Draw("same");
tl1[1]->Draw("same");
tl1[2]->Draw("same");
TLatex *lar = new TLatex();

   lar->SetTextSize(0.03);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
lar->DrawLatex(0.17, 0.3, "mass side band A");
lar->DrawLatex(0.17, 0.7, "mass side band B");
lar->DrawLatex(0.7, 0.3, "mass side band C");
lar->DrawLatex(0.7, 0.7, "mass side band D");
lar->DrawLatex(0.5, 0.7, "signal");
lar->DrawLatex(0.48, 0.3, "anti-tag");
	c->Print("cart2.pdf");
}