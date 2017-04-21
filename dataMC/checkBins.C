#include <TLegend.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TObject.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>




void setTGraph(TGraph* tg1,int i,bool setMax=0){
	
	tg1->SetLineColor(i+1);
	//tg1->SetLineColor(98-4*i);
	tg1->SetMarkerColor(i+1);
	if(i==4)tg1->SetLineColor(kOrange);
	if(i==4)tg1->SetMarkerColor(kOrange);
	//tg1->SetMarkerColor(98-4*i);
	tg1->SetMarkerStyle(20);
	
	tg1->GetXaxis()->SetTitle("m_{jj}[GeV]");
	tg1->GetYaxis()->SetTitle("");
	//limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
	tg1->SetTitle("");
	tg1->SetLineWidth(2);
	tg1->SetFillColor(0);
	//tg1->SetMaximum(1);
	//tg1->SetMinimum(0.0014);
	tg1->SetMarkerSize(1);
	//tg1->SetMarkerStyle(20);
	//tg1->SetMinimum(0.05);
	tg1->GetYaxis()->SetTitleOffset(1.5);
	//tg1->GetZaxis()->SetTitleOffset(0.65);
	// size of axis labels
	tg1->GetXaxis()->SetTitleSize(0.04);
	tg1->GetYaxis()->SetTitleSize(0.04);
	//limits->GetZaxis()->SetTitleSize(0.035);
	tg1->GetXaxis()->SetLabelSize(0.05);
	tg1->GetYaxis()->SetLabelSize(0.03); 
}

void setLeg(TLegend* leg){
	leg->SetBorderSize(1);
	//leg->SetLineColor(0);                                                                                                                     
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
}

void checkBins(){
	//TH1D* th1[3];
	//th1[0]=massPlotBase("BulkGrav_M-2000_0.root","jet1_puppi_msoftdrop");
	//th1[1]=massPlotBase("BulkGrav_M-2000_0.root","jet1_puppi_msoftdrop_TheaCorr");
	//th1[2]=massPlotBase("HCorr.root","jet1_puppi_msoftdrop_TheaCorr");
	
	TCanvas* c1;
	//setNCUStyle(true);
	c1 = new TCanvas("c1","",600,600);
	gStyle->SetNdivisions(605, "XYZ");
	TLegend *leg = new TLegend(0.50, 0.73, 0.85, 0.88);
  
	setLeg(leg);
	TFile* tf1[2];
	tf1[0]=TFile::Open("data.root");
	tf1[1]=TFile::Open("Mar7CheckRebin.root");
	
	TH1D* th1[2];
	th1[0]=(TH1D*)tf1[0]->Get("pt_j0_udcsg_antiDBT");
	th1[1]=(TH1D*)tf1[1]->Get("pt_j0_udcsg_antiDBT");
	th1[1]->Rebin(4);
	cout<<th1[0]->Integral()<<","<<th1[1]->Integral()<<endl;
	th1[0]->Draw("");
	th1[1]->SetLineColor(2);
	th1[1]->Draw("same");
	c1->Print("checkBins.pdf");
}