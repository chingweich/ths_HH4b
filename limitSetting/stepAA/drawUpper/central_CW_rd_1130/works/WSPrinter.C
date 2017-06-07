#include <cstring>
#include <cerrno>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <unistd.h>
#include <errno.h>
#include <iomanip>
// ROOT headers
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"
#include "TMath.h"

#include "Math/Math.h"
#include "Math/QuantFuncMathCore.h"
//#include "SpecFuncCephes.h"

#include "TLatex.h"
#include "TString.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

// RooFit headers
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TStyle.h"

// RooStats headers
#include "RooStats/HLFactory.h"

#include "RooAbsPdf.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsData.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooProduct.h"
#include "RooExtendPdf.h"
#include "RooBernstein.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooCmdArg.h"
#include "RooConstVar.h"
#include "RooRealVar.h"

using namespace RooFit;
using namespace RooStats ;
using namespace std ;
using namespace ROOT::Math ;

void setLeg(TLegend* leg){
	leg->SetBorderSize(1);
	//leg->SetLineColor(0);                                                                                                                     
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
}

void WSPrinter(){
	
	string DBT[2]={"TT","LL"};
	double rss0[2]={0};
	for(int j=0;j<2;j++){
		TFile * tf1 = TFile::Open(Form("w_data_%s.root",DBT[j].data()));
		RooWorkspace *ws1 = (RooWorkspace*)tf1->Get("HH4b");
		ws1->Print();
		RooDataSet* data=(RooDataSet*)ws1->data("data_obs");
		
		TFile * tf2 = TFile::Open(Form("w_data_%s_m.root",DBT[j].data()));
		RooWorkspace *ws2 = (RooWorkspace*)tf2->Get("HH4b");
		RooDataSet* data2=(RooDataSet*)ws2->data("data_obs");
		data2->SetName("data_obs2");
		ws2->Print();
		//RooAbsPdf* bkg=ws2->pdf("bgSB_");
		//RooAbsPdf* bkg2=ws2->pdf("bgSB1c_");
		RooRealVar *mass =  ws1->var("x");
		mass->setBins(10) ;
		RooPlot *plot = mass->frame(50);
		
		data->plotOn(plot,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed),Binning(50));
		data2->plotOn(plot,RooFit::LineColor(kBlue),RooFit::MarkerColor(kBlue),Binning(50));
		//bkg->paramOn(plot,Layout(0.35,0.88));
		//bkg2->plotOn(plot,RooFit::LineColor(kRed));
		plot->SetTitle("data_obs");
		TCanvas* c1;
		c1 = new TCanvas("c1","", 600, 600);
		//c1->SetLogy(1);
		plot->Draw();
		c1->SetLogy(1);
		plot->SetMinimum(0.5);
		
		TLegend *leg;
		leg = new TLegend(.56, .56, .88, .88);
		setLeg(leg);
		plot->Print();
		leg->AddEntry(plot->findObject("h_data_obs"),"Ching-Wei's");
		leg->AddEntry(plot->findObject("h_data_obs2"),"HH4b");
		leg->Draw("same");
		c1->Print(Form("%s.pdf",DBT[j].data()));
	}
	
	
}