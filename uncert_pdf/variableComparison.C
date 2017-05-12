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
#include "../untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <TProfile.h>
#include <string>
#include <sstream>
//#include "../../setNCUStyle.C"


TH1D* massPlotBaseTT(string input,int pdfIndex=10.,double down=100,double up=1000){
	TH1D* th1=new TH1D("","",100,down,up);
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			//Float_t  vari = data.GetFloat(Form("%s",variable.data()));
			//Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			Float_t  etadiff = data.GetFloat("etadiff");
			Float_t  dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
			Float_t  jet1_puppi_tau21 = data.GetFloat("jet1_puppi_tau21");
			Float_t  jet2_puppi_tau21 = data.GetFloat("jet2_puppi_tau21");
			Float_t  jet1_puppi_msoftdrop = data.GetFloat("jet1_puppi_msoftdrop");
			Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			if(jet2bbtag<0.8 || jet1bbtag<0.8 ) continue;
			if( jet1_puppi_msoftdrop<105 || jet1_puppi_msoftdrop>135)continue;
			if( jet2_puppi_msoftdrop<105 || jet2_puppi_msoftdrop>135)continue;
			
			Float_t *  pdfSF = data.GetPtrFloat("pdfSF"); 
			if(etadiff>1.3)continue;
			if(dijetmass_softdrop_corr<750)continue;
			if(jet2_puppi_tau21>0.55 || jet1_puppi_tau21>0.55)continue;
			
			
			
			th1->Fill(dijetmass_softdrop_corr,pdfSF[pdfIndex]/pdfSF[9]);
	}
	return th1;
}
TH1D* massPlotBaseLL(string input,int pdfIndex=10.,double down=100,double up=1000){
	TH1D* th1=new TH1D("","",100,down,up);
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			//Float_t  vari = data.GetFloat(Form("%s",variable.data()));
			//Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			Float_t  etadiff = data.GetFloat("etadiff");
			Float_t  dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
			Float_t  jet1_puppi_tau21 = data.GetFloat("jet1_puppi_tau21");
			Float_t  jet2_puppi_tau21 = data.GetFloat("jet2_puppi_tau21");
			Float_t  jet1_puppi_msoftdrop = data.GetFloat("jet1_puppi_msoftdrop");
			Float_t  jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			if(jet2bbtag<0.3 || jet1bbtag<0.3 ) continue;
			if(jet2bbtag>0.8 && jet1bbtag>0.8 ) continue;
			if( jet1_puppi_msoftdrop<105 || jet1_puppi_msoftdrop>135)continue;
			if( jet2_puppi_msoftdrop<105 || jet2_puppi_msoftdrop>135)continue;
			
			Float_t *  pdfSF = data.GetPtrFloat("pdfSF"); 
			if(etadiff>1.3)continue;
			if(dijetmass_softdrop_corr<750)continue;
			if(jet2_puppi_tau21>0.55 || jet1_puppi_tau21>0.55)continue;
			
			
			
			th1->Fill(dijetmass_softdrop_corr,pdfSF[pdfIndex]/pdfSF[9]);
			//th1->Fill(dijetmass_softdrop_corr,pdfSF[pdfIndex]);
	}
	//tf1->Close();
	return th1;
	
}


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

void variableComparison(){
	TH1D* th1[110];
	bool TT=0;
	gStyle->SetPadRightMargin(0.039);
gStyle->SetPadLeftMargin(0.129);
gStyle->SetNdivisions(605, "XYZ");
gStyle->SetOptStat(0);
TCanvas* c=new TCanvas("c","",0,0,600,600);
gStyle->SetPadRightMargin(0.029);
gStyle->SetPadLeftMargin(0.1509);
	
	int masspoint[6]={1200,1600,1800,2000,2500,3000};
	double masspointd[6]={1200,1600,1800,2000,2500,3000};
	double effRMS[6]={0};
	if(TT){
		for(int j=0;j<6;j++){
			TProfile * tp=new TProfile("tp","tp",100,masspoint[j]*0.7,masspoint[j]*1.3);
			TH1D* th2= new TH1D("","",100,0,0.08);
			double central=0;
			for(int i=9;i<110;i++){
				th1[i]=massPlotBaseTT(Form("dbt/BulkGrav_M-%d_0.root",masspoint[j]),i,masspoint[j]*0.7,masspoint[j]*1.3);
				if(i==9){
					th1[i]->Draw("hist");
					th1[i]->SetMaximum(th1[i]->GetMaximum()*1.3);
				}
				th1[i]->Draw("same hist");
				for(int k=0;k<th1[9]->GetNbinsX();k++){
					tp->Fill(th1[i]->GetBinCenter(k),th1[i]->GetBinContent(k));
				}
				//cout<<"int="<<;
				if(i==9){
					cout<< masspoint[j]<<",="<<th1[i]->Integral()/50000<<endl;
					central=th1[i]->Integral()/50000;
				} 
				else th2->Fill(th1[i]->Integral()/50000);
			}
			cout<< masspoint[j]<<",="<<th2->GetRMS()<<endl;
			th1[9]->SetLineColor(2);
			th1[9]->Draw("same hist");
			c->Print(Form("mjj/TT_mjj%d.pdf",masspoint[j]));
			th2->Draw();
			c->Print(Form("eff/TT_mjj%d.pdf",masspoint[j]));
			effRMS[j]=th2->GetRMS()/central;
			for(int i=9;i<110;i++){
				th1[i]->Scale(1/th1[i]->Integral());
				if(i==9){
					th1[i]->Draw("hist");
					th1[i]->SetMaximum(th1[i]->GetMaximum()*1.3);
				}
				th1[i]->Draw("same hist");
				for(int k=0;k<th1[9]->GetNbinsX();k++){
					tp->Fill(th1[i]->GetBinCenter(k),th1[i]->GetBinContent(k));
				}
			}
			th1[9]->SetLineColor(2);
			th1[9]->Draw("same hist");
			c->Print(Form("mjj/TT_mjj%d_norm.pdf",masspoint[j]));
			tp->Draw("");
			c->Print(Form("mjj/profile_TT_mjj%d.pdf",masspoint[j]));
		}
		TGraph* tg1=new TGraph(6,masspointd,effRMS);
		for(int i=0;i<6;i++)cout<<effRMS[i]<<endl;
		tg1->Draw("APL");
		c->Print("eff/uncert_pdf_eff_TT.pdf");
	}
	else {
		for(int j=0;j<6;j++){
			TProfile * tp=new TProfile("tp","tp",100,masspoint[j]*0.7,masspoint[j]*1.3);
			TH1D* th2= new TH1D("","",100,0,0.14);
			double central=0;
			for(int i=9;i<110;i++){
				th1[i]=massPlotBaseLL(Form("dbt/BulkGrav_M-%d_0.root",masspoint[j]),i,masspoint[j]*0.7,masspoint[j]*1.3);
				if(i==9){
					th1[i]->Draw("hist");
					th1[i]->SetMaximum(th1[i]->GetMaximum()*1.3);
				}
				th1[i]->Draw("same hist");
				for(int k=0;k<th1[9]->GetNbinsX();k++){
					tp->Fill(th1[i]->GetBinCenter(k),th1[i]->GetBinContent(k));
				}
				//cout<<"int="<<;
				if(i==9){
					cout<< masspoint[j]<<",="<<th1[i]->Integral()/50000<<endl;
					central=th1[i]->Integral()/50000;
				} 
				else th2->Fill(th1[i]->Integral()/50000);
			}
			cout<< masspoint[j]<<",="<<th2->GetRMS()<<endl;
			th1[9]->SetLineColor(2);
			th1[9]->Draw("same hist");
			c->Print(Form("mjj/LL_mjj%d.pdf",masspoint[j]));
			th2->Draw();
			c->Print(Form("eff/LL_mjj%d.pdf",masspoint[j]));
			effRMS[j]=th2->GetRMS()/central;
			for(int i=9;i<110;i++){
				th1[i]->Scale(1/th1[i]->Integral());
				if(i==9){
					th1[i]->Draw("hist");
					th1[i]->SetMaximum(th1[i]->GetMaximum()*1.3);
				}
				th1[i]->Draw("same hist");
				for(int k=0;k<th1[9]->GetNbinsX();k++){
					tp->Fill(th1[i]->GetBinCenter(k),th1[i]->GetBinContent(k));
				}
			}
			th1[9]->SetLineColor(2);
			th1[9]->Draw("same hist");
			c->Print(Form("mjj/LL_mjj%d_norm.pdf",masspoint[j]));
			tp->Draw("");
			c->Print(Form("mjj/profile_LL_mjj%d.pdf",masspoint[j]));
		}
		TGraph* tg1=new TGraph(6,masspointd,effRMS);
		for(int i=0;i<6;i++)cout<<effRMS[i]<<endl;
		tg1->Draw("APL");
		c->Print("eff/uncert_pdf_eff_LL.pdf");
	}
	
	
}