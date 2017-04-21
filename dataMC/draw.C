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
#include <string>
#include <sstream>
//#include "../setNCUStyle.C"

# define IntegratedLumi 35900

TCanvas* c1;

void draw(){
	//setNCUStyle();
	//c1 = new TCanvas("c1","",1360,768);
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	TFile *f[40];
	for(int i=0;i<11;i++){
		//f[i]=TFile::Open(Form("sf2/B%s.root",masspoint[i].data()));
		//f[i+11]=TFile::Open(Form("sf2/R%s.root",masspoint[i].data()));
	}
	f[22]=TFile::Open("QCD500.root");
	f[23]=TFile::Open("QCD700.root");
	f[24]=TFile::Open("QCD1000.root");
	f[25]=TFile::Open("QCD1500.root");
	f[26]=TFile::Open("QCD2000.root");
	f[27]=TFile::Open("data.root");
	
	double Xsec[5]={32100,6831,1207,119.9,25.24};
	double norm[5];
	TH1D* ct[6];
	for(int i=0;i<5;i++){
		ct[i]=(TH1D*)f[i+22]->FindObjectAny("fixScale");
		//TH1D* fixScale=(TH1D*)f[i+22]->FindObjectAny("");
		norm[i]=IntegratedLumi*Xsec[i]/ct[i]->GetBinContent(1);
		ct[i]->Scale(norm[i]);
		cout<<i<<endl;
	}
	
	for(int k=1;k<7;k++){
		double total2=0;
		for (int i=0;i<5;i++){
			//cout<<i<<","<<ct[i]->GetBinContent(i)<<endl;
			total2+=ct[i]->GetBinContent(k);
		}
		
		cout<<k<<"="<<total2<<endl;
		ct[5]=(TH1D*)f[27]->FindObjectAny("fixScale");
		cout<< ct[5]->GetBinContent(k)<<endl;
		
		
	}
	
	/*
	for(int i=0;i<11;i++){
		ct[i]=(TH1D*)f[i+11]->FindObjectAny("cutflow");
	}
	myfile<<endl;
	myfile<<""<<endl;
	myfile<<"\\hline"<<endl;
	myfile<<"Radion selections &1000 & 1200 & 1400 & 1600 & 1800 & 2000 & 2500 & 3000 & 3500 & 4000 & 4500"<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int j=0;j<14;j++){
		myfile<<ss[j]<<" & ";
		for(int i=0;i<11;i++){
		
			if (i<10 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<< ct[i]->GetBinContent(j+1)<<"\\\\";
		}
		myfile<<endl;
	}
	myfile<<endl;
	myfile<<""<<endl;
	myfile<<"\\hline"<<endl;
	myfile<<"QCD and data selections &QCD 700-1000 & QCD 1000-1500 & QCD 1500-2000 & QCD 2000-inf. &data"<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int i=0;i<5;i++){
		ct[i]=(TH1D*)f[i+22]->FindObjectAny("cutflow");
	}
	for(int j=0;j<14;j++){
		myfile<<ss[j]<<" & ";
		for(int i=0;i<5;i++){
		
			if (i<4 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<< ct[i]->GetBinContent(j+1)<<"\\\\";
		}
		myfile<<endl;
	}
	myfile<<endl;
	myfile<<"\\hline"<<endl;
	myfile<<"QCD normalized selections &QCD 700-1000 & QCD 1000-1500 & QCD 1500-2000 & QCD 2000-inf. & total"<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int i=0;i<4;i++){
		ct[i]=(TH1D*)f[i+22]->FindObjectAny("cutflowS");
	}
	for(int j=0;j<14;j++){
		myfile<<ss[j]<<" & ";
		double sum=0;
		for(int i=0;i<5;i++){
			if (i<4 )sum+=ct[i]->GetBinContent(j+1);
			if (i<4 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<< sum<<"\\\\";
		}
		myfile<<endl;
	}
	*/
}
