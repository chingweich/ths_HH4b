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
#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>


void massPlotBase(string input){
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
	
	ofstream myfile;
	myfile.open ("LL.txt");
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			if(jEntry%1000000==0)cout<<jEntry<<" out of "<<data.GetEntriesFast()<<endl;
			data.GetEntry(jEntry);
			
			Float_t  jet1_puppi_msoftdrop_TheaCorr = data.GetFloat("jet1_puppi_msoftdrop_TheaCorr");
			Float_t  jet2_puppi_msoftdrop_TheaCorr = data.GetFloat("jet2_puppi_msoftdrop_TheaCorr");
			if(jet1_puppi_msoftdrop_TheaCorr<105||jet1_puppi_msoftdrop_TheaCorr>135)continue;
			if(jet2_puppi_msoftdrop_TheaCorr<105||jet2_puppi_msoftdrop_TheaCorr>135)continue;
			Float_t  jet2pt = data.GetFloat("jet2pt");
			Float_t  jet1pt = data.GetFloat("jet1pt");
			if(jet2pt<300|| jet1pt<300)continue;
			Float_t  jet1eta = data.GetFloat("jet1eta");
			Float_t  jet2eta = data.GetFloat("jet2eta");
			
			if(fabs(jet1eta)>2.4||fabs(jet2eta)>2.4)continue;
			
			Float_t  etadiff = data.GetFloat("etadiff");
			if(etadiff>1.3)continue;
			Float_t  dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
			if(dijetmass_softdrop_corr<1100)continue;
			Float_t  jet1_puppi_tau21 = data.GetFloat("jet1_puppi_tau21");
			Float_t  jet2_puppi_tau21 = data.GetFloat("jet2_puppi_tau21");
			if(jet1_puppi_tau21>0.55||jet2_puppi_tau21>0.55)continue;
			bool HLT_PFHT800=data.GetBool("HLT_PFHT800");
			bool HLT_PFHT900=data.GetBool("HLT_PFHT900");
			bool HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=data.GetBool("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5");
			bool HLT_AK8PFJet360_TrimMass30=data.GetBool("HLT_AK8PFJet360_TrimMass30");
			bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=data.GetBool("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20");
			bool HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=data.GetBool("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50");
			bool HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=data.GetBool("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50");
			if(!(HLT_PFHT800||HLT_PFHT900||HLT_PFHT650_WideJetMJJ900DEtaJJ1p5||HLT_AK8PFJet360_TrimMass30||HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20||HLT_AK8PFHT650_TrimR0p1PT0p03Mass50||HLT_AK8PFHT700_TrimR0p1PT0p03Mass50))continue;
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			if(jet1bbtag<0.3||jet2bbtag<0.3)continue;
			if(jet1bbtag>0.8 && jet2bbtag>0.8)continue;
			
			Int_t run=data.GetInt("run");
			Int_t lumi=data.GetInt("lumi");
			Int_t event=data.GetInt("event");
			myfile<<run<<":"<<lumi<<":"<<event<<endl;
			
	}
}

void massPlotBaseTT(string input){
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
	
	ofstream myfile;
	myfile.open ("TT.txt");
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			if(jEntry%1000000==0)cout<<jEntry<<" out of "<<data.GetEntriesFast()<<endl;
			data.GetEntry(jEntry);
			Float_t  jet1_puppi_msoftdrop_TheaCorr = data.GetFloat("jet1_puppi_msoftdrop_TheaCorr");
			Float_t  jet2_puppi_msoftdrop_TheaCorr = data.GetFloat("jet2_puppi_msoftdrop_TheaCorr");
			if(jet1_puppi_msoftdrop_TheaCorr<105||jet1_puppi_msoftdrop_TheaCorr>135)continue;
			if(jet2_puppi_msoftdrop_TheaCorr<105||jet2_puppi_msoftdrop_TheaCorr>135)continue;
			Float_t  jet2pt = data.GetFloat("jet2pt");
			Float_t  jet1pt = data.GetFloat("jet1pt");
			if(jet2pt<300|| jet1pt<300)continue;
			Float_t  jet1eta = data.GetFloat("jet1eta");
			Float_t  jet2eta = data.GetFloat("jet2eta");
			
			if(fabs(jet1eta)>2.4||fabs(jet2eta)>2.4)continue;
			
			Float_t  etadiff = data.GetFloat("etadiff");
			if(etadiff>1.3)continue;
			Float_t  dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
			if(dijetmass_softdrop_corr<1100)continue;
			Float_t  jet1_puppi_tau21 = data.GetFloat("jet1_puppi_tau21");
			Float_t  jet2_puppi_tau21 = data.GetFloat("jet2_puppi_tau21");
			if(jet1_puppi_tau21>0.55||jet2_puppi_tau21>0.55)continue;
			bool HLT_PFHT800=data.GetBool("HLT_PFHT800");
			bool HLT_PFHT900=data.GetBool("HLT_PFHT900");
			bool HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=data.GetBool("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5");
			bool HLT_AK8PFJet360_TrimMass30=data.GetBool("HLT_AK8PFJet360_TrimMass30");
			bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=data.GetBool("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20");
			bool HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=data.GetBool("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50");
			bool HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=data.GetBool("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50");
			if(!(HLT_PFHT800||HLT_PFHT900||HLT_PFHT650_WideJetMJJ900DEtaJJ1p5||HLT_AK8PFJet360_TrimMass30||HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20||HLT_AK8PFHT650_TrimR0p1PT0p03Mass50||HLT_AK8PFHT700_TrimR0p1PT0p03Mass50))continue;
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			if(jet1bbtag<0.8||jet2bbtag<0.8)continue;
			//if(jet1bbtag>0.8 && jet2bbtag>0.8)continue;
			
			Int_t run=data.GetInt("run");
			Int_t lumi=data.GetInt("lumi");
			Int_t event=data.GetInt("event");
			myfile<<run<<":"<<lumi<<":"<<event<<endl;
			
	}
}

void variableRun(){
	
	
	massPlotBase("dataRun.root");
	massPlotBaseTT("dataRun.root");
	
	
	
}