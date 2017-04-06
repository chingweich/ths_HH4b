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
#include "../../HHbbbbAnalyzer/untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
//#include "../../setNCUStyle.C"


TH1D* makeEff(string input,double workingPoint=0.3){
	double bins[5]={250,350,430,840,2000};
	TH1D* th1=new TH1D("all","",4,bins);
	TH1D* th2=new TH1D("pass","",4,bins);
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			Float_t  etadiff = data.GetFloat("etadiff");
			Float_t  jet1pt = data.GetFloat("jet1pt");
			Float_t  jet2pt = data.GetFloat("jet2pt");			
			//Float_t  dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			Float_t  jet1_puppi_msoftdrop_TheaCorr = data.GetFloat("jet1_puppi_msoftdrop_TheaCorr");
			Float_t  jet2_puppi_msoftdrop_TheaCorr = data.GetFloat("jet2_puppi_msoftdrop_TheaCorr");
			//if(jet1_puppi_msoftdrop_TheaCorr<105 ||jet1_puppi_msoftdrop_TheaCorr>135)continue;
			if(jet2_puppi_msoftdrop_TheaCorr<105 ||jet2_puppi_msoftdrop_TheaCorr>135)continue;
			//if(dijetmass_softdrop_corr<mjj*0.8||dijetmass_softdrop_corr>mjj*1.2)continue;
			//if(jet1bbtag<0.3 ||jet2bbtag<0.3)continue;
			//if(jet1bbtag>0.8 &&jet2bbtag>0.8)continue;
			if(etadiff>1.3)continue;
			th1->Fill(jet1pt);
			if(jet1bbtag>workingPoint)th2->Fill(jet1pt);
			th1->Fill(jet2pt);
			if(jet2bbtag>workingPoint)th2->Fill(jet2pt);
	}
	th2->Divide(th1);
	//tf1->Close();
	return th2;
}

void skimB(string input,string output,TH1D* th1[]){
  TFile* f;
  f =  TFile::Open(Form("%s",input.data()));
 // TDirectory * dir= (TDirectory*)f->Get(Form("%sNCUGlobalTuples_%d.root:/tree",st.data(),w));
  TTree* tree;
  f->GetObject("mynewTree",tree);
  //TTree* newtree = tree->CopyTree("AK8PuppijetSDmass[0]>50","",1000000000,0);
  TFile* outputfile = new TFile(Form("%s",output.data()),"RECREATE");
 TTree* mynewTree=new TTree("mynewTree","mynewTree");
		
		//Float_t ; mynewTree->Branch("",&,"/F");
		
		Float_t jet1pt; mynewTree->Branch("jet1pt",&jet1pt,"jet1pt/F");
		Float_t jet2pt; mynewTree->Branch("jet2pt",&jet2pt,"jet2pt/F");
		Float_t jet1eta; mynewTree->Branch("jet1eta",&jet1eta,"jet1eta/F");
		Float_t jet2eta; mynewTree->Branch("jet2eta",&jet2eta,"jet2eta/F");
		Float_t jet1phi; mynewTree->Branch("jet1phi",&jet1phi,"jet1phi/F");
		Float_t jet2phi; mynewTree->Branch("jet2phi",&jet2phi,"jet2phi/F");
		Float_t jet1mass; mynewTree->Branch("jet1mass",&jet1mass,"jet1mass/F");
		Float_t jet2mass; mynewTree->Branch("jet2mass",&jet2mass,"jet2mass/F");
		Float_t etadiff; mynewTree->Branch("etadiff",&etadiff,"etadiff/F");
		Float_t dijetmass; mynewTree->Branch("dijetmass",&dijetmass,"dijetmass/F");
		Float_t dijetmass_pruned_corr; mynewTree->Branch("dijetmass_pruned_corr",&dijetmass_pruned_corr,"dijetmass_pruned_corr/F");
		Float_t dijetmass_softdrop_corr; mynewTree->Branch("dijetmass_softdrop_corr",&dijetmass_softdrop_corr,"dijetmass_softdrop_corr/F");
		Float_t dijetmass_corr_punc; mynewTree->Branch("dijetmass_corr_punc",&dijetmass_corr_punc,"dijetmass_corr_punc/F");
		Float_t dijetmass_softdrop_corr_JERup; mynewTree->Branch("dijetmass_softdrop_corr_JERup",&dijetmass_softdrop_corr_JERup,"dijetmass_softdrop_corr_JERup/F");
		Float_t dijetmass_softdrop_corr_JERdown; mynewTree->Branch("dijetmass_softdrop_corr_JERdown",&dijetmass_softdrop_corr_JERdown,"dijetmass_softdrop_corr_JERdown/F");
		Float_t dijetmass_softdrop_corr_JECup; mynewTree->Branch("dijetmass_softdrop_corr_JECup",&dijetmass_softdrop_corr_JECup,"dijetmass_softdrop_corr_JECup/F");
		Float_t dijetmass_softdrop_corr_JECdown; mynewTree->Branch("dijetmass_softdrop_corr_JECdown",&dijetmass_softdrop_corr_JECdown,"dijetmass_softdrop_corr_JECdown/F");
		
		const int nPDFMax=120;
		Int_t nPDF ;mynewTree->Branch("nPDF",&nPDF,"nPDF/I");
		Float_t pdfSF[nPDFMax];mynewTree->Branch("pdfSF",&pdfSF,"pdfSF[nPDF]/F");
		Float_t jet1pmass; mynewTree->Branch("jet1pmass",&jet1pmass,"jet1pmass/F");
		Float_t jet2pmass; mynewTree->Branch("jet2pmass",&jet2pmass,"jet2pmass/F");
		Float_t jet1bbtag; mynewTree->Branch("jet1bbtag",&jet1bbtag,"jet1bbtag/F");
		Float_t jet2bbtag; mynewTree->Branch("jet2bbtag",&jet2bbtag,"jet2bbtag/F");
		Float_t jet1_puppi_pt; mynewTree->Branch("jet1_puppi_pt",&jet1_puppi_pt,"jet1_puppi_pt/F");
		Float_t jet2_puppi_pt; mynewTree->Branch("jet2_puppi_pt",&jet2_puppi_pt,"jet2_puppi_pt/F");
		Float_t jet1_puppi_eta; mynewTree->Branch("jet1_puppi_eta",&jet1_puppi_eta,"jet1_puppi_eta/F");
		Float_t jet2_puppi_eta; mynewTree->Branch("jet2_puppi_eta",&jet2_puppi_eta,"jet2_puppi_eta/F");
		Float_t jet1_puppi_phi; mynewTree->Branch("jet1_puppi_phi",&jet1_puppi_phi,"jet1_puppi_phi/F");
		Float_t jet2_puppi_phi; mynewTree->Branch("jet2_puppi_phi",&jet2_puppi_phi,"jet2_puppi_phi/F");
		Float_t jet1_puppi_mass; mynewTree->Branch("jet1_puppi_mass",&jet1_puppi_mass,"jet1_puppi_mass/F");
		Float_t jet2_puppi_mass; mynewTree->Branch("jet2_puppi_mass",&jet2_puppi_mass,"jet2_puppi_mass/F");
		
		Float_t jet1_puppi_tau21; mynewTree->Branch("jet1_puppi_tau21",&jet1_puppi_tau21,"jet1_puppi_tau21/F");
		Float_t jet2_puppi_tau21; mynewTree->Branch("jet2_puppi_tau21",&jet2_puppi_tau21,"jet2_puppi_tau21/F");
		Float_t jet1_puppi_msoftdrop; mynewTree->Branch("jet1_puppi_msoftdrop",&jet1_puppi_msoftdrop,"jet1_puppi_msoftdrop/F");
		Float_t jet2_puppi_msoftdrop; mynewTree->Branch("jet2_puppi_msoftdrop",&jet2_puppi_msoftdrop,"jet2_puppi_msoftdrop/F");
		Float_t jet1_puppi_msoftdrop_TheaCorr; mynewTree->Branch("jet1_puppi_msoftdrop_TheaCorr",&jet1_puppi_msoftdrop_TheaCorr,"jet1_puppi_msoftdrop_TheaCorr/F");
		Float_t jet2_puppi_msoftdrop_TheaCorr; mynewTree->Branch("jet2_puppi_msoftdrop_TheaCorr",&jet2_puppi_msoftdrop_TheaCorr,"jet2_puppi_msoftdrop_TheaCorr/F");
		Float_t jet1_puppi_TheaCorr; mynewTree->Branch("jet1_puppi_TheaCorr",&jet1_puppi_TheaCorr,"jet1_puppi_TheaCorr/F");
		Float_t jet2_puppi_TheaCorr; mynewTree->Branch("jet2_puppi_TheaCorr",&jet2_puppi_TheaCorr,"jet2_puppi_TheaCorr/F");
		Float_t nTrueInt; mynewTree->Branch("nTrueInt",&nTrueInt,"nTrueInt/F");
		Float_t puWeights=1; mynewTree->Branch("puWeights",&puWeights,"puWeights/F");
		Float_t puWeightsUp=1; mynewTree->Branch("puWeightsUp",&puWeightsUp,"puWeightsUp/F");
		Float_t puWeightsDown=1; mynewTree->Branch("puWeightsDown",&puWeightsDown,"puWeightsDown/F");
		
		Float_t jet1JERup; mynewTree->Branch("jet1JERup",&jet1JERup,"jet1JERup/F");
		Float_t jet1JERcentral; mynewTree->Branch("jet1JERcentral",&jet1JERcentral,"jet1JERcentral/F");
		Float_t jet1JERdown; mynewTree->Branch("jet1JERdown",&jet1JERdown,"jet1JERdown/F");
		Float_t jet1JECup; mynewTree->Branch("jet1JECup",&jet1JECup,"jet1JECup/F");
		Float_t jet1JECdown; mynewTree->Branch("jet1JECdown",&jet1JECdown,"jet1JECdown/F");
		
		Float_t jet2JERup; mynewTree->Branch("jet2JERup",&jet2JERup,"jet2JERup/F");
		Float_t jet2JERcentral; mynewTree->Branch("jet2JERcentral",&jet2JERcentral,"jet2JERcentral/F");
		Float_t jet2JERdown; mynewTree->Branch("jet2JERdown",&jet2JERdown,"jet2JERdown/F");
		Float_t jet2JECup; mynewTree->Branch("jet2JECup",&jet2JECup,"jet2JECup/F");
		Float_t jet2JECdown; mynewTree->Branch("jet2JECdown",&jet1JECdown,"jet1JECdown/F");
		
		//Float_t isData; mynewTree->Branch("isData",&isData,"isData/F");
		
		Float_t dbtSF=1; mynewTree->Branch("dbtSF",&dbtSF,"dbtSF/F");
		Float_t dbtSFup=1; mynewTree->Branch("dbtSFup",&dbtSFup,"dbtSFup/F");
		Float_t dbtSFdown=1; mynewTree->Branch("dbtSFdown",&dbtSFdown,"dbtSFdown/F");
		//Float_t SFLoose=1; mynewTree->Branch("SFLoose",&SFLoose,"SFLoose/F");
		//Float_t SFLooseup=1; mynewTree->Branch("SFLooseup",&SFLooseup,"SFLooseup/F");
		//Float_t SFLoosedown=1; mynewTree->Branch("SFLoosedown",&SFLoosedown,"SFLoosedown/F");
		
		Float_t trigWeightUp=1; mynewTree->Branch("trigWeightUp",&trigWeightUp,"trigWeightUp/F");
		Float_t trigWeightDown=1; mynewTree->Branch("trigWeightDown",&trigWeightDown,"trigWeightDown/F");
		
		bool HLT_PFHT800=1; mynewTree->Branch("HLT_PFHT800",&HLT_PFHT800,"HLT_PFHT800/O");
		bool HLT_PFHT900=1; mynewTree->Branch("HLT_PFHT900",&HLT_PFHT900,"HLT_PFHT900/O");
		bool HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=1; mynewTree->Branch("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5",&HLT_PFHT650_WideJetMJJ900DEtaJJ1p5,"HLT_PFHT650_WideJetMJJ900DEtaJJ1p5/O");
		bool HLT_AK8PFJet360_TrimMass30=1; mynewTree->Branch("HLT_AK8PFJet360_TrimMass30",&HLT_AK8PFJet360_TrimMass30,"HLT_AK8PFJet360_TrimMass30/O");
		bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=1; mynewTree->Branch("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",&HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20,"HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20/O");
		bool HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=1; mynewTree->Branch("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50",&HLT_AK8PFHT650_TrimR0p1PT0p03Mass50,"HLT_AK8PFHT650_TrimR0p1PT0p03Mass50/O");
		bool HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=1; mynewTree->Branch("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50",&HLT_AK8PFHT700_TrimR0p1PT0p03Mass50,"HLT_AK8PFHT700_TrimR0p1PT0p03Mass50/O");
		
		
	 TreeReader data(tree);
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
	  
	  data.GetEntry(jEntry);
	  Float_t  deta = data.GetFloat("etadiff");
	  if(deta>1.3)continue;
	  Float_t  jet2SDMass = data.GetFloat("jet2_puppi_msoftdrop_TheaCorr");
	  if(jet2SDMass<105 ||jet2SDMass>135)continue;
	jet1pt = data.GetFloat("jet1pt");
	jet2pt = data.GetFloat("jet2pt");
	jet1eta = data.GetFloat("jet1eta");
	jet2eta = data.GetFloat("jet2eta");
	jet1phi = data.GetFloat("jet1phi");
	jet2phi = data.GetFloat("jet2phi");
	jet1mass = data.GetFloat("jet1mass");
	jet2mass = data.GetFloat("jet2mass");
	etadiff =deta;
	dijetmass = data.GetFloat("dijetmass");
	dijetmass_pruned_corr = data.GetFloat("dijetmass_pruned_corr");
	dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
	dijetmass_softdrop_corr_JECup= data.GetFloat("dijetmass_softdrop_corr_JECup");
	dijetmass_softdrop_corr_JECdown = data.GetFloat("dijetmass_softdrop_corr_JECdown");
	dijetmass_softdrop_corr_JERup = data.GetFloat("dijetmass_softdrop_corr_JERup");
	dijetmass_softdrop_corr_JERdown = data.GetFloat("dijetmass_softdrop_corr_JERdown");
	dijetmass_corr_punc = data.GetFloat("dijetmass_corr_punc");
	jet1pmass = data.GetFloat("jet1pmass");
	jet2pmass = data.GetFloat("jet2pmass");
	jet1bbtag = data.GetFloat("jet1bbtag");
	jet2bbtag = data.GetFloat("jet2bbtag");
	jet1_puppi_pt = data.GetFloat("jet1_puppi_pt");
	jet2_puppi_pt = data.GetFloat("jet2_puppi_pt");
	jet1_puppi_eta = data.GetFloat("jet1_puppi_eta");
	jet2_puppi_eta = data.GetFloat("jet2_puppi_eta");
	jet1_puppi_phi = data.GetFloat("jet1_puppi_phi");
	jet2_puppi_phi = data.GetFloat("jet2_puppi_phi");
	jet1_puppi_mass = data.GetFloat("jet1_puppi_mass");
	jet2_puppi_mass = data.GetFloat("jet2_puppi_mass");
	jet1_puppi_tau21 = data.GetFloat("jet1_puppi_tau21");
	jet2_puppi_tau21 = data.GetFloat("jet2_puppi_tau21");
	jet1_puppi_msoftdrop = data.GetFloat("jet1_puppi_msoftdrop");
	jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
	jet1_puppi_msoftdrop_TheaCorr = data.GetFloat("jet1_puppi_msoftdrop_TheaCorr");
	jet2_puppi_msoftdrop_TheaCorr =jet2SDMass;
	jet1_puppi_TheaCorr = data.GetFloat("jet1_puppi_TheaCorr");
	jet2_puppi_TheaCorr = data.GetFloat("jet2_puppi_TheaCorr");
	nTrueInt = data.GetFloat("nTrueInt");
	puWeights = data.GetFloat("puWeights");
	puWeightsUp = data.GetFloat("puWeightsUp");
	puWeightsDown = data.GetFloat("puWeightsDown");
	jet1JERup = data.GetFloat("jet1JERup");
	jet1JERcentral = data.GetFloat("jet1JERcentral");
	jet1JERdown = data.GetFloat("jet1JERdown");
	jet1JECup = data.GetFloat("jet1JECup");
	jet1JECdown = data.GetFloat("jet1JECdown");
	jet2JERup = data.GetFloat("jet2JERup");
	jet2JERcentral = data.GetFloat("jet2JERcentral");
	jet2JERdown = data.GetFloat("jet2JERdown");
	jet2JECup = data.GetFloat("jet2JECup");
	jet2JECdown = data.GetFloat("jet2JECdown");
	HLT_PFHT800 = data.GetBool("HLT_PFHT800");
	HLT_PFHT900 = data.GetBool("HLT_PFHT900");
	HLT_PFHT650_WideJetMJJ900DEtaJJ1p5 = data.GetBool("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5");
	HLT_AK8PFJet360_TrimMass30 = data.GetBool("HLT_AK8PFJet360_TrimMass30");
	HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20 = data.GetBool("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20");
	HLT_AK8PFHT650_TrimR0p1PT0p03Mass50 = data.GetBool("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50");
	HLT_AK8PFHT700_TrimR0p1PT0p03Mass50 = data.GetBool("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50");
		
		
	nPDF=110;
			Float_t*  pdfscaleSysWeights= data.GetPtrFloat("pdfSF");
			for (int i=0;i<nPDF;i++)pdfSF[i]=pdfscaleSysWeights[i];
	dbtSF =1;
	dbtSFup = 1;
	dbtSFdown = 1;
	
	
	for(int i=0;i<2;i++){
		double pt=0,bbtag=0;
		if(i==0){
			pt=jet1pt;
			bbtag=jet1bbtag;
		}
		else {
			pt=jet2pt;
			bbtag=jet2bbtag;
		}
		double SFL=1,SFLUp=1,SFLDown=1,effL=1;
		double SFT=1,SFTUp=1,SFTDown=1,effT=1;
		if (pt<350){
			SFL=0.96;
			SFLUp=0.99;
			SFLDown=0.94;
			effL=th1[0]->GetBinContent(1);
			
			SFT=0.92;
			SFTUp=0.95;
			SFTDown=0.89;
			effT=th1[1]->GetBinContent(1);
		}
		else if (pt <430){
			SFL=1;
			SFLUp=1.04;
			SFLDown=0.97;
			effL=th1[0]->GetBinContent(2);
			
			SFT=1.01;
			SFTUp=1.04;
			SFTDown=0.97;
			effT=th1[1]->GetBinContent(2);
		}
		else if (pt <840){
			SFL=1.01;
			SFLUp=1.03;
			SFLDown=0.97;
			effL=th1[0]->GetBinContent(3);
			
			SFT=0.92;
			SFTUp=0.95;
			SFTDown=0.87;
			effT=th1[1]->GetBinContent(3);
		}
		else {
			SFL=1.01;
			SFLUp=1.05;
			SFLDown=0.93;
			effL=th1[0]->GetBinContent(4);
			
			SFT=0.92;
			SFTUp=0.98;
			SFTDown=0.82;
			effT=th1[1]->GetBinContent(4);
		}
		
		if(bbtag>0.8){
			dbtSF *=SFT;
			dbtSFup *= SFTUp;
			dbtSFdown *= SFTDown;
		}
		else if (bbtag>0.3){
			if(effL==effT)cout<<"!pt="<<pt<<",btag="<<bbtag<<endl;
			dbtSF *=(SFL*effL-SFT*effT)/(effL-effT);
			dbtSFup *=(SFLUp*effL-SFTUp*effT)/(effL-effT);
			dbtSFdown *=(SFLDown*effL-SFTDown*effT)/(effL-effT);
		}
		else {
			dbtSF *=(1-SFL*effL)/(1-effL);
			dbtSFup *=(1-SFLUp*effL)/(1-effL);
			dbtSFdown *=(1-SFLDown*effL)/(1-effL);
		}
		
	}
	
	
	  
	  mynewTree  ->Fill();
  }
 // outputfile->cd();
  mynewTree ->Write();
  //outputfile->Write();
  TH1D* th2=(TH1D*) f->Get("CountWeighted");
  th2->Write();
    outputfile->Close();
}

void setTGraph(TH1D* tg1,int i,bool setMax=0){
	
	tg1->SetLineColor(i+1);
	//tg1->SetLineColor(98-4*i);
	tg1->SetMarkerColor(i+1);
	if(i==4)tg1->SetLineColor(kOrange);
	if(i==4)tg1->SetMarkerColor(kOrange);
	//tg1->SetMarkerColor(98-4*i);
	tg1->SetMarkerStyle(20+i);
	
	tg1->GetXaxis()->SetTitle("p_{T}[GeV]");
	tg1->GetYaxis()->SetTitle("Efficiency");
	//limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
	tg1->SetTitle("");
	tg1->SetLineWidth(2);
	tg1->SetFillColor(0);
	tg1->SetMaximum(1);
	tg1->SetMinimum(0.0014);
	tg1->SetMarkerSize(2);
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

void dbtSFv3(){
	TH1D* th1[2];
	const int nGrav=12;
	int massP[12]={750,800,900,1000,1200,1600,1800,2000,2500,3000,4000,4500};
	//int massP[nGrav]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	/*
	TCanvas * c1 = new TCanvas("c1","",600,600);
	c1->cd();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(57);
	gStyle->SetFrameLineWidth(3);
	gStyle->SetPadRightMargin(0.001);
	gStyle->SetPadLeftMargin(0.13);
	c1->SetLeftMargin(0.15);
	c1->SetRightMargin(0.01);
	gStyle->SetNdivisions(605, "XYZ");
	*/
	for(int i=0;i<nGrav;i++){
		th1[0]=makeEff(Form("minitree/BulkGrav%d.root",massP[i]),0.3);
		th1[1]=makeEff(Form("minitree/BulkGrav%d.root",massP[i]),0.8);
		skimB(Form("minitree/BulkGrav%d.root",massP[i]),Form("dbt/BulkGrav_M-%d_0.root",massP[i]),th1);
		/*
		th1[0]->Draw("");
		th1[1]->Draw("same");
		setTGraph(th1[0],0);
		setTGraph(th1[1],1);
		TLegend *leg;
		leg = new TLegend(.66, .76, .88, .88);
		setLeg(leg);
		leg->AddEntry(th1[0],"loose");
		leg->AddEntry(th1[1],"tight");
		leg->Draw("same");
		//th1[1]->SetLineColor(2);
		c1->Print(Form("plots/eff_grav_%d.pdf",massP[i]));
		*/
	}
	
	
	
	const int nRd=8;
	int massP2[nRd]={750,800,900,1000,1200,1600,3500,4500};
	//int massP2[nRd]={1000,1200,1400,1600,2500,3000,3500,4500};
	for(int i=0;i<nRd;i++){
		th1[0]=makeEff(Form("minitree/Radion%d.root",massP2[i]),0.3);
		th1[1]=makeEff(Form("minitree/Radion%d.root",massP2[i]),0.8);
		skimB(Form("minitree/Radion%d.root",massP2[i]),Form("dbt/Radion%d.root",massP2[i]),th1);
	}
	
}