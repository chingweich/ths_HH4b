#include <string>
#include <iostream>
#include <TPad.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TTree.h>
#include <TKey.h>
#include <TSystemDirectory.h>
//#include "../../setNCUStyle.C"
#include <TSystem.h>
//#include "../readHists.h"


#define  integratedLumi 35900
 //#define fixNumberGlobal 1.34975e+06/1.75969e+06
 
 
 //#define fixNumberGlobal 5.80357e+06/8.92022e+06
 #define fixNumberGlobal 5.80357e+06/8.30085e+06
 
 
 //#define fixNumberGlobal 1

double rangeUserUp=0,rangeUserDown=0;
int isSetRange=0;

void myPlot(vector< TH1D*> h_Z,
           vector< TH1D*> v_data,
	    TH1D* h_data, TH1D* h_bkg,int option=0,string xtitle=""){

  h_data->Reset();
  for(unsigned int i=0;i<v_data.size();i++)h_data->Add(v_data[i]);
  
  TLegend *leg = new TLegend(0.8, 0.68, 0.92, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
int colorNumQCD[4]={46,15,38,30};
   h_bkg->Reset();
    THStack *h_stack = new THStack("h_stack", "");
  for(unsigned int i=0;i<h_Z.size();i++){
	  h_Z[i]->SetFillColor(colorNumQCD[i]);
	  h_Z[i]->SetLineColor(kBlack);
	  h_bkg->Add(h_Z[i]);
	   h_stack->Add(h_Z[i]);
	  
  }
   
   
   if(option==1){
	   leg->Clear();
	   leg->AddEntry(h_Z[3], "light", "f");
	leg->AddEntry(h_Z[2], "cc,c", "f");
	leg->AddEntry(h_Z[1], "b", "f");
	leg->AddEntry(h_Z[0], "bb", "f");
	}
	else{
		leg->AddEntry(h_Z[0], "QCD700", "f");
   leg->AddEntry(h_Z[1], "QCD1000", "f");
   leg->AddEntry(h_Z[2], "QCD1500", "f");
   leg->AddEntry(h_Z[3], "QCD2000", "f");
	}
  /*
leg->AddEntry(h_Zjets, "Z+Jets", "f");
  leg->AddEntry(h_TT, "t#bar{t}", "f");
  leg->AddEntry(h_WW, "WW", "f");
  leg->AddEntry(h_WZ, "WZ", "f");
  leg->AddEntry(h_ZZ, "ZZ", "f");
  leg->AddEntry(h_ZH, "ZH", "f");
*/

  h_data->SetLineColor(kBlack);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1);
  h_data->GetYaxis()->SetTitleOffset(1.3);
  h_data->GetXaxis()->SetTitle("");
 // h_stack->GetXaxis()->SetTitle("?");
  h_data->GetXaxis()->SetLabelOffset(999);
  h_data->GetXaxis()->SetLabelSize(0);
  
  int bmin=0,bmax=0;

		for (int k=1;k<25001;k++){
			bmin=k;
			if (h_data->GetBinContent(k)/h_data->GetMaximum()>0.02) break;
	}

		for (int k=25000;k>0;k--){
			bmax=k;
			if (h_data->GetBinContent(k)/h_data->GetMaximum()>0.02) break;
	}
	double width=h_data->GetBinWidth(1);
	//rangeUserUp=(bmax-0.5)*width+h_data->GetBinCenter(1);
	//rangeUserDown=(bmin-0.5)*width+h_data->GetBinCenter(1);
	//rangeUserDown=0;
	if(isSetRange)h_data->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
	//if(isSetRange)h_stack->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
	//h_stack->GetXaxis()->SetRangeUser((bmin-0.5)*width+h_data->GetBinCenter(1),(bmax-0.5)*width+h_data->GetBinCenter(1));

  
	  h_stack->SetTitle(h_data->GetTitle());
	if(xtitle.find("Mass")!= std::string::npos &&xtitle.find("total")== std::string::npos
	&& xtitle.find("no")== std::string::npos)h_stack->SetMaximum(h_stack->GetMaximum()*1.3);
	if(xtitle.find("h_nvtx")!= std::string::npos )h_stack->SetMaximum(h_stack->GetMaximum()*1.3);
	h_stack->SetMaximum(h_stack->GetMaximum()>h_data->GetMaximum()?h_stack->GetMaximum()*1.2:h_data->GetMaximum()*1.2);
	
    h_stack->Draw("histe");
	if(xtitle.find("double")!= std::string::npos)h_stack->SetMaximum(h_stack->GetMaximum()*1.5);
	if(h_data->GetBinWidth(1)>1)h_stack->GetHistogram()->GetYaxis()->SetTitle(Form("Events / %2.0f units",h_data->GetBinWidth(1)));
    else h_stack->GetHistogram()->GetYaxis()->SetTitle(Form("Events / %0.2f units",h_data->GetBinWidth(1)));
		if(xtitle.find("pt")!= std::string::npos|| xtitle.find("SD")!= std::string::npos||xtitle.find("jj")!= std::string::npos)h_stack->GetHistogram()->GetYaxis()->SetTitle(Form("Events / %2.0f GeV",h_data->GetBinWidth(1)));
    h_stack->GetHistogram()->GetXaxis()->SetTitle(Form("%s",xtitle.data()));
	
   // h_stack->GetHistogram()->GetYaxis()->SetTitle("Event Numbers");
    //h_stack->GetHistogram()->GetYaxis()->SetTitle("");
    //h_stack->GetHistogram()->GetXaxis()->SetTitle(Form("%s",xtitle.data()));
    h_stack->GetHistogram()->GetYaxis()->SetTitleSize(h_data->GetYaxis()->GetTitleSize());
    h_stack->GetHistogram()->GetYaxis()->SetLabelSize(h_data->GetYaxis()->GetLabelSize());
    h_stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.35);
     h_stack->GetHistogram()->GetXaxis()->SetTitleSize(0.045);
     h_stack->GetHistogram()->GetXaxis()->SetTitleFont(62);
	 h_stack->GetHistogram()->SetMarkerSize(1.6);

 h_stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.7);
 h_stack->GetHistogram()->GetZaxis()->SetTitleOffset(0.75);
 // size of axis labels
 h_stack->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
 h_stack->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
 h_stack->GetHistogram()->GetZaxis()->SetTitleSize(0.035);
 h_stack->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
 h_stack->GetHistogram()->GetYaxis()->SetLabelSize(0.04); 
 h_stack->GetHistogram()->GetZaxis()->SetLabelSize(0.025);
     if(isSetRange)h_stack->GetHistogram()->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
       if(isSetRange)cout<<"set"<<rangeUserUp<<","<<rangeUserDown<<endl;
   // if(xtitle.find("M_{jj}")!= std::string::npos )h_stack->GetHistogram()->GetXaxis()->SetRangeUser(1000,3000);
    //h_stack->GetHistogram()->GetXaxis()->SetTickLength(0);
    //h_stack->GetHistogram()->GetXaxis()->SetLabelOffset(999);
h_data->Draw("elsame");
  
  
   

    
  
  leg->AddEntry(h_data, "Data", "lp");
  leg->Draw();

  TLatex *lar = new TLatex();

  lar->SetTextSize(0.044);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.7, 0.94, " 35.9 fb^{-1} (13 TeV)");
  
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
		lar->DrawLatex(0.2, 0.83, "CMS");
		//lar2->DrawLatex(0.31, 0.77, "#it{#bf{Simulation}} ");
		lar2->DrawLatex(0.31, 0.8, "#it{#bf{Preliminary}} ");

}

void myRatio(TH1D* h_data, TH1D *h_bkg,string xtitle="",int option=0){
cout<<h_bkg->Integral()/h_data->Integral()<<endl;
  TH1D* h_ratio = (TH1D*)h_bkg->Clone("h_ratio");
	
  h_ratio->Reset();

  Int_t nbin = h_ratio->GetNbinsX();
  Float_t ratio[nbin];
  Float_t error[nbin];
  Float_t numer_nbincontent[nbin];
  Float_t denom_nbincontent[nbin];
  Float_t numer_binerror[nbin];
  Float_t denom_binerror[nbin];

  for(Int_t i = 1; i <= nbin; i++){
	
	
    numer_nbincontent[i] = h_data->GetBinContent(i);
    denom_nbincontent[i] = h_bkg ->GetBinContent(i);
    numer_binerror[i]    = h_data->GetBinError(i);
    denom_binerror[i]    = h_bkg ->GetBinError(i);

    if( denom_nbincontent[i] <= 0 || numer_nbincontent[i] <= 0 ) continue;
    if( denom_binerror[i] <= 0 || numer_binerror[i] <= 0 ) continue;

    ratio[i] = (Float_t)numer_nbincontent[i]/denom_nbincontent[i];
    error[i] = (ratio[i])*sqrt(pow(numer_binerror[i]/numer_nbincontent[i],2)+pow(denom_binerror[i]/denom_nbincontent[i],2));

    h_ratio->SetBinContent(i,ratio[i]);
    h_ratio->SetBinError(i,error[i]);

  }

  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerStyle(8);
  h_ratio->SetMarkerSize(1);
  h_ratio->SetTitle("");
  h_ratio->SetXTitle(xtitle.data());
  if(xtitle.find("eta")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_ratio->SetXTitle("#eta (j1)");
	if(xtitle.find("eta")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_ratio->SetXTitle("#eta (j2)");
	if(xtitle.find("pt")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_ratio->SetXTitle("p_{T} (j1) [GeV]");
	if(xtitle.find("pt")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_ratio->SetXTitle("p_{T} (j2) [GeV]");
	if(xtitle.find("SD")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_ratio->SetXTitle("PUPPI soft-drop mass(j1)");
	if(xtitle.find("Tau")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_ratio->SetXTitle("#tau_{21} (j1)");
	if(xtitle.find("Tau")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_ratio->SetXTitle("#tau_{21} (j2)");
	if(xtitle.find("SD")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_ratio->SetXTitle("PUPPI soft-drop mass(j2)");
	if(xtitle.find("delta")!= std::string::npos)h_ratio->SetXTitle("#Delta#eta (j1,j2)");
	if(xtitle.find("double")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_ratio->SetXTitle("double-b tagger (j1)");
	if(xtitle.find("double")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_ratio->SetXTitle("double-b tagger (j2)");
	if(xtitle.find("nVtxWeighted")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_ratio->SetXTitle("number of verteces (weighted)");
	else if(xtitle.find("nVtx")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_ratio->SetXTitle("number of verteces");
  h_ratio->GetYaxis()->SetTitle("Data/MC");
  h_ratio->GetYaxis()->SetTitleOffset(0.45);
  h_ratio->GetXaxis()->SetLabelSize(0.125);
  h_ratio->GetXaxis()->SetLabelOffset(0.005);
  h_ratio->GetXaxis()->SetTitleSize(0.125);
  h_ratio->GetXaxis()->SetTitleOffset(0.92);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.12);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetRangeUser(0,2);
  if(option==1){
	   h_ratio->GetXaxis()->SetLabelSize(0.18);
	  h_ratio->GetXaxis()->SetBinLabel(1,"0b");
	  h_ratio->GetXaxis()->SetBinLabel(2,"1b");
	  h_ratio->GetXaxis()->SetBinLabel(3,"2b");
	  h_ratio->GetXaxis()->SetBinLabel(4,"3b");
	  h_ratio->GetXaxis()->SetBinLabel(5,"4b");
	  
		//h_ratio->GetXaxis()->SetBinLabel(6,"3bHPHP");
  }
  h_ratio->Draw();
  
  h_ratio->SetTitle("");
 //h_data->SetMaximum(3000);

  Float_t x0 = h_bkg->GetXaxis()->GetXmin();
  Float_t x1 = h_bkg->GetXaxis()->GetXmax();
  Float_t y0 = 1.;
  Float_t y1 = 1.;
if(isSetRange){
	x0=rangeUserDown;
	x1=rangeUserUp;
	     if(xtitle.find("M_{jj}")!= std::string::npos ){
		    // x0=1000;
		    // x1=3000;
	     }
		
	
}
  TLine* one = new TLine(x0,y0,x1,y1);

  one->SetLineColor(2);
  one->SetLineStyle(1);
  one->SetLineWidth(2);
  one->Draw("same");
 if(isSetRange)h_ratio->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
 // if(isSetRange)h_stack->GetHistogram()->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
     
   //  if(xtitle.find("M_{jj}")!= std::string::npos )h_ratio->GetXaxis()->SetRangeUser(1000,3000);
  h_ratio->Draw("same");

}

void dataMCplots_v2(string antiTau21=""){

  //setNCUStyle(true);
 
  Float_t up_height     = 0.8;
  Float_t dw_correction = 1.355;
  Float_t dw_height= (1-up_height)*dw_correction;

  
  // To get the name of histograms
  cout<<"reading root file"<<endl;
  TFile* tf1[10];
  tf1[0]=TFile::Open("QCD500.root");
  tf1[1]=TFile::Open("QCD700.root");
  tf1[2]=TFile::Open("QCD1000.root");
  tf1[3]=TFile::Open("QCD1500.root");
  tf1[4]=TFile::Open("QCD2000.root");
  tf1[5]=TFile::Open("data.root");
  cout<<"finishing reading root file"<<endl;
  
  
  TFile* tf2[10];
   //tf2[0]=TFile::Open("sf2/B1400.root");
   //tf2[1]=TFile::Open("sf2/B1800.root");
   //tf2[2]=TFile::Open("sf2/B2500.root");
   //tf2[3]=TFile::Open("sf2/R1800.root");
   
  vector<std::string> h_name;
 
  
 h_name.push_back("pt_j0");
 h_name.push_back("eta_j0");
 h_name.push_back("puppiSDMassThea_j0");
 h_name.push_back("doubleSV_j0");
 h_name.push_back("puppiTau21_j0");
 h_name.push_back("nMu_j0");
 h_name.push_back("nEle_j0");
 h_name.push_back("nCh_j0");
 h_name.push_back("MuEnFr_j0");
 h_name.push_back("EmEnFr_j0");
 h_name.push_back("nSv_j0");
 h_name.push_back("svMass0_j0");
 h_name.push_back("svMass1_j0");
	
h_name.push_back("pt_j1");
h_name.push_back("eta_j1");
h_name.push_back("puppiSDMassThea_j1");
 h_name.push_back("doubleSV_j1");
 h_name.push_back("puppiTau21_j1");
 h_name.push_back("nMu_j1");
 h_name.push_back("nEle_j1");
 h_name.push_back("nCh_j1");
 h_name.push_back("MuEnFr_j1");
 h_name.push_back("EmEnFr_j1");
 h_name.push_back("nSv_j1");
 h_name.push_back("svMass0_j1");
 h_name.push_back("svMass1_j1");
 
 h_name.push_back("deltaEta");
 h_name.push_back("totalMass");
 h_name.push_back("totalMassRed");
 h_name.push_back("nVtx");
 h_name.push_back("nVtxWeighted");
 
 //string antiTau21="_antiDBT";
	//string antiTau21="_antiTau21";
	//;
 TFile* outFile= new TFile(Form("dataMC%s.root",antiTau21.data()),"recreate");
	
  for(unsigned int i = 0; i < h_name.size(); i++){
	  TCanvas c("c","",0,0,600,600);
	  gStyle->SetOptStat(0);
 gStyle->SetPaintTextFormat("2.1f");
 gStyle->SetPalette(57);
 gStyle->SetFrameLineWidth(3);
 gStyle->SetPadRightMargin(0.059);
 gStyle->SetPadLeftMargin(0.159);
	  
  c.Divide(1,2);

  //TPad* c_up = (TPad*) c.GetListOfPrimitives()->FindObject("c_1");
  TPad* c_up = (TPad*) c.cd(1);
  //c_up->SetLogy(1);
  TPad* c_dw = (TPad*) c.GetListOfPrimitives()->FindObject("c_2"); 

  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);
//	 for(unsigned int i = 0; i < 1; i++){
	  
	cout<<h_name[i]<<endl;
	bool drawSignal=0;
	isSetRange=0;
	
	//---------
	 TH1D* thh[6][4];
  double fixNumber=fixNumberGlobal;//22295/30584.3;//7306/11857.3;
  //if(antiTau21.find("_antiTau21")!= std::string::npos)fixNumber=36587/29119.8;
  //else if(antiTau21.find("_antiDBT")!= std::string::npos)fixNumber=1.30839e+06/1.49543e+06;
  //else if(h_name[i].find("noPr")!= std::string::npos)fixNumber=(1.34975e+06*2.8913e+06/2.7067e+07)/(1.75969e+06*2.0524e+06/2.31872e+07);
  
  if(antiTau21.find("_antiTau21")!= std::string::npos)fixNumber=55008/40530.4;
  else if(antiTau21.find("_antiDBT")!= std::string::npos)fixNumber=1.02528e+06/1.04613e+06;
  else if(antiTau21.find("_mjj_antiDBT")!= std::string::npos)fixNumber=1.30839e+06/1.39172e+06;
  else if(antiTau21.find("_mjj_antiTau21")!= std::string::npos)fixNumber=1.30839e+06/1.39172e+06;
  else if(antiTau21.find("_mjj")!= std::string::npos)fixNumber=1.30839e+06/1.39172e+06;
  
 //fixNumber*=(36813.0/35900);
double Xsec[5]={32100,6831,1207,119.9,25.24};
  //string categories[4]={"0b","1b","2b","DSV"};
  string hadflv[4]={"bb","b","cc","udcsg"};
  
  cout<<"fixNumber="<<fixNumber<<endl;
  
  TString endfix;
	endfix=gSystem->GetFromPipe(Form("file=%s; test=${file%%*_}; echo \"${test}\"",h_name[i].data()));
	
  	
  
  for(int k=0;k<5;k++){
	  for(int j=0;j<4;j++){
		  
		  TH1D *th2=(TH1D* )tf1[k]->FindObjectAny("fixScale");
		  thh[k][j]=(TH1D*)tf1[k]->FindObjectAny(Form("%s_%s%s",h_name[i].data(),hadflv[j].data(),antiTau21.data()));
		  thh[k][j]->Scale(fixNumber* integratedLumi*Xsec[k]/th2->GetBinContent(1));
		 cout<<k<<","<<j<<","<<Form("%s_%s",h_name[i].data(),hadflv[j].data())<<","<<thh[k][j]->Integral()<<","<<thh[k][j]->GetNbinsX()<<endl;
	  }
	  
	  if(k==4){
		  for(int j=0;j<4;j++)thh[k][j]->Add(thh[0][j]);
		  for(int j=0;j<4;j++)thh[k][j]->Add(thh[1][j]);
		  for(int j=0;j<4;j++)thh[k][j]->Add(thh[2][j]);
		  for(int j=0;j<4;j++)thh[k][j]->Add(thh[3][j]);
	  }
  }
  cout<<"here"<<endl;
	vector<TH1D* > v2;
	vector<TH1D* > vd;
	string ntvxtemp=h_name[i];
	
	
	thh[5][0]=(TH1D*)tf1[5]->FindObjectAny(Form("%s_udcsg%s",h_name[i].data(),antiTau21.data()));
	//cout<<thh[4][0]->Integral();
	
	h_name[i]=ntvxtemp;
	double temp_scale=0;
	for(int i=0;i<4;i++){
		temp_scale+=thh[4][i]->Integral();
	}
	cout<< temp_scale<<",data="<<thh[5][0]->Integral()<<endl;
	temp_scale=thh[5][0]->Integral()/temp_scale;
	cout<<"facor="<<temp_scale<<endl;
	if(h_name[i].find("svMass0")!= std::string::npos)for(int i=0;i<4;i++)thh[4][i]->Scale(temp_scale);
	//else if(antiTau21.find("_antiTau21")!= std::string::npos)for(int i=0;i<4;i++)thh[4][i]->Scale(temp_scale);
	if(h_name[i].find("pt_")!= std::string::npos||h_name[i].find("eta_")!= std::string::npos){
		for(int g=0;g<4;g++)thh[4][g]->Rebin(2);
	}
	v2.push_back(thh[4][0]);
	v2.push_back(thh[4][1]);
	v2.push_back(thh[4][2]);
	v2.push_back(thh[4][3]);
	
	
	vd.push_back(thh[5][0]);
	
	//--------
	
	TH1D *h_bkg  = (TH1D* )thh[3][0]->Clone("h_bkg");
	h_bkg->Clear();
	//TH1D *temp = (TH1D* )th1[0]->Clone("h_bkg");
	
	
    TH1D *h_data = (TH1D* )thh[4][0]->Clone("h_data");
    h_data->Clear();
	 h_data->SetTitle("");
	 h_bkg->SetTitle("");
	 c_up->SetTitle("");
	 c_up->cd();
	
	 if(h_name[i].find("pt")!= std::string::npos){
		 rangeUserUp=800;
		 rangeUserDown=200;
		 isSetRange=1;	
	 }	 
	 if(h_name[i].find("total")!= std::string::npos){
		 isSetRange=1;
rangeUserUp=1700;
		 rangeUserDown=500;		 
	 }
	

	
	//vd.push_back(th1[4]);
	cout<<h_name[i]<<endl;
	
	
	
	string tempName=h_name[i];
	if(h_name[i].find("eta_")!= std::string::npos){
		TString endfix2;
		 endfix2=gSystem->GetFromPipe(Form("file=%s; test=${file##Eta}; echo \"${test}\"",h_name[i].data()));
		 h_name[i]=Form("#eta%s",endfix2.Data());
	}
	
	if(h_name[i].find("tau")!= std::string::npos){
		TString endfix2;
		 endfix2=gSystem->GetFromPipe(Form("file=%s; test=${file##tau}; echo \"${test}\"",h_name[i].data()));
		 h_name[i]=Form("#tau%s",endfix2.Data());
	}
	
	if (h_name[i].find("total")!= std::string::npos)h_name[i]="M_{jj}[GeV]";
	
	else if (h_name[i].find("_j0")!= std::string::npos){
		TString endfix2;
		 endfix2=gSystem->GetFromPipe(Form("file=%s; test=${file%%*_j0*}; echo \"${test}\"",h_name[i].data()));
		 //cout<<endfix2<<endl;
		 if (h_name[i].find("_sj0")!= std::string::npos) h_name[i]=Form("%s^{Jet0}_{subjet0}",endfix2.Data());
		 else  if (h_name[i].find("_sj1")!= std::string::npos) h_name[i]=Form("%s^{Jet0}_{subjet1}",endfix2.Data());
		 else h_name[i]=Form("%s^{Jet0}",endfix2.Data());
	}
	
	else if (h_name[i].find("_j1")!= std::string::npos){
		TString endfix2;
		 endfix2=gSystem->GetFromPipe(Form("file=%s; test=${file%%*_j1*}; echo \"${test}\"",h_name[i].data()));
		 //cout<<endfix2<<endl;
		 if (h_name[i].find("_sj0")!= std::string::npos) h_name[i]=Form("%s^{Jet1}_{subjet0}",endfix2.Data());
		 else  if (h_name[i].find("_sj1")!= std::string::npos) h_name[i]=Form("%s^{Jet1}_{subjet1}",endfix2.Data());
		 else h_name[i]=Form("%s^{Jet1}",endfix2.Data());
	}
	
    else if(h_name[i].find("pt")!= std::string::npos||
	//h_name[i].find("total")!= std::string::npos||
	h_name[i].find("ass")!= std::string::npos)	h_name[i]=Form("%s[GeV]",endfix.Data());
	
	else h_name[i]=Form("%s",endfix.Data());
	
	
	
    myPlot(v2,
vd,
	   h_data, h_bkg,1,h_name[i]);

    //c_up->RedrawAxis();
   
    h_name[i]=tempName;
    
    
    if(drawSignal){
	    int colorNum[4]={kRed,kBlue,kOrange+2,kViolet};
	    double xsec2[4]={1.9,0.155,0.0158,17.3};
		TH1D* th_signal[4];
		
		TLegend *leg = new TLegend(0.64, 0.68, 0.79, 0.87);
  string legendS[4]={"M_{G}=1.4TeV","M_{G}=1.8TeV","M_{G}=2.5TeV","M_{R}=1.8TeV"};
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
		for(int k=0;k<3;k++){
			
	    th_signal[k]=(TH1D*)tf2[k]->FindObjectAny(Form("%s",h_name[i].data()));
	    leg->AddEntry(th_signal[k],Form("%s",legendS[k].data()));
	    th_signal[k]->SetLineColor(colorNum[k]);
	    th_signal[k]->SetLineWidth(2);
	   // th_signal[k]->SetLineStyle(k+1);
	    
	    
	     TString endfix;
		  if(h_name[i].find("4b")!= std::string::npos){
		  endfix=gSystem->GetFromPipe(Form("file=%s; test=${file%%*_4b}; echo \"${test}\"",h_name[i].data()));
			  TH1D *th2;
			
			  th2=(TH1D* )tf2[k]->FindObjectAny(Form("%s_0b",endfix.Data()));
 			   th_signal[k]->Add(th2);
			  th2=(TH1D* )tf2[k]->FindObjectAny(Form("%s_1b",endfix.Data()));
 			   th_signal[k]->Add(th2);
			  th2=(TH1D* )tf2[k]->FindObjectAny(Form("%s_2b",endfix.Data()));
 			   th_signal[k]->Add(th2);
			   th2=(TH1D* )tf2[k]->FindObjectAny(Form("%s_3b",endfix.Data()));
 			   th_signal[k]->Add(th2);
		  }
	    
	     TH1D *th3=(TH1D* )tf2[k]->FindObjectAny("cutflow");
	     //cout<<"k="<<k<<","<<th2->GetBinContent(1)<<","<<12.883846147301*xsec2[k]*200/th2->GetBinContent(1)<<endl;
	     //th_signal[k]->Scale(thh[4][0]->Integral()/(th_signal[k]->Integral()*2));
	     th_signal[k]->Scale(integratedLumi* 50/th3->GetBinContent(1));
	     /*
	     if(h_name[i].find("deltaR")!= std::string::npos)th_signal[k]->Scale(0.5);
	     if(h_name[i].find("noPr")!= std::string::npos &&
	     h_name[i].find("noPr_")== std::string::npos&&
	     h_name[i].find("tau21")== std::string::npos)th_signal[k]->Scale(0.4);
	    */
	    //if(k<3)th_signal[k]->Scale(10);
	     cout<<th_signal[k]->Integral()<<",";
	     th_signal[k]->Draw("same hist ");
		}
		leg->Draw("same");
		
    }
    
    
    c_dw->cd();
 tempName=h_name[i];
	if(h_name[i].find("Eta")!= std::string::npos){
		TString endfix2;
		 endfix2=gSystem->GetFromPipe(Form("file=%s; test=${file##Eta}; echo \"${test}\"",h_name[i].data()));
		 h_name[i]=Form("#eta%s",endfix2.Data());
	}
	
	if(h_name[i].find("tau")!= std::string::npos){
		TString endfix2;
		 endfix2=gSystem->GetFromPipe(Form("file=%s; test=${file##tau}; echo \"${test}\"",h_name[i].data()));
		 h_name[i]=Form("#tau%s",endfix2.Data());
	}
	
	if (h_name[i].find("total")!= std::string::npos)myRatio(h_data, h_bkg,"M_{jj}[GeV]");
	
	else if (h_name[i].find("_j0")!= std::string::npos){
		TString endfix2;
		 endfix2=gSystem->GetFromPipe(Form("file=%s; test=${file%%*_j0*}; echo \"${test}\"",h_name[i].data()));
		 //cout<<endfix2<<endl;
		 if (h_name[i].find("_sj0")!= std::string::npos) myRatio(h_data, h_bkg,Form("%s^{Jet0}_{subjet0}",endfix2.Data()));
		 else  if (h_name[i].find("_sj1")!= std::string::npos) myRatio(h_data, h_bkg,Form("%s^{Jet0}_{subjet1}",endfix2.Data()));
		 else myRatio(h_data, h_bkg,Form("%s^{Jet0}",endfix2.Data()));
	}
	
	else if (h_name[i].find("_j1")!= std::string::npos){
		TString endfix2;
		 endfix2=gSystem->GetFromPipe(Form("file=%s; test=${file%%*_j1*}; echo \"${test}\"",h_name[i].data()));
		 //cout<<endfix2<<endl;
		 if (h_name[i].find("_sj0")!= std::string::npos) myRatio(h_data, h_bkg,Form("%s^{Jet1}_{subjet0}",endfix2.Data()));
		 else  if (h_name[i].find("_sj1")!= std::string::npos) myRatio(h_data, h_bkg,Form("%s^{Jet1}_{subjet1}",endfix2.Data()));
		 else myRatio(h_data, h_bkg,Form("%s^{Jet1}",endfix2.Data()));
	}
	
    else if(h_name[i].find("Pt")!= std::string::npos||
	//h_name[i].find("total")!= std::string::npos||
	h_name[i].find("ass")!= std::string::npos)myRatio(h_data, h_bkg,Form("%s[GeV]",endfix.Data()));
	else if (h_name[i].find("Nbtagje")!= std::string::npos)myRatio(h_data, h_bkg,"Nbtagjet",1);
	
	else myRatio(h_data, h_bkg,Form("%s",endfix.Data()));
	
	
	h_name[i]=tempName;
	
	
   // c.Draw();
  
	c.SaveAs(Form("dataMC_trig%s/%s.png",antiTau21.data(),h_name[i].data()));
	c.SaveAs(Form("dataMC_trig%s/%s.pdf",antiTau21.data(),h_name[i].data()));
	c.SetName(Form("%s",h_name[i].data()));
	c.SetTitle(Form("%s",h_name[i].data()));
	c.Write();
  }
  
  outFile->Close();
}