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
 #define fixNumberGlobal 4.49681/5.4227

double rangeUserUp=0,rangeUserDown=0;
int isSetRange=0;

void myPlot(vector< TH1D*> h_Z,
           vector< TH1D*> v_data,
	    TH1D* h_data, TH1D* h_bkg,int option=0,string xtitle=""){

	    gStyle->SetNdivisions(605, "XYZ");
  h_data->Reset();
  for(unsigned int i=0;i<v_data.size();i++)h_data->Add(v_data[i]);
  
  TLegend *leg = new TLegend(0.84, 0.68, 0.92, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

   h_bkg->Reset();
    THStack *h_stack = new THStack("h_stack", "");
	int colorNumQCD[4]={46,15,38,30};
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
  h_data->SetMarkerSize(1.5);
  h_data->GetYaxis()->SetTitleOffset(2);
  h_data->GetXaxis()->SetTitle("");
 // h_stack->GetXaxis()->SetTitle("?");
  h_data->GetXaxis()->SetLabelOffset(999);
  h_data->GetXaxis()->SetLabelSize(0);
  
  h_data->SetTitle("");
 //h_data->SetMaximum(3000);
 h_data->SetMarkerSize(1.6);

 h_data->GetYaxis()->SetTitleOffset(1.3);
 h_data->GetZaxis()->SetTitleOffset(0.75);
 // size of axis labels
 h_data->GetXaxis()->SetTitleSize(0.04);
 h_data->GetYaxis()->SetTitleSize(0.04);
 h_data->GetZaxis()->SetTitleSize(0.035);
 h_data->GetXaxis()->SetLabelSize(0.04);
 h_data->GetYaxis()->SetLabelSize(0.04); 
 h_data->GetZaxis()->SetLabelSize(0.025);
  
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
	rangeUserUp=(bmax-0.5)*width+h_data->GetBinCenter(1);
	//rangeUserDown=(bmin-0.5)*width+h_data->GetBinCenter(1);
	rangeUserDown=0;
	if(isSetRange)h_data->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
	//if(isSetRange)h_stack->GetXaxis()->SetRangeUser(rangeUserDown,rangeUserUp);
	//h_stack->GetXaxis()->SetRangeUser((bmin-0.5)*width+h_data->GetBinCenter(1),(bmax-0.5)*width+h_data->GetBinCenter(1));

  
	  h_stack->SetTitle(h_data->GetTitle());
	if(xtitle.find("Mass")!= std::string::npos &&xtitle.find("total")== std::string::npos
	&& xtitle.find("no")== std::string::npos)h_stack->SetMaximum(h_stack->GetMaximum()*1.4);
	//if(xtitle.find("double")!= std::string::npos)h_stack->GetHistogram()SetMaximum(h_stack->GetHistogram()GetMaximum()*1.3);
    h_stack->Draw("histe");
	 
     h_stack->SetMaximum( h_stack->GetMaximum()*1.3);
	
   // h_stack->GetHistogram()->GetYaxis()->SetTitle("Event Numbers");
    if(h_data->GetBinWidth(1)>1)h_stack->GetHistogram()->GetYaxis()->SetTitle(Form("Events / %2.0f units",h_data->GetBinWidth(1)));
    else h_stack->GetHistogram()->GetYaxis()->SetTitle(Form("Events / %0.2f units",h_data->GetBinWidth(1)));
		if(xtitle.find("pt")!= std::string::npos|| xtitle.find("SD")!= std::string::npos||xtitle.find("jj")!= std::string::npos)h_stack->GetHistogram()->GetYaxis()->SetTitle(Form("Events / %2.0f GeV",h_data->GetBinWidth(1)));
    h_stack->GetHistogram()->GetXaxis()->SetTitle(Form("%s",xtitle.data()));
	if(xtitle.find("eta")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("#eta (j1)");
	if(xtitle.find("eta")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("#eta (j2)");
	if(xtitle.find("pt")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("p_{T} (j1) [GeV]");
	if(xtitle.find("pt")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("p_{T} (j2) [GeV]");
	if(xtitle.find("SD")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("PUPPI soft-drop mass(j1)");
	if(xtitle.find("Tau")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("#tau_{21} (j1)");
	if(xtitle.find("Tau")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("#tau_{21} (j2)");
	if(xtitle.find("SD")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("PUPPI soft-drop mass(j2)");
	if(xtitle.find("delta")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("#Delta#eta (j1,j2)");
	if(xtitle.find("double")!= std::string::npos&& xtitle.find("0")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("double-b tagger (j1)");
	if(xtitle.find("double")!= std::string::npos&& xtitle.find("1")!= std::string::npos)h_stack->GetHistogram()->GetXaxis()->SetTitle("double-b tagger (j2)");
	
    h_stack->GetHistogram()->GetYaxis()->SetTitleSize(h_data->GetYaxis()->GetTitleSize());
    h_stack->GetHistogram()->GetYaxis()->SetLabelSize(h_data->GetYaxis()->GetLabelSize());
    h_stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.35);
     h_stack->GetHistogram()->GetXaxis()->SetTitleSize(0.045);
     h_stack->GetHistogram()->GetXaxis()->SetTitleFont(62);
	 
	 h_stack->GetHistogram()->SetTitle("");
 //h_data->SetMaximum(3000);
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
    // if(xtitle.find("double")!= std::string::npos)h_stack->GetHistogram()->GetYaxis()->SetRangeUser(0,h_stack->GetHistogram()GetMaximum()*1.3);
     if(xtitle.find("{jj}")!= std::string::npos )h_stack->GetHistogram()->GetXaxis()->SetRangeUser(500,3000);
    //h_stack->GetHistogram()->GetXaxis()->SetTickLength(0);
    //h_stack->GetHistogram()->GetXaxis()->SetLabelOffset(999);
//h_data->Draw("elsame");
  
  
   

    
  
  //leg->AddEntry(h_data, "Data", "lp");
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
		lar->DrawLatex(0.2, 0.8, "CMS");
		lar2->DrawLatex(0.31, 0.77, "#it{#bf{Simulation Preliminary}} ");

}


void MCplots_v2(){

  //setNCUStyle(true);
 
  Float_t up_height     = 0.8;
  Float_t dw_correction = 1.455;
  Float_t dw_height     = 0;

  
 // c.Divide(1,2);

 // TPad* c_up = (TPad*) c.GetListOfPrimitives()->FindObject("c_1");
  //TPad* c_dw = (TPad*) c.GetListOfPrimitives()->FindObject("c_2"); 

  //c_up->SetPad(0,1-up_height,1,1);
  //c_dw->SetPad(0,0,1,dw_height);
//c_dw->SetBottomMargin(0.25);

  // To get the name of histograms
  
  TFile* tf1[10];
  tf1[0]=TFile::Open("QCD700.root");
  tf1[1]=TFile::Open("QCD1000.root");
  tf1[2]=TFile::Open("QCD1500.root");
  tf1[3]=TFile::Open("QCD2000.root");
  tf1[4]=TFile::Open("QCD500.root");
//  tf1[4]=TFile::Open("sf2/data.root");
  
  TFile* tf2[10];
   tf2[0]=TFile::Open("1.root");
   tf2[1]=TFile::Open("2.root");
   tf2[2]=TFile::Open("3.root");
   tf2[3]=TFile::Open("rd1400.root");
   tf2[4]=TFile::Open("rd1800.root");
   tf2[5]=TFile::Open("rd2500.root");
//   tf2[3]=TFile::Open("sf2/R1800.root");
   
  vector<std::string> h_name;
  h_name.push_back("nVtx");
 h_name.push_back("pt_j0");
 h_name.push_back("eta_j0");
 h_name.push_back("puppiSDMassThea_j0");
 h_name.push_back("puppiTau21_j0");
 h_name.push_back("doubleSV_j0");
h_name.push_back("pt_j1");
h_name.push_back("eta_j1");
h_name.push_back("puppiSDMassThea_j1");
 h_name.push_back("puppiTau21_j1");
 h_name.push_back("doubleSV_j1");
 h_name.push_back("deltaEta");
 h_name.push_back("totalMass");
 h_name.push_back("totalMassRed");
 
 TFile* outFile= new TFile("MCcompare.root","recreate");

  for(unsigned int i = 0; i < h_name.size(); i++){
	  TCanvas c("c","",0,0,600*1.3,600);
	  
	  gStyle->SetOptStat(0);
 gStyle->SetPaintTextFormat("2.1f");
 gStyle->SetPalette(57);
 gStyle->SetFrameLineWidth(3);
 gStyle->SetPadRightMargin(0.059);
 gStyle->SetPadLeftMargin(0.159);
	  
//	 for(unsigned int i = 0; i < 1; i++){
	  
	//cout<<h_name[i]<<endl;
	bool drawSignal=1;
	
	//---------
	 TH1D* thh[6][4];
  double fixNumber=22295/30584.3;//7306/11857.3;
double Xsec[5]={6831,1207,119.9,25.24,32100};
  //string categories[4]={"0b","1b","2b","DSV"};
  string hadflv[4]={"bb","b","cc","udcsg"};
  
  fixNumber*=(36813.0/35900);
  
  TString endfix;
	endfix=gSystem->GetFromPipe(Form("file=%s; test=${file%%*_}; echo \"${test}\"",h_name[i].data()));
	
  
  
  for(int k=0;k<5;k++){
	  for(int j=0;j<4;j++){
		  
		  TH1D *th2=(TH1D* )tf1[k]->FindObjectAny("fixScale");
		  
	
		  thh[k][j]=(TH1D*)tf1[k]->FindObjectAny(Form("%s_%s",h_name[i].data(),hadflv[j].data()));
		   TString endfix;
		
		  thh[k][j]->Scale(fixNumber* integratedLumi*Xsec[k]/th2->GetBinContent(1));
		 // cout<<k<<","<<j<<","<<Form("%s_%s",h_name[i].data(),hadflv[j].data())<<","<<thh[k][j]->Integral()<<endl;
	  }
	  
	  if(k==4){
		  for(int j=0;j<4;j++)thh[k][j]->Add(thh[0][j]);
		  for(int j=0;j<4;j++)thh[k][j]->Add(thh[1][j]);
		  for(int j=0;j<4;j++)thh[k][j]->Add(thh[2][j]);
		  for(int j=0;j<4;j++)thh[k][j]->Add(thh[3][j]);
	  }
  }
vector<TH1D* > v2;
	vector<TH1D* > vd;
	
	thh[5][0]=(TH1D*)tf1[3]->FindObjectAny(Form("%s_udcsg",h_name[i].data()));
	
	cout<<thh[4][0]->Integral();
	thh[4][0]->SetName("qcd_bb");
	thh[4][1]->SetName("qcd_b");
	thh[4][2]->SetName("qcd_cc");
	thh[4][3]->SetName("qcd_udcsg");
	
	v2.push_back(thh[4][0]);
	v2.push_back(thh[4][1]);
	v2.push_back(thh[4][2]);
	v2.push_back(thh[4][3]);
	
	
	
	double temp_scale=0;
	for(int i=0;i<4;i++){
		temp_scale+=thh[3][i]->Integral();
		
	}
	
	for(int i=0;i<4;i++){
		if(h_name[i].find("HT")!= std::string::npos)continue;
		//thh[3][i]->Scale(thh[4][0]->Integral()/temp_scale);
	}
	
	vd.push_back(thh[5][0]);
	
	//--------
	
	TH1D *h_bkg  = (TH1D* )thh[3][0]->Clone("h_bkg");
	h_bkg->Clear();
	//TH1D *temp = (TH1D* )th1[0]->Clone("h_bkg");
	
	
    TH1D *h_data = (TH1D* )thh[4][0]->Clone("h_data");
    h_data->Clear();
    /*
	h_data->SetTitle(Form("%s",endfix.Data()));
	h_bkg->SetTitle(Form("%s",endfix.Data()));
	c_up->SetTitle(Form("%s",endfix.Data()));
	c.SetTitle(Form("%s",endfix.Data()));
	 */
	 h_data->SetTitle("");
	 h_bkg->SetTitle("");
	// c_up->SetTitle("");
	// c_up->cd();
	 
	 //if (h_name[i].find("Nbtagjet")!= std::string::npos)c_up->SetLogy();	
	// if(h_name[i].find("cutflow")!= std::string::npos)c_up->SetLogy();	
	 if(h_name[i].find("Pt")!= std::string::npos)isSetRange=1;	
	 if(h_name[i].find("logPt")!= std::string::npos){
		   isSetRange=0;	
		  // c_up->SetLogy(1);	
	 }
	 
	 
	  if(h_name[i].find("total")!= std::string::npos)isSetRange=1;	
	
	
	double scaleTemp[4];
	
	
	
	
	//vd.push_back(th1[4]);
	cout<<h_name[i]<<endl;
	
	
	
	string tempName=h_name[i];
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
	if (h_name[i].find("Red")!= std::string::npos)h_name[i]="M^{reduced}_{jj}[GeV]";
	else if (h_name[i].find("total")!= std::string::npos)h_name[i]="M_{jj}[GeV]";
	
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
	
    else if(h_name[i].find("Pt")!= std::string::npos||
	//h_name[i].find("total")!= std::string::npos||
	h_name[i].find("ass")!= std::string::npos)	h_name[i]=Form("%s[GeV]",endfix.Data());
	
	else h_name[i]=Form("%s",endfix.Data());
	
	
	
	
	
	
    myPlot(v2,
vd,
	   h_data, h_bkg,1,h_name[i]);

    //c_up->RedrawAxis();
    
    h_name[i]=tempName;
    
    
    if(drawSignal){
	    int colorNum[4]={kRed,kBlue,kOrange+2,};
	    int colorNumRd[4]={kViolet,kGreen+3,kOrange+2,kViolet};
	    double xsec2[4]={ 0.057*0.8+0.559*0.2, 0.018, 1.86E-03,17.3};
		TH1D* th_signal[4][4];
		
		TLegend *leg = new TLegend(0.63, 0.68, 0.78, 0.87);
  string legendS[4]={"M_{G}=1.4TeV","M_{G}=1.8TeV","M_{G}=2.5TeV","M_{R}=1.8TeV"};
  string legendSRd[4]={"M_{R}=1.4TeV","M_{R}=1.8TeV","M_{R}=2.5TeV","M_{R}=1.8TeV"};
  string legendSName[4]={"bg_1400","bg_1800","bg_2500","M_{R}=1.8TeV"};
  string legendSNameRd[4]={"rd_1400","rd_1800","rd_2500","M_{R}=1.8TeV"};
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
		for(int k=0;k<3;k++){
			if(k==1)continue;
	       for(int j=0;j<4;j++)  th_signal[k][j]=(TH1D*)tf2[k]->FindObjectAny(Form("%s_%s",h_name[i].data(),hadflv[j].data()));
	  th_signal[k][3]->Add(th_signal[k][2]);
	  th_signal[k][3]->Add(th_signal[k][1]);
	  th_signal[k][3]->Add(th_signal[k][0]);
	    th_signal[k][3]->SetName(Form("%s",legendSName[k].data()));
	    
	    leg->AddEntry(th_signal[k][3],Form("%s",legendS[k].data()));
	    th_signal[k][3]->SetLineColor(colorNumRd[k]);
	    th_signal[k][3]->SetLineWidth(2);
	   // th_signal[k]->SetLineStyle(k+1);
	   
	     TH1D *th3=(TH1D* )tf2[k]->FindObjectAny("fixScale");
	     //cout<<"k="<<k<<","<<th2->GetBinContent(1)<<","<<12.883846147301*xsec2[k]*200/th2->GetBinContent(1)<<endl;
	     //th_signal[k]->Scale(thh[4][0]->Integral()/(th_signal[k]->Integral()*2));
	     th_signal[k][3]->Scale(integratedLumi*20/th3->GetBinContent(1));
	  
	     cout<<th_signal[k][3]->Integral()<<",";
	     th_signal[k][3]->Draw("same hist ");
		}
		
		for(int k=0;k<3;k++){
			if(k==1)continue;
	       for(int j=0;j<4;j++)  th_signal[k][j]=(TH1D*)tf2[k+3]->FindObjectAny(Form("%s_%s",h_name[i].data(),hadflv[j].data()));
	  th_signal[k][3]->Add(th_signal[k][2]);
	  th_signal[k][3]->Add(th_signal[k][1]);
	  th_signal[k][3]->Add(th_signal[k][0]);
	    th_signal[k][3]->SetName(Form("%s",legendSNameRd[k].data()));
	    
	    leg->AddEntry(th_signal[k][3],Form("%s",legendSRd[k].data()));
	    th_signal[k][3]->SetLineColor(colorNumRd[k]);
	    th_signal[k][3]->SetLineStyle(7);
	    th_signal[k][3]->SetLineWidth(2);
	 //   th_signal[k]->SetLineStyle(k+1);
	   
	     TH1D *th3=(TH1D* )tf2[k]->FindObjectAny("fixScale");
	     //cout<<"k="<<k<<","<<th2->GetBinContent(1)<<","<<12.883846147301*xsec2[k]*200/th2->GetBinContent(1)<<endl;
	     //th_signal[k]->Scale(thh[4][0]->Integral()/(th_signal[k]->Integral()*2));
	     th_signal[k][3]->Scale(integratedLumi*20/th3->GetBinContent(1));
	  
	     cout<<th_signal[k][3]->Integral()<<",";
	     th_signal[k][3]->Draw("same hist ");
		}
		leg->Draw("same");
		
    }
    
    
    
    
     
		  
	
		
/*
    
    c_dw->cd();
	string tempName=h_name[i];
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
	*/
	
	if(h_name[i].find("Pt")!= std::string::npos)isSetRange=0;	
	 if(h_name[i].find("total")!= std::string::npos)isSetRange=0;	
   // c.Draw();
  
	if(h_name[i].find("0b")!= std::string::npos)c.SaveAs(Form("dataMC_v2/0b/%s.png",h_name[i].data()));
	else if(h_name[i].find("1b")!= std::string::npos)c.SaveAs(Form("dataMC_v2/1b/%s.png",h_name[i].data()));
	else if(h_name[i].find("2b")!= std::string::npos)c.SaveAs(Form("dataMC_v2/2b/%s.png",h_name[i].data()));
    else if(h_name[i].find("4b")!= std::string::npos)c.SaveAs(Form("dataMC_v2/all/%s.png",endfix.Data()));
	else c.SaveAs(Form("dataMC_v2/all/%s.png",h_name[i].data()));
	
	if(h_name[i].find("0b")!= std::string::npos)c.Print(Form("dataMC_v2/0b/%s.pdf",h_name[i].data()));
	else if(h_name[i].find("1b")!= std::string::npos)c.Print(Form("dataMC_v2/1b/%s.pdf",h_name[i].data()));
	else if(h_name[i].find("2b")!= std::string::npos)c.Print(Form("dataMC_v2/2b/%s.pdf",h_name[i].data()));
    else if(h_name[i].find("4b")!= std::string::npos)c.Print(Form("dataMC_v2/all/%s.pdf",endfix.Data()));
	else c.Print(Form("MC_N1/%s.pdf",h_name[i].data()));
	
	c.SetName(Form("%s",h_name[i].data()));
	c.SetTitle(Form("%s",h_name[i].data()));
	c.Write();
	
	if(h_name[i].find("logPt")!= std::string::npos){
		   isSetRange=0;	
		   //c_up->SetLogy(0);	
	 }
  }
  
  outFile->Close();
}