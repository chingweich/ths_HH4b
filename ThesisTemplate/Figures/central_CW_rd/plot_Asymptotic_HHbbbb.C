#include <Riostream.h>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <cstdlib>
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "setNCUStyle.C"
//#define nXm 6
//#define nXm_thr 6
#define nXm 7                                                                                         
#define nXm_thr 6

const float intLumi = 27.6;
const string dirXSect = "./";

// 10 = radion_2cat, 11 = radion_4btag_cat0, 12 = radion_3btag_HPHP_cat1 
// 20 = graviton_2cat, 21 = graviton_4btag_cat0, 22 = graviton_3btag_HPHP_cat1 
string xTitle;
string sCat;
string sXsec;
string sThTitle;
double BFXHH = 0;
double BFHH4b = 0;

void plot_Asymptotic_HHbbbb();
void setFPStyle();
void scaleGraph(TGraphAsymmErrors* g, double factor)
{
  int npoints = g->GetN();
  for(int i=0; i!=npoints; ++i) {
    double x = g->GetX()[i];
    double y = g->GetY()[i];
    double eyh = g->GetEYhigh()[i];
    double eyl = g->GetEYlow()[i];
    y = (y*factor);
    eyh = (eyh*factor);
    eyl = (eyl*factor);
    g->SetPoint(i,x,y);
    g->SetPointEYhigh(i, eyh);
    g->SetPointEYlow(i, eyl);
  }

}

double expo_interp(double s2, double s1,  double newM, double m2, double m1)
{
  if (m1 > m2) {
    double tmps = s1;
    double tmpm = m1;
    s1 = s2;
    m1 = m2;
    s2 = tmps;
    m2 = tmpm;
  }
  double deltaM = m2 - m1;
  double alpha = (log(s2) - log(s1)) / deltaM;
  double newS = s1 * pow(exp(newM - m1), alpha);
  return newS;
}



double linear_interp(double s2, double s1, double mass, double m2, double m1)
{
  if (m1 > m2) {
    double tmps = s1;
    double tmpm = m1;
    s1 = s2;
    m1 = m2;
    s2 = tmps;
    m2 = tmpm;
  }
  return (s1 + (s2 - s1) * (mass - m1) / (m2 - m1));
}



void plot_Asymptotic_HHbbbb()
{
	
	TStyle* ts=setNCUStyle(0);
ts->SetTitleSize(0.045,"XYZ");
ts->SetTitleOffset(1.3, "Y");
ts->SetLabelSize(0.04, "XYZ");

ts->SetPadGridX(false);
  ts->SetPadGridY(false);
  ts->SetGridColor(0);
  ts->SetGridStyle(3);
  ts->SetGridWidth(1);
   ts->SetNdivisions(605, "XYZ");
  int sigHyp=10,  subtr=0;

  if (sigHyp == 10 || sigHyp == 11 || sigHyp == 12) {
    xTitle = string("M_{R} [GeV]");
    sXsec = string("radion_toHH_toBB_lambda1.txt");
    sThTitle = string("Radion (#Lambda_{R} = 3TeV)");
    BFXHH = 1;

    BFHH4b = 1;
    //    sXsec = "radion_lambda1";
  }
  else if (sigHyp == 20 || sigHyp == 21 || sigHyp == 22) {
    xTitle = string("M_{Gkk} [GeV]");
    sXsec = string("bulk_graviton_toHH_toBB_kmpl05.txt");
    sThTitle = string("G_{Bulk} (#kappa/\bar{M_{Pl}} = 0.5)");
    BFXHH = 1;
    BFHH4b = 1;
  }

  if (subtr == 0){
    if (sigHyp == 10) sCat = string("Radion");
    else if (sigHyp == 11) sCat = string("Radion_4btag_cat0");
    else if (sigHyp == 12) sCat = string("Radion_3btag_cat1");
    
    else if (sigHyp == 20) sCat = string("Graviton");
    else if (sigHyp == 21) sCat = string("Graviton_4btag_cat0");
    else if (sigHyp == 22) sCat = string("Graviton_3btag_cat1");
  } else {
    if (sigHyp == 10) sCat = string("Radion_subtr");
    else if (sigHyp == 11) sCat = string("Radion_subtr_4btag_cat0");
    else if (sigHyp == 12) sCat = string("Radion_subtr_3btag_cat1");
    
    else if (sigHyp == 20) sCat = string("Graviton_subtr");
    else if (sigHyp == 21) sCat = string("Graviton_subtr_4btag_cat0");
    else if (sigHyp == 22) sCat = string("Graviton_subtr_3btag_cat1");

  }

  //TString outfilename = TString(outputdir.c_str())+".root";
  //  TFile *fout = new TFile(outfilename,"RECREATE");
  TFile *fout = new TFile("limit.root","RECREATE");  
bool useNewStyle = true;
  if (useNewStyle)  setFPStyle();
//  gROOT->LoadMacro("CMS_lumi.C");
 
  TFile *fFREQ[nXm];
  TTree *t[nXm];
  //  int Xmass[nXm]={1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500};  
  //int Xmass_thr[nXm_thr]={1200,1400,1600,1800,2000,2500};
  //int Xmass[nXm]={1200,1400,1600,1800,2000,2500};
  int Xmass_thr[nXm_thr]={1000,1500,1800,2000,2500,3000};  
  int Xmass[nXm]={1200,1400,1600,1800,2000,2500,3000};  
  vector<double> v_mh, v_median, v_68l, v_68h, v_95l, v_95h, v_obs;
 
 
  for(int n=0;n<nXm;n++)
  {
    char limitfilename[100];

    string sDatacard("higgsCombine");
    sDatacard = sDatacard + "" + sCat + Form("_bumphunt.Asymptotic.mH%d.root",Xmass[n]);

    TString limitfile = sDatacard.c_str();
    cout<<" Read limit file: "<<limitfile << " and Xmass = " << Xmass[n]<<endl;
    //fFREQ[n] = new TFile(limitfile, "READ");
    fFREQ[n] = new TFile(Form("CMS_%d_HH4b_13TeV_asymptoticCLs.root",Xmass[n]), "READ");
    //higgsCombineTest.Asymptotic.mH1000.root

    t[n] = (TTree*)fFREQ[n]->Get("limit");
  
    double mh, limit;
    float quant;
    t[n]->SetBranchAddress("mh", &mh);
    t[n]->SetBranchAddress("limit", &limit);
    t[n]->SetBranchAddress("quantileExpected", &quant);
  
    
    
    //int iMH = 0;
    //while (iMH < n) {
 
      for (int i = 0; i < t[n]->GetEntries(); i++) {

        t[n]->GetEntry(i);

	//        cout<<" quant : "<<quant<<" limit : " <<limit<<endl;
        /// Map: mh --> observed, 95low, 68low, expected, 68hi, 95hi, xsec
        if (quant > -1.01 && quant < -0.99) {
        v_obs.push_back(limit*10);
        } 
        else if (quant > 0.02 && quant < 0.03) {
	v_95l.push_back(limit*10);
        }
        else if (quant > 0.15 && quant < 0.17) {
	v_68l.push_back(limit*10);
        }
        else if (quant > 0.49 && quant < 0.51) {
	v_median.push_back(limit*10);
        v_mh.push_back(mh);
        }
        else if (quant > 0.83 && quant < 0.85) {
	v_68h.push_back(limit*10);
        }
        else if (quant > 0.965 && quant < 0.98) {
	v_95h.push_back(limit*10);
        }
        else {
        cout << "Error! Quantile =  " << quant << endl;
        }
     }
      //iMH++;
	//     }//end while loop

  }//file loop

  string xsect_file_th =  sXsec.c_str(); 


  ifstream xsect_file(xsect_file_th.c_str(), ios::in);
  if (! xsect_file.is_open()) {
    cout << "Failed to open file with xsections: " << xsect_file_th << endl;
  }

  float mH, CS;
  vector<float> v_mhxs, v_xs, v_toterrh, v_toterrl;
  while (xsect_file.good()) {
    xsect_file >> mH >> CS;
  
    v_mhxs.push_back(mH);
    v_xs.push_back(CS*BFXHH*BFHH4b);//*BRZZ2l2q (multyply by BRZZ2l2q only if exp rates in cards are for process X->ZZ->2l2q !)
    
    //unavailable theory errors for graviton

    float tot_err_p = 0.0;
    float tot_err_m = 0.0;

    v_toterrh.push_back(1.0 + (tot_err_p));
    v_toterrl.push_back(1.0 - (tot_err_m));
  }
  cout << "Size of theory xsects vector" << v_mhxs.size() << endl;
  xsect_file.close();
  ///////////////////////////
  // END THEORY INPUT PART //
  ///////////////////////////


  /// Here we multiply the limits in terms of signal strength by the cross-section.
  /// There are also some hooks to exclude sick mass points.
  
  double mass[nXm], obs_lim_cls[nXm];
  double medianD[nXm];
  double up68err[nXm], down68err[nXm], up95err[nXm], down95err[nXm];
  int nMassEff = 0;
  
  string sOutput = string("BrazilianFlags/Limits_") + sCat + "_HH.txt";

  ofstream myfile;
  myfile.open (sOutput.c_str());
  myfile << "mX\tobserved\texpected\n";

  for (int im = 0; im < nXm; im++) {

    double fl_xs = 1;
    double fl_xs10 = 0;

    mass[nMassEff] = Xmass[im];

    /// This is the part where we multiply the limits in terms of signal strength
    /// by the cross-section, in order to have limits in picobarns.
    //std::cerr << mass[nMassEff] << ":" << v_obs.at(im) << std::endl;
      obs_lim_cls[nMassEff] = v_obs.at(im) * fl_xs;
 
      medianD[nMassEff] = v_median.at(im) * fl_xs;
      up68err[nMassEff] = (v_68h.at(im) - v_median.at(im)) * fl_xs;
      down68err[nMassEff] = (v_median.at(im) - v_68l.at(im)) * fl_xs;

      up95err[nMassEff] = (v_95h.at(im) - v_median.at(im)) * fl_xs;
      down95err[nMassEff] = (v_median.at(im) - v_95l.at(im)) * fl_xs;
    
      cout<<"i " << im <<" median_lim_cls : " <<medianD[nMassEff] << " obs: " << obs_lim_cls[nMassEff] <<" mass : "<<mass[nMassEff]<<endl;
 
      myfile << Form("%.1f\t%.2f\t%.2f\n", mass[nMassEff],  obs_lim_cls[nMassEff], medianD[nMassEff]);

      nMassEff++;
    
    
  }//end loop over im (mass points)


  myfile.close();

  double xs[nXm_thr], xs_uperr[nXm_thr], xs_downerr[nXm_thr];
  double xs10[nXm_thr], xs10_uperr[nXm_thr], xs10_downerr[nXm_thr];
  double mass_thr[nXm_thr];
  int nMassEff_thr = 0;
  
 for (int im = 0; im < nXm_thr; im++) {
  //scale factor 100 for making the xsect visible

    double fl_xs = double(v_xs.at(im)); //*1000.0
    double fl_xs10 = double(v_xs.at(im))*10;//double(v_xs10.at(ind)); //*1000.0

   mass_thr[nMassEff_thr] = Xmass_thr[im]; 
   xs[nMassEff_thr] =  double(v_xs.at(im)); //*100.0;    
 
  xs_uperr[nMassEff_thr] = double(v_toterrh.at(im)) * xs[nMassEff_thr] - xs[nMassEff_thr];
  xs_downerr[nMassEff_thr] =  xs[nMassEff_thr] - double(v_toterrl.at(im)) * xs[nMassEff_thr];

  xs10[nMassEff_thr] = fl_xs10; //*100.0;
  xs10_uperr[nMassEff_thr] = double(v_toterrh.at(im)) * xs10[nMassEff_thr] - xs10[nMassEff_thr];
  xs10_downerr[nMassEff_thr] =  xs10[nMassEff_thr] - double(v_toterrl.at(im)) * xs10[nMassEff_thr];

  nMassEff_thr++;

 }
  

  /// The TGraphs themselves.

  //cout<<"Working on TGraph"<<endl;
  TGraphAsymmErrors *grobslim_cls = new TGraphAsymmErrors(nMassEff, mass, obs_lim_cls);
  grobslim_cls->SetName("LimitObservedCLs");
  TGraphAsymmErrors *grmedian_cls = new TGraphAsymmErrors(nMassEff, mass, medianD);
  grmedian_cls->SetName("LimitExpectedCLs");
  TGraphAsymmErrors *gr68_cls = new TGraphAsymmErrors(nMassEff, mass, medianD, 0, 0, down68err, up68err);
  gr68_cls->SetName("Limit68CLs");
  TGraphAsymmErrors *gr95_cls = new TGraphAsymmErrors(nMassEff, mass, medianD, 0, 0, down95err, up95err);
  gr95_cls->SetName("Limit95CLs");

  // TGraphAsymmErrors *grthSM=new TGraphAsymmErrors(nMassEff1,mass1,xs,0,0,0,0);//xs_downerr,xs_uperr);
  TGraph *grthSM=new TGraph(nMassEff_thr,mass_thr,xs);//xs_downerr,xs_uperr);
  grthSM->SetName("SMXSection");


  // TGraphAsymmErrors *grthSM10=new TGraphAsymmErrors(nMassEff1,mass1,xs10,0,0,0,0);
  TGraph *grthSM10=new TGraph(nMassEff_thr,mass_thr,xs10);
  grthSM10->SetName("SMXSection_2nd");

  // double fr_left = 590.0, fr_down = 1E-5, fr_right = 2000.0, fr_up = 0.5; 
  // double fr_left = 1190.0, fr_down = 5E-5, fr_right = 2500.0, fr_up = 5;
  double fr_left = 1190.0, fr_down = 5E-3, fr_right = 3100.0, fr_up = 1E4;

  TCanvas *cMCMC = new TCanvas("c_lim_Asymptotic", "canvas with limits for Asymptotic CLs", 630, 600);
  cMCMC->cd();
  //cMCMC->SetGridx(1);
  //cMCMC->SetGridy(1);
  // draw a frame to define the range

  TH1F *hr = cMCMC->DrawFrame(fr_left, fr_down, fr_right, fr_up, "");
  TString VV = "ZH";
  
  hr->SetXTitle(xTitle.c_str());
  hr->SetYTitle("#sigma(pp#rightarrowX)#timesBR(X#rightarrowHH#rightarrowb#bar{b}b#bar{b})[fb]"); // #rightarrow 2l2q
  hr->SetMinimum(0.06);
  hr->SetMaximum(1000);

  gr95_cls->SetFillColor(kOrange);
  gr95_cls->SetFillStyle(1001);//solid
  gr95_cls->SetLineStyle(kDashed);
  gr95_cls->SetLineWidth(3);
  gr95_cls->GetXaxis()->SetTitle(xTitle.c_str());
  gr95_cls->GetYaxis()->SetTitle("95% CLs on #sigma(Z`#rightarrow#chi#bar{@chi}H)#timesBR(H#rightarrowb#bar{b})[fb] "); // #rightarrow 2l2q
  //gr95_cls->GetYaxis()->SetTitle("95% CLs on #sigma(Z`#rightarrow#chi#bar{#chi}H)[pb] "); // #rightarrow 2l2q
  gr95_cls->GetXaxis()->SetRangeUser(fr_left, fr_right);
  //gr95_cls->GetYaxis()->SetRangeUser(1E-2, 1E4);
  gr95_cls->Draw("3");
  //  gr95_cls->SetMinimum(0.00001);
  //gr95_cls->SetMaximum(1000.0);
  
  //grmedian_cls->SetMinimum(0.00001);
  //grmedian_cls->SetMaximum(1000.0);
  
  gr68_cls->SetFillColor(kGreen+1);
  gr68_cls->SetFillStyle(1001);//solid
  gr68_cls->SetLineStyle(kDashed);
  gr68_cls->SetLineWidth(3);
  gr68_cls->Draw("3same");
  grmedian_cls->GetXaxis()->SetTitle(xTitle.c_str());
  grmedian_cls->GetYaxis()->SetTitle("95% CLs on #sigma(Z`#rightarrow#chi#bar{#chi}H)[fb]"); // #rightarrow 2l2q
  grmedian_cls->SetMarkerStyle(24);//25=hollow squre
  grmedian_cls->SetMarkerColor(kBlack);
  grmedian_cls->SetLineStyle(2);
  grmedian_cls->SetLineWidth(3);


  grobslim_cls->SetMarkerColor(kBlack);
  grobslim_cls->SetMarkerStyle(21);//24=hollow circle
  grobslim_cls->SetMarkerSize(1.0);
  grobslim_cls->SetLineStyle(1);
  grobslim_cls->SetLineWidth(3);

  grthSM->SetLineColor(kRed);
  grthSM->SetLineWidth(2);
  grthSM->SetLineStyle(kSolid);
  grthSM->SetFillColor(kRed);
  grthSM->SetFillStyle(3344);

  grthSM10->SetLineColor(kRed);
  grthSM10->SetLineWidth(2);
  grthSM10->SetLineStyle(1);
  grthSM10->SetLineStyle(kDashed);
  grthSM10->SetFillColor(kRed);
  grthSM10->SetFillStyle(3344);


  grthSM->Draw("L3");
  grmedian_cls->Draw("L");
  // observed limit
  grobslim_cls->Draw("LP");

  /*
  TFile *fUnMPlus=new TFile("AsymptoticCLs_UnmatchedPlus_TGraph.root","READ");
  TGraph *grobs_ump=(TGraph*)fUnMPlus->Get("LimitObservedCLs");
  TGraph *grmedian_ump=(TGraph*)fUnMPlus->Get("LimitExpectedCLs");
  grobs_ump->SetName("LimitObs_UnmatchedPlus");
  grmedian_ump->SetName("LimitExp_UnmatchedPlus");
  grobs_ump->SetMarkerColor(kBlue);
  grobs_ump->SetLineColor(kBlue);
  grobs_ump->SetMarkerStyle(25);
  grmedian_ump->SetMarkerColor(kBlue);
  grmedian_ump->SetLineColor(kBlue);
  grmedian_ump->SetMarkerStyle(25);
  grobs_ump->Draw("P");
  grmedian_ump->Draw("L");
  */

  //draw grid on top of limits
  gStyle->SetOptStat(0);
  TH1D* postGrid = new TH1D("postGrid", "", 1, fr_left, fr_right);
  postGrid->GetYaxis()->SetRangeUser(fr_down, fr_up);
  postGrid->Draw("AXIGSAME");

  //more graphics

  TLegend *leg = new TLegend(.50, .65, .85, .90);
  //   TLegend *leg = new TLegend(.35,.71,.90,.90);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  //   leg->SetBorderMode(0);
  
  //leg->AddEntry(grobslim_cls, "CL_{S} Observed", "LP");
  //leg->SetHeader(sCat.c_str());
  leg->AddEntry(grmedian_cls, "CL_{S} Expected", "L");
  leg->AddEntry(gr68_cls, "CL_{S} Expected #pm 1#sigma", "LF");
  leg->AddEntry(gr95_cls, "CL_{S} Expected #pm 2#sigma", "LF");
  leg->AddEntry(grthSM, sThTitle.c_str(), "L");
//    leg->AddEntry(grthSM, "#sigma_{TH} x BR(Z' #rightarrow " + VV + "), #tilde{k}=0.50", "L"); // #rightarrow 2l2q
//    leg->AddEntry(grthSM10, "#sigma_{TH} x BR(Z' #rightarrow " + VV + "), #tilde{k}=0.20", "L"); // #rightarrow 2l2q
  leg->Draw();

    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->SetTextAlign(31);
    latex->SetTextAlign(11); // align left
    //latex->DrawLatex(0.18, 0.96, "CMS preliminary 2015");
   // latex->DrawLatex(0.60, 0.96, Form("%.2f fb^{-1} at #sqrt{s} = 13 TeV", intLumi));
  TLatex *lar = new TLatex();

  lar->SetTextSize(0.044);
 lar->SetTextAlign(12);
 lar->SetNDC(kTRUE);
 lar->SetTextFont(42);
  //lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.66, 0.975, " 35.9 fb^{-1} (13 TeV)");
  
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

  // cMCMC->RedrawAxis("");
  gPad->RedrawAxis("");
  // hr->GetYaxis()->DrawClone();
  cMCMC->Update();
  char fnam[50];
  //string outputname="shape2d";
  //string outputHHbbbbpe1d";
  //string outputname="counting";

  //  sprintf(fnam, "HHtobbbb_%s_Asymptotic.png", outputdir.data());
  //  sprintf(fnam, "HHtobbbb_%s_Asymptotic.pdf", outputdir.data());
  //  cMCMC->SaveAs(fnam);
  gPad->SetLogy(1);

  sOutput = string("BrazilianFlags/Limits_") + sCat + "_HH_log.png";
  cMCMC->SaveAs(sOutput.c_str());

  sOutput =  string("BrazilianFlags/Limits_") + sCat + "_HH_log.pdf";
  cMCMC->SaveAs(sOutput.c_str());

	cMCMC->Print("limit.pdf");
	cMCMC->Print("limit.png");

  /*
    sprintf(fnam, "XZHllbb_%s_Asymptotic.root",outputdir.data() );
    cMCMC->SaveAs(fnam);
    sprintf(fnam, "XZHllbb_%s_Asymptotic.eps", outputdir.data());
    cMCMC->SaveAs(fnam);
    
    sprintf(fnam, "XZHllbb_%s_Asymptotic.pdf", outputdir.data());
    cMCMC->SaveAs(fnam);
    gPad->SetLogy();
    sprintf(fnam, "XZHllbb_%s_Asymptotic_log.eps", outputdir.data());
    cMCMC->SaveAs(fnam);
    sprintf(fnam, "XZHllbb_%s_Asymptotic_log.pdf", outputdir.data());
    cMCMC->SaveAs(fnam);
  */

  cMCMC->Draw();

  fout->cd();
  grmedian_cls->Write();
  fout->Close();

}//end main

void setFPStyle()
{
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);//0.13);
  gStyle->SetPadLeftMargin(0.15);//0.16);
  gStyle->SetPadRightMargin(0.05);//0.02);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the Frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(605, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);


  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // gStyle->SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(1.15);//0.9);
  gStyle->SetTitleYOffset(1.3); // => 1.15 if exponents
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
}


