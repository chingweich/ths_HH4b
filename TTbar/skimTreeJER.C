//Histogram 
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>

//vector ,string ,stream
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

//root feature
#include <TLegend.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <TSystem.h>
#include <TF1.h>
//math 
#include <cmath>
#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "../untuplizer.h"
//#include "ApplyJER.h"
//#include "jetEnergyScale.h"

#include "standalone_LumiReWeighting.cc"
#include "standalone_LumiReWeighting.h"
//#include "BTagCalibrationStandalone.h"
//#include "BTagCalibrationStandalone.cpp"

float getPUPPIweight(float puppipt, float puppieta ){

   TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(1.00626, -1.06161, 0.07999,1.20454 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_cen->SetParameters(1.09302,-0.000150068,3.44866e-07,-2.68100e-10,8.67440e-14,-1.00114e-17);
  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters( 1.27212,-0.000571640,8.37289e-07,-5.20433e-10,1.45375e-13,-1.50389e-17);
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  totalWeight = genCorr * recoCorr;
  return totalWeight;
}

double ApplyJERp4(double eta, int jerShift) {

  eta = abs(eta) ;
  if (eta >= 4.7) eta = 4.699 ;
  double jerscale(0) ; 

  if (eta >= 3.2 && eta < 5.0) {
    if (jerShift == 1) jerscale = 1.16  ;  
    else if (jerShift == 2) jerscale = 1.16 + 0.029 ;  
    else if (jerShift == -1) jerscale = 1.16 - 0.029 ;  
  }
  else if (eta >= 3.0 && eta < 3.2) {
    if (jerShift == 1) jerscale = 1.328  ;  
    else if (jerShift == 2) jerscale = 1.328 + 0.022;  
    else if (jerShift == -1) jerscale = 1.328 - 0.022;  
  }
  else if (eta >= 2.8 && eta < 3.0) {
    if (jerShift == 1) jerscale =  1.857 ;  
    else if (jerShift == 2) jerscale =  1.857 + 0.071;  
    else if (jerShift == -1) jerscale =  1.857 - 0.071;  
  }
  else if (eta >= 2.5 && eta < 2.8) {
    if (jerShift == 1) jerscale = 1.364  ;  
    else if (jerShift == 2) jerscale = 1.364 + 0.039;  
    else if (jerShift == -1) jerscale = 1.364 - 0.039;  
  }
  else if (eta >= 2.3 && eta < 2.5) {
    if (jerShift == 1) jerscale = 1.177 ;  
    else if (jerShift == 2) jerscale = 1.177 + 0.041;  
    else if (jerShift == -1) jerscale = 1.177 - 0.041;  
  }
  else if (eta >= 2.1 && eta < 2.3) {
    if (jerShift == 1) jerscale = 1.067 ;  
    else if (jerShift == 2) jerscale = 1.067 + 0.053;  
    else if (jerShift == -1) jerscale = 1.067 - 0.053 ;  
  }
  else if (eta >= 1.9 && eta < 2.1) {
    if (jerShift == 1) jerscale = 1.140  ;  
    else if (jerShift == 2) jerscale = 1.140 + 0.047;  
    else if (jerShift == -1) jerscale = 1.140 - 0.047;  
  }
  else if (eta >= 1.7 && eta < 1.9) {
    if (jerShift == 1) jerscale = 1.082  ;  
    else if (jerShift == 2) jerscale = 1.082 + 0.035;  
    else if (jerShift == -1) jerscale = 1.082 - 0.035 ;  
  }
  else if (eta >= 1.3 && eta < 1.7) {
    if (jerShift == 1) jerscale = 1.084  ;  
    else if (jerShift == 2) jerscale = 1.084 + 0.011;  
    else if (jerShift == -1) jerscale = 1.084 - 0.011 ;  
  }
  else if (eta >= 1.1 && eta < 1.3) {
    if (jerShift == 1) jerscale = 1.123  ;  
    else if (jerShift == 2) jerscale = 1.123 + 0.024;  
    else if (jerShift == -1) jerscale = 1.123 - 0.024;  
  }
  else if (eta >= 0.8 && eta < 1.1) {
    if (jerShift == 1) jerscale = 1.114  ;  
    else if (jerShift == 2) jerscale = 1.114 + 0.013;  
    else if (jerShift == -1) jerscale = 1.114 - 0.013;  
  }
  else if (eta >= 0.5 && eta < 0.8) {
    if (jerShift == 1) jerscale = 1.138  ;  
    else if (jerShift == 2) jerscale = 1.138 + 0.013;  
    else if (jerShift == -1) jerscale = 1.138 - 0.013;  
  }
  else if (eta >= 0.0 && eta < 0.5) {
    if (jerShift == 1) jerscale = 1.109  ;  
    else if (jerShift == 2) jerscale = 1.109 + 0.008;  
    else if (jerShift == -1) jerscale = 1.109 - 0.008;  
  }

  return jerscale ; 
}


void skimTreeJERBase(string w , string st){
	
	
	//option-----------------------------------------------------------
	
	
	standalone_LumiReWeighting LumiWeights_central(0),LumiWeights_up(1),LumiWeights_down(-1);
		
		TFile* tf1=TFile::Open(st.data());
		TTree* tree=(TTree*)tf1->Get("treeMaker");
		TreeReader 	data(tree);
		//data.Print();
		//cout<<"here"<<endl;
		TFile* outFile = new TFile(Form("minitree/%s.root",w.data()),"recreate");
		//TTree* hh4b=new TTree("hh4b","hh4b");
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
		Float_t dijetmass_softdrop_corr_JERup; mynewTree->Branch("dijetmass_softdrop_corr_JERup",&dijetmass_softdrop_corr_JERup,"dijetmass_softdrop_corr_JERup/F");
		Float_t dijetmass_softdrop_corr_JERdown; mynewTree->Branch("dijetmass_softdrop_corr_JERdown",&dijetmass_softdrop_corr_JERdown,"dijetmass_softdrop_corr_JERdown/F");
		Float_t dijetmass_softdrop_corr_JECup; mynewTree->Branch("dijetmass_softdrop_corr_JECup",&dijetmass_softdrop_corr_JECup,"dijetmass_softdrop_corr_JECup/F");
		Float_t dijetmass_softdrop_corr_JECdown; mynewTree->Branch("dijetmass_softdrop_corr_JECdown",&dijetmass_softdrop_corr_JECdown,"dijetmass_softdrop_corr_JECdown/F");
		
		Float_t dijetmass_corr_punc; mynewTree->Branch("dijetmass_corr_punc",&dijetmass_corr_punc,"dijetmass_corr_punc/F");
		
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
		Float_t gen1Pt; mynewTree->Branch("gen1Pt",&gen1Pt,"gen1Pt/F");
		Float_t gen1phi; mynewTree->Branch("gen1phi",&gen1phi,"gen1phi/F");
		Float_t gen1Eta; mynewTree->Branch("gen1Eta",&gen1Eta,"gen1Eta/F");
		Float_t gen1Mass; mynewTree->Branch("gen1Mass",&gen1Mass,"gen1Mass/F");
		Float_t gen1Id; mynewTree->Branch("gen1Id",&gen1Id,"gen1Id/F");
		Float_t gen2Pt; mynewTree->Branch("gen2Pt",&gen2Pt,"gen2Pt/F");
		Float_t gen2phi; mynewTree->Branch("gen2phi",&gen2phi,"gen2phi/F");
		Float_t gen2Eta; mynewTree->Branch("gen2Eta",&gen2Eta,"gen2Eta/F");
		Float_t gen2Mass; mynewTree->Branch("gen2Mass",&gen2Mass,"gen2Mass/F");
		
		
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
		
		//Float_t isData; mynewTree->Branch("isData",&isData,"isData/F");
		
		Float_t jet1JERup; mynewTree->Branch("jet1JERup",&jet1JERup,"jet1JERup/F");
		Float_t jet1JERcentral; mynewTree->Branch("jet1JERcentral",&jet1JERcentral,"jet1JERcentral/F");
		Float_t jet1JERdown; mynewTree->Branch("jet1JERdown",&jet1JERdown,"jet1JERdown/F");
		Float_t jet1JECup; mynewTree->Branch("jet1JECup",&jet1JECup,"jet1JECup/F");
		Float_t jet1JECdown; mynewTree->Branch("jet1JECdown",&jet1JECdown,"jet1JECdown/F");
		
		Float_t jet2JERup; mynewTree->Branch("jet2JERup",&jet2JERup,"jet2JERup/F");
		Float_t jet2JERcentral; mynewTree->Branch("jet2JERcentral",&jet2JERcentral,"jet2JERcentral/F");
		Float_t jet2JERdown; mynewTree->Branch("jet2JERdown",&jet2JERdown,"jet2JERdown/F");
		Float_t jet2JECup; mynewTree->Branch("jet2JECup",&jet2JECup,"jet2JECup/F");
		Float_t jet2JECdown; mynewTree->Branch("jet2JECdown",&jet2JECdown,"jet2JECdown/F");
		
		Float_t SFTight=1; mynewTree->Branch("SFTight",&SFTight,"SFTight/F");
		Float_t SFTightup=1; mynewTree->Branch("SFTightup",&SFTightup,"SFTightup/F");
		Float_t SFTightdown=1; mynewTree->Branch("SFTightdown",&SFTightdown,"SFTightdown/F");
		Float_t SFLoose=1; mynewTree->Branch("SFLoose",&SFLoose,"SFLoose/F");
		Float_t SFLooseup=1; mynewTree->Branch("SFLooseup",&SFLooseup,"SFLooseup/F");
		Float_t SFLoosedown=1; mynewTree->Branch("SFLoosedown",&SFLoosedown,"SFLoosedown/F");
		
		Float_t trigWeightUp=1; mynewTree->Branch("trigWeightUp",&trigWeightUp,"trigWeightUp/F");
		Float_t trigWeightDown=1; mynewTree->Branch("trigWeightDown",&trigWeightDown,"trigWeightDown/F");
		
		
		bool HLT_PFHT800=1; mynewTree->Branch("HLT_PFHT800",&HLT_PFHT800,"HLT_PFHT800/O");
		bool HLT_PFHT900=1; mynewTree->Branch("HLT_PFHT900",&HLT_PFHT900,"HLT_PFHT900/O");
		bool HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=1; mynewTree->Branch("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5",&HLT_PFHT650_WideJetMJJ900DEtaJJ1p5,"HLT_PFHT650_WideJetMJJ900DEtaJJ1p5/O");
		bool HLT_AK8PFJet360_TrimMass30=1; mynewTree->Branch("HLT_AK8PFJet360_TrimMass30",&HLT_AK8PFJet360_TrimMass30,"HLT_AK8PFJet360_TrimMass30/O");
		bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=1; mynewTree->Branch("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",&HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20,"HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20/O");
		bool HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=1; mynewTree->Branch("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50",&HLT_AK8PFHT650_TrimR0p1PT0p03Mass50,"HLT_AK8PFHT650_TrimR0p1PT0p03Mass50/O");
		bool HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=1; mynewTree->Branch("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50",&HLT_AK8PFHT700_TrimR0p1PT0p03Mass50,"HLT_AK8PFHT700_TrimR0p1PT0p03Mass50/O");
		
		int nPass[20]={0},total=0;
		total+=data.GetEntriesFast();
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
			data.GetEntry(jEntry);
			
			bool isData  = data.GetBool("isData"); // only apply trigger if it's data

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    HLT_PFHT800=0;
    HLT_PFHT900=0;
    HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=0;
    HLT_AK8PFJet360_TrimMass30=0;
    HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=0;
    HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=0;
    HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=0;
    
    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
        bool results = trigResult[it];
		if(thisTrig.find("HLT_PFHT800_v")!= std::string::npos )HLT_PFHT800=results;
		if(thisTrig.find("HLT_PFHT900_v")!= std::string::npos )HLT_PFHT900=results;
		if(thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos )HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=results;
		if(thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos )HLT_AK8PFJet360_TrimMass30=results;
		if(thisTrig.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v")!= std::string::npos )HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=results;
		if(thisTrig.find("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")!= std::string::npos )HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=results;
		if(thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos )HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=results;
      }

    if(!(HLT_PFHT800||HLT_PFHT900||HLT_PFHT650_WideJetMJJ900DEtaJJ1p5||HLT_AK8PFJet360_TrimMass30||HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20||
	HLT_AK8PFHT650_TrimR0p1PT0p03Mass50||HLT_AK8PFHT700_TrimR0p1PT0p03Mass50))continue;
    nPass[0]++;

    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[1]++;

    //1. veto events containing 1 or 2 leptons
	
    Int_t nEle         = data.GetInt("nEle");
    TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");
    Float_t* eleSCEta                = data.GetPtrFloat("eleScEta");
    Float_t* eleSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
    Float_t* eleHoverE               = data.GetPtrFloat("eleHoverE");
    Float_t* eleEcalPFClusterIso     = data.GetPtrFloat("eleEcalPFClusterIso");
    Float_t* eleHcalPFClusterIso     = data.GetPtrFloat("eleHcalPFClusterIso");
    Float_t* eleDr03TkSumPt          = data.GetPtrFloat("eleDr03TkSumPt");
    Float_t* eledEtaAtVtx            = data.GetPtrFloat("eledEtaAtVtx");
    Float_t* eledPhiAtVtx            = data.GetPtrFloat("eledPhiAtVtx");
    Int_t*   eleCharge               = data.GetPtrInt("eleCharge");
    vector<bool>    &isMediumEle= *((vector<bool>*) data.GetPtr("eleIsPassMVAMedium"));
    vector<bool>     &isTightEle= *((vector<bool>*) data.GetPtr("eleIsPassMVATight"));

    // check if there is any single tight electron
    
    bool hasSingleTightElectron=false;


    for(int ie=0; ie < nEle; ie++){

      TLorentzVector* thisEle= (TLorentzVector*)eleP4->At(ie);
      if( fabs(thisEle->Eta())>2.4 )continue;
      if( thisEle->Pt() < 25 )continue;
      bool preSelect= 

	( fabs(eleSCEta[ie]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[ie] < 0.012 && 
	  eleHoverE[ie] < 0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.37 && ( eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.25 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 
	  && fabs(eledEtaAtVtx[ie]) < 0.0095 && fabs(eledPhiAtVtx[ie]) < 0.065 ) ||

	( fabs(eleSCEta[ie]) > 1.5660 && fabs(eleSCEta[ie]) < 2.5 && eleSigmaIEtaIEtaFull5x5[ie] < 0.033 && 
	  eleHoverE[ie] <0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.45 && (eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.28 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 );
    
      if(!preSelect)continue;

      if( !isTightEle[ie] )continue;

	
      hasSingleTightElectron=true;
      break;
      
    }
    if( hasSingleTightElectron )continue;
    nPass[2]++;


    // check if there is any double loose electron
    bool hasLooseElectronPair=false;
    for(int ie=0; ie < nEle; ie++){
      for(int je=0; je < ie; je++){

	if(eleCharge[ie]*eleCharge[je]>0)continue;

   	TLorentzVector* thisEle= (TLorentzVector*)eleP4->At(ie);

    	if( fabs(thisEle->Eta())>2.4 )continue;
    	if( thisEle->Pt() < 15 )continue;

   	TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(je);
    	if( fabs(thatEle->Eta())>2.4 )continue;
    	if( thatEle->Pt() < 15 )continue;

	bool preSelect_ie= 

	( fabs(eleSCEta[ie]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[ie] < 0.012 && 
	  eleHoverE[ie] < 0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.37 && ( eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.25 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 
	  && fabs(eledEtaAtVtx[ie]) < 0.0095 && fabs(eledPhiAtVtx[ie]) < 0.065 ) ||

	( fabs(eleSCEta[ie]) > 1.5660 && fabs(eleSCEta[ie]) < 2.5 && eleSigmaIEtaIEtaFull5x5[ie] < 0.033 && 
	  eleHoverE[ie] <0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.45 && (eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.28 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 );
	
	if(!preSelect_ie)continue;


	bool preSelect_je= 

	( fabs(eleSCEta[je]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[je] < 0.012 && 
	  eleHoverE[je] < 0.09 && (eleEcalPFClusterIso[je]/thatEle->Pt()) < 0.37 && ( eleHcalPFClusterIso[je]/thatEle->Pt()) < 0.25 && (eleDr03TkSumPt[je]/thatEle->Pt()) < 0.18 
	  && fabs(eledEtaAtVtx[je]) < 0.0095 && fabs(eledPhiAtVtx[je]) < 0.065 ) ||

	( fabs(eleSCEta[je]) > 1.5660 && fabs(eleSCEta[je]) < 2.5 && eleSigmaIEtaIEtaFull5x5[je] < 0.033 && 
	  eleHoverE[je] <0.09 && (eleEcalPFClusterIso[je]/thatEle->Pt()) < 0.45 && (eleHcalPFClusterIso[je]/thatEle->Pt()) < 0.28 && (eleDr03TkSumPt[je]/thatEle->Pt()) < 0.18 );
	
	if(!preSelect_je)continue;


	if( !isMediumEle[ie] )continue;
	if( !isMediumEle[je] )continue;
	
	hasLooseElectronPair=true;
	break;
      
      }
    }

    if( hasLooseElectronPair )continue;
    nPass[3]++;
    
    
    // veto of single tight muon or two oppositely-charged muon pair
    Int_t nMu          = data.GetInt("nMu");
    TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t*   muCharge        = data.GetPtrInt("muCharge");
    vector<bool>    &isLooseMuon= *((vector<bool>*) data.GetPtr("isLooseMuon"));
    vector<bool>    &isTightMuon= *((vector<bool>*) data.GetPtr("isTightMuon"));
    Float_t*          muIso1  = data.GetPtrFloat("muChHadIso");
    Float_t*          muIso2  = data.GetPtrFloat("muNeHadIso");
    Float_t*          muIso3  = data.GetPtrFloat("muGamIso");
    Float_t*          muIso4  = data.GetPtrFloat("muPUPt");
    
    // check if there is any single tight muon
    bool hasSingleTightMuon=false;

    for(int im=0; im < nMu; im++){

      TLorentzVector* thisMu = (TLorentzVector*)muP4->At(im);

      if( fabs(thisMu->Eta())>2.4 )continue;
      if( thisMu->Pt() < 25 )continue;
      if( !isTightMuon[im] )continue;

      float iso =  (muIso1[im] + max(0., muIso2[im]+muIso3[im] - 0.5*muIso4[im]))/thisMu->Pt();
      if( iso > 0.15)continue;
	
      hasSingleTightMuon=true;
      break;
      
    }
    if( hasSingleTightMuon )continue;
    nPass[4]++;


    // check if there is any double loose muons
    bool hasLooseMuonPair=false;
    for(int im=0; im < nMu; im++){
      for(int jm=0; jm < im; jm++){

	if(muCharge[im]*muCharge[jm]>0)continue;
	if( !isLooseMuon[im] )continue;
	if( !isLooseMuon[jm] )continue;

   	TLorentzVector* thisMu = (TLorentzVector*)muP4->At(im);

    	if( fabs(thisMu->Eta())>2.4 )continue;
    	if( thisMu->Pt() < 15 )continue;

   	TLorentzVector* thatMu = (TLorentzVector*)muP4->At(jm);
    	if( fabs(thatMu->Eta())>2.4 )continue;
    	if( thatMu->Pt() < 15 )continue;

	float iso1=  (muIso1[im] + max(0., muIso2[im]+muIso3[im] - 0.5*muIso4[im]))/thisMu->Pt();
	if( iso1> 0.25)continue;

	float iso2=  (muIso1[jm] + max(0., muIso2[jm]+muIso3[jm] - 0.5*muIso4[jm]))/thatMu->Pt();
	if( iso2> 0.25)continue;
	
	hasLooseMuonPair=true;
	break;
      
      }
    }
    if( hasLooseMuonPair )continue;

    nPass[5]++;

    // apply jet requirement

    int nJets         = data.GetInt("FATnJet");
    TClonesArray* fatjetP4   = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    TClonesArray* puppijetP4 = (TClonesArray*) data.GetPtrTObject("FATjetPuppiP4");

    Float_t*  FATjetPuppiTau1 = data.GetPtrFloat("FATjetPuppiTau1");
    Float_t*  FATjetPuppiTau2 = data.GetPtrFloat("FATjetPuppiTau2");
    Float_t*  FATjetPuppiTau3 = data.GetPtrFloat("FATjetPuppiTau3");
    Float_t*  fatjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
   
    vector<bool>    &passFatJetTightID = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));

    Float_t*  fatjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
    Float_t*  fatjetMuoEF = data.GetPtrFloat("FATjetMuoEF");

    int nADDJets         = data.GetInt("ADDnJet");
    TClonesArray* addjetP4   = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doublesv = data.GetPtrFloat("ADDjet_DoubleSV");

    if(nJets<1)continue;
    nPass[6]++;

    if(nADDJets<1)continue;
    nPass[7]++;

	//float leading_thea_mass=0;
	float thea_mass[2]={0};
	bool vetoBoost=0;
    if(nJets>1){
		int nGoodJets=0;
		
		int addJetIndex[2]={-1,-1};
		for(int ij=0; ij<2; ij++){
			if( !passFatJetTightID[ij] )continue;
			TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
			if( thisJet->Pt()<300 )continue;
			if( fabs(thisJet->Eta())>2.4 )continue;
			float tau21_puppi = FATjetPuppiTau2[ij]/FATjetPuppiTau1[ij];
			if(tau21_puppi>0.55)continue;
			TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(ij);
			float thea_corr = getPUPPIweight(thisJetPuppi->Pt(),thisJetPuppi->Eta());
			float raw_mass = fatjetPuppiSDmass[ij];
			thea_mass[ij] = raw_mass*thea_corr;
			//if(i==0)leading_thea_mass=thea_mass[ij];
			if(thea_mass[ij] < 105 || thea_mass[ij] > 135)continue;
			for(int k=0; k < nADDJets; k++){
			TLorentzVector* thatJet = (TLorentzVector*)addjetP4->At(k);
				if(thisJet->DeltaR(*thatJet)<0.8 && addJetIndex[ij]<0)
				{
				addJetIndex[ij]=k;
				break;
				}
			}
			nGoodJets++;
		} 
		TLorentzVector* higgsJet[2];
		for(int i=0;i<2;i++)higgsJet[i] = (TLorentzVector*)fatjetP4->At(i);
		float dEta = fabs(higgsJet[0]->Eta()-higgsJet[1]->Eta());
		float mjj = (*higgsJet[0]+*higgsJet[1]).M();
		float mjjred = mjj + 250 - thea_mass[0]-thea_mass[1];
		
		vetoBoost= nGoodJets && (dEta<1.3) && (mjjred>750) && (addjet_doublesv[addJetIndex[0]]>0.3) && (addjet_doublesv[addJetIndex[1]]>0.3);
	}
	if(vetoBoost)continue;
	nPass[8]++;
	
	int isAK8jet=-1;
	for(int i=0;i<nJets;i++){
		if( !passFatJetTightID[i] )continue;
		TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(i);
		if( thisJet->Pt()<250 )continue;
		if( fabs(thisJet->Eta())>2.4 )continue;
		bool isdoubleb=0;
		for(int k=0; k < nADDJets; k++){
		TLorentzVector* thatJet = (TLorentzVector*)addjetP4->At(k);
			if(thisJet->DeltaR(*thatJet)<0.8 ){
				isdoubleb=(addjet_doublesv[k]>0.8);
				jet1bbtag=addjet_doublesv[k];
				break;
			}
		}
		if(!isdoubleb)continue;
		if(thea_mass[i] < 105 || thea_mass[i] > 135)continue;
		float tau21_puppi = FATjetPuppiTau2[i]/FATjetPuppiTau1[i];
		if(tau21_puppi>0.55)continue;
		isAK8jet=i;
		jet1_puppi_tau21=tau21_puppi;
		jet1_puppi_msoftdrop_TheaCorr=thea_mass[i];
		
		break;
	}
	
	if(isAK8jet<0)continue;
	//cout<<isAK8jet<<endl;
	nPass[9]++;
	TLorentzVector* AK8Jet = (TLorentzVector*)fatjetP4->At(isAK8jet);
	
    int THINnJet         = data.GetInt("THINnJet");
	int nGoodThinJets=0;
	
	vector<int> nGoodThinJetsIndex;
	TClonesArray* THINjetP4   = (TClonesArray*) data.GetPtrTObject("THINjetP4");
	TLorentzVector AK4s;
	for(int i=0;i<THINnJet;i++){
		TLorentzVector* thisJet = (TLorentzVector*)THINjetP4->At(i);
		if(thisJet->DeltaR(*AK8Jet)<0.8)continue;
		if( thisJet->Pt()<30 )continue;
		if( fabs(thisJet->Eta())>2.4 )continue;
		nGoodThinJets++;
		/*
		for(int j=i+1;j<THINnJet;j++){
			TLorentzVector* thatJet = (TLorentzVector*)THINjetP4->At(j);
			if(thisJet->DeltaR(*thatJet)>1.5)continue;
			pair=1;
			break;
		}*/
		nGoodThinJetsIndex.push_back(i);
	}
	if(nGoodThinJets<1)continue;
	nPass[10]++;
	
	bool pair=0;
	int AK4Index[2]={-1,-1};
	float AK4AddMass=0;
	for(unsigned int i=0;i<nGoodThinJetsIndex.size();i++){
		TLorentzVector* thisJet = (TLorentzVector*)THINjetP4->At(i);
		for(unsigned int j=i+1;j<nGoodThinJetsIndex.size();j++){
			TLorentzVector* thatJet = (TLorentzVector*)THINjetP4->At(j);
			if(thisJet->DeltaR(*thatJet)>1.5)continue;
			AK4AddMass=((*thisJet)+(*thatJet)).M();
			if(AK4AddMass<105 || AK4AddMass>135)continue;
			
			float dR=2;
			int nearestIndex=0;
			for(unsigned int k=0;k<nGoodThinJetsIndex.size();k++){
				if(k==i ||k==j)continue;
				TLorentzVector* nearestJet = (TLorentzVector*)THINjetP4->At(k);
				if(nearestJet->DeltaR(((*thisJet)+(*thatJet)))<dR){
					dR=nearestJet->DeltaR(((*thisJet)+(*thatJet)));
					nearestIndex=k;
				}
			}
			TLorentzVector* nearestJet = (TLorentzVector*)THINjetP4->At(nearestIndex);
			pair=1;
			AK4Index[0]=i;
			AK4Index[1]=j;
			if(((*thisJet)+(*thatJet)+(*nearestJet)).M()<200)continue;
			AK4s=((*thisJet)+(*thatJet));
			break;
		}
	}
	if(!pair)continue;
	//cout<<AK4Index[0]<<","<<AK4Index[1]<<endl;
	nPass[11]++;
	
	
			jet1pt=AK8Jet->Pt();
			jet1eta=AK8Jet->Eta();
			jet1phi=AK8Jet->Phi();
			jet1mass=AK8Jet->M();
			
			Int_t nGenPar        = data.GetInt("nGenPar");
			Int_t* genParId      = data.GetPtrInt("genParId");
			TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
			
			for(int i=0;i<nGenPar;i++){
				TLorentzVector* genJet = (TLorentzVector*)genParP4->At(i);
				if(genJet->DeltaR(*AK8Jet)<0.4){
					gen1Id=genParId[i];
					gen1Pt=genJet->Pt();
					gen1Eta=genJet->Eta();
					gen1phi=genJet->Phi();
					gen1Mass=genJet->M();
					
				}
			}
			
			etadiff=fabs(AK8Jet->Eta()-AK4s.Eta());
			dijetmass=(AK4s+*AK8Jet).M();
			
	/*		
			Long64_t runId        = data.GetLong64("runId");
			Long64_t lumiSection        = data.GetLong64("lumiSection");
			Long64_t eventId        = data.GetLong64("eventId");
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nJets);
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nJets);
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nJets);
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nJets);
			//vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			Float_t*  FATjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
			Int_t*  FATjetHadronFlavor = data.GetPtrInt("FATjetHadronFlavor");
			Int_t*  FATnSubSDJet = data.GetPtrInt("FATnSubSDJet");
			//Float_t*  FATjetPuppiSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr");
			Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  FATjetSDmass = data.GetPtrFloat("FATjetSDmass");
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  FATjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
			Float_t*  FATjetMuoEF = data.GetPtrFloat("FATjetMuoEF");
			
			Float_t  ptHat = data.GetFloat("ptHat");

			Float_t ntrue= data.GetFloat("pu_nTrueInt");
			double PU_weight[3]={1,1,1};
			
			
			if(ntrue<51){
					PU_weight[0] = LumiWeights_central.weight(ntrue);
					PU_weight[1]= LumiWeights_up.weight(ntrue);
					PU_weight[2] = LumiWeights_down.weight(ntrue);
				}
			else {
					PU_weight[0] = LumiWeights_central.weight(50);
					PU_weight[1] = LumiWeights_up.weight(50);
					PU_weight[2]= LumiWeights_down.weight(50);
			}
			puWeights=PU_weight[0] ;
			puWeightsUp=PU_weight[1];
			puWeightsDown=PU_weight[2];
			
			nTrueInt=ntrue;
			//isData=isData;
			
			//jet2pt=((*higgsJet[1])*ptsmearGlobal[0][1]).Pt();
			//jet2eta=((*higgsJet[1])*ptsmearGlobal[0][1]).Eta();
			//jet2phi=((*higgsJet[1])*ptsmearGlobal[0][1]).Phi();
			//jet2mass=((*higgsJet[1])*ptsmearGlobal[0][1]).M();
			
			jet1JERup=ptsmearGlobal[1][0];
			jet1JERcentral=ptsmearGlobal[0][0];
			jet1JERdown=ptsmearGlobal[2][0];
			jet1JECup=FATjetCorrUncUp[0];
			jet1JECdown=FATjetCorrUncDown[0];
			
			//jet2JERup=ptsmearGlobal[1][1];
			//jet2JERcentral=ptsmearGlobal[0][1];
			//jet2JERdown=ptsmearGlobal[2][1];
			//jet2JECup=FATjetCorrUncUp[1];
			//jet2JECdown=FATjetCorrUncDown[1];
			
			TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(0);
			TLorentzVector* thatJetPuppi = (TLorentzVector*)puppijetP4->At(1);
			jet1_puppi_pt=thisJetPuppi->Pt();
			jet1_puppi_eta=thisJetPuppi->Eta();
			jet1_puppi_phi=thisJetPuppi->Phi();
			jet1_puppi_mass=thisJetPuppi->M();
			//jet2_puppi_pt=thatJetPuppi->Pt();
			//jet2_puppi_eta=thatJetPuppi->Eta();
			//jet2_puppi_phi=thatJetPuppi->Phi();
			//jet2_puppi_mass=thatJetPuppi->M();
			
			
			
			dijetmass_softdrop_corr=mjjred;
			dijetmass_softdrop_corr_JERup=mjjredJERUp;
			dijetmass_softdrop_corr_JERdown=mjjredJERDown;
			
			dijetmass_softdrop_corr_JECup=mjjredJECUp;
			dijetmass_softdrop_corr_JECdown=mjjredJECDown;
			
			dijetmass_pruned_corr=mjj;
			dijetmass_corr_punc=mjj;
			
			jet1pmass=fatjetPRmassL2L3Corr[0];
			jet2pmass=fatjetPRmassL2L3Corr[1];
			addjet_doublesv[addJetIndex[0]];
			jet2bbtag=addjet_doublesv[addJetIndex[1]];
		
			
			jet2_puppi_tau21=FATjetPuppiTau2[1]/FATjetPuppiTau1[1];
			
			jet1_puppi_msoftdrop=fatjetPuppiSDmass[0];
			jet2_puppi_msoftdrop=fatjetPuppiSDmass[1];
		
			
			jet2_puppi_msoftdrop_TheaCorr=thea_mass[1];
			
			jet1_puppi_TheaCorr=thea_mass[0]/fatjetPuppiSDmass[0];
			jet2_puppi_TheaCorr=thea_mass[1]/fatjetPuppiSDmass[1];
			
			
			*/
			mynewTree ->Fill();
		}//end event loop----------------------------------------------------------------------------------------
		mynewTree->Write();
		
		TH1F* th1f= new  TH1F("CountWeighted","",1,0,1);
		th1f->SetBinContent(1,total);
		th1f->Write();
		
		outFile->Close();
		cout<<"total="<<total<<endl;
		for(int i=0;i<15;i++)cout<<"npass["<<i<<"]="<<nPass[i]<<endl;
	
	
	
	
}

void skimTreeJER(){
	skimTreeJERBase("TTbar_100events_v2","skimmed/NCUGlobalTuples_1.root");
	
}
