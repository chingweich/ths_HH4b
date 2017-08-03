#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TColor.h>
#include <TLegend.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFrame.h>
#include <TLatex.h>
//#include "../macro/mkPlotsLivia/CMS_lumi.C"
#include <iostream>
#include <vector>
//Histogram 
#include <TH1.h>
//#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
//#include <TGraph.h>
//#include <TGraphErrors.h>

//vector ,string ,stream
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
//#include <sstream>

//root feature
//#include <TLegend.h>
//#include <TRandom.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TFile.h>
//#include <TCanvas.h>
//#include "TSystem.h"
//#include "TStyle.h"
#include <TClonesArray.h>
#include <TF1.h>
//math 
#include <cmath>
#include "../../untuplizer.h"


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

void skimTreeJERBase(string outputStr,string st){
	double nPass[30]={0},total=0,fixScaleNum[8]={0};

	TreeReader data(st.data());
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
		if(jEntry%10000==0)cout<<jEntry<<" out of "<<data.GetEntriesFast() <<" events are processed"<<endl;
		data.GetEntry(jEntry);
		nPass[0]++;
		
		//0. has a good vertex
		Int_t nVtx        = data.GetInt("nVtx");
		if(nVtx<1)continue;
		bool isData  = data.GetBool("isData"); // only apply trigger if it's data
		
		std::string* trigName = data.GetPtrString("hlt_trigName");
		vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
		const Int_t nsize = data.GetPtrStringSize();
		
		int nJets         = data.GetInt("FATnJet");
		if(nJets<2)continue;
		bool passTrigger=false;
		for(int it=0; it< nsize; it++){
			std::string thisTrig= trigName[it];
			bool results = trigResult[it];
			if( (thisTrig.find("HLT_PFHT800_v")!= std::string::npos && results==1) ||
				(thisTrig.find("HLT_PFHT900_v")!= std::string::npos && results==1) || // for runH
				(thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos && results==1) ||
				(thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos && results==1) ||
				(thisTrig.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v")!= std::string::npos && results==1) ||
				(thisTrig.find("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")!= std::string::npos && results==1)){
					passTrigger=true;
					break;
			}
		}

		if(!passTrigger )continue;
		nPass[1]++;

		TClonesArray* fatjetP4   = (TClonesArray*) data.GetPtrTObject("FATjetP4");
		TLorentzVector* higgsJet[2];
		for(int i=0;i<2;i++)higgsJet[i] = (TLorentzVector*)fatjetP4->At(i);
		if( higgsJet[0]->Pt()<300 || higgsJet[1]->Pt()<300)continue;
		if( fabs(higgsJet[0]->Eta())>2.4 || fabs(higgsJet[1]->Eta())>2.4)continue;
		nPass[2]++;
		
		vector<bool>    &passFatJetTightID = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
		if( !passFatJetTightID[0] ||  !passFatJetTightID[1] )continue;
		nPass[3]++;
		
		float dEta = fabs(higgsJet[0]->Eta()-higgsJet[1]->Eta());
		if(dEta>1.3)continue;
		nPass[4]++;
		
		
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
			bool preSelect= ( fabs(eleSCEta[ie]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[ie] < 0.012 && 
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

				bool preSelect_ie= ( fabs(eleSCEta[ie]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[ie] < 0.012 && 
									eleHoverE[ie] < 0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.37 && ( eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.25 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 
									&& fabs(eledEtaAtVtx[ie]) < 0.0095 && fabs(eledPhiAtVtx[ie]) < 0.065 ) ||
									( fabs(eleSCEta[ie]) > 1.5660 && fabs(eleSCEta[ie]) < 2.5 && eleSigmaIEtaIEtaFull5x5[ie] < 0.033 && 
									eleHoverE[ie] <0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.45 && (eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.28 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 );
	
				if(!preSelect_ie)continue;
				bool preSelect_je= ( fabs(eleSCEta[je]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[je] < 0.012 && 
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
		Float_t*  FATjetPuppiTau1 = data.GetPtrFloat("FATjetPuppiTau1");
		Float_t*  FATjetPuppiTau2 = data.GetPtrFloat("FATjetPuppiTau2");
		if(FATjetPuppiTau2[0]/FATjetPuppiTau1[0]>0.55 || FATjetPuppiTau2[1]/FATjetPuppiTau1[1]>0.55)continue;
		nPass[6]++;
		
		TClonesArray* puppijetP4 = (TClonesArray*) data.GetPtrTObject("FATjetPuppiP4");

		
		Float_t*  fatjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
   
		
		Float_t*  fatjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
		Float_t*  fatjetMuoEF = data.GetPtrFloat("FATjetMuoEF");

		int nADDJets         = data.GetInt("ADDnJet");
		TClonesArray* addjetP4   = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
		Float_t*  addjet_doublesv = data.GetPtrFloat("ADDjet_DoubleSV");
		
		

		if(nADDJets<2)continue;
		

		
		
		
		
		Int_t nGoodJets=0;
		const float dRMax=0.8;
		int addJetIndex[2]={-1,-1};
		float thea_mass[2]={-1,-1};
		for(int ij=0; ij<2; ij++){
			TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
			//float tau21_puppi = FATjetPuppiTau2[ij]/FATjetPuppiTau1[ij];
			//if( tau21_puppi > 0.6 )continue;
			TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(ij);
			float thea_corr = getPUPPIweight(thisJetPuppi->Pt(),thisJetPuppi->Eta());
			float raw_mass = fatjetPuppiSDmass[ij];
			thea_mass[ij] = raw_mass*thea_corr;

			if(thea_mass[ij] < 105 || thea_mass[ij] > 135)continue;
			if(thea_mass[ij] < 50)continue;
			// now look for add jets
			for(int k=0; k < nADDJets; k++){
				TLorentzVector* thatJet = (TLorentzVector*)addjetP4->At(k);
				if(thisJet->DeltaR(*thatJet)<dRMax && addJetIndex[ij]<0){
					addJetIndex[ij]=k;
					break;
				}
			}
			nGoodJets++;
		} // end of loop over two leading jets
		if(nGoodJets<2)continue;
		nPass[7]++;
		float mjj = (*higgsJet[0]+*higgsJet[1]).M();

		float mjjred = mjj + 250 - thea_mass[0]-thea_mass[1];
		if(mjjred<750)continue;

		if(nGoodJets<2)continue;
		nPass[8]++;

		if(addJetIndex[0]==addJetIndex[1])continue;
		if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;
		

		// check if any jet fail loose double b-tag

		const float DBLoose = 0.3;
		const float DBTight = 0.8;

		if(addjet_doublesv[addJetIndex[0]]< DBLoose)continue;
		if(addjet_doublesv[addJetIndex[1]]< DBLoose)continue;
		nPass[9]++;
		
		if(addjet_doublesv[addJetIndex[0]]< DBTight)continue;
		if(addjet_doublesv[addJetIndex[1]]< DBTight)continue;
		nPass[10]++;
		
	}
	TFile* output=new TFile(Form("%s.root",outputStr.data()),"recreate");
	TH1D* thf=new TH1D("cut","cut",11,-0.5,11.5);
	for(int i=0;i<11;i++)thf->SetBinContent(i+1,nPass[i]);
	thf->Write();
	output->Close();
}

void radionCutflow(){
	skimTreeJERBase("BulkGrav_M-750_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/BulkGravTohhTohbbhbb/GluGluToBulkGravitonToHHTo4B_M-750_narrow_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-800_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/BulkGravTohhTohbbhbb/GluGluToBulkGravitonToHHTo4B_M-800_narrow_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-900_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/BulkGravTohhTohbbhbb/GluGluToBulkGravitonToHHTo4B_M-900_narrow_13TeV-madgraph.root");

	skimTreeJERBase("BulkGrav_M-1000_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-1200_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-1400_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-1600_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-1800_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-2000_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-2500_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-3000_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-4000_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root");
	skimTreeJERBase("BulkGrav_M-4500_0","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root");

	
	skimTreeJERBase("Radion750","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/RadionTohhTohbbhbb/GluGluToRadionToHHTo4B_M-750_narrow_13TeV-madgraph.root");
	skimTreeJERBase("Radion800","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/RadionTohhTohbbhbb/GluGluToRadionToHHTo4B_M-800_narrow_13TeV-madgraph.root");
	skimTreeJERBase("Radion900","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/RadionTohhTohbbhbb/GluGluToRadionToHHTo4B_M-900_narrow_13TeV-madgraph.root");
	
	skimTreeJERBase("Radion1000","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root");
	skimTreeJERBase("Radion1200","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root");
	skimTreeJERBase("Radion1400","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root");
	skimTreeJERBase("Radion1600","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root");
	skimTreeJERBase("Radion2500","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root");
	skimTreeJERBase("Radion3000","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root");
	skimTreeJERBase("Radion3500","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root");
	skimTreeJERBase("Radion4500","/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/doublebtagv4_JECv3/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root");
	skimTreeJERBase("Radion1800","/data7/chchen/Mar222017/R1800.root");
	skimTreeJERBase("Radion2000","/data7/chchen/Mar222017/R2000.root");
}