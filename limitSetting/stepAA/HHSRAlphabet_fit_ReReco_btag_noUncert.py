import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import scipy
import pdb

# Our functions:
import Alphabet_Header
from Alphabet_Header import *
import Plotting_Header
from Plotting_Header import *
import Converters
from Converters import *
import Distribution_Header
from Distribution_Header import *
import Alphabet
from Alphabet import *

lowBin=1100
highBin=3000

def GetNom(file_string):
	tempFile = TFile(file_string)
	tempHist = tempFile.Get("CountWeighted")
	norm = tempHist.GetBinContent(1)
	tempFile.Close()
	return norm

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--B', '--binsize', metavar='Bin', type='string', dest='bin', default="15")

parser.add_option('--T2', '--Selection', metavar='T32', type='string', dest='tightpre')
parser.add_option('--T1', '--Cut', metavar='T13', type='float', dest='passcut', default = 0.8)
parser.add_option('--Fail', '--FailCut', metavar='T33', type='float', dest='failcut', default = 0.3)

parser.add_option('--N', '--name', metavar='Name', type='string', dest='name', default="test")
parser.add_option('--L', '--lumi', metavar='Name', type='float', dest='lumi', default="35900")

parser.add_option("--data", action="store_true", dest="isData", default=True)
parser.add_option("--qcd", action="store_false", dest="isData")

parser.add_option("--quad", action="store_false", dest="Linear", default=False)
parser.add_option("--lin", action="store_true", dest="Linear")

parser.add_option("--blind", action="store_false", dest="Truth", default=False)
parser.add_option("--unblind", action="store_true", dest="Truth")

parser.add_option("--finebins", action="store_false", dest="finebins", default=False)
parser.add_option("--dijetbins", action="store_true", dest="finsbines")

parser.add_option("--log", action="store_true", dest="log", default=False)
parser.add_option("--nolog", action="store_false", dest="log")

parser.add_option("--sig", action="store_true", dest="Sig", default=True)
parser.add_option("--nosig", action="store_false", dest="Sig")

parser.add_option('-I', '--inject', metavar='Inj', type='string', dest='inject', default="none")

parser.add_option('--LL', action="store_true", dest='LL_DoubleB_Region', default=True)
parser.add_option('--TT', action="store_false", dest='LL_DoubleB_Region')


parser.add_option('--workspace', metavar='WSPC', type='string', dest='workspace', default="alphabet")
(Options, args) = parser.parse_args()

preselection    =       "&jet2pt>300&jet1pt>300&abs(jet1eta-jet2eta)<1.3 & dijetmass_softdrop_corr>750&abs(jet1eta)<2.4&abs(jet2eta)<2.4 " 
tauselection = "&jet1_puppi_tau21<0.55&jet2_puppi_tau21<0.55" 
#triggerselection = "&(HLT_PFHT900_v==1||HLT_PFHT800_v==1||HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v==1||HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v==1||HLT_AK8PFJet360_V==1||HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v==1||HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v==1)"
triggerselection = "&(HLT_PFHT900||HLT_PFHT800||HLT_AK8PFHT650_TrimR0p1PT0p03Mass50||HLT_AK8PFHT700_TrimR0p1PT0p03Mass50||HLT_AK8PFJet360_TrimMass30||HLT_PFHT650_WideJetMJJ900DEtaJJ1p5||HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20)"
#triggerselection = "&1"
#triggerselection = "&(HLT_PFHT800_v==1||HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v==1||HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v==1||HLT_AK8PFJet360_V==1||HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v==1)"

TightPre 		=	Options.tightpre + preselection + tauselection 
if Options.isData : TightPre = TightPre+triggerselection
if Options.LL_DoubleB_Region:
  TightPreAT              =       TightPre + "&jet2bbtag<0.8"
  TightAT                 =       TightPreAT + "&jet1_puppi_msoftdrop_TheaCorr>105&jet1_puppi_msoftdrop_TheaCorr<135&(jet1bbtag<0.3)"
  TightT          =       TightPre + "&jet1_puppi_msoftdrop_TheaCorr>105&jet1_puppi_msoftdrop_TheaCorr<135&(jet1bbtag>"+str(Options.passcut)+")"
else:
  TightPreAT              =       TightPre
  TightAT                 =       TightPreAT + "&jet1_puppi_msoftdrop_TheaCorr>105&jet1_puppi_msoftdrop_TheaCorr<135&(jet1bbtag<0.3)"
  TightT          =       TightPre + "&jet1_puppi_msoftdrop_TheaCorr>105&jet1_puppi_msoftdrop_TheaCorr<135&(jet1bbtag>"+str(Options.passcut)+")"

if Options.LL_DoubleB_Region:
  TightT2         = "1"+preselection + tauselection +triggerselection +" & jet2_puppi_msoftdrop_TheaCorr > 105 & jet2_puppi_msoftdrop_TheaCorr < 135  & (!( jet1bbtag > 0.8 & jet2bbtag > 0.8))& jet2bbtag > 0.3 & jet1_puppi_msoftdrop_TheaCorr>105&jet1_puppi_msoftdrop_TheaCorr<135 & jet1bbtag > 0.3"  #orthogonality with TT
else:
  TightT2         = "1"+preselection + tauselection +triggerselection +" & jet2_puppi_msoftdrop_TheaCorr > 105 & jet2_puppi_msoftdrop_TheaCorr < 135  & jet2bbtag > 0.8 & jet1_puppi_msoftdrop_TheaCorr>105&jet1_puppi_msoftdrop_TheaCorr<135 & jet1bbtag > 0.8"

#if Options.LL_DoubleB_Region:
#  highBin = 2843
#else:
#  highBin = 2671

Options.finebins = True

if Options.finebins:
	binBoundaries=[]
	for i in range(0,highBin-lowBin):	
		binBoundaries.append(lowBin+i*1)
else:
	binBoundaries =[750, 775, 800, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4681, 4853, 5025]

variable = "dijetmass_softdrop_corr"
variable2 = "dijetmass_softdrop_corr"
#variable = "dijetmass_corr"


sigpath = "/afs/cern.ch/work/c/chchen/public/dbtSFTriggerBits/"


if Options.workspace == "alphabet":
	print "creating workspace and datacard: ALPHABET"

	mass=[1200,1400,1600,1800,2000, 2500, 3000]
	for m in mass:
		print str(m)
		SF_tau21 = 1.03*1.03
		UD = ['Up','Down']

		output_file = TFile("outputs/datacards/HH_mX_"+Options.name+"_%s"%(m)+"_13TeV.root", "RECREATE")
		vh=output_file.mkdir("vh")
		vh.cd()

		Signal_mX = TH1F("Signal_mX_%s_"%(m)+Options.name, "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_antitag = TH1F("Signal_mX_antitag_%s"%(m)+Options.name, "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_trig_up = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_trigUp", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_trig_down = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_trigDown", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_btag_up = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_btagUp", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_btag_down = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_btagDown", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_pu_up = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_puUp", "", len(binBoundaries)-1, array('d',binBoundaries))
	 	Signal_mX_pu_down = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_puDown", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_FJEC_Up = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_JECUp", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_FJEC_Down = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_JECDown", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_FJER_Up = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_JERUp", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_FJER_Down = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_JERDown", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_MJEC_Up = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_massJECUp", "", len(binBoundaries)-1, array('d',binBoundaries))
		Signal_mX_MJEC_Down = TH1F("Signal_mX_%s_"%(m)+Options.name+"_CMS_eff_massJECDown", "", len(binBoundaries)-1, array('d',binBoundaries))
		
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX, variable2, TightT,"puWeights*dbtSFup/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_antitag, variable2, TightAT,"puWeights*dbtSFup/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_btag_up, variable2, TightT,"puWeights*dbtSFup/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_btag_down, variable2, TightT,"puWeights*dbtSFdown/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_trig_up, variable2, TightT,"trigWeightUp_Update*puWeights*dbtSF/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_trig_down, variable2, TightT,"trigWeightDown_Update*puWeights*dbtSF/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_pu_up, variable2, TightT,"puWeightsUp*dbtSF/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_pu_down, variable2, TightT,"puWeightsDown*dbtSF/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_FJEC_Up, variable2, TightT,"puWeights*dbtSF*(1+jet1JECup)*(1+jet2JECup)/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_FJEC_Down, variable2, TightT,"puWeights*dbtSF*(1-jet1JECdown)*(1-jet2JECdown)/1.")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_FJER_Up, variable2, TightT,"puWeights*dbtSF*jet1JERup*jet2JERup/(1.*jet1JERcentral*jet2JERcentral)")
                quickplot(sigpath+"BulkGrav_M-%s_0.root"%(m), "mynewTree", Signal_mX_FJER_Down, variable2, TightT,"puWeights*dbtSF*jet1JERdown*jet2JERdown/(1.*jet1JERcentral*jet2JERcentral)")
		

		norm = GetNom(sigpath+"BulkGrav_M-%s_0.root"%(m))

		btaglnNup= 1. + abs(Signal_mX_btag_up.GetSumOfWeights()-Signal_mX_btag_down.GetSumOfWeights())/(2.*Signal_mX_btag_up.GetSumOfWeights())
                btaglnNdown= 1. - abs(Signal_mX_btag_up.GetSumOfWeights()-Signal_mX_btag_down.GetSumOfWeights())/(2.*Signal_mX_btag_up.GetSumOfWeights())
		PUlnN= 1. + abs(Signal_mX_pu_up.GetSumOfWeights()-Signal_mX_pu_down.GetSumOfWeights())/(2.*Signal_mX.GetSumOfWeights())

		Signal_mX.Scale(SF_tau21*Options.lumi*0.01/norm)
		Signal_mX_antitag.Scale(SF_tau21*Options.lumi*0.01/norm)
		Signal_mX_btag_up.Scale(SF_tau21*Options.lumi*0.01/norm)
		Signal_mX_btag_down.Scale(SF_tau21*Options.lumi*0.01/norm)
		Signal_mX_trig_up.Scale(SF_tau21*0.01*Options.lumi/norm)
		Signal_mX_trig_down.Scale(SF_tau21*0.01*Options.lumi/norm)
		Signal_mX_pu_up.Scale(Options.lumi*SF_tau21*0.01/norm)
		Signal_mX_pu_down.Scale(Options.lumi*SF_tau21*0.01/norm)

		if Options.LL_DoubleB_Region:
		  if str(m) == "1200":
		    PDFup = 0.996
		    PDFdown = 0.997
		  elif str(m) == "1400":
                    PDFup = 0.998
                    PDFdown = 0.998
                  elif str(m) == "1600":
                    PDFup = 0.998
                    PDFdown = 0.998
                  elif str(m) == "1800":
                    PDFup = 0.998
                    PDFdown = 0.998
                  elif str(m) == "2000":
                    PDFup = 0.998
                    PDFdown = 0.998
                  elif str(m) == "2500":
                    PDFup = 1.004
                    PDFdown = 1.003
                  elif str(m) == "3000":
                    PDFup = 0.987
                    PDFdown = 0.988
		else:
                  if str(m) == "1200":
                    PDFup = 0.996
                    PDFdown = 0.996
                  elif str(m) == "1400":
                    PDFup = 0.997
                    PDFdown = 0.998
                  elif str(m) == "1600":
                    PDFup = 0.996
                    PDFdown = 0.997
                  elif str(m) == "1800":
                    PDFup = 0.993
                    PDFdown = 0.994
                  elif str(m) == "2000":
                    PDFup = 0.993
                    PDFdown = 0.994
                  elif str(m) == "2500":
                    PDFup = 0.990
                    PDFdown = 0.991 
                  elif str(m) == "3000":
                    PDFup = 1.003
                    PDFdown = 1.003


                HTaggingUnc = (1. - math.exp(-0.125052 + 32.5054/(float(m)/2)))*2+ 1.
		MJEClnN= 1.02 ## add variation from ntuples
		print "Signal_mX_FJER_Up.GetSumOfWeights()"
		print Signal_mX_FJER_Up.GetSumOfWeights()
		print "Signal_mX_FJER_Down.GetSumOfWeights()"
		print Signal_mX_FJER_Down.GetSumOfWeights()

		FJEClnN= 1. + abs(Signal_mX_FJEC_Up.GetSumOfWeights()-Signal_mX_FJEC_Down.GetSumOfWeights())/(2.*Signal_mX_FJEC_Up.GetSumOfWeights())
		FJERlnN= 1. + abs(Signal_mX_FJER_Up.GetSumOfWeights()-Signal_mX_FJER_Down.GetSumOfWeights())/(2.*Signal_mX_FJER_Up.GetSumOfWeights())
		#TRIGlnN= 1. +abs(Signal_mX_trig_up.GetSumOfWeights()-Signal_mX_trig_down.GetSumOfWeights())/(2.*Signal_mX_trig_up.GetSumOfWeights())
                print"FJEClnN"
                print FJEClnN
		print"FJERlnN"
		print FJERlnN

		#signal_integral = Signal_mX.Integral()
		signal_integral = Signal_mX.Integral(Signal_mX.FindBin(lowBin),Signal_mX.FindBin(highBin))
		signal_integral_anti = Signal_mX_antitag.Integral(Signal_mX_antitag.FindBin(lowBin),Signal_mX_antitag.FindBin(highBin))
		vh.cd()
		Signal_mX.Write()
		Signal_mX_antitag.Write()
		output_file.Close()

	

		text_file = open("outputs/datacards/HH_mX_%s_"%(m)+Options.name+"_13TeV.txt", "w")

		data_integral = -1
		text_file.write("max    1     number of categories\n")
		text_file.write("jmax   1     number of samples minus one\n")
		text_file.write("kmax    *     number of nuisance parameters\n")
		text_file.write("-------------------------------------------------------------------------------\n")
		#text_file.write("shapes * * HH_mX_%s_"%(m)+Options.name+"_13TeV.root vh/$PROCESS vh/$PROCESS_$SYSTEMATIC\n")
		if Options.LL_DoubleB_Region:
		  text_file.write("shapes Signal_mX_%s_"%(m)+Options.name+"      HH4b w_signal_LL_%s.root      HH4b:signal_fixed_ \n"%(m))
                  text_file.write("shapes "+Options.name+"EST HH4b w_background_LL.root HH4b:bg_\n")
                  text_file.write("shapes data_obs   HH4b w_data_LL.root                HH4b:data_obs\n")
		else:
                  text_file.write("shapes Signal_mX_%s_"%(m)+Options.name+"      HH4b w_signal_TT_%s.root      HH4b:signal_fixed_ \n"%(m))
                  text_file.write("shapes "+Options.name+"EST HH4b w_background_TT.root HH4b:bg_\n")
                  text_file.write("shapes data_obs   HH4b w_data_TT.root                HH4b:data_obs\n")

		text_file.write("-------------------------------------------------------------------------------\n")
		text_file.write("bin                                            HH4b\n")
		text_file.write("observation                                    %f\n"%(data_integral))
		text_file.write("-------------------------------------------------------------------------------\n")
		text_file.write("bin                                             HH4b            HH4b\n")
		text_file.write("process                                          0      1\n")
		text_file.write("process                                         Signal_mX_%s_"%(m)+Options.name+"  "+Options.name+"EST\n")
		text_file.write("rate                                            %f  1.00\n"%(signal_integral))
		text_file.write("-------------------------------------------------------------------------------\n")
		#text_file.write("lumi_13TeV lnN                          1.025       -\n")	
	
		# text_file.write("CMS_eff_tau21_sf lnN                    1.30/0.74        -\n") #(0.028/0.979)
		# text_file.write("CMS_eff_Htag lnN                    %f       -\n"%(HTaggingUnc))
		# text_file.write("CMS_JEC lnN 		     %f        -\n"%(FJEClnN)) 	
		# text_file.write("CMS_massJEC lnN                 %f        -\n"%(MJEClnN))
		# text_file.write("CMS_eff_bbtag_sf lnN                    %f/%f       -\n"%(btaglnNup,btaglnNdown))
		# text_file.write("CMS_JER lnN                    %f        -\n"%(FJERlnN))
		# text_file.write("CMS_PU lnN                    %f        -\n"%(PUlnN))
#                text_file.write("CMS_eff_trig shapeN2           1.000   -\n")
#               text_file.write("CMS_eff_trig lnN           %f   -\n"%(TRIGlnN))	 	
		#text_file.write("CMS_scale"+Options.name+"_13TeV shapeN2                           -       1.000\n")
		# text_file.write("CMS_PDF_Scales lnN   %f/%f        -\n"%(PDFup,PDFdown))

#		for bin in range(0,len(binBoundaries)-1):
#			text_file.write("CMS_stat"+Options.name+"_13TeV_bin%s shapeN2                           -       1.000\n"%(bin))

                if Options.LL_DoubleB_Region:
                  text_file.write("R_LL param 0.208185541036 0.00503830220307 \n")
                  text_file.write("n_exp_binHH4b_proc_EST_LL_  rateParam HH4b "+Options.name+"EST @0*@1 bgSB_norm_LL,R_LL\n")
		else:
                  text_file.write("R_TT param 0.0505660156451 0.00320347783357 \n")
                  text_file.write("n_exp_binHH4b_proc_EST_TT_  rateParam HH4b "+Options.name+"EST @0*@1 bgSB_norm_TT,R_TT\n")

                text_file.close()


                text_filea = open("outputs/datacards/HH_mX_%s_"%(m)+Options.name+"_13TeV_fail.txt", "w")
                text_filea.write("imax    1     number of categories\n")
                text_filea.write("jmax    1     number of samples minus one\n")
                text_filea.write("kmax    *     number of nuisance parameters\n")
                text_filea.write("-------------------------------------------------------------------------------\n")
                if Options.LL_DoubleB_Region:
                  text_filea.write("shapes Signal_mX_antitag_%s_"%(m)+Options.name+"      HH4b w_signal_antitag_LL_%s.root      HH4b:signal_fixed_antitag_ \n"%(m))
                  text_filea.write("shapes "+Options.name+"EST_antitag HH4b w_background_LL.root HH4b:bgSB_\n")
                  text_filea.write("shapes data_obs   HH4b w_data_LL.root                HH4b:data_obs_sb\n")
		else:
                  text_filea.write("shapes Signal_mX_antitag_%s_"%(m)+Options.name+"      HH4b w_signal_antitag_TT_%s.root      HH4b:signal_fixed_antitag_ \n"%(m))
                  text_filea.write("shapes "+Options.name+"EST_antitag HH4b w_background_TT.root HH4b:bgSB_\n")
                  text_filea.write("shapes data_obs   HH4b w_data_TT.root                HH4b:data_obs_sb\n")

                text_filea.write("-------------------------------------------------------------------------------\n")
                text_filea.write("bin                                            HH4b                   \n")
                text_filea.write("observation                                    -1.0                           \n")
                text_filea.write("-------------------------------------------------------------------------------\n")
                text_filea.write("bin                                             HH4b            HH4b  \n")
                text_filea.write("process                                         Signal_mX_antitag_%s_"%(m)+Options.name+"  "+Options.name+"EST_antitag\n")
                text_filea.write("process                                          0      1     \n")
                text_filea.write("rate                                            %f    1.0000  \n"%(signal_integral_anti))
                text_filea.write("-------------------------------------------------------------------------------\n")

	        # text_filea.write("lumi_13TeV lnN                          1.025       -\n")

                # text_filea.write("CMS_eff_tau21_sf lnN                    1.30/0.74       -\n") #(0.028/0.979)
                # text_filea.write("CMS_eff_Htag lnN                    %f       -\n"%(HTaggingUnc))
                # text_filea.write("CMS_JEC lnN                 %f        -\n"%(FJEClnN))
                # text_filea.write("CMS_massJEC lnN                 %f        -\n"%(MJEClnN))
                # text_filea.write("CMS_eff_bbtag_sf lnN                    %f/%f       -\n"%(btaglnNdown,btaglnNup))
                # text_filea.write("CMS_JER lnN                    %f        -\n"%(FJERlnN))
                # text_filea.write("CMS_PU lnN                    %f        -\n"%(PUlnN))
#                text_filea.write("CMS_eff_trig shapeN2           1.000   -\n")
#               text_filea.write("CMS_eff_trig lnN           %f   -\n"%(TRIGlnN))
                if Options.LL_DoubleB_Region:	
		  text_filea.write("bgSB_norm_LL rateParam HH4b "+Options.name+"EST_antitag 3104.27800006  \n")
		else:
		  text_filea.write("bgSB_norm_TT rateParam HH4b "+Options.name+"EST_antitag 1171.45600007 \n")
	
		

                text_filea.close()

if Options.workspace == "fit":
	print "creating workspace and datacard: ALPHABET ASSISTED FIT"








