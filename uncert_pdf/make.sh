#!/bin/bash
a=(750 800 900 1000 1200 1600 1800 2000 2500 3000 4000 4500)
for ((k=0; k<3; k=k+1 ))
do
	root -q -l pdfSkimTree.C\(\"BulkGrav${a[$k]}\",\"/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/BulkGravTohhTohbbhbb/GluGluToBulkGravitonToHHTo4B_M-${a[$k]}_narrow_13TeV-madgraph.root\"\)
done
for ((k=3; k<12; k=k+1 ))
do
	root -q -l pdfSkimTree.C\(\"BulkGrav${a[$k]}\",\"/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-${a[$k]}_13TeV-madgraph.root\"\)
done

b=(750 800 900 1000 1200 1600 3500 4500)
for ((k=0; k<3; k=k+1 ))
do
	root -q -l pdfSkimTree.C\(\"Radion${b[$k]}\",\"/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/RadionTohhTohbbhbb/GluGluToRadionToHHTo4B_M-${b[$k]}_narrow_13TeV-madgraph.root\"\)
done
for ((k=3; k<12; k=k+1 ))
do
	root -q -l pdfSkimTree.C\(\"Radion${b[$k]}\",\"/data7/syu/NCUGlobalTuples/80X_Moriond/80X_puppi/pdfscaleWeights/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-${b[$k]}_13TeV-madgraph.root\"\)
done