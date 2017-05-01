which=$1
sigHypIn=$2
subtr=$3

sigHyp=$(expr $((sigHypIn+0)))
# 0 : run them all
# 10 = radion_2cat, 11 = radion_4btag_cat0, 12 = radion_3btag_HPHP_cat1 
# 20 = graviton_2cat, 21 = graviton_4btag_cat0, 22 = graviton_3btag_HPHP_cat1 
if [ $sigHyp.eq.0 ]
then
    echo "run all"
elif [ $sigHyp.eq.10 ]
then
    echo "radion 2 categories together"
elif [ $sigHyp.eq.11 ]
then
    echo "radion 4btag category"
elif [ $sigHyp.eq.12 ]
then
    echo "radion 3btag category with HPHP jets"
elif [ $sigHyp.eq.20 ]
then
    echo "graviton 2 categories together"
elif [ $sigHyp.eq.21 ]
then
    echo "graviton 4btag category"
elif [ $sigHyp.eq.22 ]
then
    echo "graviton 3btag category with HPHP jets"
fi

compar=$((sigHyp-10 > -1))

if  [ $compar -eq 1 ]
then
    echo "here"
    root -l -b -q plot_Asymptotic_HHbbbb.C\(\"$which\"\,${sigHyp}\,${subtr}\)
else
    root -l -b -q plot_Asymptotic_HHbbbb.C\(\"$which\"\,10\,${subtr}\)
    root -l -b -q plot_Asymptotic_HHbbbb.C\(\"$which\"\,11\,${subtr}\)
    root -l -b -q plot_Asymptotic_HHbbbb.C\(\"$which\"\,12\,${subtr}\)
    root -l -b -q plot_Asymptotic_HHbbbb.C\(\"$which\"\,20\,${subtr}\)
    root -l -b -q plot_Asymptotic_HHbbbb.C\(\"$which\"\,21\,${subtr}\)
    root -l -b -q plot_Asymptotic_HHbbbb.C\(\"$which\"\,22\,${subtr}\)
fi