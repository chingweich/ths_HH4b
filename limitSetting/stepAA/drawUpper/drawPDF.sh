#!/bin/bash

for ((k=10; k<=10; k=k+1 ))
do
    root -l plot_Asymptotic_HHbbbb.C++\(\"pdf/pdf_${k}\",\"pdf/pdf_${k}\",20,1\)
done
