#!/bin/bash
cd ~/Programs/nptool/Projects/e793s;
cmake ./;
make -j6;

#====================================================
#rfile='Reaction/47Kdp_08Nov.reaction'
#rfile='Reaction/47Kdd_08Nov.reaction'
#rfile='Reaction/47Kpp_08Nov.reaction'
#rfile='Reaction/47Kdt_08Nov.reaction'
#----------------------------------------------------
#rfile='Reaction/47Kdp_11Jul22.reaction'
#rfile='Reaction/47Kdt_11Jul22.reaction'
rfile='Reaction/47Kdd_11Jul22.reaction'
#rfile='Reaction/47Kpp_11Jul22.reaction'
#====================================================
#dfile='Detector/mugast_08Nov.detector'
#----------------------------------------------------
dfile='Detector/mugast_11Jul22.detector'
#====================================================


ptI='_PartI'
ptII='_PartII'

npanalysis --definition Exp -R RunToTreat_PartI.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptI;
npanalysis --definition Exp -R RunToTreat_PartII.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptII;
