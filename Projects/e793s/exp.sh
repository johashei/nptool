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
#rfile='Reaction/47Kdd_11Jul22.reaction'
#rfile='Reaction/47Kpp_11Jul22.reaction'
#----------------------------------------------------
#rfile='Reaction/47Kdp_30Aug22_2p06um.reaction'
#----------------------------------------------------
#rfile='Reaction/47Kdp_14Oct22.reaction'
#rfile='Reaction/47Kdp_14Oct22_2.reaction'
#----------------------------------------------------
rfile='Reaction/47Kdp_18Oct22.reaction'
#rfile='Reaction/47Kdt_18Oct22.reaction'
#====================================================
#dfile='Detector/mugast_08Nov.detector'
#----------------------------------------------------
#dfile='Detector/mugast_11Jul22.detector'
#dfile='Detector/mugast_Test2um.detector'
#dfile='Detector/mugast_17Jul22_MM+200mm.detector'
#dfile='Detector/mugast_17Jul22_MM+203mm.detector'
#dfile='Detector/mugast_17Jul22_MM+210mm.detector'
#dfile='Detector/mugast_17Jul22_MM+190mm.detector' #DO ME!!!
#----------------------------------------------------
#dfile='Detector/mugast_01Sep22.detector'
#----------------------------------------------------
#dfile='Detector/mugast_14Oct22.detector'
#dfile='Detector/mugast_14Oct22_2.detector'
dfile='Detector/mugast_18Oct22.detector'
#====================================================


smallTest='_Test'
ptI='_PartI'
ptII='_PartII'
ptA='_PartA'
ptB='_PartB'
ptC='_PartC'

echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

#npanalysis --definition Exp -R RunToTreat_PartI.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptI;
#npanalysis --definition Exp -R RunToTreat_PartI_2.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptI;
#npanalysis --definition Exp -R RunToTreat_PartII.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptII;
npanalysis --definition Exp -R RunToTreat_PartA.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptA;
npanalysis --definition Exp -R RunToTreat_PartB.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptB;
npanalysis --definition Exp -R RunToTreat_PartC.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptC;

echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 'CURRENTLY EXCLUDING RUNS 51 AND 52!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#npanalysis --definition Exp --definition ExcludeThePoor -R RunToTreat_PartI.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptI;
#npanalysis --definition Exp --definition ExcludeThePoor -R RunToTreat_PartII.txt -E $rfile -D $dfile -C Calibration.txt -O $1$ptII;
