#!/bin/bash
cd ~/Programs/nptool/Projects/e793s;
cmake ./;
make -j6;
#npanalysis --definition Exp -R RunToTreat_PartII.txt -E Reaction/47Kdp_08Nov.reaction -D Detector/mugast_08Nov.detector -C Calibration.txt -O $1;
npanalysis --definition Exp -R RunToTreat_PartI.txt -E Reaction/47Kdp_08Nov.reaction -D Detector/mugast_18May22_TestingThickerTarget.detector -C Calibration.txt -O $1;
