#!/bin/bash 
for i in {418..441}
do
  npanalysis -D s034.detector -C calibration.txt -T root/mrdc/gamma/run$i.root  RawTree -O gamma/test_run$i
done
