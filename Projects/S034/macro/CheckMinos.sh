# convert physics run
ridf2nptool raw/run0582.ridf.gz root/mrdc/run582.root db/s034_list.txt 
# convert an empty target run
ridf2nptool raw/run0824.ridf.gz root/mrdc/run824.root db/s034_list.txt 
# analysis of both run
npcompilation -l && npanalysis -D s034.detector -C calibration.txt -T root/mrdc/run582.root RawTree -O run582CheckMinos -L 100000
npcompilation -l && npanalysis -D s034.detector -C calibration.txt -T root/mrdc/run824.root RawTree -O run824CheckMinos -L 100000
# Plot the results
root macro/CheckMinos.cxx
