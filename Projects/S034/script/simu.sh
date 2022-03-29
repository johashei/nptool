

for i in {1..30}
do
  nohup npsimulation -D minos_dali_short.detector -E reaction/pp.reaction -B script/run.mac -O simu_tpad_short_$i --random-seed $i & 
done
