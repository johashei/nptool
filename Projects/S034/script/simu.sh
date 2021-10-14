

for i in {1..20}
do
  nohup npsimulation -D minos.detector -E reaction/pp.reaction -B run.mac -O simu_tpad_$i --random-seed $i & 
done
