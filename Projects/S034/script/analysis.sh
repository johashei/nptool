

for i in {1..30}
do
  nohup npanalysis -D minos_dali.detector -E reaction/pp.reaction -T root/simulation/simu_tpad_$i.root SimulatedTree -O simu_tpad_$i &
done
