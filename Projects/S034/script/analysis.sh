

for i in {1..20}
do
  nohup npanalysis -D minos.detector -E reaction/pp.reaction -T root/simulation/simu_tpad_$i.root SimulatedTree -O simu_tpad_$i &
done
