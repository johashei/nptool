

for i in {1..30}
do
  nohup npanalysis  -T root/simulation/simu_tpad_d5mm$i.root SimulatedTree -O simu_tpad_d5mm$i &
done
