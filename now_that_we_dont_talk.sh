rm run_file
for i in $(seq 1 1 3); do
seed=$RANDOM
echo "./dym-mod-metro ./groups/S1080ct 4 32 32 0 8.407 0 -1.65 0.00 0 0 0 0 0 0 0 0.00 0.00 0 0 0 $seed > data_s1080/out_bf8.407_ba-1.65_nt32_nx32_s${seed}.log" >> run_file
done
