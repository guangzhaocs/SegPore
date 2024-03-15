cd ../code/HierHmmCuda
for((idx=0;idx<=3999;idx++));  
do 
   ./hmm_one_read ../../demo/4_hhmm/hhmm_input/HEK293T_WT_rep1_FAK27249_demo_0 ../../demo/4_hhmm/hhmm_output/HEK293T_WT_rep1_FAK27249_demo_0  $idx
done
