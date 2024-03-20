#!/bin/bash -l

cd ..
python code/2_0_generate_dataset_from_hmm_final.py

cd code/Resquiggle
chmod 777 resquiggle_2D 
./resquiggle_2D ../../demo/5_align/HEK293T_WT_rep1_FAK27249_demo_0 segpore_eventalign_2D.txt ../../demo/0_reference/model_kmer_m6A_without_header.csv

cd ../..
python code/2_1_combine_segpore_eventalign.py
cat demo/5_align/HEK293T_WT_rep1_FAK27249_demo_0/segpore_eventalign_2D_combined.txt | awk -F "\t" '{if($4=="GGACT" ){print $0}}' > demo/5_align/HEK293T_WT_rep1_FAK27249_demo_0/segpore_eventalign_2D_combined_GGACT.txt
