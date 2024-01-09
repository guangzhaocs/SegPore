#!/bin/bash -l

cd ../demo
rm Hek293T_config.yml

echo "creat folders ..."
mkdir 0_reference
mkdir 1_fast5
mkdir 2_fastq
mkdir 3_nanopolish
mkdir 4_hhmm
mkdir 5_align

mv demo.fa 0_reference/
mv demo.gtf 0_reference/

echo "combine single fast5 files to multi_fast5 ..."
single_to_multi_fast5 --input_path data/HEK293T-WT-rep1/fast5/ --save_path 1_fast5/ --filename_base HEK293T_WT_rep1_FAK27249_demo --batch_size 4000 --recursive

mkdir 2_fastq/base_calling_res
mv data 0_origin_data

echo "0_data_proc.sh done ..."

