#!/usr/bin/env bash

cd ..

/scratch/work/chengg1/ont-guppy-cpu/bin/guppy_basecaller -c /scratch/work/chengg1/ont-guppy-cpu/data/rna_r9.4.1_70bps_hac.cfg --num_callers 20 --cpu_threads_per_caller 20 -i demo/1_fast5 -s demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0
echo "base calling done ..."

mkdir demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0/all
cat demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0/pass/*.fastq > demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0/all/pass.fastq
cat demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0/fail/*.fastq > demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0/all/fail.fastq
cat demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0/all/*.fastq > demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0.fastq
rm -rf demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0






