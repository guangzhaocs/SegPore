#!/usr/bin/env bash

cd ..

guppy_basecaller -c rna_r9.4.1_70bps_hac.cfg --num_callers 20 --cpu_threads_per_caller 20 -i demo/1_fast5 -s demo/2_fastq/base_calling_res
echo "base calling done ..."

mkdir demo/2_fastq/base_calling_res/all
cat demo/2_fastq/base_calling_res/pass/*.fastq > demo/2_fastq/base_calling_res/all/pass.fastq
cat demo/2_fastq/base_calling_res/fail/*.fastq > demo/2_fastq/base_calling_res/all/fail.fastq
cat demo/2_fastq/base_calling_res/all/*.fastq > demo/2_fastq/HEK293T_WT_rep1_FAK27249_demo_0.fastq
rm -rf demo/2_fastq/base_calling_res
echo "post-base-calling processing done ..."
