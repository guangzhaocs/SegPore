#!/usr/bin/env bash

cd ../demo

nanopolish index --directory=1_fast5 2_fastq/HEK293T_WT_rep1_FAK27249_demo_0.fastq
echo " *** nanopolish index is over."

cd 2_fastq
minimap2 -ax map-ont --MD -t 3 --secondary=no ../0_reference/demo.fa HEK293T_WT_rep1_FAK27249_demo_0.fastq | samtools sort -o reads-ref.sorted.bam -T reads.tmp
echo " *** minimap2 is over."

samtools index reads-ref.sorted.bam
samtools quickcheck reads-ref.sorted.bam
cd ..
echo " *** samtools is over."

nanopolish eventalign --reads 2_fastq/HEK293T_WT_rep1_FAK27249_demo_0.fastq \
--bam 2_fastq/reads-ref.sorted.bam \
--genome 0_reference/demo.fa \
--signal-index \
--scale-events \
--summary 3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_summary.txt \
--threads 32 > 3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_eventalign.txt
echo " *** nanopolish eventalign is over."

nanopolish polya --threads=10 --reads=2_fastq/HEK293T_WT_rep1_FAK27249_demo_0.fastq --bam=2_fastq/reads-ref.sorted.bam --genome=0_reference/demo.fa > 3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya.tsv
echo " *** polyA is over."


grep -E 'PASS|readname' 3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya.tsv > 3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya_pass.tsv
echo " *** filter polyA is over (grep PASS)."


cd ..
python code/0_1_add_strand_and_order.py
python code/0_2_combine_nanopolish_eventalign.py
cp demo/1_fast5/HEK293T_WT_rep1_FAK27249_demo_0.fast5 demo/1_fast5/HEK293T_WT_rep1_FAK27249_demo_0.standardized.fast5
python code/0_3_standardize_multi_fast5.py 
echo " *** Post-processing of nanopolish is over."