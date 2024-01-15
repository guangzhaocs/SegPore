# SegPore Tutorials
```
git clone https://github.com/guangzhaocs/SegPore.git
cd SegPore
```
## Environment setup

```
conda env create -f environment.yml
```
The environment will get installed in their default conda environment path. If you want to specify a different install path than the default for your system, add `-p`:
```
conda env create -f environment.yml -p /home/user/anaconda3/envs/env_name
```
Finally, activate the environment:
```
conda activate segpore_env
```

Download demo data from [xPore](https://xpore.readthedocs.io/en/latest/index.html). Our SegPore is the single-mode method, so we only use the WT data in the demo data.
```
cd SegPore
wget https://zenodo.org/record/5162402/files/demo.tar.gz
tar -xvf demo.tar.gz
cd scripts
sh 0_data_proc.sh
```
After running `0_data_proc.sh`, your folders will like this: 
```
SegPore
  | -- environment.yml
  | -- Segpore
          | -- scripts
          | -- code
          | -- demo
                 | -- 0_origin_data
                 | -- 0_reference
                 | -- 1_fast5
                 | -- 2_fastq
                 | -- 3_nanopolish
                 | -- 4_hhmm
                         | -- hhmm_init
                         | -- hhmm_input
                         | -- hhmm_output
                         | -- hhmm_final
                 | -- 5_align
```
### Step 1: Basecalling, mapping and preprocessing
```
sh 1_basecalling.sh
sh 1_nanopolish.sh
```
Next, SegPore will use the PloyA to standardize the reads, so we only keep the `PASS` reads. 
```
wc -l ../demo/3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya.tsv
# 165 ../demo/3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya.tsv
wc -l ../demo/3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya_pass.tsv
# 104 ../demo/3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya_pass.tsv
```
The nanopolish eventalign also maps one or multiple events to one kmer, so here we combine the eventalign results:
```
contig	position	reference_kmer	read_index	strand	event_index	event_level_mean	event_stdv	event_length	model_kmer	model_mean	model_stdv	standardized_level	start_idx	end_idx
ENST00000273480.3	14	TAGGC	0	t	4	79.04	0.758	0.00299	TAGGC	93.67	7.84	-1.67	53396	53405
ENST00000273480.3	14	TAGGC	0	t	5	95.28	3.837	0.00730	TAGGC	93.67	7.84	0.18	53374	53396
```
After combining:
```
read_idx	contig	pos	kmer	mean	start_idx	end_idx
0	ENST00000273480.3	16	TAGGC	90.565	53374	53405
```
The mean after combining is the weighted average of the length:
```
53405 - 53396 = 9
53396 - 53374 = 22
79.04 * 9 / (9+ 22) + 95.28* 22 / (9+ 22) = 90.565
```
In nanopolish eventalign and summary results, one read may be processed multiple times (one read_name may correspond to multi read_idx). 
```
6cd0cd47-db93-4e73-99cc-e91c68f45268 [91, 157]
a7aa6921-712f-481e-8496-4a963618b786 [134, 163]
```
So when standardizing the fast5 file, one read_name only choose one read_idx. The example of a standardized multi-fast5 file is as follows:

<div align="center">
<img src=../pics/standardized_multi_fast5_example.jpg width=60% />
</div>

### Step 2: Hierarchical hidden Markov model (HHMM) for signal segmentation
Firstly, prepare the input of HHMM.
```
sh 2_hhmm_prepare.sh
```
If the `code/HierHmmCuda/hmm_one_read` can not run on your cluster, you can compile it as follows:
```
module load gcc
module load cuda
cd ../code/HierHmmCuda
srun --mem=1G --time=00:10:00 --gres=gpu nvcc -o hmm_one_read hmm_one_read.cu
cd ../../scripts
```
Next, run HHMM on CUDA:
```
sh 2_hhmm_GPU.sh
```
The output illustration of HHMM:
<div align="center">
<img src=../pics/hhmm_output.jpg width=80% />
</div>

Finally, generate the final output:

```
sh 2_hhmm_post_proc.sh
```

### Step 3: Alignment of signal segments with reference sequence
If the `code/Resquiggle/resquiggle_2D` can not run on your cluster, you can compile it as follows:
```
module load gcc
cd ../code/Resquiggle
g++ -O3 -Werror -Wall --pedantic -std=c++17 -march=native -fopenmp -o resquiggle_2D main_2D.cpp
cd ../../scripts
```
Run alignment algorithm:
```
sh 3_alignment.sh
```

The output file `segpore_eventalign_2D.txt` is as follows:
```
0	da848fa9-1322-4fea-b550-7efb32b014b6	ENST00000273480.3	TAGGC	16	80.5237	1.3826	53397	53423	26	
0	da848fa9-1322-4fea-b550-7efb32b014b6	ENST00000273480.3	TAGGC	16	95.779	2.7758	53376	53393	17	
0	da848fa9-1322-4fea-b550-7efb32b014b6	ENST00000273480.3	TAGGC	16	102.227	2.8678	53339	53372	7	
0	da848fa9-1322-4fea-b550-7efb32b014b6	ENST00000273480.3	AGGCA	17	111.711	4.7163	53278	53332	36	
0	da848fa9-1322-4fea-b550-7efb32b014b6	ENST00000273480.3	AGGCA	17	103.782	5.6484	53260	53274	13	
0	da848fa9-1322-4fea-b550-7efb32b014b6	ENST00000273480.3	AGGCA	17	116.527	5.3756	53229	53256	25	
0	da848fa9-1322-4fea-b550-7efb32b014b6	ENST00000273480.3	AGGCA	17	100.959	5.2579	53216	53225	4	
0	da848fa9-1322-4fea-b550-7efb32b014b6	ENST00000273480.3	AGGCA	17	111.676	3.778	53201	53210	7	
```
After combining:
```
0	ENST00000273480.3	16	TAGGC	810	88.749	53339	53423	50
0	ENST00000273480.3	17	AGGCA	164	110.284	53135	53332	108
```
The density of all `GGACT` is as follows:
<div align="center">
<img src=../pics/SegPore_GGACT.jpg width=40% />
</div>

### Step 4: GMM to update 5mer parameter table
Fix the mean of the first component of GMM. For GGACT, the fixed mean is 123.83.
```
sh 4_gmm.sh
```
Use the results of GMM to update the 5mer parameter table and iteratively run Step 3 and Step 4.

In this demo experiment, the 5mer parameter table is `demo/0_reference/model_kmer_m6A_without_header.csv`. In the first round, the fixed mean is from the kmer_model (https://github.com/nanoporetech/kmer_models) of ONT.   

