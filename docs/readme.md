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
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       GGCAC   18      105.39  2.6755  49149   49181   19
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       GGCAC   18      110.209 1.2636  49128   49145   17
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       GGCAC   18      112.569 2.5778  49115   49124   9
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       GCACC   19      71.4983 1.6224  49084   49110   14
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       CACCA   20      79.83   1.9821  49073   49080   7
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       CACCA   20      78.9216 0.8394  49055   49069   14
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       ACCAC   21      76.8632 1.1586  49042   49051   9
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       CCACG   22      77.0395 1.334   48989   49038   49
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       CACGG   23      79.8887 1.1274  48978   48985   7
1       00831d10-b2ea-4388-89ee-6fc63448e87f    ENST00000273480.3       CACGG   23      77.3341 1.2753  48965   48974   9
```

### Step 4: GMM to update 5mer parameter table
```
sh todo.sh
```

