# SegPore Tutorials
```
git clone https://github.com/guangzhaocs/SegPore.git
cd SegPore
```
## Environment setup

```
conda env create -f environment.yml
```
The environment will get installed in their default conda environment path. If you want to specify a different install path than the default for your system, add `-p`: `conda env create -f environment.yml -p /home/user/anaconda3/envs/env_name`

Finally, activate the environment, install `ont2cram` and `Guppy`:
```
conda activate segpore_env
pip3 install git+https://github.com/EGA-archive/ont2cram
```
For Guppy installation, some references: 
- https://help.nanoporetech.com/en/articles/6628042-how-do-i-install-stand-alone-guppy
- https://ontpipeline2.readthedocs.io/en/latest/GetStarted.html

## Environment test

```
nanopolish --version             #  nanopolish version 0.13.2  
minimap2 --version               #  2.24-r1122
guppy_basecaller --version       #  ... Version 6.0.1+652ffd1
single_to_multi_fast5 --version  #  4.0.2
gcc --version                    #  gcc (Spack GCC) 9.3.0
nvcc --version                   #  ... Cuda compilation tools, release 11.0, V11.0.194 ...
```

## Download Demo Data

We use the demo data from [xPore](https://xpore.readthedocs.io/en/latest/index.html). Our SegPore is the single-mode method, so we only use the WT data in the demo data.
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
## Example demo output
Maybe the SegPore workflow is complex, and we also provide the example demo output ([example_demo_output.tar.gz](https://drive.google.com/file/d/1y0IhL0zABeB0Wl1pmH92OYxeJZfgkeI8/view?usp=drive_link)). Following the SegPore workflow, if you can get the same outputs as the `example_demo_output.tar.gz`, the SegPore runs successfully.

## SegPore Workflow

### Step 1: Basecalling, mapping and preprocessing
```
sh 1_basecalling.sh
sh 1_nanopolish.sh
```
Below are some explanations about what happened in the script `1_nanopolish.sh`. If you are interested in this, you can read the following. Otherwise, ship to the next step directly!

SegPore uses the `PloyA` tail to standardize the reads, so only the `PASS` reads are kept. 
```
wc -l ../demo/3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya.tsv
# 165 ../demo/3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya.tsv
wc -l ../demo/3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya_pass.tsv
# 104 ../demo/3_nanopolish/HEK293T_WT_rep1_FAK27249_demo_0_polya_pass.tsv
```
The nanopolish eventalign also maps one or multiple events to one kmer, so here we also combine the eventalign results:
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
<img src=media/standardized_multi_fast5_example.jpg width=60% />
</div>

### Step 2: Hierarchical hidden Markov model (HHMM) for signal segmentation
#### Step 2.1: Firstly, prepare the input of HHMM.
```diff
! ATTENTION: This step may take a lot of time, and we will update the CUDA version soon.
```
```
sh 2_hhmm_prepare.sh
```
**Output**: If the script runs successfully, the folder `demo/4_hhmm/hhmm_init` will have four files: `HEK293T_WT_rep1_FAK27249_demo_0-border.csv`, `HEK293T_WT_rep1_FAK27249_demo_0-peaks.csv`, `HEK293T_WT_rep1_FAK27249_demo_0-readname.csv` and `HEK293T_WT_rep1_FAK27249_demo_0-signal.csv`, and the folder `demo/4_hhmm/hhmm_input/HEK293T_WT_rep1_FAK27249_demo_0` will have two files: `init_border.csv` and `signal.csv`.

#### Step 2.2: Next, run HHMM on CUDA:
```
sh 2_hhmm_GPU.sh
```
**Output**: If the script runs successfully, the folder `demo/4_hhmm/hhmm_output/HEK293T_WT_rep1_FAK27249_demo_0` will have two files: `res_border.csv` and `res_state.csv`.

If the above script has no errors, you can run next Step 2.3 direactly. If the above script has errors or the `code/HierHmmCuda/hmm_one_read` can not run on your cluster, you can re-compile it as follows, and then run `sh 2_hhmm_GPU.sh`.
```
cd ../code/HierHmmCuda
nvcc -o hmm_one_read hmm_one_read.cu
cd ../../scripts
```
The output illustration of HHMM:
<div align="center">
<img src=media/hhmm_output.jpg width=80% />
</div>

#### Step 2.3: Finally, generate the final output:
```
sh 2_hhmm_post_proc.sh
```
**Output**: If the script runs successfully, the folder `demo/4_hhmm/hhmm_final/HEK293T_WT_rep1_FAK27249_demo_0` will have the `mu`, `sigma` and `len` resluts files for `curr`, `prev` and `next`.

### Step 3: Alignment of signal segments with reference sequence
Run alignment algorithm:
```
sh 3_alignment.sh
```
**Output**: If the script runs successfully, the folder `demo/5_align/HEK293T_WT_rep1_FAK27249_demo_0` will contain `segpore_eventalign_2D.txt` and `segpore_eventalign_2D_combined.txt`. files.

If the above script has no errors, you can run next step direactly. If the above script has errors or the `code/Resquiggle/resquiggle_2D` can not run on your cluster, you can re-compile it as follows, and then run `sh 3_alignment.sh`.
```
cd ../code/Resquiggle
g++ -O3 -Werror -Wall --pedantic -std=c++17 -march=native -fopenmp -o resquiggle_2D main_2D.cpp
cd ../../scripts
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
After combining the results are as following (demo/5_align/HEK293T_WT_rep1_FAK27249_demo_0/segpore_eventalign_2D_combined.txt). The last column `mod` represents the modification state (0 is for unmodified, 1 is for modified).
```
read_idx	contig	pos	kmer	kmer_idx	mean	start_idx	end_idx	event_len	mod
1	ENST00000273480.3	18	GGCAC	657	108.646	49115	49181	45	0
1	ENST00000273480.3	19	GCACC	581	71.498	49084	49110	14	0
1	ENST00000273480.3	20	CACCA	276	79.224	49055	49080	21	0
1	ENST00000273480.3	21	ACCAC	81	76.863	49042	49051	9	0
```
The density of all `GGACT` is as follows:
<div align="center">
<img src=media/SegPore_GGACT.jpg width=40% />
</div>

### Step 4: GMM to update 5mer parameter table
Fix the mean of the first component of GMM. For GGACT, the fixed mean is 123.83.
```
sh 4_gmm.sh
```
You will get the output:
```
 * New GGACT estimated paras : mean_1 = 123.83, sigma_1 = 3.27, w_1 = 0.63, mean_2 = 117.81, sigma_2 = 2.64, w_2 = 0.37.
 * This is only the simple demo for GGACT.
```

Here we only estimate the kmer `GGACT`.

Use the results of GMM to update the 5mer parameter table `(demo/0_reference/model_kmer_m6A_without_header.csv)` manually and iteratively run Step 3 and Step 4.

In this demo experiment, the 5mer parameter table is `demo/0_reference/model_kmer_m6A_without_header.csv`. In the first round, the fixed mean is from the kmer_model (https://github.com/nanoporetech/kmer_models) of ONT. Each round, the 5mer parameter table will be updated. And after training, the 5mer parameter table is fixed and used for testing.

<div align="center">
<img src=media/segpore_workflow.jpg width=80% />
</div>

