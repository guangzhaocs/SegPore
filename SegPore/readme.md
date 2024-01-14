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

In nanopolish eventalign and summary results, one read may be processed multiple times (one read_name may correspond to multi read_idx).
```
6cd0cd47-db93-4e73-99cc-e91c68f45268 [91, 157]
a7aa6921-712f-481e-8496-4a963618b786 [134, 163]
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


### Step 2: Hierarchical hidden Markov model (HHMM) for signal segmentation
```
sh todo.sh
```

### Step 3: Alignment of signal segments with reference sequence
```
sh todo.sh
```

### Step 4: GMM to update 5mer parameter table
```
sh todo.sh
```

