# SegPore: raw signal segmentation for estimating RNA modifications and structures from Nanopore direct RNA sequencing data

### Preprint: https://www.biorxiv.org/content/10.1101/2024.01.11.575207v1

<div align="center">
<img src=pics/SegPore_anim.gif width=80% />
<img src=pics/github_pic.png width=95% />
</div>

## SegPore Workflow

More details on [SegPore tutorials](docs/).

### Environment setup
```
git clone https://github.com/guangzhaocs/SegPore.git
cd SegPore
conda env create -f environment.yml -p /home/user/anaconda3/envs/env_name
conda activate segpore_env
```
Here `/home/user/anaconda3/envs/env_name` is your prefix path if you want to specify a different install path than the default for your system. 

### Download demo data

Here we use the WT demo data of [xPore](https://xpore.readthedocs.io/en/latest/index.html).
```
cd SegPore
wget https://zenodo.org/record/5162402/files/demo.tar.gz
tar -xvf demo.tar.gz
cd scripts
sh 0_data_proc.sh
```

### Step 1: Basecalling, mapping and preprocessing
```
sh 1_basecalling.sh
sh 1_nanopolish.sh
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

### Reference
1. Pratanwanich, P.N., Yao, F., Chen, Y. et al. Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. Nat Biotechnol 39, 1394–1402 (2021). https://doi.org/10.1038/s41587-021-00949-w
2. Zhong, ZD., Xie, YY., Chen, HX. et al. Systematic comparison of tools used for m6A mapping from nanopore direct RNA sequencing. Nat Commun 14, 1906 (2023). https://doi.org/10.1038/s41467-023-37596-5
