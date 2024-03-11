## SegPore: Raw Signal Segmentation for Estimating RNA Modifications and Structures from Nanopore Direct RNA Sequencing Data

<div align="center">
<img src=docs/media/SegPore_anim.gif width=80% />
<img src=docs/media/github_pic.png width=95% />
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
sh 2_hhmm_post_proc.sh
```

### Step 3: Alignment of signal segments with reference sequence
```
sh 3_alignment.sh
```

### Step 4: GMM to update 5mer parameter table
Fix the mean of the first component of GMM.
```
sh 4_gmm.sh
```
Use the results of GMM to update the 5mer parameter table and iteratively run Step 3 and Step 4.


## Citation
Cheng, Guangzhao, Aki Vehtari, and Lu Cheng. "Raw signal segmentation for estimating RNA modifications and structures from Nanopore direct RNA sequencing data." bioRxiv (2024): 2024-01. https://doi.org/10.1101/2024.01.11.575207

## Reference
1. Pratanwanich, P.N., Yao, F., Chen, Y. et al. Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. Nat Biotechnol 39, 1394–1402 (2021). https://doi.org/10.1038/s41587-021-00949-w
2. Zhong, ZD., Xie, YY., Chen, HX. et al. Systematic comparison of tools used for m6A mapping from nanopore direct RNA sequencing. Nat Commun 14, 1906 (2023). https://doi.org/10.1038/s41467-023-37596-5
