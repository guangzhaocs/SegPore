## SegPore: Raw Signal Segmentation for Estimating RNA Modification from Nanopore Direct RNA Sequencing Data

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
conda env create -f environment.yml
conda activate segpore_env
pip3 install git+https://github.com/EGA-archive/ont2cram
# install Guppy
```

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
2.1 Firstly, prepare the input of HHMM.
```
sh 2_hhmm_prepare.sh
```
2.2 Next, run HHMM on CUDA:
```
sh 2_hhmm_GPU.sh
```
2.3 Finally, generate the final output:
```
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

## Cite SegPore

Guangzhao Cheng, Aki Vehtari, Lu Cheng (2025) **Raw signal segmentation for estimating RNA modification from Nanopore direct RNA sequencing data**, *eLife*, **14**:RP104618
https://doi.org/10.7554/eLife.104618.1

```
@article{Cheng_2025,
  title={Raw signal segmentation for estimating RNA modification from Nanopore direct RNA sequencing data},
  url={http://dx.doi.org/10.7554/eLife.104618.1}, DOI={10.7554/elife.104618.1},
  publisher={eLife Sciences Publications, Ltd},
  author={Cheng, Guangzhao and Vehtari, Aki and Cheng, Lu},
  year={2025},
  month=mar
}
```

## Reference
1. Pratanwanich, P.N., Yao, F., Chen, Y. et al. Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. Nat Biotechnol 39, 1394–1402 (2021). https://doi.org/10.1038/s41587-021-00949-w
2. Zhong, ZD., Xie, YY., Chen, HX. et al. Systematic comparison of tools used for m6A mapping from nanopore direct RNA sequencing. Nat Commun 14, 1906 (2023). https://doi.org/10.1038/s41467-023-37596-5
