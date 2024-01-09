# SegPore: raw signal segmentation for estimating RNA modifications and structures from Nanopore direct RNA sequencing data

<div align="center">
<img src=pics/SegPore_anim.gif width=80% />
<img src=pics/github_pic.png width=95% />
</div>

## SegPore Workflow

### Environment setup
```
git clone https://github.com/guangzhaocs/SegPore.git
cd SegPore
conda env create -f environment.yml
conda activate segpore_env
```

### Download demo data

Here we use the WT demo data of [xPore](https://xpore.readthedocs.io/en/latest/index.html).
```
cd SegPore
wget https://zenodo.org/record/5162402/files/demo.tar.gz
tar -xvf demo.tar.gz
```

### Step 1: Basecalling, mapping and preprocessing
```
sh basecalling.sh
```

### Step 2: Hierarchical hidden Markov model for signal segmentation
```
sh todo.sh
```

### Step 3: Alignment of signal segments with reference sequence
```
sh todo.sh
```
