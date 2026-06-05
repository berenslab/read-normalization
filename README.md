# Compound Pearson residuals for single-cell RNA-seq data without UMIs

This repository holds the code needed to reproduce the analyses and figures presented in our preprint [Lause et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.08.02.551637v2), including additional analysis that was requested during peer review for a journal submission. The code for the earlier version of the preprint can be found under earlier releases ([v1](https://github.com/berenslab/read-normalization/releases/tag/v1.0), [v2](https://github.com/berenslab/read-normalization/releases/tag/v2.0)). Code release v1 corresponds to v1 of the preprint, code release v2 to the current version of the preprint. Release v3.0 contains the most recent revision analyses.

# Code

Some of the notebooks and R scripts depend on each other and are best run in the order indicated. For the R scripts, use our separate R environment (see below for setup instructions). Notebooks/Scripts 1-18 use the Tasic 2018 dataset. Notebooks 19-24 are based on the reads-per-UMI tables from the Ziegenhain/Hagemann-Jensen datasets.


# Datasets

- Download the reads-per-UMI tables from [zenodo](https://zenodo.org/record/8172702) and save them to `.data/reads_per_umi_tables/`. R code to obtain the same tables from the public raw data is available in `data/reads_per_umi_tables/prepare_data.R`.
- Download the Tasic raw count data from [brain-map.org](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq) via the `Gene-level (exonic and intronic) read count values for all samples (zip)` link. From these `*.zip` files, extract the `mouse_ALM_2018-06-14_exon-matrix.csv` and `mouse_VISp_2018-06-14_exon-matrix.csv` to `.data/tasic/`.
- All required metadata tables are contained in this repository for convenience.

# Compute environment

We ran all notebooks in Python 3.8.10 on an Ubuntu machine with 40 CPUs and 440 GB RAM. The following package versions were used:

- `scanpy 1.9.0`
- `anndata 0.8.0`
- `sklearn 1.0.2`
- `numpy 1.21.5`
- `matplotlib 3.5.1`
- `openTSNE 0.6.0`
- `pandas 1.4.1`
- `seaborn 0.11.2`
- `mygene 3.2.2.`
- `scipy 1.8.0`

Census and qUMI where run in a separate R conda environment specified in `r41_env.yml`. To install it, create the environment from that file with

```
conda env create -f r41_env.yml
```

Then, to install qUMI, activate the environment with `conda activate r41_env_full`, start R and run 

```
remotes::install_github("willtownes/quminorm")
```
