# Compound models and Pearson residuals for single-cell RNA-seq data without UMIs

This repository holds the code needed to reproduce the analyses and figures presented in [Lause et al. (2024)](https://www.biorxiv.org/content/10.1101/2023.08.02.551637v2). The code for the earlier version of the preprint from August 2023 can be found under the [v1 release](https://github.com/berenslab/read-normalization/releases/tag/v1.0).

# Code

Some of the notebooks depend on each other.

All plots based on the Tasic 2018 dataset require `01_prepare_tasic` to run first. Then,

- to reproduce Figures for the homogeneous within-cluster data (Figure 2, S1, S2, S4), run notebook `02`
- to reproduce Figures for the full Tasic data (Figure 3, S5, S6, S8), run notebook `03` to compute t-SNEs etc. and notebook `04` to make the figures
- to reproduce Figure S7 for the Census/qUMI comparison, run Census and qUMI with `05_compute_tasic_qumis_census.R` (using our separate R environment, see below for setup instructions), and then run notebooks `06`-`08` to load the data, process and plot
- to reproduce Figure S9 for the Tasic-like simulated data, use notebooks `09`-`11` to simulate, process and plot
- to reproduce Figure 6, run notebook `16` to prepare simulated data and notebook `17` to plot

All plots based on the reads-per-UMI tables from the Ziegenhain/Hagemann-Jensen datasets requires `12_prepare_ziegenhain` to run first. Then,

- to reproduce the main Figures 4 and 5, run notebook `13`
- to reproduce Figure S3 on Pseudogenes, run notebook `14`
- to reproduce Figure S10 and S11 on per-cell amplification parameter estimates, run notebook `15`


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



