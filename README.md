# Compound models and Pearson residuals for normalization of single-cell RNA-seq data without UMIs

This repository holds the code needed to reproduce the analyses and figures presented in Lause et al. (2023).

# Code

Some of the notebooks depend on each other.

All plots based on the Tasic datasets require `01_prepare_tasic` to run first. Then,

- to reproduce Figure 2, S1, S2, S3, run `02`
- to reproduce Figure 3, S4, S5, S6, run `03` to compute t-SNEs and `04` to make the figures
- to reproduce Figure 6, run `09` to prepare simulated data and `10` to make the figure

All plots based on the reads-per-UMI tables requires `05_prepare_umi_datasets` to run first. Then,

- to reproduce Figures 4 and 5, run `06`
- to reproduce Figure S7, run `07`
- to reproduce Figure S8, run `08`

# Datasets

- Download the reads-per-UMI tables from [zenodo](https://zenodo.org/record/8172702) and save them to `.data/reads_per_umi_tables/`. R code to obtain the same tables from the public raw data is available in `data/reads_per_umi_tables/prepare_data.R`.
- Download the Tasic raw count data from [brain-map.org](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq) via the `Gene-level (exonic and intronic) read count values for all samples (zip)` link. From these `*.zip` files, extract the `mouse_ALM_2018-06-14_exon-matrix.csv` and `mouse_VISp_2018-06-14_exon-matrix.csv` to `.data/tasic/`.
- All required metadata tables are contained in this repository for convenience.

# Compute environment

We ran the code in Python 3.8.10 on an Ubuntu machine with 40 CPUs and 440 GB RAM. The following package versions were used:

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


