# Replication Origin Architecture analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23353039.svg)](https://doi.org/10.5281/zenodo.23353039)

This repository contains the R code used for the quantitative analyses in Figures 2 and 3 of the manuscript *Replication Origin Architecture: A Multidimensional Framework for Genome Duplication*.

The script performs three linked analyses:

1. It compares inter-origin spacing in nine circular archaeal chromosomes with uniformly randomized origin arrangements.
2. It examines the association between origin density and spacing variability across the 18 datasets summarized in manuscript Table 1C.
3. It groups 421 *Saccharomyces cerevisiae* origins by mean untreated MCM ChIP-seq signal and compares the spatial organization of the upper and lower licensing cohorts with randomized cohort assignments.

## Repository structure

```text
analysis/
  data_analysis.R          Complete analysis and figure-generation script
data/
  README.md                Source, citation and placement of the yeast dataset
results/
  .gitkeep                 Generated files are written here and are not tracked
LICENSE                    MIT License for the analysis code
CITATION.cff               Citation metadata for GitHub and reference managers
.zenodo.json               Metadata used when archiving releases with Zenodo
```

## Software

The manuscript analysis used R 4.2.0 with:

- readxl 1.5.0
- dplyr 1.2.1
- ggplot2 4.0.3
- patchwork 1.3.2

The script checks that these packages are installed and records the complete R environment in `results/session_info.txt` after each run.

## Input data

The yeast spreadsheet is not redistributed because it is the unchanged supplementary dataset published by Dukaj and Rhind (2021), cited as reference 14 in the manuscript. Download the source spreadsheet and save it as:

```text
data/OriDep.xlsx
```

Full citation and file-identification details are provided in `data/README.md`.

The archaeal coordinates and the 18-value cross-system dataset are the values reported in manuscript Table 1 and are written directly in the analysis script so that their use is transparent.

## Run the analysis

From the repository root:

```bash
Rscript analysis/data_analysis.R
```

The default analysis uses 10,000 randomizations and seed 123, as reported in the manuscript. For a short code check, the number of randomizations can be reduced without editing the script:

```bash
NSIM=100 Rscript analysis/data_analysis.R
```

An input file or output directory outside the default repository structure can also be supplied through environment variables:

```bash
DATA_FILE=/path/to/OriDep.xlsx OUTPUT_DIR=/path/to/results Rscript analysis/data_analysis.R
```

## Generated outputs

The script writes:

- `archaeal_empirical_results.csv`
- `yeast_licensing_empirical_results.csv`
- `DukajRhind2021_ranked_origins.csv`
- `origin_density_statistics.txt`
- `figure_2.pdf`
- `figure_3.pdf`
- `session_info.txt`

Only PDF figures are generated. The analysis does not require `svglite` and does not create SVG files.

## Statistical design

For each archaeal chromosome, random arrangements preserve chromosome length and origin number. Adjacent distances include the interval across the circular coordinate boundary. One-sided empirical tests assess whether observed spacing CV and the largest gap are lower than expected under random placement.

For budding yeast, the two untreated 0 micromolar auxin MCM ChIP-seq replicates are averaged. Origins are ranked by this mean signal and divided into four quartiles for Figure 3B. The strongest and strong quartiles form the upper cohort, while the weak and weakest quartiles form the lower cohort. Figures 3C and 3D analyse only these two cohorts. Randomization reassigns cohort membership among the 421 fixed origin positions while preserving cohort sizes.

Empirical P values are calculated as `(b + 1) / (N + 1)`, and Benjamini-Hochberg correction is applied to the test families described in the manuscript methods.

## License

The analysis code is released under the MIT License. See `LICENSE` for details.

## Citation

Citation metadata for the software and its authors is provided in `CITATION.cff`.

To cite the archived v1.0.0 release:

> Al Mamun, M., Calzada, A., & Bermejo, R. (2026). *Replication Origin
> Architecture analysis* (v1.0.0). Zenodo.
> https://doi.org/10.5281/zenodo.23353040

The version-independent DOI for all releases is
[`10.5281/zenodo.23353039`](https://doi.org/10.5281/zenodo.23353039).
