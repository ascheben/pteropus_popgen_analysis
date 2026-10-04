# Analysis of the population genetics of four flying fox species (*Pteropus*)

We investigate the population genomics of the highly mobile flying foxes *Pteropus alecto*, *P. conspicillatus*, *P. poliocephalus* and *P. scapulatus* across large parts of their overlapping ranges in Australia, Indonesia and New Guinea. Using reduced-representation resequencing of 242 individuals and a range of population genetic analyses, we examine the extent to which panmixia, isolation-by-distance and hybridization shaped these populations. This repository supplements our manuscript with the datasets and scripts used to generate our results and figures. Raw reads are in the SRA (BioProject [PRJNA1230740](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1230740)).

## Repository contents

Each numbered directory is one part of the analysis. It has a `README.md` with the commands and parameters used, a `data/` directory of inputs and intermediate results, and, where figures are produced, `scripts/` and `figures/` directories.

| Directory | Content |
|---|---|
| `0_Metadata` | Sample metadata, popmaps and population code definitions |
| `1_SNP_Calling` | Alignment, variant calling and filtering; the three filtered SNP sets used throughout |
| `2_SNP_Stats` | Heterozygosity, nucleotide diversity, F<sub>ST</sub> and kinship |
| `3_SNP_Structure` | PCA, fastSTRUCTURE, RAxML/CASTER phylogenies, triangle plots, isolation-by-distance |
| `4_Treemix` | Population tree with migration edges (TreeMix) |
| `5_Fastsimcoal` | Demographic model comparison and parameter estimation for population pairs (fastsimcoal2) |
| `6_BEAST` | Species divergence times (BEAST2/SNAPP) |
| `7_Map_Damage` | DNA damage in historical museum samples (mapDamage) |

## Reproducing the figures

Each R script is run from its own `scripts/` directory and writes to the module's `figures/` directory:

```
cd <module>/scripts && Rscript <script>.R
```

| Manuscript | Script | Output in `figures/` |
|---|---|---|
| Figure 3A, 3B, 3C, 3D; Figure S1 | `3_SNP_Structure/scripts/pca_and_phylo_plot.R` | `Geography_*.pdf`, `PCA_*.pdf` |
| Figure 3E | `3_SNP_Structure/scripts/structure_plot.R` | `Admixture_barplot_faststructure_bff_sff_k4.pdf` |
| Figures 4C, 4D; Table S17 | `3_SNP_Structure/scripts/pca_and_phylo_plot.R` | `Distance_from_Sawu_by_PC1.pdf`, `GeneticDistance_by_GeoDistance.pdf`, `mantel_test_*.txt` |
| Figure S9 | `3_SNP_Structure/scripts/pca_and_phylo_plot.R` | `PCA_PC1_PC2_bff_sff_species_regions_nomiss.pdf` |
| Figure S2 | `7_Map_Damage/scripts/plot_mapdamage.R` | `CtoT_extant_vs_museum.pdf` |
| Figures S3–S6; Tables S3–S6 | `2_SNP_Stats/scripts/plot_stats.R` | `Heterozygosity_*`, `Nucleotide_diversity_*`, `Fst_*` |
| Figures S10, S11 | `3_SNP_Structure/scripts/triangle_hybridization_plot.R` | `Triangle_plot_BFF_SFF*.pdf` |
| Figures S12–S15 | `3_SNP_Structure/scripts/pca_and_phylo_plot.R` | `Phylogram_*.pdf` (RAxML and CASTER-site) |
| Figure S16A | `4_Treemix/scripts/treemix_plot.R` | `Treemix_m2_*.pdf` |
| Figures S17, S18 | `5_Fastsimcoal/scripts/obs_vs_exp.R` | `2D_SFS_*.pdf` |
| Table 1, Tables S14–S16 | `5_Fastsimcoal/scripts/aic.R`, `calc_CI.R` | `5_Fastsimcoal/results/` |

The scripts produce the data panels. For the manuscript, labels, legends, colours of tip labels and panel layout were finalised using Inkscape.

The following are not reproduced by scripts in this repository:
- Figure 1 (photographs and IUCN ranges) and Figure 2 (model schematic)
- Figure 4A (admixture map) and Figure 4B (FEEMSmix)
- Figures S7B and S8 (alternative PCAs)
- Figure S16B (migration map)
- Figure S19 (BEAST tree; tree file in `6_BEAST/results/`)

See the module READMEs for details.

## Software

Command-line tools (versions as in the manuscript):
- stacks 2.52
- BWA-MEM 0.7.17
- GATK 4.6.2.0
- BCFtools 1.19
- VCFtools 0.1.16
- PLINK2 2.00a2.3LM
- mapDamage 2.2.2
- pixy 1.2.11.beta1
- RAxML 8.2.12
- CASTER-site 1.23.2.6
- fastSTRUCTURE 1.0
- TreeMix 1.13
- easySFS
- fastsimcoal2 2.8
- BEAST2 with SNAPP

The R scripts were checked with R 4.5.0 and these package versions: adegenet 2.1.11, ape 5.8.1, colorspace 2.1.2, cowplot 1.2.0, dplyr 1.2.0, gdsfmt 1.46.0, geosphere 1.6.8, ggplot2 4.0.2, ggpubr 1.0.0, ggrepel 0.9.7, ggtree 4.0.5, ggpointdensity 0.2.1, gridExtra 2.3, pophelper 2.3.1 (GitHub royfrancis/pophelper), poppr 2.9.8, R.utils 2.13.0, RColorBrewer 1.1.3, reshape2 1.4.5, rnaturalearth 1.2.0, rnaturalearthdata 1.0.0, rnaturalearthhires 1.0.0.9000, sf 1.1.3, SNPRelate 1.40.0 (Bioconductor), tibble 3.3.1, tidyr 1.3.2, triangulaR 0.0.1 (GitHub omys-omics/triangulaR), vcfR 1.16.0, vegan 2.7.3, viridis 0.6.5. So that newer SNPRelate version reproduce the same plots, `pca_and_phylo_plot.R` sets `missing.rate = NaN` explicitly.

`4_Treemix/scripts/plotting_funcs.R` contains the TreeMix plotting functions by Joe Pickrell, taken from [joepickrell/pophistory-tutorial](https://github.com/joepickrell/pophistory-tutorial/blob/master/example2/plotting_funcs.R).

## Citation

>Scheben, A., McKeown, A., Walsh, T., Westcott, D.A., Metcalfe, S.S., Vanderduys, E.P., and Webber, B.L. Extreme mobility creates contrasting patterns of panmixia, isolation-by-distance and hybridization in four flying-fox species (*Pteropus*). 2026. *bioRxiv*. doi: 10.1101/2025.04.06.647431.
