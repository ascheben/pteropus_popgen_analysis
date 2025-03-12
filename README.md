# Analysis of the population genetics of four flying fox species (*Pteropus*)

We investigate the population genomics of the highly mobile flying foxes *Pteropus alecto*, *P. conspicillatus*, *P. poliocephalus* and *P. scapulatus* across large parts of their overlapping ranges in Australia, Indonesia and New Guinea. Using reduced-representation resequencing and an array of population genetics analyses, we examine the extent to which panmixia, isolation-by-distance and hybridization shaped these populations. This repository supplements our manuscript by providing additional datasets and scripts used to generate our results and figures.

## Repository contents

This repository is organized into subdirectories corresponding to independent parts of the population genetics analysis. Each subdirectory contains a `README.md` file with a description of the work and can also include `data` and `scripts` directories. 

0. `Metadata`: Sample metadata including sampling location and phenotype 
1. `SNP_Calling`: Variant calling pipeline and results
2. `SNP_Stats`: Population genetic summary statistics calculated from the variant set 
3. `SNP_Structure`: PCA, phylogenetic analysis, and admixture inference
4. `Treemix`: Analysis of migration using Treemix
5. `Fastsimcoal`: Inference of demographic parameters using population pairs and coalescent modelling with fastsimcoal2 
6. `BEAST`: Estimation of divergence time between species using BEAST2 and SNAPP

## Citation

Pending.
