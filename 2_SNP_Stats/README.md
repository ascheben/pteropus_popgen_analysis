## Genetic diversity, differentiation and relatedness

We calculated genome-wide heterozygosity, nucleotide diversity (π), pairwise F<sub>ST</sub> and kinship for the four flying-fox species and for the regional populations of *P. alecto* and *P. conspicillatus*. These support Figures S3–S6 and Tables S3–S10.

Group labels: `BFF` *P. alecto*, `SFF` *P. conspicillatus*, `GHFF` *P. poliocephalus*, `LRFF` *P. scapulatus*, `alecto_alecto`/`IBFF`/`INDO_outgroup` *P. alecto alecto*. The regions are `INDO` Indonesia, `N_AUS` North West Australia, `NQ` North East Australia, `E_COAST` East Australia, `Wet_tropics` Wet Tropics Australia and `PNG` New Guinea. Popmaps are in `0_Metadata/data/`.

### Heterozygosity

Per-sample heterozygosity was calculated from the VCF that includes invariant sites. It was filtered only for depth (≥ 5), genotype quality (≥ 20) and missingness (≤ 20%); see `1_SNP_Calling`. Heterozygosity is the proportion of callable (non-missing) sites that are heterozygous, `Het_frac = Het / (Total - Missing)`, from per-sample genotype counts produced with BCFtools. The counts in `data/pteropus_heterozygosity.txt` (`Hom_Ref`, `Hom_Alt`, `Het`, `Missing`) are equivalent to the per-sample `PSC` statistics:

```
bcftools stats -s - <species>.withinvariant.vcf.gz | grep '^PSC' | \
  awk 'BEGIN{OFS="\t"} {print $3, $4, $5, $6, $14, $6/($4+$5+$6)}'   # sample, Hom_Ref, Hom_Alt, Het, Missing, Het_frac
```

### Nucleotide diversity

π was computed with [pixy](https://pixy.readthedocs.io/) 1.2.11.beta1 from the same variant + invariant VCF, in non-overlapping 100 kb windows. Windows with fewer than 50 sites are removed in `plot_stats.R`.

```
pixy --populations pteropus.popmap --vcf pteropus.withinvariant.vcf.gz --stats pi --n_cores 16 --output_prefix pteropus_pi --window_size 100000
```

### F<sub>ST</sub>

F<sub>ST</sub> was calculated with stacks `populations` from the filtered unlinked SNPs (`1_SNP_Calling/data/`). It was run for the species ("all" set), for the regional populations ("alecto + conspicillatus" set), and for two genetic-population groupings used for TreeMix and fastsimcoal2 (Tables S9, S10). Significance was assessed by permuting population labels and recalculating F<sub>ST</sub> 1,000 times.

```
populations -V pteropus.vcf.gz -M pteropus.popmap -t 12 --fstats -O pteropus
populations -V palecto_pconspicillatus_broad.vcf.gz -M palecto_pconspicillatus.popmap -t 12 --fstats -O palecto_pconspicillatus
```

### Kinship

KING kinship coefficients were calculated with PLINK2. Pairs with kinship > 0.2 (duplicates and first-degree relatives) were used to exclude samples. The *P. alecto*–*P. conspicillatus* table is also used for Figure 4D in `3_SNP_Structure`.

```
plink2 --vcf pteropus.vcf.gz --allow-extra-chr --make-king-table --out pteropus
```

### Data

| File | Content |
|---|---|
| `data/pteropus_heterozygosity.txt` | per-sample genotype counts and heterozygosity (`Pop`, `Sample`, `Het_frac`, `Hom_Ref`, `Hom_Alt`, `Het`, `Missing`, `Total`) |
| `data/pteropus_nucleotide_diversity.txt.gz` | pixy π per 100 kb window, species and regions |
| `data/pteropus.fst_summary.tsv` | F<sub>ST</sub> between species (Figure S5, Table S7) |
| `data/palecto_pconspicillatus.fst_summary.tsv` | F<sub>ST</sub> between regions (Figure S6, Table S8) |
| `data/palecto_pconspicillatus_genetic_scenario_a.fst_summary.tsv` | F<sub>ST</sub> between genetic populations INDO (INDO+3NWA), BFFAUS (NWA+EA+NEA), SFF and IBFF |
| `data/palecto_pconspicillatus_genetic_scenario_b.fst_summary.tsv` | F<sub>ST</sub> between genetic populations N_AUS (NWA+BURK), BFFEAST (EA+NEA), INDO (INDO+3NWA), SFF and IBFF |
| `data/pteropus.kinship.tsv`, `data/palecto_pconspicillatus.kinship.tsv` | PLINK2 KING kinship tables |

### Reproducing Figures S3–S6

```
cd scripts && Rscript plot_stats.R
```

Outputs in `figures/`:
- `Heterozygosity_species_populations.pdf` (S3)
- `Nucleotide_diversity_all_species.pdf` (S4)
- `Fst_species.pdf` (S5)
- `Fst_bff_sff_regional.pdf` (S6)
- `*_summary.tsv`: means and SDs per group
- `*_wilcoxon.tsv`: two-sided pairwise Wilcoxon rank-sum tests, Bonferroni-adjusted

