## Inferring migration between populations using fastsimcoal2

We used the coalescent simulator [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2) (v2.8) to compare four demographic models for pairs of *P. alecto* and *P. conspicillatus* genetic populations: strict isolation, constant migration, ancient migration and recent migration (manuscript Figure 2). Migration could be asymmetric. These analyses support Table 1, Tables S14–S16, Figure S16B, Figure S17 and Figure S18.

### Directory layout

```
5_Fastsimcoal/
├── pairs.tsv        # population pairs and haploid sample sizes (n0, n1) used in the .tpl files
├── templates/       # model templates m1–m4 (.tpl with @N0@, @N1@, @MU@ placeholders)
├── est/             # parameter search ranges (m3 and m4 share m3_m4_migration.est)
├── sfs/             # observed joint SFS per pair (<pair>_jointMAFpop1_0.obs, from easySFS)
├── scripts/
│   ├── get_best_from_preview.sh  # choose easySFS projection
│   ├── make_fsc_inputs.sh        # build .tpl/.est/.obs for all pairs x models for a mutation rate
│   ├── aic.R                     # AIC model selection from results tables
│   ├── calc_CI.R                 # bootstrap confidence intervals
│   └── obs_vs_exp.R              # observed vs simulated 2D SFS (Figures S17, S18)
├── data/<rate>_rate/best_replicates/   # fastsimcoal2 output of the best replicate per pair
├── results/<rate>_rate/                # likelihoods, AIC tables and bootstrap results
└── figures/                            # output of obs_vs_exp.R
```

`<rate>` is `human` (2.5×10⁻⁸ per site per generation, the main results) or `mammal` (2.2×10⁻⁹ per site per year × 4-year generation time = 8.8×10⁻⁹ per site per generation, Table S16).

### Population codes

| File code | Manuscript code | Description |
|---|---|---|
| `BFFINDOplusWAmuseum` | INDO+3NWA | Indonesian *P. alecto* plus the three North West Australian museum samples (M53754, M54659, M54660) that cluster with them |
| `BFFINDO` | INDO | Indonesian *P. alecto* only |
| `BFFNAUS` | NWA+BURK | remaining North West Australian and all Burketown *P. alecto* |
| `BFFEAST` | EA+NEA | East and North East Australian *P. alecto* |
| `BFFAUS` | NWA+EA+NEA | all Australian *P. alecto* except the three samples above |
| `SFF` | WTA+NG | *P. conspicillatus* (Wet Tropics Australia and New Guinea) |

*P. alecto alecto* is not included. Pairs are named `<pop0>_<pop1>`. In the `.obs` files, columns (`d0_*`) are pop0 and rows (`d1_*`) are pop1; in the `.tpl`, the first sample size is pop0.

### 1. Preparing the site frequency spectrum

The SFS was built from SNPs filtered for depth, genotype quality, missingness and biallelic sites, with **no** MAF or LD filter (see `1_SNP_Calling`). Missing data were handled by down-projection with [easySFS](https://github.com/isaacovercast/easySFS). The max-missing filter was reapplied to each pair-specific VCF so that no SNP was entirely missing in either population. `popmap.txt` assigns samples to the two populations of the pair.

```
easySFS.py -i pop1_pop2.vcf -p popmap.txt --preview -a > pop1_pop2.preview
scripts/get_best_from_preview.sh pop1_pop2.preview   # prints the projection maximizing segregating sites per population
```

The total number of sites (variant + invariant) was counted from the pair-specific VCF that includes invariant sites:

```
zcat pop1_pop2_snps_withinvariant.vcf.gz | grep -v '^#' | wc -l
easySFS.py -i pop1_pop2.vcf -p popmap.txt --proj <n_pop1>,<n_pop2> -a --total-length <n_sites> -o pop1_pop2_out
```

The resulting `*_jointMAFpop1_0.obs` files are in `sfs/`.

### 2. Generating the fastsimcoal2 input files

fastsimcoal2 needs a `.tpl`, a `.est` and an `.obs` file with the same prefix for each pair and model. All 48 input sets per mutation rate differ only in sample sizes and mutation rate, so they are generated from `templates/`, `est/`, `sfs/` and `pairs.tsv`:

```
cd 5_Fastsimcoal
scripts/make_fsc_inputs.sh 2.5e-8 runs/human_rate
scripts/make_fsc_inputs.sh 8.8e-9 runs/mammal_rate
```

The generated files are identical to the input files used for the published runs.

### 3. Running fastsimcoal2 and selecting the best model

Each pair × model was run as 100 independent replicates (note that replicate numbering can go over 100 because some replicates were rerun after failing due to excess runtime or memory usage). Each replicate was run in its own directory, with the `.tpl`, `.est` and `.obs` files renamed with a `_rep<N>` suffix:

```
fsc28 -t <pair>_<model>_rep<N>.tpl -e <pair>_<model>_rep<N>.est -m -n 50000 -c 1 -B 1 -L 30 -s 0 -M -C 10
```

For each pair × model, the replicate with the highest maximum composite likelihood (`MaxEstLhood` in `.bestlhoods`) was retained. These replicates are tabulated in `results/<rate>_rate/fastsimcoal2_bestlhoods_per_scenario_<rate>_mutation_rate.txt`. Models were then compared with AIC:

```
cd scripts && Rscript aic.R   # writes results/<rate>_rate/model_selection_AIC_<rate>_mutation_rate.tsv
```

The fastsimcoal2 output for the best replicate of the best model per pair is in `data/<rate>_rate/best_replicates/`. For the human rate this includes `<rep>.bestlhoods`, `.pv`, `_maxL.par` and the expected SFS `_jointMAFpop1_0.txt`. For the mammal rate, only the initial `<rep>.par` is included.

### 4. Parametric bootstrap confidence intervals

For each best replicate (`BASENAME`, e.g. `BFFINDO_SFF_m4_recent_migration_rep102`), 100 SFS were simulated from the maximum-likelihood parameters, and parameters were re-estimated for each one:

```
BASENAME=BFFINDO_SFF_m4_recent_migration_rep102
REP=data/human_rate/best_replicates/${BASENAME}
cp ${REP}/${BASENAME}.pv .
sed 's/^1 0$/10000 0/; s/^FREQ 1/DNA 100/' ${REP}/${BASENAME}_maxL.par > ${BASENAME}.par
fsc28 -i ${BASENAME}.par -n 100 -j -m -s0 -x -I -q     # simulates ${BASENAME}/${BASENAME}_{1..100}

# .tpl/.est for the best model, e.g. from scripts/make_fsc_inputs.sh, renamed to ${BASENAME}
for m in $(seq 1 100); do
    cp ${BASENAME}.est ${BASENAME}.tpl ${BASENAME}.pv ${BASENAME}/${BASENAME}_$m
    cd ${BASENAME}/${BASENAME}_$m
    fsc28 -t ${BASENAME}.tpl -e ${BASENAME}.est --initValues ${BASENAME}.pv -m -n 50000 -c 1 -B 1 -L 30 -s 0 -M -C 10
    cd ../..
done

find . -name "${BASENAME}*bestlhoods" | head -1 | xargs head -1 > ${BASENAME}.bootstraps.txt
find . -name "${BASENAME}*bestlhoods" | while read f; do tail -1 $f; done >> ${BASENAME}.bootstraps.txt
Rscript scripts/calc_CI.R ${BASENAME}.bootstraps.txt   # writes ${BASENAME}.confidence_intervals.txt
```

The bootstrap estimates and confidence intervals for both mutation rates are in `results/<rate>_rate/bootstrap/`. The manuscript reports the quantile-based intervals (`lower_quantile_CI`, `upper_quantile_CI`). Times are in generations; multiply by 4 for years.

### 5. Observed vs simulated SFS (Figures S17, S18)

```
cd scripts && Rscript obs_vs_exp.R
```

This writes `figures/2D_SFS_basic_plot_small.pdf` (Figure S17: the six Table 1 pairs) and `figures/2D_SFS_plot_big_poor_fit.pdf` (Figure S18: BFFAUS vs BFFINDOplusWAmuseum and vs BFFINDO). These figures are untidy and were edited with Inkscape for the manuscript.

### Notes

* An earlier version of this analysis used 15 pairs of the six geographic populations (PNG, INDO, NAUS, NQ, ECOAST, WETTROPICS). It has been replaced by the genetic-population pairs above, and its files remain in the git history.
* The parameter counts in the `params` column are the values used for the published AIC comparison.
