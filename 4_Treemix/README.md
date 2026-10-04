## Inferring migration edges with TreeMix

We used [TreeMix](https://bitbucket.org/nygcresearch/treemix) 1.13 to infer a maximum likelihood population tree with migration edges for genetic populations of *P. alecto* and *P. conspicillatus*. The Indonesian *P. alecto alecto* (IBFF) was used as the outgroup. This supports Figure S16A.

### Populations

| TreeMix label | Manuscript code | Description |
|---|---|---|
| `INDO` | INDO+3NWA | Indonesian *P. alecto* plus the three North West Australian museum samples (M53754, M54659, M54660) clustering with them |
| `N_AUS` | NWA+BURK | remaining North West Australian and all Burketown *P. alecto* |
| `BFFEAST` | EA+NEA | East and North East Australian *P. alecto* |
| `BFFAUS` | NWA+EA+NEA | all Australian *P. alecto* except the three samples above (alternative grouping) |
| `SFF` | WTA+NG | *P. conspicillatus* |
| `IBFF` | – | *P. alecto alecto* (outgroup) |

### Running TreeMix

The filtered unlinked SNPs of the "alecto + conspicillatus" set (`1_SNP_Calling/data/palecto_pconspicillatus_broad.vcf.gz`) were converted to TreeMix format using [stacks](https://catchenlab.life.illinois.edu/stacks/) `populations`, with a popmap that assigns samples to the populations above:

```
populations -V palecto_pconspicillatus_broad.vcf.gz -M treemix.popmap --treemix -O .
gzip -c palecto_pconspicillatus_broad.p.treemix > bffsff_<grouping>.treemix.gz   # remove the stacks header line if present or treemix will throw an error
```

TreeMix was run with 0–7 migration edges (m), 10 replicates each, 500 SNPs per block, global rearrangements, a random seed and bootstrap resampling. The input file is the first argument:

```
for m in {0..7}; do
    for i in {1..10}; do
        s=$RANDOM
        treemix -bootstrap -n_warn 1 -i $1 -o ${1%%.treemix.gz}.${i}.${m} -global -m ${m} -k 500 -seed ${s} -root IBFF
    done
done
```

The optimal number of migration edges was chosen with the Evanno method in the R package [OptM](https://cran.r-project.org/package=OptM) (`optM(<dir_with_m0-7_outputs>)`). We then ran 100 further replicates at the optimal m (same command, `i` in 1..100). We used m = 2 for the primary grouping (INDO, N_AUS, BFFEAST, SFF, IBFF) and m = 1 for the alternative grouping with all Australian *P. alecto* combined (BFFAUS; mentioned but not shown in the manuscript).

### Data

* `data/m2_reps_NAUS_NQECOAST_INDOplusWAmuseum_SFF_IBFF/`: 100 replicates with m = 2, primary grouping.
* `data/m1_reps_AUS_INDOplusWAmuseum_SFF_IBFF/`: 100 replicates with m = 1, alternative grouping.

Each directory contains the log-likelihoods of all replicates (`<prefix>.<rep>.<m>.llik`). The full TreeMix outputs (`.cov.gz`, `.covse.gz`, `.edges.gz`, `.modelcov.gz`, `.treeout.gz`, `.vertices.gz`) are in `replicate_outputs.tar.gz` (extract with `tar xzf replicate_outputs.tar.gz`), except for the plotted replicate (rep 99, m = 2), whose outputs are kept unpacked. Rep 99 is one of 28 replicates tied for the highest log-likelihood (112.664); all of them have the same topology and migration edges.

### Reproducing Figure S16A

```
cd scripts && Rscript treemix_plot.R
```

This writes `figures/Treemix_m2_NAUS_NQECOAST_INDOplusWAmuseum_SFF_IBFF.pdf`. Population labels were renamed to the manuscript codes in a graphics editor. `scripts/plotting_funcs.R` contains the TreeMix plotting functions by Joe Pickrell, taken from [joepickrell/pophistory-tutorial](https://github.com/joepickrell/pophistory-tutorial/blob/master/example2/plotting_funcs.R).

The TreeMix input files and the m = 0–7 runs used for OptM are not included in this repository.
