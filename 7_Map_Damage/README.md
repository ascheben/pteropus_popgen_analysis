## Assessing DNA degradation in historical samples with mapDamage

To test whether the historical *P. alecto* specimens from the Western Australian Museum (collected 1989–2003) show more DNA degradation than extant samples, we compared cytosine deamination (C to T substitutions at 5' read ends) using [mapDamage](https://ginolhac.github.io/mapDamage/) 2.2.2. This supports Figure S2 and the comparison reported in Results section 3.1.

### Running mapDamage

mapDamage was run on each read-group-tagged alignment of *P. alecto* to the *P. alecto* reference genome (see `1_SNP_Calling`). The alignments had the 5 bp low-quality read flanks clipped, so position 6 in the mapDamage output is the first retained 5' base:

```
mapDamage -i <sample>_rg.bam -r GCF_000325575.1_ASM32557v1_genomic.fna --merge-reference-sequences
```

The per-position C to T frequencies (`5pCtoT_freq.txt` in each mapDamage output directory) were combined into one table, with the sample name added as a column. For example, with the default output directories `results_<sample>_rg`:

```
for d in results_*_rg; do
    s=${d#results_}; s=${s%_rg}
    tail -n +2 $d/5pCtoT_freq.txt | awk -v s=$s 'BEGIN{OFS="\t"}{print $1,$2,s}'
done | sed '1i pos\tCtoT_fraction\tsample' > 5pCtoT_freq_clipped_all.txt
```

### Data

* `data/5pCtoT_freq_clipped_all.txt`: C to T substitution fraction for 5' positions 1–25 for each sample (columns `pos`, `CtoT_fraction`, `sample`; `sample` matches `Sample Identifier` in `0_Metadata/data/pteropus_metadata.txt`).

### Reproducing Figure S2

```
cd scripts && Rscript plot_mapdamage.R
```

The script compares the C to T fraction at the first retained base between extant (`Source` = CSIRO, n = 66) and historical (`Source` = WA museum, n = 15) *P. alecto* samples with an unpaired *t*-test (p = 0.45). It writes `figures/CtoT_extant_vs_museum.pdf`. In the manuscript the groups were relabelled "Extant" and "Historical".
