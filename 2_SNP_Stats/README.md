## Calculating Fst, pi and heterozygosity statistics

We calculated Fst with [stacks populations]() using filtered unlinked SNPs for 1) all species and 2) *P. alecto* and *P. conspicillatus* populations.

```
populations -V pteropus.vcf.gz -M pteropus.popmap -t 12 --fstats -O pteropus
populations -V palecto_pconspicillatus_broad.vcf.gz -M palecto_pconspicillatus.popmap -t 12 --fstats -O palecto_pconspicillatus
```

Using the variant and invariant SNPs, we then use [pixy](https://pixy.readthedocs.io/en/latest/) to calculate nucleotide diversity in 100kb windows.

```
pixy --populations pteropus.popmap --vcf pteropus.withinvariant.vcf.gz --stats pi --n_cores 16 --output_prefix pteropus_pi --window_size 100000
```

Heterozygosity was calculated per sample using [vcftools](https://vcftools.github.io). 

```
vcftools --gzvcf pteropus.withinvariant.vcf.gz --het --out pteropus.withinvariant 
```


The statistics were visualised using a custom R script:

```
Rscript plot_stats.R
```

### Results files 

* `data/palecto_pconspicillatus.fst_summary.tsv`
* `data/pteropus.fst_summary.tsv`
* `data/pteropus_heterozygosity.txt`
* `data/pteropus_nucleotide_diversity.txt.gz`
