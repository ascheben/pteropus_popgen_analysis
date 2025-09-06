## Inferring migration between populations using fastsimcoal2

The coalescent simulator [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2) was used to test four demographic models (isolation,constant migration,ancient migration,recent migration) for pairs of *P. alecto* and *P. conspicillatus* populations (PNG, INDO, NAUS, NQ, ECOAST, WET).

## Preparing SNP data for fastsimcoal2

The site frequency spectrum cannot easily be inferred from a SNP matrix with missing values. Because ddRAD-seq is prone to substantial levels of missing data, we maximized the number of segregating sites by using a down-projection approach implemented in [EasySFS](https://github.com/isaacovercast/easySFS). We used a a SNP set that had no MAF filter applied. All other filters were the same as for the SNP sets used for other analyses. However, the max-missing filter and the plink LD filter was reapplied to the population pair specific *P. alecto* and *P. conspicillatus* SNP sets that were used for constructing the site frequency spectrum. This ensures that there are no entirely missing SNPs for each population pair.

The next step is to select projection values as follows:

```
easySFS.py -i pop1_pop2.vcf -p popmap.txt --preview -a > pop1_pop2.preview
```

A simple script can then parse the output files to select the downprojection that maximizes segregating sites. 
```
cat get_best_from_preview.sh
grep -B 1 "(2" $1 | grep -v '(2'| grep -v "\-\-"
grep "(2," $1| head -1| tr '\t' '\n'| sed 's/(//'| sed 's/)//'| sed 's/ //'| sed '/^$/d'| tr ',' '\t'| sort -n -k2,2|tail -1
grep "(2," $1| tail -1| tr '\t' '\n'| sed 's/(//'| sed 's/)//'| sed 's/ //'| sed '/^$/d'| tr ',' '\t'| sort -n -k2,2| tail -1
```

This script can take the preview file as input and will print the recommended projection number as well as the number of sites to stdout.

```
get_best_from_preview.sh pop1_pop2.preview
pop1
pop2
26      33514
18      34838
```

Next, we need to calculate the total number of sites assessed to detect the SNPs used for inference of the site frequency spectrum, which easySFS refers to as the "total length". Knowing the number of monomorphic ('zero bin') sites is useful when fastsimcoal2 is going to be calibrated with just a mutation rate, which is what we will do here. For this reason, monomorphic sites were called together with polymorphic sites using bcftools in our initial variant calling.

The total number of sites can thus be easily counted for ach SNP set.

```
zcat pop1_pop2_snps_withinvariant.vcf.gz| grep -v '^#'| wc -l 
3162569
```

Finally, the projection sample number and the total length parameters can be used as inputs to calculate the joint site frequency spectrum.

```
easySFS.py -i pop1_pop2.vcf -p popmap.txt --proj 26,18 -a --total-length 3162569 -o pop1_pop2_out
```

The output file with the suffix `jointMAFpop1_0.obs` is the site frequency spectrum that will be used for model fitting with fastsimcoal2.

## Running fastsimcoal2 to find the best-fitting demographic model

We require three input files for fastsimcoal2: 

1. The joint site frequency spectrum file `.obs`
2. A model template file `.tpl`
3. An estimation file `.est`

We provide these input files for each population pair in the directory `5_Fastsimcoal/data/`. Importantly, the template file contains the custom mutation rate used for calibrating the estimated parameters, in this case we use `2.5*e-8`. Now we execute fastsimcoal2 for 100 replicates per scenario. Note that for the command below the file `pop1_pop2_m4_recent_migration_<rep>_jointMAFpop1_0.obs` is expected in the same directory where the command is executed.

```
fastsimcoal2 -t pop1_pop2_m4_recent_migration_<rep_n>.tpl -e pop1_pop2_m4_recent_migration_<rep>.est -m -n 50000 -c 1 -B 1 -L 30 -s 0 -M -C 10
```

By generating 100 replicates per scenario we ensure that the model fitting has converged to a robust set of parameters. For each model we then select the replicate with the highest maximum likelihood value. Because the models we are comparing contain different numbers of parameters and more parameters generally allows for a better fit, we must account for this when comparing the fits. Here we apply the widely used Akike Information Criterion (AIC) to compare model fits.

We calculate AIC in R as follows for each row in a dataframe `fsc2_df` of fastsimcoal2 results.
```
fsc2_df$AIC <- 2*fsc2_df$params-2*(fsc2_df$MaxEstLhood/log10(exp(1)))
```

## Estimating demographic parameter confidence intervals using parametric bootstrapping

The point estimates for migration rates, effective population size, divergence time and other demographic parameters are more useful when combined with confidence intervals. For fastsimcoal2, it is recommended to calculate these using parametric bootstrapping. Here, we run 100 parametric bootstraps per sample and then calculate the 95% confidence interval.

To simplify this step, we define a basename to run the bootstraps on, in this case the replicate of the best demographic model with the highest max likelihood.  

```
BASENAME="<pop1_pop2>"
```

Next, we prepare 100 replicate runs for 10,000 loci of length 100 each using the commands below.

```
cp ../reps/${BASENAME}/${BASENAME}.pv ${BASENAME}.pv
cat ../reps/${BASENAME}/${BASENAME}_maxL.par| sed 's/^1 0$/10000 0/'| sed 's/^FREQ 1/DNA 100/' > ${BASENAME}.par
cp ../reps/${BASENAME}.est ${BASENAME}.est
cp ../reps/${BASENAME}.tpl ${BASENAME}.tpl
cp ../reps/${BASENAME}_jointMAFpop1_0.obs ${BASENAME}_jointMAFpop1_0.obs

fastsimcoal2 -i ${BASENAME}.par -n 100 -j -m -s0 -x -I -q

```

Now we have to copy over our estimate, template and parameter values file to each replicate directory and rerun fastsimcoal2 passing the parameter values to `--initValues`. 

```
seq 1 100| while read m; do 
    cp ${BASENAME}.est ${BASENAME}/${BASENAME}_$m 
    cp ${BASENAME}.tpl ${BASENAME}/${BASENAME}_$m
    cp ${BASENAME}.pv ${BASENAME}/${BASENAME}_$m
    cd ${BASENAME}/${BASENAME}_$m
    fastsimcoal2 -t ${BASENAME}.tpl -e ${BASENAME}.est --initValues ${BASENAME}.pv -m -n 50000 -c 1 -B 1 -L 30 -s 0 -M -C 10
    cd ../..
;done
```
Once this step is complete, we concatenate all of the best parameter estimates in the `.bestlhoods` output files and then calculate the 95% confidence intervals using a custom R script `scripts/calc_CI.R`.

```
# Add the header
find . -name "${BASENAME}*bestlhoods"| head -1| while read m; do head -1 $m  > ${BASENAME}.bootstraps.txt;done
# Add the point estimates for parameters
find . -name "${BASENAME}*bestlhoods"| while read m; do tail -1 $m  >> ${BASENAME}.bootstraps.txt;done
# Calculate confidence intervals
Rscript calc_CI.R ${BASENAME}.bootstraps.txt
```

The R script generates an output file with the suffix `.confidence_intervals.txt`.

