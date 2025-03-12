## Divergence time estimation with BEAST2

We used genome-wide SNPs generated from the alignment of the ddRAD-seq reads in step 1 to estimate the divergence times of the four flying fox species we sampled. SNP-based estimates relied on the BEAST2 package [SNAPP](https://www.beast2.org/snapp/). We used unlinked SNPs from the previously generated file `../1_SNP_Calling/data/pteropus.vcf.gz`, selecting three representative non-admixed individuals per species. For convenience, the input filtered SNP file is also provided here: `data/pteropus_beast2.vcf.gz`. Our divergence time estimation approach is based on work by Michael Matschiner, specifically the [snapp_prep](https://github.com/mmatschiner/snapp_prep/) script and the [snapp tutorial](https://github.com/ForBioPhylogenomics/tutorials/blob/main/divergence_time_estimation_with_snp_data/README.md).

After installing BEAST2, we installed the SNAPP using the BEAST2 packagemanager.

```
packagemanager -add SNAPP
```

We next downloaded the ruby script `snapp_prep.rb` to generate the appropriate xml input file for BEAST2 and SNAPP.

```
wget https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/snapp_prep.rb
```

Besides the VCF file, we need to provide at least one time constraint for the analysis. Here we constrained the crown group containing all species but *P. scapulatus* based on the time tree in [Almeida et al. (2014)](https://www.sciencedirect.com/science/article/pii/S1055790314001092). 

```
cat data/pteropus_constraints.txt
lognormal(0,3.54,0.3)   crown   BFF,GHFF,SFF,IBFF
monophyletic    NA      BFF,GHFF,SFF,IBFF
```

*Warning*: The lognormal constraint mean (in this case 3.54) is in real space, not in log space, by default due to the flag `meanInRealSpace="true"` in the xml file produced by snapp_prep.rb. The offset is not in real space.

We also provided a table of species assignments per sample `data/pteropus_samples_snapp.txt` and a starting tree `data/pteropus_starting_tree.nwk` based on the prior RAxML analysis and known supported topologies from the literature. We could then execute `snapp_prep.rb` to produce the xml input file for BEAST2, selecting 50m MCMC generations to ensure convergence. 

```
ruby snapp_prep.rb --vcf data/pteropus_beast2.vcf --table data/pteropus_samples_snapp.txt --xml pteropus_beast2_snapp --constraints data/pteropus_constraints.txt --length 50000000 --starting-tree data/pteropus_starting_tree.nwk
```

We then ran BEAST2 using 32 threads.

```
beast -threads 32 pteropus_beast2_snapp.xml
```

Due to the high number of MCMC generations, this step took >2 weeks and was then interrupted at ~40m generations based on convergence criteria having been reached. The output logs `results/pteropus_beast2_snapp.log.gz` and output tree files `results/pteropus_beast2_snapp.trees.gz` were assessed using Tracer and BEAST-associated GUIs. TreeAnnotator was then used to produce a maximum clade credibility tree with mean node heights `results/pteropus_beast2_snapp.mean_height.tre`.    

