## Inference of admixture and population structure

We first investigated population structure and admixture using conventional approaches including fastSTRUCTURE analysis, maximum likelihood phylogenetic inference, and principal component analysis (PCA).


### Admixture analysis with fastSTRUCTURE 

To run fastSTRUCTURE, we converted the filtered unlinked SNPs to PLINK format.

```
plink2 --set-all-var-ids @:# --allow-extra-chr --make-bed --out palecto_pconspicillatus_broad --vcf palecto_pconspicillatus_broad.vcf 
```

We then executed the analysis using default paramters and a range of values for the number of clusters *k* (1-15).
 
```
seq 1 15 | while read n; do 
structure.py -K $n --input=palecto_pconspicillatus_broad --output=palecto_pconspicillatus_broad_k{n}_out
done
```

We visualised the fastSTRUCTURE results using R and the package pophelper:

```
Rscript scripts/structure_plot.R
```

### Phylogenetic inference with RAxML

SNPs in VCF format were converted to phylip format using [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip). Phylogenetic inference was then conducted using RAxML's rapid bootstrap approach and 100 replicates.

```
raxmlHPC-PTHREADS -T 24 -f a -m GTRCAT -p 12345 -x 12345 -# 100 -s pteropus.phy -n pteropus
raxmlHPC-PTHREADS -T 24 -f a -m GTRCAT -p 12345 -x 12345 -# 100 -s palecto_pconspicillatus.phy -n palecto_pconspicillatus
```

### PCA of SNP matrix

PCA was run using the R package SNPRelate and the results of both phylogenetic inference and PCA were visualised using a custom R script:

```
Rscript scripts/plot_fst_snp_struc.R
```

### Investigating hybridization with a triangular plot 

Although fastSTRUCTURE, RAxML and PCA can all help identify hybridization, the most robust approach may be a triangle plot analysis. We therefore conducted this analysis using the R package triangulaR.

```
Rscript scripts/triangulaR.R
```
