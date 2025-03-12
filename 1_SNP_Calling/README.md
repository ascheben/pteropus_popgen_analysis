## Alignment and variant calling

Reads were aligned to the Pteropus alecto reference genome `GCF_000325575.1_ASM32557v1_genomic.fna` from GenBank using bwa mem 0.7.17-r1198-dirty with default settings. Initial exploratory analyses using the stacks pipeline with and without a reference showed that a bwa-bcftools approach produced the most high-quality SNPs after filtering. 

After adapter-trimming and quality control, paired reads were aligned to the reference.

```
bwa mem -t 1 GCF_000325575.1_ASM32557v1_genomic.fna <sample>.1.fq.gz  <sample>.2.fq.gz | samtools sort --threads 1 > <sample>.bam
```

SNPs were then called with bcftools and invariant sites were included to inform demographic modelling.
```
bcftools mpileup --threads 48 -q 20 -Q 10 -f GCF_000325575.1_ASM32557v1_genomic.fna --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -b rg_bam.list | bcftools call --threads 48 -mO z -f GQ -o snps_withinvariant.vcf.gz
```

SNPs were filtered using VCFtools. First multiallelic and monomorphic SNPs were removed and all SNPs within the pre-calculated 5bp flanks of aligned reads were also removed. 

```
vcftools --gzvcf snps_withinvariant.vcf.gz --exclude-bed cluster_start_end_blacklist_regions.bed --remove-indels --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --stdout |gzip -c > snps_no_5bp_flanks.biallelic.vcf.gz
```

Next we filtered alleles with a genotype depth <5 and then we removed SNPs with >50% missing alleles or a minor allele frequency below 5%.

```
myvcf='snps_no_5bp_flanks.biallelic.vcf.gz'
# Remove samples with high missingness and apply standard filters
#vcftools --gzvcf $myvcf --minDP 5 --recode --recode-INFO-all --stdout|gzip -c > ${myvcf%%.vcf.gz}_minDP5.vcf.gz
# Filter rare alleles
vcftools --gzvcf ${myvcf%%.vcf.gz}_mim09_biallelic_minDP5.vcf.gz  --maf 0.05 --max-missing 0.5 --recode --recode-INFO-all --stdout | gzip -c > ${myvcf%%.vcf.gz}_minD5_mm05_maf005.vcf.gz
```

Finally linked SNPs were pruned using plink.
```
plink2 --set-all-var-ids @# --allow-extra-chr --indep-pairwise 50 5 0.3 --vcf ${myvcf}  --out ${myvcf%%.vcf.gz}_unlinked
```

The following VCF files of filtered SNPs (provided here in the `data/` directory) were used as the basis of all analyses in this repository.
 
* `pteropus.vcf.gz`: Filtered unlinked SNPs for all *Pteropus* species.
* `palecto_pconspicillatus_broad.vcf.gz`: Filtered unlinked SNPs for *P. conspicillatus* and *P. alecto*, including *P. alecto alecto*.
* `palecto_pconspicillatus_narrow.vcf.gz`: Filtered unlinked SNPs for *P. conspicillatus* and *P. alecto*, excluding *P. alecto alecto*.


