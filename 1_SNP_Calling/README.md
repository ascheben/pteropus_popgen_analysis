## Alignment and variant calling

Reads were aligned to the Pteropus alecto reference genome `GCF_000325575.1_ASM32557v1_genomic.fna` from GenBank using bwa mem 0.7.17-r1198-dirty with default settings. Initial exploratory analyses using the stacks pipeline with and without a reference showed that a bwa-bcftools approach produced the most high-quality SNPs after filtering. 

After adapter-trimming and quality control, paired reads were aligned to the reference.

```
bwa mem -t 1 GCF_000325575.1_ASM32557v1_genomic.fna <sample>.1.fq.gz  <sample>.2.fq.gz | samtools sort --threads 1 > <sample>.bam
```

The 5bp flanks were then trimmed from the alignments.
```
gatk ClipReads \
    -I <infile>.bam \
    -O <infile>.clip.bam \
    --cycles-to-trim "1-5,39-43" \
    --clip-representation WRITE_NS
```

SNPs were next called with bcftools and invariant sites were included to inform demographic modelling.
```
bcftools mpileup --threads 48 -q 20 -Q 10 -f GCF_000325575.1_ASM32557v1_genomic.fna --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -b rg_bam.list | bcftools call --threads 48 -G groups.txt -mO z -f GQ -o calls_withinvariant.vcf.gz
```

Genotype calls were filtered using VCFtools. First we removed calls with low genotype depth and/or low genotype quality as well as indels.

```
vcftools --gzvcf calls_withinvariant.vcf.gz --remove-indels --recode --recode-INFO-all --remove-indels --minDP 5 --minGQ 20 --stdout |gzip -c > snps.minDP5.minGQ20.vcf.gz
```

Next we filtered rare alleles with a minor allele frequency below 5% and we also removed SNPs with over 20% missing alleles.

```
vcftools --gzvcf snps.minDP5.minGQ20.vcf.gz --maf 0.05 --max-missing 0.8 --recode --recode-INFO-all --stdout | gzip -c > snps.minDP5.minGQ20.mm02.maf005.vcf.gz
```

We then remove multi-allelic sites from the call set.
```
bcftools view --types snps -m 2 -M 2 -Oz -o snps.minDP5.minGQ20.mm02.maf005.biallelic.vcf.gz snps.minDP5.minGQ20.mm02.maf005.vcf.gz
```

Finally linked SNPs were pruned using plink.
```
plink2 --set-all-var-ids @# --allow-extra-chr --indep-pairwise 50 5 0.3 --vcf snps.minDP5.minGQ20.mm02.maf005.biallelic.vcf.gz --out snps.minDP5.minGQ20.mm02.maf005.biallelic.unlinked
```

The following VCF files of filtered SNPs (provided here in the `data/` directory) were used as the basis of all analyses in this repository.
 
* `pteropus.vcf.gz`: Filtered unlinked SNPs for all *Pteropus* species.
* `palecto_pconspicillatus_broad.vcf.gz`: Filtered unlinked SNPs for *P. conspicillatus* and *P. alecto*, including *P. alecto alecto*.
* `palecto_pconspicillatus_narrow.vcf.gz`: Filtered unlinked SNPs for *P. conspicillatus* and *P. alecto*, excluding *P. alecto alecto*.


