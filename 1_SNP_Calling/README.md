## Alignment and variant calling

Raw reads are available from the SRA (BioProject [PRJNA1230740](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1230740)). Libraries were prepared with PstI/EcoRI double digestion and sequenced with 75 bp paired-end reads on an Illumina NextSeq. Reads were demultiplexed and trimmed to 43 bp with stacks 2.52 `process_radtags`, and quality-checked with FastQC.

Exploratory analyses with the stacks pipeline, with and without a reference, showed that a BWA–BCFtools approach produced the most high-quality SNPs after filtering. Reads were aligned to the *P. alecto* reference genome `GCF_000325575.1_ASM32557v1_genomic.fna` (GenBank) using BWA-MEM 0.7.17-r1198-dirty with default settings:

```
bwa mem -t 1 GCF_000325575.1_ASM32557v1_genomic.fna <sample>.1.fq.gz <sample>.2.fq.gz | samtools sort --threads 1 > <sample>.bam
```

The low-quality 5 bp flanks were clipped from the alignments with GATK 4.6.2.0:

```
gatk ClipReads \
    -I <sample>.bam \
    -O <sample>.clip.bam \
    --cycles-to-trim "1-5,39-43" \
    --clip-representation WRITE_NS
```

Variants were called with BCFtools 1.19, including invariant sites to inform demographic modelling. Hardy-Weinberg equilibrium was assumed within populations (`-G groups.txt`, which assigns each sample to its population). `rg_bam.list` lists the read-group-tagged BAM files.

```
bcftools mpileup --threads 48 -q 20 -Q 10 -f GCF_000325575.1_ASM32557v1_genomic.fna --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -b rg_bam.list | bcftools call --threads 48 -G groups.txt -mO z -f GQ -o calls_withinvariant.vcf.gz
```

After an initial round of calling, quality control removed 46 of 288 samples (alignment metrics, preliminary phylogenetic species assignment, duplicates and first-degree relatives; Text S2, Table S1). Duplicates and first-degree relatives were identified as pairs with KING kinship > 0.2 (`plink2 --make-king-table`, see `2_SNP_Stats`). Variant calling was then repeated with the same parameters for the final 242 samples.

### Filtering

The calls were split into three datasets, and each was filtered separately to minimize the number of SNPs lost to missingness in species-specific analyses (Table S2). Using VCFtools 0.1.16, we first removed indels and genotype calls with low depth or genotype quality:

```
vcftools --gzvcf calls_withinvariant.vcf.gz --remove-indels --minDP 5 --minGQ 20 --recode --recode-INFO-all --stdout | gzip -c > snps.minDP5.minGQ20.vcf.gz
```

Next, we removed rare alleles (MAF < 0.05) and SNPs with more than 20% missing genotypes:

```
vcftools --gzvcf snps.minDP5.minGQ20.vcf.gz --maf 0.05 --max-missing 0.8 --recode --recode-INFO-all --stdout | gzip -c > snps.minDP5.minGQ20.mm02.maf005.vcf.gz
```

Multi-allelic sites were removed:

```
bcftools view --types snps -m 2 -M 2 -Oz -o snps.minDP5.minGQ20.mm02.maf005.biallelic.vcf.gz snps.minDP5.minGQ20.mm02.maf005.vcf.gz
```

Linked SNPs were pruned with PLINK2 2.00a2.3LM. The retained variants in `<prefix>.prune.in` were then extracted from the VCF:

```
plink2 --set-all-var-ids @# --allow-extra-chr --indep-pairwise 50 5 0.3 --vcf snps.minDP5.minGQ20.mm02.maf005.biallelic.vcf.gz --out snps.minDP5.minGQ20.mm02.maf005.biallelic.unlinked
```

### Modified filters for specific analyses

* **Heterozygosity and nucleotide diversity** (`2_SNP_Stats`): only the depth, genotype quality and missingness filters were applied. Invariant, low-MAF, multiallelic and linked sites were retained.
* **fastsimcoal2** (`5_Fastsimcoal`): depth, genotype quality, missingness and biallelic filters were applied, with no MAF or LD filter. The missingness filter was reapplied to each population-pair-specific VCF.

### Data

These filtered unlinked SNP sets are the basis of the analyses in this repository:

| File | Dataset (manuscript) | Samples | SNPs |
|---|---|---|---|
| `data/pteropus.vcf.gz` | "all": all four species | 242 | 11,818 |
| `data/palecto_pconspicillatus_broad.vcf.gz` | "alecto + conspicillatus", including *P. alecto alecto* | 150 | 17,708 |
| `data/palecto_pconspicillatus_narrow.vcf.gz` | "alecto ex alecto alecto + conspicillatus" | 141 | 18,285 |
