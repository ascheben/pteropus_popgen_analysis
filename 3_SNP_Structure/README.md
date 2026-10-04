## Population structure, phylogenies and hybridization

We characterized population structure with principal component analysis (PCA), fastSTRUCTURE, maximum likelihood (RAxML) and coalescent (CASTER-site) phylogenies, and triangle plots of hybrid index and interclass heterozygosity. These support Figures 3, 4C, 4D, S1, S9–S15 and Tables S11–S13 and S17.

The filtered unlinked SNP sets from `1_SNP_Calling/data/` are used:
- **"all"** (`pteropus.vcf.gz`, 242 samples, 11,818 SNPs)
- **"alecto + conspicillatus"** (`palecto_pconspicillatus_broad.vcf.gz`, 150 samples, 17,708 SNPs)
- **"alecto ex alecto alecto + conspicillatus"** (`palecto_pconspicillatus_narrow.vcf.gz`, 141 samples, 18,285 SNPs)

### PCA

PCA was run with SNPRelate 1.40.0 inside `scripts/pca_and_phylo_plot.R`, using the "all" set for the species PCA (Figure 3C) and the "alecto ex alecto alecto + conspicillatus" set for *P. alecto* and *P. conspicillatus* (Figure 3D). Figure S9 repeats the latter PCA with only the 216 SNPs that have no missing genotypes (equivalent to `vcftools --max-missing 1`). The VCFs are converted to GDS files on the first run (`data/*.gds`, not tracked).

### Isolation by distance

For the 86 *P. alecto* samples (`Species` = BFF in the metadata), the genetic distance between each pair of samples is the proportion of allelic differences (`poppr::diss.dist`), calculated from the 216 SNPs without missing genotypes. It was compared with the great-circle geographic distance using a Mantel test (Pearson, 999 permutations; r = 0.26, p = 0.001; Figure 4D, Table S17).

### Admixture with fastSTRUCTURE

The "alecto + conspicillatus" SNPs were converted to PLINK format, and fastSTRUCTURE 1.0 was run with the default prior for K = 1–15. K = 4 was selected with `chooseK.py`.

```
plink2 --vcf palecto_pconspicillatus_broad.vcf.gz --set-all-var-ids @:# --allow-extra-chr --make-bed --out palecto_pconspicillatus_broad
for K in $(seq 1 15); do
    structure.py -K $K --input=palecto_pconspicillatus_broad --output=structure_palecto_pconspicillatus_out
done
chooseK.py --input=structure_palecto_pconspicillatus_out
```

### Phylogenetic inference

SNPs were converted to PHYLIP format with [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip). Trees were inferred with RAxML 8.2.12 (GTRCAT, rapid bootstrap, 100 replicates) and with CASTER-site 1.23.2.6:

```
python vcf2phylip.py -i pteropus.vcf.gz # output is phylip format file pteropus.phy
raxmlHPC-PTHREADS -T 24 -f a -m GTRCAT -p 12345 -x 12345 -# 100 -s pteropus.phy -n pteropus
raxmlHPC-PTHREADS -T 24 -f a -m GTRCAT -p 12345 -x 12345 -# 100 -s palecto_pconspicillatus.phy -n palecto_pconspicillatus
```

### Hybrid index and triangle plots

Ancestry-informative markers (AIMs) with an allele frequency difference ≥ 0.85 between *P. alecto* and *P. conspicillatus* were selected with [triangulaR](https://github.com/omys-omics/triangulaR), using the "alecto ex alecto alecto + conspicillatus" set (15 AIMs). Individuals missing one third or more of the AIMs were excluded.

### Data

| File | Content |
|---|---|
| `data/structure/structure_palecto_pconspicillatus_out.<K>.{meanQ,meanP,log}` | fastSTRUCTURE output for K = 1–15 |
| `data/structure/structure_palecto_pconspicillatus_labels.txt` | sample labels in `.meanQ` row order |
| `data/structure/*.fam` | PLINK `.fam` of the fastSTRUCTURE input; sample IDs in `.meanQ` row order |
| `data/raxml/RAxML_bipartitions.{pteropus,palecto_pconspicillatus}` | RAxML best trees with bootstrap support |
| `data/raxml/CASTER.{pteropus,palecto_pconspicillatus}` | CASTER-site trees (Figures S14, S15) |

### Reproducing the figures

```
cd scripts
Rscript pca_and_phylo_plot.R           # Figures 3A-D, 4C, 4D, S1, S9, S12-S15; Tables S11, S17
Rscript structure_plot.R               # Figure 3E
Rscript triangle_hybridization_plot.R  # Figures S10, S11
```

Outputs are written to `figures/` (see `figures/README.md`). The species trees are rooted with *P. scapulatus*, and the *P. alecto* / *P. conspicillatus* (and, for CASTER, *P. alecto alecto*) clades are collapsed for Figures S13 and S15. The CASTER *P. alecto* / *P. conspicillatus* tree is rooted with *P. alecto alecto* (Figure S14).
