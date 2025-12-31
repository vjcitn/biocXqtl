library(biocXqtl)
library(VariantAnnotation)
data(mageSE_19) # produced from MAGE "DESeq2" quantifications
cc = colnames(mageSE_19)
data(tiles38)
library(GenomeInfoDb)
seqlevelsStyle(tiles38) = "Ensembl"
t19 = tiles38[which(seqnames(tiles38)=="19")]
bigvcf = filterVcfByGenotypeCounts("~/MAGE/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", samples=cc, min_per_genotype=10, regions=range(t19))
