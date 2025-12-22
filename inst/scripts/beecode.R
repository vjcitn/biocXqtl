
library(biocXqtl)
data(geuv19)
mo = function(x) as.numeric(assay(geuv19[x,]))
sn = function(x) as.numeric(colData(geuv19)[,x])
beeswarm::beeswarm(jitter(mo("ENST00000529442"))~sn("snp_19_19631444"), ylab="a transcript of ABCA7 (jittered)", xlab="minor allele count at GRCh38 19:19631444", main="GEUVADIS RNA-seq")

