
#' produce data.frame with minor allele counts in a region, from VCF
#' @import VariantAnnotation
#' @import GenomicRanges
#' @import Seqinfo
#' @param vcfpath character path to VCF
#' @param region a GRanges instance
#' @param genome character(1)
#' @return a GRanges with mcols giving the count of minor alleles for each sample at each locus
#' @examples
#' vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
#' mp = minorAlleleCounts(vp)
#' mp[,1:5]
#' dim(mcols(mp))
#' @export
minorAlleleCounts = function(vcfpath, region = GRanges(c("19:1-100000","19:100001-200000")), genome="GRCh38") {
 rcnt = function (x) 
 {
     z = strsplit(x, "")
     sapply(z, function(x) sum(x=="1"))
 }
 param = ScanVcfParam(which=region)
 dat = scanVcf(vcfpath, param=param)
 rr = lapply(dat, function(x) x$rowRanges)
 mc = unlist(unname(GRangesList(rr)))
 ans = do.call(rbind, lapply(dat, function(x) x$GENO$GT))
 nans = matrix(numeric(1), nrow=nrow(ans), ncol=ncol(ans))
 nans[] = rcnt(as.character(ans))
 colnames(nans) = colnames(ans)
 mcols(mc) = nans
 Seqinfo::genome(mc) = genome
 mc
}
  
   
