#' add minor allele counts from a VCF to SummarizedExperiment with common samples
#' @param se SummarizedExperiment
#' @param vcfpath path to VCF with genotypes for samples
#' @param variantRanges range to filter variants before adding to SummarizedExperiment
#' @return an updated SummarizedExperiment with metadata elements snpaddrs and snpchrs and colData
#' augmented with genotype calls.
#' @note If a metadata element nonCallVars is present, the function will exit with error.
#' @examples
#' vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
#' data("mageSE_19", package="biocXqtl")
#' nn = add_calls_from_vcf(mageSE_19, vp, GRanges("19:1-150000"))
#' @export
add_calls_from_vcf = function(se, vcfpath, variantRanges) {
   nonc = metadata(se)$nonCallVars
   if (length(nonc)>0) stop("calls already present?  nonCallVars is set in metadata.")
   m = minorAlleleCounts(vcfpath, variantRanges)
   w = width(m)
   dro = which(w>1)
   if (length(dro)>0) {
     message("non-SNV variants found, dropping")
     m = m[-dro,]
     }
   saddrs = start(m)
   calls = t(as.matrix(mcols(m)))
   oksamp = intersect(colnames(se), rownames(calls))
   se = se[, oksamp]
   ncv = colnames(colData(se))
   metadata(se)$nonCallVars = ncv
   colData(se) = cbind(DataFrame(calls[oksamp,]), colData(se))
   metadata(se)$snpaddrs = saddrs
   metadata(se)$snpchrs = as.character(seqnames(m))
   se
}
