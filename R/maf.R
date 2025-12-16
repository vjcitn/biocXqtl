#' simple function for MAF calculation
#' @param snpmat an entity that answers `apply(...,2,f)` with columns corresponding to SNPs,
#' rows corresponding to minor allele count for individuals
#' @examples
#' data("geuv19", package="biocXqtl")
#' m = maf(colData(geuv19))
#' summary(m)
#' @export
maf = function(snpmat) apply(snpmat,2, function(x) sum(x)/(2*length(x)))
