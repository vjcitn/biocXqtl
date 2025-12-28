#' provide path to a demonstration VCF file
#' @note derived from output of
#' `aws s3 cp --no-sign-request s3://1000genomes/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz .` 26 Dec 2025.
#' @note Variants from positions 1-5000000 were filtered by restricting to samples in
#' `mageSE_19` and variants with at least 10 individuals in each genotype class (00, 01, 11).
#' @examples
#' v = demo_vcf()
#' VariantAnnotation::scanVcfHeader(v)
#' @export
demo_vcf = function() {
 system.file("vcf", "filtered_19.vcf.gz", package="biocXqtl")
}
