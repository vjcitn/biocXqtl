
#' extend RangedSummarizedExperiment to include genotype calls in a GRanges
#' @rawNamespace importFrom("methods", "callNextMethod", "new", "slot", "slot<-")
#' @note We use RangedSummarizedExperiment to ensure we can identify genomic
#' distance between molecular features and variants
#' @examples
#' data(mageSE_19)
#' vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
#' mp = minorAlleleCounts(vp)
#' mp[,1:5]
#' mp = mp[, colnames(mageSE_19)]
#' new("XqtlExperiment", mageSE_19, calls=mp)
#' @export
setClass("XqtlExperiment", contains="RangedSummarizedExperiment",
  slots=c("calls"="GRanges"))

#' present concise view of XqtlExperiment
#' @rawNamespace import(BiocGenerics, except=c(Position))
#' @param object instance of XqtlExperiment
#' @export
setMethod("show", "XqtlExperiment", function(object) {
 callNextMethod()
 cat(sprintf("  %d genotype calls present.\n", length(slot(object, "calls"))))
 cat(sprintf("  use getCalls() to see them with addresses.\n"))
})

#' XqtlExperiment constructor
#' @param se SummarizedExperiment instance
#' @param calls GRanges instance
#' @examples
#' data(mageSE_19)
#' vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
#' mp = minorAlleleCounts(vp)
#' XqtlExperiment(mageSE_19, mp)
#' @export
XqtlExperiment = function(se, calls) {
 okids = intersect(colnames(mcols(calls)), colnames(se))
 se = se[,okids]
 calls = calls[, okids]
 new("XqtlExperiment", se, calls=calls)
}

#' getter for genotype calls
#' @param xse instance of XqtlExperiment
#' @examples
#' data(mageSE_19)
#' vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
#' mp = minorAlleleCounts(vp)
#' xmage19 = XqtlExperiment(mageSE_19, mp)
#' getCalls(xmage19)[1:4,1:4]
#' @export
getCalls = function(xse) slot(xse, "calls")

#' compute putative minor allele frequency for XqtlExperiment
#' @param xse XqtlExperiment instance
#' @export
maf = function(xse) apply(mcols(getCalls(xse)),1, function(x) sum(x)/(2*length(x)))

#' filter the call component of an XqtlExperiment
#' @param xse XqtlExperiment instance
#' @param keep vector of indices or logical selections
#' @export
filterCalls = function(xse, keep) { slot(xse, "calls") = slot(xse, "calls")[keep,]; xse }
