#' a RangedSummarizedExperiment with data from GeuvadisTranscriptExpr
#' @docType data
#' @usage data(geuv19)
#' @format RangedSummarizedExperiment
#' @note Minor allele counts are in colData.  Some entries are -1.  This
#' is not explained in GeuvadisTranscriptExpr.
#' @examples
#' data(geuv19)
#' rowRanges(geuv19[1:3,])
"geuv19"

#' data frame with sample information for geuv19
#' @docType data
#' @usage data(geuv19_samples)
#' @format data.frame
"geuv19_samples"

#' SummarizedExperiment from MAGE
#' @docType data
#' @usage data(mageSE_19)
#' @format SummarizedExperiment
"mageSE_19"

#' a RangedSummarizedExperiment with data from GeuvadisTranscriptExpr, in XqtlExperiment format
#' @docType data
#' @usage data(geuv19xse)
#' @format RangedSummarizedExperiment
#' @note Minor allele counts are in calls slot.  Some entries are -1.  This
#' is not explained in GeuvadisTranscriptExpr.
#' @examples
#' data(geuv19xse)
#' rowRanges(geuv19xse[1:3,])
"geuv19xse"

#' a GRanges tiling GRCh38
#' @docType data
#' @note Created with `GenomicRanges::tileGenome(BSgenome.Hsapiens.UCSC.hg38)`
#' @usage data(tiles38)
#' @format GRanges
"tiles38"
