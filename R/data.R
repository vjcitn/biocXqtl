#' a RangedSummarizedExperiment with data from GeuvadisTranscriptExpr
#' @docType data
#' @format RangedSummarizedExperiment
#' @note Minor allele counts are in colData.  Some entries are -1.  This
#' is not explained in GeuvadisTranscriptExpr.
#' @examples
#' data(geuv19)
#' rowRanges(geuv19[1:3,])
"geuv19"
