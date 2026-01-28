
#' operate with BiocFileCache to retrieve SummarizedExperiment banovich 450k data
#' @import SummarizedExperiment
#' @param ca BiocFileCache instance
#' @examples
#' bano = get_banovich_SE()
#' bano
#' @export
get_banovich_SE = function(ca = BiocFileCache::BiocFileCache()) {
  url = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocXQTLsupport/banoMethSE2.rda"
  inf = bfcquery(ca, "banoMethSE2.rda")
  nr = nrow(inf)
  if (nr == 0) {
    bfcadd(ca, fpath=url, rname="banoMethSE2.rda", download=TRUE)
  }
  inf = bfcquery(ca, "banoMethSE2.rda")$rpath
  get(load(inf))
}
