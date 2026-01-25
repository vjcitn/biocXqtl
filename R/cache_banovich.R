
#' operate with BiocFileCache to retrieve SummarizedExperiment banovich 450k data
#' @import SummarizedExperiment
#' @param ca BiocFileCache instance
#' @examples
#' bano = get_banovich_SE()
#' bano
#' @export
get_banovich_SE = function(ca = BiocFileCache::BiocFileCache()) {
  url = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocXQTLsupport/banoMethSE.rda"
  inf = bfcquery(ca, "banoMethSE.rda")
  nr = nrow(inf)
  if (nr == 0) {
    bfcadd(ca, fpath=url, rname="banoMethSE.rda", download=TRUE)
  }
  inf = bfcquery(ca, "banoMethSE.rda")$rpath
  get(load(inf))
}
