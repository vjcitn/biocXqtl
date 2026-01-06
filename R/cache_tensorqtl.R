#' operate with BiocFileCache to retrieve zip file of tensorQTL output
#' @import BiocFileCache
#' @param ca BiocFileCache instance
#' @examples
#' get_tensor_example()
#' @export
get_tensor_example_path = function(ca = BiocFileCache::BiocFileCache()) {
  url = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocXQTLsupport/GEUVADIS.445_samples.cis_qtl_pairs.chr18.parquet"
  inf = bfcquery(ca, "GEUVADIS.445_samples.cis_qtl_pairs.chr18.parquet")
  nr = nrow(inf)
  if (nr == 0) {
    bfcadd(ca, fpath=url, rname="GEUVADIS.445_samples.cis_qtl_pairs.chr18.parquet", download=TRUE)
  }
  inf = bfcquery(ca, "GEUVADIS.445_samples.cis_qtl_pairs.chr18.parquet")
  inf$rpath
}

#' produce arrow dataset with tensorQTL cis example output
#' @import arrow
#' @examples
#' example_tensorQTL()
#' @export
example_tensorQTL = function() {
  pa = get_tensor_example_path()
  arrow::open_dataset(pa)
}

