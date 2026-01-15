#' utility to add annotation to trans results
#' @importFrom dplyr left_join mutate select filter
#' @param xexp XqtlExperiment instance
#' @param trout output of trans.map_trans
#' @return data.frame instance with variant and phenotype feature annotation as
#' new columns.  Note that there may be some excess conversion to character class.
#' @export
add_loc_trans = function(xexp, trout) {
  varlocdf = as.data.frame(granges(getCalls(xexp)))  # is getCalls costly?
  varlocdf$variant_id = rownames(varlocdf)
#  phlocdf = as.data.frame(rowRanges(xexp)) -> fails because gene-id is used as rownames
  rrx = rowRanges(xexp)
  rrcols = colnames(mcols(rrx))
  mrr = data.frame(phenotype_id=names(rrx), seqnames=seqnames(rrx), start=start(rrx), end=end(rrx),
         width=width(rrx), strand=strand(rrx))
  extr = setdiff(rrcols, colnames(mrr))
  ex = lapply(mcols(rrx)[,extr], as.character)  # bad for numerics for now
  phlocdf = cbind(mrr, data.frame(ex))
  lj = left_join(trout, mutate(varlocdf, var_chr=seqnames, var_start=start, var_end=end,
           var_width=width, var_strand=strand, .keep="unused"),
          by="variant_id")
  #lj = left_join(mutate(lj, gene_id=phenotype_id), 
  lj = left_join(lj,
                   mutate(phlocdf, ph_chr=seqnames, ph_start=start, ph_end=end,
                   ph_width=width, ph_strand=strand, .keep="unused"),
          by="phenotype_id")
  lj
}
