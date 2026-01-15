#' utility to add annotation to trans results
#' @importFrom dplyr left_join mutate select filter
#' @param xexp XqtlExperiment instance
#' @param trout output of trans.map_trans
#' @return data.frame instance with variant and phenotype feature annotation as
#' new columns
#' @export
add_loc_trans = function(xexp, trout) {
  varlocdf = as.data.frame(granges(getCalls(xexp)))  # is getCalls costly?
  varlocdf$variant_id = rownames(varlocdf)
  phlocdf = as.data.frame(rowRanges(xexp))
  lj = left_join(trout, mutate(varlocdf, var_chr=seqnames, var_start=start, var_end=end,
           var_width=width, var_strand=strand, .keep="unused"),
          by="variant_id")
  lj = left_join(mutate(lj, gene_id=phenotype_id), 
                   mutate(phlocdf, ph_chr=seqnames, ph_start=start, ph_end=end,
                   ph_width=width, ph_strand=strand, .keep="unused"),
          by="gene_id")
  lj
}
