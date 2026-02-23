transtab2gr = function(transtab) {
  gr = GRanges(transtab$var_chr, IRanges::IRanges(transtab$var_start, width=1))
  mcols(gr) = DataFrame(transtab)
  gr
  }
 
gwdf2gr = function(gwascatdf) {
  pos = gwascatdf$CHR_POS
  isna = which(is.na(as.numeric(pos)))
  if (length(isna)>0) gwascatdf = gwascatdf[-isna,]
  gwascatdf$CHR_POS = as.numeric(gwascatdf$CHR_POS)
  gr = GRanges(gwascatdf$CHR_ID, IRanges(gwascatdf$CHR_POS, width=1))
  mcols(gr) = DataFrame(gwascatdf)
  gr
  }

#' combine an output of tq_xexp_trans with a GWAS catalog-like resource
#' @importFrom S4Vectors queryHits subjectHits
#' @import IRanges
#' @param transtab output of tq_xexp_trans
#' @param gwascatdf data.frame with CHR_POS, CHR_ID, SNPS (if using linkby='rsid')
#' @param verbose logical(1)
#' @param phtag character(1) name of variable in transtab regarded as molecular phenotype,
#' typically 'symbol' (default)
#' @param linkby character(1), "rsid" or "pos"
#' @export
join_xqtr_gwas = function (transtab, gwascatdf, verbose = TRUE, 
    phtag = "symbol", linkby="rsid") {
    ac = as.character
    if (linkby == "rsid") {
      uuu = inner_join(transtab, mutate(gwascatdf, variant_id = SNPS), 
        by = "variant_id", relationship = "many-to-many")
      }
    else if (linkby == "pos") {
      ttgr = transtab2gr(transtab)     
      gwgr = gwdf2gr(gwascatdf)
      ov = findOverlaps(ttgr, gwgr)
      uuu = cbind(mcols(ttgr)[queryHits(ov),], mcols(gwgr)[subjectHits(ov),])
      uuu = as.data.frame(uuu)  # FIXME -- don't you need to check crossChr here?
    }
    else stop("linkby must be rsid or pos")
    NInputGenes = length(unique(transtab[[phtag]]))
    NInputSNPs = length(unique(transtab[["variant_id"]]))
    NOutputSNPs = length(unique(uuu[["variant_id"]]))
    NOutputGenes = length(unique(uuu[[phtag]]))
    if (verbose) {
        cat(sprintf("%d snps in, %d snps out, %d genes in, %d genes out\n", 
            NInputSNPs, NOutputSNPs, NInputGenes, NOutputGenes))
    }
    class(uuu) = c("xqgw", "data.frame")
    attr(uuu, "parms") = list(phtag = phtag)
    uuu
}

