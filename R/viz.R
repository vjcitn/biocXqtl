#' produce a plotly display of statistics of xQTL association
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @import SummarizedExperiment
#' @importFrom tidyr pivot_longer
#' @param se SummarizedExperiment post "bind_Zs"
#' @param jitter_fac plot employs jittering of SNP address by this factor
#' @param ptcolor point color
#' @param midchop stats with absolute value less than this value are omitted from display
#' @param xlabel axis label
#' @param ylabel axis label
#' @examples
#' data(geuv19xse)
#' sds = MatrixGenerics::rowSds(assay(geuv19xse), na.rm=TRUE)
#' qq = quantile(sds, .8)
#' ok = which(sds > qq)
#' lk = geuv19xse[ok,]
#' mafs = maf(lk) 
#' mins = apply(data.matrix(mcols(getCalls(lk))), 1, min, na.rm=TRUE)
#' lk = filterCalls(lk, which(mafs>.25 & mins > -1))
#' lk <- bind_Zs(lk)
#' viz_stats(lk)
#' @export
viz_stats = function(se, jitter_fac=500, ptcolor="blue", midchop=2,
   xlabel="SNP addr", ylabel="xQTL association Z") {
#  md = metadata(se)
#  stopifnot("snpaddrs" %in% names(md))
#  uchr = unique(md$snpchrs)
#  if (length(uchr)>1) message("snpchrs metadata indicate more than one chromosome, SNPs from different chrs may be superimposed")
  feat <- name <- value <- NULL
#  tt = rowData(se)
#  sn = colnames(tt)[-c(1:6)]
#  addr = md$snpaddrs
  sndf = as.data.frame(mcols(rowRanges(se)[,-c(1:6)]))
  ii = intersect(colnames(sndf), names(getCalls(se)))  # not all calls are tested
  addr = start(getCalls(se)[ii,]) # get starts for those tested
  sndf = sndf[,ii] # make compatible
#  sndf$feat = rownames(sndf)  # adds column
  pp = pivot_longer(sndf, cols=1:ncol(sndf))
  pp$feat = rep(rownames(sndf), ncol(sndf))
  longaddr = rep(addr, nrow(sndf)) # for feat also
  pp$addr = longaddr
  pp = pp[abs(pp$value)>midchop,]
  #pp$addr = as.numeric(sapply(strsplit(pp$name, "_"), "[", 3))
  pl <- ggplot(pp, aes(x=jitter(addr,jitter_fac), y=value, 
       text=sprintf("%s<br>%s",feat, name))) + geom_point(colour=ptcolor) +
     xlab(xlabel) + ylab(ylabel)
  ggplotly(pl)
}
