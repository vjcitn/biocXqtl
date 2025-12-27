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
#' data(geuv19)
#' sds = MatrixGenerics::rowSds(assay(geuv19), na.rm=TRUE)
#' qq = quantile(sds, .8)
#' ok = which(sds > qq)
#' lk = geuv19[ok,]
#' mafs = maf(colData(lk)) # only snps here
#' mins = apply(data.matrix(as.data.frame(colData(lk))), 2, min, na.rm=TRUE) # some -1 values
#' colData(lk) = colData(lk)[,which(mafs>.25 & mins > -1)]
#' lk <- bind_Zs(lk, colselector = function(se) colnames(colData(se)))
#' viz_stats(lk)
#' @export
viz_stats = function(se, jitter_fac=500, ptcolor="blue", midchop=2,
   xlabel="SNP addr", ylabel="xQTL association Z") {
  md = metadata(se)
  stopifnot("snpaddrs" %in% names(md))
  uchr = unique(md$snpchrs)
  if (length(uchr)>1) message("snpchrs metadata indicate more than one chromosome, SNPs from different chrs may be superimposed")
  feat <- name <- value <- NULL
  tt = rowData(se)
  sn = colnames(tt)[-c(1:6)]
  addr = md$snpaddrs
  sndf = as.data.frame(tt[,-c(1:6)])
  nc = ncol(sndf)
  sndf$feat = rownames(sndf)
  pp = pivot_longer(sndf, cols=1:nc)
  longaddr = rep(addr, nc)
  pp$addr = longaddr
  pp = pp[abs(pp$value)>midchop,]
  #pp$addr = as.numeric(sapply(strsplit(pp$name, "_"), "[", 3))
  pl <- ggplot(pp, aes(x=jitter(addr,jitter_fac), y=value, 
       text=sprintf("%s<br>%s",feat, name))) + geom_point(colour=ptcolor) +
     xlab(xlabel) + ylab(ylabel)
  ggplotly(pl)
}
