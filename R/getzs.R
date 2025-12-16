#' create a matrix of Z statistics for an additive effect of rare allele count in calls on molecular
#' phenotype in molec
#' @import RcppEigen
#' @param molec numeric matrix of molecular phenotype values, columns are samples
#' @param calls numeric matrix of rare allele counts, rows are samples, columns are loci
#' @note For a SNP with MAF 0, NA is returned.
#' @return A matrix with one row per
#' @examples
#' data(geuv19)
#' mol = assay(geuv19)
#' calls = data.matrix(as.data.frame(colData(geuv19)))
#' m = maf(calls)
#' allz = getzs(mol[1:100,], calls[,m>.3])
#' summary(as.numeric(allz))
#' @export
getzs = function(molec, calls) {
   stopifnot(ncol(molec) == nrow(calls))
   nmolec = nrow(molec)
   ncalls = ncol(calls)
   get1 = function(y) function(x) {
     fit = RcppEigen::fastLmPure(cbind(1.,x),y)
     if (any(is.na(fit$coefficients))) return(NA) 
     fit$coefficients[2]/fit$se[2]
     }
   prep = apply(molec,1,get1)
   ans = apply(calls,2,function(x) sapply(prep, function(z) z(x)))
   rownames(ans) = rownames(molec)  # get rid of modifications
   ans
}

#' bind a matrix of Z statistics created with getzs to the rowRanges of a SummarizedExperiment
#' @param se SummarizedExperiment assumed to have molecular phenotype data in assay
#' @param colselector function with argument "se" returning indices of SNP genotypes in colData(se)
#' @examples
#' data(geuv19)
#' lk = geuv19[1:20,]
#' mafs = maf(colData(lk)) # only snps here
#' colData(lk) = colData(lk)[,which(mafs>.25)]
#' lk = bind_Zs(lk)
#' head(rowRanges(lk)[,7])
#' plpp2tx = as.numeric(assay(lk["ENST00000434325",]))
#' hiz = colData(lk)[,"snp_19_5694231"]
#' plot(plpp2tx~jitter(hiz), xlab="snp at 19:5694231", ylab="counts for a transcript of PLPP2")
#' summary(lm(plpp2tx~hiz))
#' @export
bind_Zs = function(se, colselector = function(se) 1:10) {
  molec = as.matrix(assay(se))
  calls = data.matrix(colData(se)[, colselector(se)])
  stats = getzs(molec, calls)
  stopifnot(all(rownames(stats) == rownames(se)))
  mcols(rowRanges(se)) = cbind(mcols(rowRanges(se)), stats)
  se
}
