#' create a matrix of Z statistics for an additive effect of rare allele count in calls on molecular
#' phenotype in molec
#' @import RcppEigen
#' @param molec numeric matrix of molecular phenotype values, columns are samples
#' @param calls numeric matrix of rare allele counts, rows are samples, columns are loci
#' @param covdf data.frame of additional covariates, to be modeled along with genotype
#' @param statfun a function with arguments x, y, with the intention that x is a design matrix
#' lacking a column of 1s and y is a response vector with nrow(x) elements.  `statfun` must
#' return a list with elements `coefficients` and `se`.  `getzs` processes the second element
#' of each of these to form a Z-score, assuming that the first element corresponds to an intercept.
#' @note For a SNP with MAF 0, NA is returned.
#' @return A matrix with rows corresponding to molecular phenotype features and columns corresponding to SNPs
#' @examples
#' data(geuv19)
#' mol = assay(geuv19)
#' sds = MatrixGenerics::rowSds(mol)
#' mol = mol[which(sds>0), ] # drop constant features
#' calls = data.matrix(as.data.frame(colData(geuv19)))
#' csds = apply(calls,2,sd, na.rm=TRUE)
#' mins = apply(calls,2,min, na.rm=TRUE)  # some snps include -1 values
#' m = maf(calls)
#' allz = getzs(mol[1:100,], calls[,m>.3 & csds>0 & mins > -1])
#' summary(as.numeric(allz))
#' \donttest{  # 90 sec?
#'  if (requireNamespace("MASS")) {
#'  nbz = function(x,y) { 
#'  # error messages will be thrown from glm ... maybe condition to allow warning/error
#'    f = tryCatch(MASS::glm.nb(y~x, data=data.frame(y=y, x=x[,-1])), # switch to -1 here to allow covariates
#'             error=function(e) return(list(coefficients=NA, se=NA) ))  # getzs adds column of 1s
#'    #if (inherits(f, "try-error")) return(list(coefficients=NA, se=NA))
#'    dat = summary(f)$coefficients
#'    list(coefficients=dat[,1], se=dat[,2])
#'  }  
#'  # for now NA is noisily returned
#'  allz2 = suppressWarnings(getzs(mol[1:100,], calls[,m>.3 & csds>0 & mins > -1], statfun = nbz))
#'  # do NB results differ substantially from OLS?
#'  plot(allz2[1,], allz[1,], xlab="NB", ylab="OLS", main=sprintf("txQTL Zs for %s", rownames(mol)[1]))
#'  }
#' }
#' @export
getzs = function(molec, calls, covdf = data.frame(), statfun = function(x,y) RcppEigen::fastLmPure(X=x, y=y)) {
   stopifnot(ncol(molec) == nrow(calls))
   nmolec = nrow(molec)
   ncalls = ncol(calls)
   get1 = function(y, covdf) function(x) {
     if (nrow(covdf)>0) x = cbind(x, data.matrix(covdf)) # genotype will always be second component of coeff
     fit = try(statfun(cbind(1,x), y), silent=TRUE) 
     if (inherits(fit, "try-error")) return(NA) 
     if (any(is.na(fit$coefficients))) return(NA)
     fit$coefficients[2]/fit$se[2]
     }
   prep = apply(molec,1,get1,covdf=covdf)
   ans = apply(calls,2,function(x) sapply(prep, function(z) try(z(x))))
   rownames(ans) = rownames(molec)  # get rid of modifications
   ans
}

#' bind a matrix of Z statistics created with getzs to the rowRanges of a SummarizedExperiment
#' @import SummarizedExperiment
#' @param se SummarizedExperiment assumed to have molecular phenotype data in assay.  A metadata
#' component (list element) named nonCallVars will be checked and associated colData elements
#' will be used as covariates in models for effect of dosage of minor allele
#' @param colselector function with argument "se" returning indices of SNP genotypes in colData(se)
#' @examples
#' data(geuv19)
#' lk = geuv19[1:20,]
#' mafs = maf(colData(lk)) # only snps here
#' mins = apply(data.matrix(as.data.frame(colData(lk))), 2, min, na.rm=TRUE) # some -1 values
#' colData(lk) = colData(lk)[,which(mafs>.25 & mins > -1)]
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
  noncallvars = metadata(se)$nonCallVars
  covdf = data.frame()
  if (length(noncallvars)>0) covdf = as.data.frame(colData(se)[,noncallvars])
  stats = getzs(molec, calls, covdf)
  stopifnot(all(rownames(stats) == rownames(se)))
  mcols(rowRanges(se)) = cbind(mcols(rowRanges(se)), stats)
  se
}
