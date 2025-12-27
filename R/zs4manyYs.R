#' use procedure tailored to multiple responses for a fixed design matrix
#' @import BiocParallel
#' @import MatrixGenerics
#' @param se RangedSummarizedExperiment
#' @param seq4mol quantile of rowSds(assay(se)) at which to filter away features with low variability over samples
#' @param minmaf lower bound on MAF for genotype calls
#' @param BPPARAM instance of BiocParallelParam to control parallel execution over genotype calls
#' @return matrix with rows corresponding to molecular features and columns corresponding to genotypes
#' @examples
#' data(geuv19)
#' sds = MatrixGenerics::rowSds(assay(geuv19))
#' BiocParallel::register(BiocParallel::SerialParam())
#' mafs = maf(colData(geuv19)) # only snps here
#' mins = apply(data.matrix(as.data.frame(colData(geuv19))), 2, min, na.rm=TRUE) # some -1 values
#' print(quantile(mins))
#' colData(geuv19) = colData(geuv19)[,which(mafs>.25 & mins > -1)]
#' selz = zs4manyYs(geuv19, minmaf = median(mafs))
#' print(dim(selz))
#' print(selz[1:5,1:5])
#' # include adjustment for sex
#' data(geuv19)
#' data(geuv19_samples)
#' sds = rowSds(assay(geuv19), na.rm=TRUE)
#' qq = quantile(sds, .8)
#' ok = which(sds > qq)
#' lk = geuv19[ok,]
#' mafs = maf(colData(lk)) # only snps here
#' mins = apply(data.matrix(as.data.frame(colData(lk))), 2, min, na.rm=TRUE) # some -1 values
#' colData(lk) = colData(lk)[,which(mafs>.25 & mins > -1)]
#' namedSex = geuv19_samples$Sex
#' names(namedSex) = geuv19_samples[["Sample name"]]
#' snpn = colnames(colData(lk))
#' lk$Sex = namedSex[colnames(lk)]
#' table(lk$Sex)
#' metadata(lk) = list(nonCallVars="Sex")
#' run2 <- zs4manyYs(lk, seq4mol = 0.5, minmaf = .3)
#' @export
zs4manyYs = function (se, seq4mol = 0.5, minmaf = 0.3, BPPARAM=BiocParallel::bpparam())
{
    ses = MatrixGenerics::rowSds(assay(se))
    ok = which(ses > seq4mol)
    se = se[ok, ]
    ncv = metadata(se)$nonCallVars
    cdn = colnames(colData(se))
    covmat = NULL
    if (length(ncv)>0) {
      snpn = setdiff(cdn, ncv)
      covmat = data.matrix(as.data.frame(colData(se)[,ncv]))
      } else snpn = cdn
    calls = data.matrix(as.data.frame(colData(se)[,snpn]))
    colnames(calls) = colnames(colData(se)[,snpn]) #names get munged otherwise
    m = maf(calls)
    if (any(m > 0.5)) message("some variants have MAF > 0.5, omitting")
    mins = apply(calls,2,min,na.rm=TRUE)  
    maxs = apply(calls,2,max,na.rm=TRUE)  
    if (any(mins < 0)) message("negative minor allele count found")
    if (any(maxs > 2)) message("minor allele count(s) > 2 found, omitting")
    okc = which(m > minmaf & mins >=0 & m <= .5)
    calls = calls[, okc]
    Y = t(as.matrix(assay(se)))
    test_bld = function(cmat) function(x) {
        if (!is.null(cmat)) X = cbind(calls[,x], cmat)
        else X = calls[,x]
        tmp = fastLmMany_R(X, Y)
        tmp$coef[2, ]/tmp$se[2, ]
    }
    tester = test_bld(covmat)
    ans = BiocParallel::bplapply(colnames(calls), function(x) {
        tester(x)
    }, BPPARAM=BPPARAM)
    ans = do.call(cbind, ans)
    colnames(ans) = colnames(colData(se)[,okc])
    ans
}
