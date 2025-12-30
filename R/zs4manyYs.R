#' use procedure tailored to multiple responses for a fixed design matrix
#' @import BiocParallel
#' @import MatrixGenerics
#' @param se RangedSummarizedExperiment
#' @param omit.hi.maf logical(1) if TRUE, variants for which MAF exceeds 0.5 are omitted
#' @param BPPARAM instance of BiocParallelParam to control parallel execution over genotype calls
#' @return matrix with rows corresponding to molecular features and columns corresponding to genotypes
#' @examples
#' data(geuv19xse)
#' sds = MatrixGenerics::rowSds(assay(geuv19xse))
#' print(quantile(sds))
#' BiocParallel::register(BiocParallel::SerialParam())
#' mafs = maf(geuv19xse)
#' mins = apply(data.matrix(mcols(getCalls(geuv19xse))), 1, min, na.rm=TRUE) # some -1 values
#' print(quantile(mins))
#' lk = filterCalls(geuv19xse, which(mafs>.25 & mins > -1))
#' lk = lk[which(sds>median(sds, na.rm=TRUE)),]
#' chk1 = zs4manyYs(lk)
#' data(geuv19_samples)
#' namedSex = geuv19_samples$Sex
#' names(namedSex) = geuv19_samples[["Sample name"]]
#' lk$Sex = namedSex[colnames(lk)]
#' table(lk$Sex)
#' chk2 <- zs4manyYs(lk) # use covariate
#' plot(as.numeric(chk2), as.numeric(chk1))
#' @export
zs4manyYs = function (se, omit.hi.maf=FALSE, BPPARAM=BiocParallel::bpparam())
{
    covmat = NULL
    calls = data.matrix(mcols(getCalls(se)))
    m = maf(se)
    omtag = "; omitting"
    if (!omit.hi.maf) omtag = ""
    if (any(m > 0.5)) message(sprintf("some variants have MAF > 0.5 %s", omtag))
    mins = apply(calls,1,min,na.rm=TRUE)  
    maxs = apply(calls,1,max,na.rm=TRUE)  
    if (any(mins < 0)) message("negative minor allele count found, dropping SNP")
    if (any(maxs > 2)) message("minor allele count(s) > 2 found, omitting")
    okc = which(mins >=0)
    if (omit.hi.maf) okc = which(mins >=0 & m <= .5)
    se = filterCalls(se, okc)
    Y = t(as.matrix(assay(se)))
    calls = t(data.matrix(mcols(getCalls(se))))
    covmat = NULL
    if (ncol(colData(se))>0) covmat = data.matrix(as.data.frame(colData(se)))
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
    colnames(ans) = colnames(calls)
    ans
}
