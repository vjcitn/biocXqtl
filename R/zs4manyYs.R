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
#' m = maf(colData(geuv19)) # only genotypes
#' sds = MatrixGenerics::rowSds(assay(geuv19))
#' print(quantile(m)) # note negatives!
#' print(quantile(sds))
#' selz = zs4manyYs(geuv19, minmaf = median(m))
#' print(dim(selz))
#' print(selz[1:5,1:5])
#' @export
zs4manyYs = function (se, seq4mol = 0.5, minmaf = 0.3, BPPARAM=BiocParallel::bpparam())
{
    ses = MatrixGenerics::rowSds(assay(se))
    ok = which(ses > seq4mol)
    se = se[ok, ]
    calls = data.matrix(as.data.frame(colData(se)))
    colnames(calls) = colnames(colData(se))
    m = maf(calls)
    okc = which(m > minmaf)
    calls = calls[, okc]
    Y = t(as.matrix(assay(se)))
    ans = BiocParallel::bplapply(colnames(calls), function(x) {
        tmp = fastLmMany_R(calls[, x], Y)
        tmp$coef[2, ]/tmp$se[2, ]
    }, BPPARAM=BPPARAM)
    ans = do.call(cbind, ans)
    colnames(ans) = colnames(colData(se)[,okc])
    ans
}
