
#' internal function to combine new xQTL test results with existing XqtlExperiment instance
#' @param xse XqtlExperiment instance
#' @param mat a matrix of test results conformant with the XqtlExperiment: rownames agree,
#' @param check.cols logical(1) will message if redundant variants are found when new scores
#' are bound in
#' columns will be added to the mcols(rowRanges(xse))
add_scores = function(xse, mat, check.cols=TRUE) {
 stopifnot(inherits(xse, "XqtlExperiment"))
 stopifnot(inherits(mat, "matrix"))
 stopifnot(all.equal(rownames(xse), rownames(mat)))
 if (check.cols) {
   ov = intersect(colnames(mat), colnames(mcols(xse)))
   if (length(ov)>0) message(sprintf("%d of new variants found in mcols(xse)", length(ov)))
   }
 mcols(xse) = cbind(mcols(xse), mat)
 xse
}

#' compute Z-statistics for xQTL association and bind them to the input x-ome+genome XqtlExperiment
#' @param xse XqtlExperiment instance
#' @param omit.hi.maf logical(1), passed to zs4manyYs, defaults to FALSE
#' @param BPPARAM instance of BiocParallelParam
#' @examples
#' data(geuv19xse)
#' lk = geuv19xse[1:20,]
#' mafs = maf(lk) # only snps here
#' mins = apply(data.matrix(mcols(getCalls(lk))), 1, min, na.rm=TRUE) # some -1 values
#' #filterCalls = function(xse, keep) { slot(xse, "calls") = slot(xse, "calls")[keep,]; xse }
#' lk = filterCalls(lk, which(mafs>.25 & mins > -1))
#' lk = bind_Zs(lk)
#' head(rowRanges(lk)[,7])
#' plpp2tx = as.numeric(assay(lk["ENST00000434325",]))
#' hiz = as.numeric(data.matrix(mcols(getCalls(lk))["snp_19_5694231",]))
#' plot(plpp2tx~jitter(hiz), xlab="snp at 19:5694231", ylab="counts for a transcript of PLPP2")
#' summary(lm(plpp2tx~hiz))
#' @export 
bind_Zs = function(xse, omit.hi.maf=FALSE, BPPARAM=BiocParallel::bpparam()) {
  ans = zs4manyYs(xse, omit.hi.maf, BPPARAM=BPPARAM)
  add_scores(xse, ans, check.cols=TRUE)
}
