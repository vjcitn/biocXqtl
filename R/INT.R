#' inverse normal transformation
#' @param x matrix for which rows are to be transformed
#' @examples
#' m = matrix(rbeta(2000,2,1),4,500)
#' par(mfrow=c(2,2))
#' mt = INT(m)
#' hist(m[1,]); hist(mt[1,])
#' hist(m[2,]); hist(mt[2,])
#' @export
INT = function(x) {
 stopifnot(length(dim(x))==2)
 t(apply(x,1,function(y) {
  r = rank(y);
  qnorm((r-.5)/length(y))}))
}
