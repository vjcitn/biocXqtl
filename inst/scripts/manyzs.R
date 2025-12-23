 
manyzs = function(se, seq4mol=.5, minmaf=.3) {
  mol = t(assay(se))
  ses = colSds(mol)
  ok = which(ses > seq4mol)
  se = se[ok,]
  calls = data.matrix(as.data.frame(colData(se))) # only calls
  m = maf(calls)
  okc = which(m > minmaf)
  calls = calls[,okc]
  Y = t(assay(se))
  ans = bplapply(colnames(calls), function(x) {
   tmp = fastLmMany_R(calls[,x], Y) 
   tmp$coef[2,]/tmp$se[2,]
   })
  names(ans) = colnames(calls)
  do.call(cbind, ans)
  }

