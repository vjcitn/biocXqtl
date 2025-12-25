
test_that("getzs with one covariate works", {
# build RSE with nonCallVars in metadata
   if (!exists("geuv19")) data(geuv19)
   data(geuv19_samples)
   namedSex = geuv19_samples$Sex
   names(namedSex) = geuv19_samples[["Sample name"]]
   geuv19$Sex = namedSex[colnames(geuv19)]
   table(geuv19$Sex)
   metadata(geuv19) = list(nonCallVars="Sex")
# process the metadata to develop calls and mol
   cd = colData(geuv19)
   kp = setdiff(colnames(cd), metadata(geuv19)$nonCallVars)
   calls = cd[,kp]
   m = maf(calls)
   calls = calls[,which(m>.3)]
   covar = as.data.frame(cd[,metadata(geuv19)$nonCallVars])
   names(covar) = metadata(geuv19)$nonCallVars
   mol = assay(geuv19)
# use lm for a simple model
   litdf = data.frame(y=as.numeric(mol[1,]), snp=calls[,1], Sex=covar[,1])
   m1 = summary(lm(y~snp+Sex, data=litdf))
   b = m1$coefficients[2,1]
   se = m1$coefficients[2,2]
# use getzs for the same model
   full = getzs(mol[1:10,], calls, covar)
   expect_equal(b/se, full[1,1])
})


