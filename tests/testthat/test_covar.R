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
   litdf = data.frame(y=as.numeric(mol[1,]), snp=calls[,1], Sex=covar[,1])

test_that("getzs with one covariate works", {
# use lm for a simple model
   m1 = summary(lm(y~snp+Sex, data=litdf))
   b = m1$coefficients[2,1]
   se = m1$coefficients[2,2]
# use getzs for the same model
   full = getzs(mol[1:10,], calls, covar)
   expect_equal(b/se, full[1,1])
})


# note that the GeuvadisTranscriptExpr/geuv19 "counts" are expected counts from FluxCapacitor...
# thus rounding is used below
test_that("negative binomial fits work with covariate", { 
  nbz = function(x,y) {
  # error messages will be thrown from glm ... maybe condition to allow warning/error
    f = tryCatch(MASS::glm.nb(y~x, data=data.frame(y=y, x=x[,-1])), # switch to -1 here to allow covariates
             error=function(e) return(list(coefficients=NA, se=NA) ))  # getzs adds column of 1s
    #if (inherits(f, "try-error")) return(list(coefficients=NA, se=NA))
    dat = summary(f)$coefficients
    list(coefficients=dat[,1], se=dat[,2])
  } 
   expect_warning({
      full = getzs(round(mol[1:10,],0), calls, covar, statfun=nbz) # includes sex adjustment
   })
   lit = summary(MASS::glm.nb(round(y,0)~snp+Sex, data=litdf))
   expect_true(abs(full[1,1]-lit$coefficients[2,3])<1e-12)
})


