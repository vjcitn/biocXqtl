# build RSE with nonCallVars in metadata
   if (!exists("geuv19xse")) data(geuv19xse)
   data(geuv19_samples)
   namedSex = geuv19_samples$Sex
   names(namedSex) = geuv19_samples[["Sample name"]]
   geuv19xse$Sex = namedSex[colnames(geuv19xse)]
   table(geuv19xse$Sex)
# process the metadata to develop calls and mol
   m = maf(geuv19xse)
   calls = t(data.matrix(mcols(getCalls(geuv19xse)[which(m>.3),])))
   covar = as.data.frame(geuv19xse$Sex)
   mol = assay(geuv19xse)
   litdf = data.frame(y=as.numeric(mol[1,]), snp=calls[,1], Sex=covar[,1])

test_that("getzs with one covariate works", {
   source(system.file("old", "getzs.R", package="biocXqtl"))
# use lm for a simple model
   m1 = summary(lm(y~snp+Sex, data=litdf))
   b = m1$coefficients[2,1]
   se = m1$coefficients[2,2]
# use getzs for the same model
   full = getzs(mol[1:10,], calls, covar)
   expect_equal(b/se, full[1,1])
})


# note that the GeuvadisTranscriptExpr/geuv19xse "counts" are expected counts from FluxCapacitor...
# thus rounding is used below
test_that("negative binomial fits work with covariate", { 
  source(system.file("old", "getzs.R", package="biocXqtl"))
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


