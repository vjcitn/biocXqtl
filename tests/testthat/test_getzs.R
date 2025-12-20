#> dput(tt[1:2,1:2])
kn = structure(c(0.469670093243737, -0.302763172921657, 0.107161454143419, 
0.0767932518195293), dim = c(2L, 2L), dimnames = list(c("ENST00000429344", 
"ENST00000248420"), c("snp_19_1400679", "snp_19_1400766")))


test_that("getzs works", {
   if (!exists("geuv19")) data(geuv19)
   mol = assay(geuv19)
   calls = data.matrix(as.data.frame(colData(geuv19)))
   m = maf(calls)
   tt = getzs(mol[300:301,], calls[,m>.3])
   expect_true(all.equal(tt[1:2,1:2], kn))
})

