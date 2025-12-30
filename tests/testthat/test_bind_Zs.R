library(S4Vectors)

targ_dim = c(20,259)

targ_mcolnames = c("tx_id", "tx_biotype", "tx_cds_seq_start", "tx_cds_seq_end", 
"gene_id", "tx_name", "snp_19_1400679", "snp_19_1400766", "snp_19_5694231", 
"snp_19_5694630", "snp_19_5696245", "snp_19_11451163", "snp_19_11451648", 
"snp_19_11459260", "snp_19_11464003", "snp_19_11465129")

targ_mcrow1 = new("DFrame", rownames = "ENST00000545779", nrows = 1L, elementType = "ANY", 
    elementMetadata = NULL, metadata = list(), listData = list(
        tx_id = "ENST00000545779", tx_biotype = "unprocessed_pseudogene", 
        tx_cds_seq_start = NA_integer_, tx_cds_seq_end = NA_integer_, 
        gene_id = "ENSG00000266818", tx_name = "ENST00000545779", 
        snp_19_1400679 = 0.0361026244593439, snp_19_1400766 = 0.664277177180487, 
        snp_19_5694231 = -1.26656099928513, snp_19_5694630 = -0.958176164275918, 
        snp_19_5696245 = -1.53183698742795, snp_19_11451163 = -0.390475458383216, 
        snp_19_11451648 = -1.44345290966637, snp_19_11459260 = -0.950967711073628, 
        snp_19_11464003 = -0.390475458383216, snp_19_11465129 = -0.950967711073628))

test_that("bind_Zs result as expected", {
  data(geuv19xse)
  lk = geuv19xse[1:20,]
  mafs = maf(lk) 
  mins = apply(data.matrix(mcols(getCalls(lk))), 1, min, na.rm=TRUE) # some -1 values
  lk = filterCalls(lk,which(mafs>.25 & mins > -1))
  lk = bind_Zs(lk)
  expect_equal(dim(mcols(rowRanges(lk))), targ_dim)
  expect_equal(colnames(mcols(rowRanges(lk)))[1:16], targ_mcolnames)
  expect_equal(mcols(rowRanges(lk))[1,1:16], targ_mcrow1)
})

newdim = c(20,142)
test_that("omit.hi.maf works", {
  data(geuv19xse)
  lk = geuv19xse[1:20,]
  mafs = maf(lk) 
  mins = apply(data.matrix(mcols(getCalls(lk))), 1, min, na.rm=TRUE) # some -1 values
  lk = filterCalls(lk,which(mafs>.25 & mins > -1))
  lk = bind_Zs(lk, omit.hi.maf=TRUE)
  expect_equal(newdim, dim(mcols(lk)))
})
