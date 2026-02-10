# utility for conversion of XqtlExperiment contents to tensorqtl inputs

xexp2dfs = function(xexp, useINT=FALSE) {
  stopifnot(inherits(xexp, "XqtlExperiment"))
  xmat = as.matrix(assay(xexp))
  if (useINT) xmat = INT(xmat)
  phenotype_df = as.data.frame(xmat)
  rownames(phenotype_df) = names(rowRanges(xexp))
  colnames(phenotype_df) = colnames(xexp)
  phenotype_pos_df = data.frame(chr=seqnames(rowRanges(xexp)), pos=start(rowRanges(xexp)))
  rownames(phenotype_pos_df) = rownames(xexp)
  cal = getCalls(xexp)
  genotype_df = as.data.frame(data.matrix(mcols(cal)))
  rownames(genotype_df) = names(cal) # colnames should be good
  variant_df = data.frame(chrom=seqnames(cal), pos=IRanges::start(cal)) #, index=seq_len(length(cal))-1L)
  rownames(variant_df) = rownames(genotype_df)
  covariates_df = as.data.frame(colData(xexp))
  r2p = function(x) reticulate::r_to_py(x, convert=TRUE)
  list(phenotype_df = r2p(phenotype_df), phenotype_pos_df = r2p(phenotype_pos_df),
    genotype_df = r2p(genotype_df), variant_df = r2p(variant_df), covariates_df = r2p(covariates_df) )
}
  
#genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix,
#                covariates_df=None, paired_covariate_df=None, maf_threshold=0, interaction_df=None, maf_threshold_interaction=0.05,               
#                group_s=None, window=1000000, run_eigenmt=False, logp=False,
#                output_dir='.', write_top=True, write_stats=True, logger=None, verbose=True)

#' tq_xexp_cis runs tensorqtl cis.map_nominal on components from an XqtlExperiment instance
#' @importFrom reticulate import py_require
#' @param xexp XqtlExperiment instance
#' @param window integer(1) region around start to search
#' @param maf_threshold numeric(1) defaults to 0.
#' @param prefix character(1) used as stem for output, defaults to 'tqrun'
#' @param write_parquet_path path to folder used for output, defaults to `tempdir()`
#' @param useINT logical, if TRUE, compute inverse normal transform for each feature
#' @examples
#' example(XqtlExperiment) # makes nn
#' td = tempdir()
#' sexn = ifelse(colData(nn)$sex == "XY", 0., 1.)
#' colData(nn) = NULL
#' nn$sex = sexn
#' tq_xexp_cis(nn, prefix="example_tq", write_parquet_path=td)
#' dir(td, full=TRUE)
#' @export
tq_xexp_cis = function(xexp, window=1000000L, maf_threshold=0, prefix = "tqrun",
    write_parquet_path = tempdir(), useINT=FALSE) {
 stopifnot(inherits(xexp, "XqtlExperiment"))
 lkcd = colData(xexp)
 dcd = sapply(lkcd, class)
 if (any(dcd == "character")) stop("can't use character variables in colData(xexp)")
 reticulate::py_require("tensorqtl")
 reticulate::py_require("rpy2")
 reticulate::py_require("pandas==2.3.3")  # discovered on linux
 reticulate::py_require("pandas_plink==2.3.1")
 reticulate::import("pandas_plink") # needed?
 reticulate::py_require("torch")
 reticulate::py_require("pyarrow")
 reticulate::py_require("fastparquet")
 fp = reticulate::import("fastparquet")  # needed?
 reticulate::py_require("pyarrow")
 pd = reticulate::import("pandas")
 tor = reticulate::import("torch")
 tq = reticulate::import("tensorqtl")
 dev = tor$device("cuda")
 cis = reticulate::import("tensorqtl.cis")

 if (isTRUE(options()$verbose))  print(reticulate::py_config())

 conv = xexp2dfs(xexp, useINT=useINT)

 cis$map_nominal(conv$genotype_df, conv$variant_df, conv$phenotype_df, 
    conv$phenotype_pos_df, prefix, covariates_df=conv$covariates_df,
    window=window, maf_threshold=maf_threshold,
    write_stats=TRUE, write_top=TRUE, output_dir=write_parquet_path)
}


#   map_trans(genotype_df, phenotype_df, covariates_df=None, interaction_s=None, 
#      return_sparse=True, pval_threshold=1e-05, maf_threshold=0.05, 
#      alleles=2, return_r2=False, batch_size=20000, logp=False, 
#      logger=None, verbose=True)
#        Run trans-QTL mapping

#' tq_xexp_trans runs tensorqtl map_trans on components from an XqtlExperiment instance
#' @importFrom reticulate import py_require
#' @param xexp XqtlExperiment instance
#' @param pval_threshold numeric(1) defaults to 1e-05,
#' @param maf_threshold numeric(1) defaults to 0.05
#' @param batch_size numeric(1) defaults to 20000
#' @param useINT logical, if TRUE, apply inverse normal transformation to all features
#' @examples
#' example(XqtlExperiment) # makes nn
#' td = tempdir()
#' sexn = ifelse(colData(nn)$sex == "XY", 0., 1.)
#' colData(nn) = NULL
#' nn$sex = sexn
#' Sys.setenv(OMP_NUM_THREADS = 1)  # these seem important for avoiding crash in macos
#' Sys.setenv(OPENBLAS_NUM_THREADS = 1)
#' Sys.setenv(MKL_NUM_THREADS = 1)
#' lk = tq_xexp_trans(nn, maf_threshold=.01, pval_threshold=1e-4)
#' head(lk)
#' lkint = tq_xexp_trans(nn, maf_threshold=.01, pval_threshold=1e-4, useINT=TRUE)
#' head(lkint)
#' @export
tq_xexp_trans = function(xexp, pval_threshold=1e-5, maf_threshold=0.05, batch_size=20000L,
   useINT=FALSE ) {
 stopifnot(inherits(xexp, "XqtlExperiment"))
 lkcd = colData(xexp)
 dcd = sapply(lkcd, class)
 if (any(dcd == "character")) stop("can't use character variables in colData(xexp)")
 reticulate::py_require("tensorqtl")
 reticulate::py_require("rpy2")
 reticulate::py_require("pandas==2.3.3")
 reticulate::py_require("pandas_plink==2.3.1")
 reticulate::import("pandas_plink") # needed?
 reticulate::py_require("torch")
 reticulate::py_require("pyarrow")
 reticulate::py_require("fastparquet")
 fp = reticulate::import("fastparquet")  # needed?
 reticulate::py_require("pyarrow")
 pd = reticulate::import("pandas")
 tor = reticulate::import("torch")
 tq = reticulate::import("tensorqtl")
 dev = tor$device("cuda")
 trans = reticulate::import("tensorqtl.trans")

 if (isTRUE(options()$verbose))  print(reticulate::py_config())

 conv = xexp2dfs(xexp, useINT=useINT)

#   map_trans(genotype_df, phenotype_df, covariates_df=None, interaction_s=None, 
#      return_sparse=True, pval_threshold=1e-05, maf_threshold=0.05, 
#      alleles=2, return_r2=False, batch_size=20000, logp=False, 
#      logger=None, verbose=True)
#        Run trans-QTL mapping

 transtab = trans$map_trans(conv$genotype_df, conv$phenotype_df, conv$covariates_df,
    pval_threshold=pval_threshold, maf_threshold=maf_threshold,
    batch_size=batch_size)
 transtab = add_loc_trans(xexp=xexp, trout=transtab)
 transtab
}
