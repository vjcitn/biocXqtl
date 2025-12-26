<div id="main" class="col-md-9" role="main">

# use procedure tailored to multiple responses for a fixed design matrix

<div class="ref-description section level2">

use procedure tailored to multiple responses for a fixed design matrix

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
zs4manyYs(se, seq4mol = 0.5, minmaf = 0.3, BPPARAM = BiocParallel::bpparam())
```

</div>

</div>

<div class="section level2">

## Arguments

-   se:

    RangedSummarizedExperiment

-   seq4mol:

    quantile of rowSds(assay(se)) at which to filter away features with
    low variability over samples

-   minmaf:

    lower bound on MAF for genotype calls

-   BPPARAM:

    instance of BiocParallelParam to control parallel execution over
    genotype calls

</div>

<div class="section level2">

## Value

matrix with rows corresponding to molecular features and columns
corresponding to genotypes

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
data(geuv19)
sds = MatrixGenerics::rowSds(assay(geuv19))
BiocParallel::register(BiocParallel::SerialParam())
mafs = maf(colData(geuv19)) # only snps here
mins = apply(data.matrix(as.data.frame(colData(geuv19))), 2, min, na.rm=TRUE) # some -1 values
print(quantile(mins))
#>   0%  25%  50%  75% 100% 
#>   -1    0    0    0    2 
colData(geuv19) = colData(geuv19)[,which(mafs>.25 & mins > -1)]
selz = zs4manyYs(geuv19, minmaf = median(mafs))
print(dim(selz))
#> [1] 3629  253
print(selz[1:5,1:5])
#>                 snp_19_1400679 snp_19_1400766 snp_19_5694231 snp_19_5694630
#> ENST00000545779     0.03610262      0.6642772     -1.2665610    -0.95817616
#> ENST00000318050     1.19105119      0.7522053     -0.2852815    -0.03788855
#> ENST00000327790     1.06226054      0.6297555      1.5138651    -1.59259101
#> ENST00000434325     0.21061342      0.6310024      2.9397710    -1.76515522
#> ENST00000269812     0.68286599      0.6720288     -0.4720674    -0.84038324
#>                 snp_19_5696245
#> ENST00000545779     -1.5318370
#> ENST00000318050      0.1532092
#> ENST00000327790     -0.4960955
#> ENST00000434325      0.8775969
#> ENST00000269812     -1.7454723
# include adjustment for sex
data(geuv19)
data(geuv19_samples)
sds = rowSds(assay(geuv19), na.rm=TRUE)
qq = quantile(sds, .8)
ok = which(sds > qq)
lk = geuv19[ok,]
mafs = maf(colData(lk)) # only snps here
mins = apply(data.matrix(as.data.frame(colData(lk))), 2, min, na.rm=TRUE) # some -1 values
colData(lk) = colData(lk)[,which(mafs>.25 & mins > -1)]
namedSex = geuv19_samples$Sex
names(namedSex) = geuv19_samples[["Sample name"]]
snpn = colnames(colData(lk))
lk$Sex = namedSex[colnames(lk)]
table(lk$Sex)
#> 
#> female   male 
#>     46     45 
metadata(lk) = list(nonCallVars="Sex")
run2 <- zs4manyYs(lk, seq4mol = 0.5, minmaf = .3)
```

</div>

</div>

</div>
