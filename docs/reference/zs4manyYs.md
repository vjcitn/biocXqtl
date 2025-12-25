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
m = maf(colData(geuv19)) # only genotypes
sds = MatrixGenerics::rowSds(assay(geuv19))
print(quantile(m)) # note negatives!
#>          0%         25%         50%         75%        100% 
#> -0.07142857  0.00000000  0.04945055  0.18131868  1.00000000 
print(quantile(sds))
#>           0%          25%          50%          75%         100% 
#>     0.000000     2.220802    24.784532   152.732556 51865.584335 
selz = zs4manyYs(geuv19, minmaf = median(m))
print(dim(selz))
#> [1] 3629  647
print(selz[1:5,1:5])
#>                 snp_19_1397443 snp_19_1398143 snp_19_1399056 snp_19_1400679
#> ENST00000545779      0.3803266     0.40013841      2.1452455     0.03610262
#> ENST00000318050      0.6361947    -0.95286456     -1.2390239     1.19105119
#> ENST00000327790     -0.6058847    -0.53425465     -0.4662258     1.06226054
#> ENST00000434325      0.2122657    -0.07532599     -0.1338932     0.21061342
#> ENST00000269812     -0.4329748    -0.47618967     -0.8334847     0.68286599
#>                 snp_19_1400766
#> ENST00000545779      0.6642772
#> ENST00000318050      0.7522053
#> ENST00000327790      0.6297555
#> ENST00000434325      0.6310024
#> ENST00000269812      0.6720288
```

</div>

</div>

</div>
