<div id="main" class="col-md-9" role="main">

# simple function for MAF calculation

<div class="ref-description section level2">

simple function for MAF calculation

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
maf(snpmat)
```

</div>

</div>

<div class="section level2">

## Arguments

-   snpmat:

    an entity that answers \`apply(...,2,f)\` with columns corresponding
    to SNPs, rows corresponding to minor allele count for individuals

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
data("geuv19", package="biocXqtl")
m = maf(colData(geuv19))
summary(m)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.07143  0.00000  0.04945  0.14139  0.18132  1.00000 
```

</div>

</div>

</div>
