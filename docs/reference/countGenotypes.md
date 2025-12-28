<div id="main" class="col-md-9" role="main">

# Count genotypes across samples for each variant

<div class="ref-description section level2">

Count genotypes across samples for each variant

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
countGenotypes(gt)
```

</div>

</div>

<div class="section level2">

## Arguments

-   gt:

    matrix of genotype calls (variants x samples) with form "0\|0" or
    "0/0", as produced by VariantAnnotion::geno applied to a VCF object

</div>

<div class="section level2">

## Value

DataFrame with counts for each genotype category

</div>

</div>
