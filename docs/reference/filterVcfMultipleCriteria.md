<div id="main" class="col-md-9" role="main">

# Filter VCF with multiple criteria

<div class="ref-description section level2">

Filter VCF with multiple criteria

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
filterVcfMultipleCriteria(
  vcf,
  min_per_genotype = NULL,
  min_maf = NULL,
  min_call_rate = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   vcf:

    VCF object

-   min\_per\_genotype:

    Minimum samples per genotype group

-   min\_maf:

    Minimum minor allele frequency

-   min\_call\_rate:

    Minimum call rate

</div>

<div class="section level2">

## Value

Filtered VCF object

</div>

</div>
