<div id="main" class="col-md-9" role="main">

# add minor allele counts from a VCF to SummarizedExperiment with common samples

<div class="ref-description section level2">

add minor allele counts from a VCF to SummarizedExperiment with common
samples

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
add_calls_from_vcf(se, vcfpath, variantRanges)
```

</div>

</div>

<div class="section level2">

## Arguments

-   se:

    SummarizedExperiment

-   vcfpath:

    path to VCF with genotypes for samples

-   variantRanges:

    range to filter variants before adding to SummarizedExperiment

</div>

<div class="section level2">

## Value

an updated SummarizedExperiment with metadata elements snpaddrs and
snpchrs and colData augmented with genotype calls.

</div>

<div class="section level2">

## Note

If a metadata element nonCallVars is present, the function will exit
with error.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
data("mageSE_19", package="biocXqtl")
nn = add_calls_from_vcf(mageSE_19, vp, GRanges("19:1-150000"))
#> non-SNV variants found, dropping
```

</div>

</div>

</div>
