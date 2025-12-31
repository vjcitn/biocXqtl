<div id="main" class="col-md-9" role="main">

# XqtlExperiment constructor

<div class="ref-description section level2">

XqtlExperiment constructor

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
XqtlExperiment(se, calls)
```

</div>

</div>

<div class="section level2">

## Arguments

-   se:

    SummarizedExperiment instance

-   calls:

    GRanges instance

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
data(mageSE_19)
vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
mp = minorAlleleCounts(vp)
nn = XqtlExperiment(mageSE_19, mp)
nn
#> class: XqtlExperiment 
#> dim: 1377 731 
#> metadata(0):
#> assays(1): logcounts
#> rownames(1377): ENSG00000176695 ENSG00000141934 ... ENSG00000213753
#>   ENSG00000099326
#> rowData names(6): gene_id gene_name ... symbol entrezid
#> colnames(731): HG00096 HG00100 ... NA21129 NA21130
#> colData names(13): SRA_accession internal_libraryID ...
#>   RNAQubitTotalAmount_ng RIN
#>   272 genotype calls present.
#>   use getCalls() to see them with addresses.
```

</div>

</div>

</div>
