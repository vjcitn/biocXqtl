<div id="main" class="col-md-9" role="main">

# compute putative minor allele frequency for XqtlExperiment

<div class="ref-description section level2">

compute putative minor allele frequency for XqtlExperiment

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
maf(xse)
```

</div>

</div>

<div class="section level2">

## Arguments

-   xse:

    XqtlExperiment instance

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
example(XqtlExperiment)  # makes nn
#> 
#> XqtlEx> data(mageSE_19)
#> 
#> XqtlEx> vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
#> 
#> XqtlEx> mp = minorAlleleCounts(vp)
#> 
#> XqtlEx> nn = XqtlExperiment(mageSE_19, mp)
#> 
#> XqtlEx> nn
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
summary(maf(nn))
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.000000 0.000000 0.000000 0.017784 0.000684 0.999316 
```

</div>

</div>

</div>
