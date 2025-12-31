<div id="main" class="col-md-9" role="main">

# getter for genotype calls

<div class="ref-description section level2">

getter for genotype calls

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
getCalls(xse)
```

</div>

</div>

<div class="section level2">

## Arguments

-   xse:

    instance of XqtlExperiment

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
data(mageSE_19)
vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
mp = minorAlleleCounts(vp)
xmage19 = XqtlExperiment(mageSE_19, mp)
getCalls(xmage19)[1:4,1:4]
#> GRanges object with 4 ranges and 4 metadata columns:
#>               seqnames      ranges strand |   HG00096   HG00100   HG00105
#>                  <Rle>   <IRanges>  <Rle> | <numeric> <numeric> <numeric>
#>   rs541392352       19       60842      * |         0         0         0
#>   rs534193774       19 62935-62937      * |         0         0         0
#>   rs559839262       19       64705      * |         0         0         0
#>   rs372156287       19       69984      * |         0         0         1
#>                 HG00108
#>               <numeric>
#>   rs541392352         0
#>   rs534193774         0
#>   rs559839262         0
#>   rs372156287         0
#>   -------
#>   seqinfo: 1 sequence from GRCh38 genome; no seqlengths
```

</div>

</div>

</div>
