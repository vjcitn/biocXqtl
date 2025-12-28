<div id="main" class="col-md-9" role="main">

# produce data.frame with minor allele counts in a region, from VCF

<div class="ref-description section level2">

produce data.frame with minor allele counts in a region, from VCF

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
minorAlleleCounts(
  vcfpath,
  region = GRanges(c("19:1-100000", "19:100001-200000")),
  genome = "GRCh37"
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   vcfpath:

    character path to VCF

-   region:

    a GRanges instance

-   genome:

    character(1)

</div>

<div class="section level2">

## Value

a GRanges with mcols giving the count of minor alleles for each sample
at each locus

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
vp = system.file("vcf", "chr19_50k.vcf.gz", package="biocXqtl")
mp = minorAlleleCounts(vp)
mp[,1:5]
#> GRanges object with 272 ranges and 5 metadata columns:
#>               seqnames      ranges strand |   HG00096   HG00097   HG00099
#>                  <Rle>   <IRanges>  <Rle> | <numeric> <numeric> <numeric>
#>   rs541392352       19       60842      * |         0         0         0
#>   rs534193774       19 62935-62937      * |         0         0         0
#>   rs559839262       19       64705      * |         0         0         0
#>   rs372156287       19       69984      * |         0         0         0
#>   rs201816663       19 80840-80842      * |         1         0         0
#>           ...      ...         ...    ... .       ...       ...       ...
#>   rs577974513       19      199649      * |         0         0         0
#>   rs545361133       19      199690      * |         0         0         0
#>   rs557202556       19      199765      * |         0         0         0
#>   rs575319884       19      199777      * |         0         0         0
#>   rs543254795       19      199879      * |         0         0         0
#>                 HG00100   HG00101
#>               <numeric> <numeric>
#>   rs541392352         0         0
#>   rs534193774         0         0
#>   rs559839262         0         0
#>   rs372156287         0         0
#>   rs201816663         1         1
#>           ...       ...       ...
#>   rs577974513         0         0
#>   rs545361133         0         0
#>   rs557202556         0         0
#>   rs575319884         0         0
#>   rs543254795         0         0
#>   -------
#>   seqinfo: 1 sequence from GRCh37 genome; no seqlengths
dim(mcols(mp))
#> [1]  272 2504
```

</div>

</div>

</div>
