<div id="main" class="col-md-9" role="main">

# extend RangedSummarizedExperiment to include genotype calls in a GRanges

<div class="ref-description section level2">

extend RangedSummarizedExperiment to include genotype calls in a GRanges

</div>

<div class="section level2">

## Note

We use RangedSummarizedExperiment to ensure we can identify genomic
distance between molecular features and variants

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
data(mageSE_19)
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
mp = mp[, colnames(mageSE_19)]
new("XqtlExperiment", mageSE_19, calls=mp)
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
