<div id="main" class="col-md-9" role="main">

# a RangedSummarizedExperiment with data from GeuvadisTranscriptExpr, in XqtlExperiment format

<div class="ref-description section level2">

a RangedSummarizedExperiment with data from GeuvadisTranscriptExpr, in
XqtlExperiment format

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
data(geuv19xse)
```

</div>

</div>

<div class="section level2">

## Format

RangedSummarizedExperiment

</div>

<div class="section level2">

## Note

Minor allele counts are in calls slot. Some entries are -1. This is not
explained in GeuvadisTranscriptExpr.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
data(geuv19xse)
rowRanges(geuv19xse[1:3,])
#> GRanges object with 3 ranges and 6 metadata columns:
#>                   seqnames            ranges strand |           tx_id
#>                      <Rle>         <IRanges>  <Rle> |     <character>
#>   ENST00000545779       18 15254418-15271744      + | ENST00000545779
#>   ENST00000408051       19       71973-72110      + | ENST00000408051
#>   ENST00000318050       19     110643-111696      + | ENST00000318050
#>                               tx_biotype tx_cds_seq_start tx_cds_seq_end
#>                              <character>        <integer>      <integer>
#>   ENST00000545779 unprocessed_pseudogene             <NA>           <NA>
#>   ENST00000408051                  miRNA             <NA>           <NA>
#>   ENST00000318050         protein_coding           110679         111596
#>                           gene_id         tx_name
#>                       <character>     <character>
#>   ENST00000545779 ENSG00000266818 ENST00000545779
#>   ENST00000408051 ENSG00000275604 ENST00000408051
#>   ENST00000318050 ENSG00000176695 ENST00000318050
#>   -------
#>   seqinfo: 319 sequences (1 circular) from GRCh38 genome
```

</div>

</div>

</div>
