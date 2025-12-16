<div id="main" class="col-md-9" role="main">

# biocXqtl: flexible molecular QTL assessment

<div class="section level2">

## Introduction

eQTL analysis generally involves large data volumes of genotype calls.
No standard representation for genotypes has emerged, though VCFs and
plink2 files are frequently encountered.

In this package we provide a flexible approach to estimating
associations between genotype and molecular phenotypes such as gene
expression or DNA methylation. We will use the
RangedSummarizedExperiment structure as the primary element for
structuring the inputs to association test procedures.

</div>

<div class="section level2">

## Preliminary illustration

`geuv19` is a RangedSummarizedExperiment derived from the
[GeuvadisTranscriptExpr](https://bioconductor.org/package/GeuvadisTranscriptExpr)
package. The assay component provides transcript level counts. Genotypes
are available in the colData.

<div id="cb1" class="sourceCode">

``` r
library(biocXqtl)
data(geuv19)
geuv19
```

</div>

    ## class: RangedSummarizedExperiment 
    ## dim: 4278 91 
    ## metadata(1): snplocs
    ## assays(1): counts
    ## rownames(4278): ENST00000545779 ENST00000408051 ... ENST00000313038
    ##   ENST00000401283
    ## rowData names(6): tx_id tx_biotype ... gene_id tx_name
    ## colnames(91): NA06984 NA06985 ... NA12889 NA12890
    ## colData names(1330): snp_19_1392636 snp_19_1393723 ... snp_19_58273441
    ##   snp_19_58274425

<div id="cb3" class="sourceCode">

``` r
assay(geuv19[1:3,1:4])
```

</div>

    ##                    NA06984 NA06985 NA06986 NA06989
    ## ENST00000545779 180.000000      40      58     172
    ## ENST00000408051   0.000000       0       0       0
    ## ENST00000318050   0.805323       0       0       0

<div id="cb5" class="sourceCode">

``` r
rowRanges(geuv19[1:3,])
```

</div>

    ## GRanges object with 3 ranges and 6 metadata columns:
    ##                   seqnames            ranges strand |           tx_id
    ##                      <Rle>         <IRanges>  <Rle> |     <character>
    ##   ENST00000545779       18 15254418-15271744      + | ENST00000545779
    ##   ENST00000408051       19       71973-72110      + | ENST00000408051
    ##   ENST00000318050       19     110643-111696      + | ENST00000318050
    ##                               tx_biotype tx_cds_seq_start tx_cds_seq_end
    ##                              <character>        <integer>      <integer>
    ##   ENST00000545779 unprocessed_pseudogene             <NA>           <NA>
    ##   ENST00000408051                  miRNA             <NA>           <NA>
    ##   ENST00000318050         protein_coding           110679         111596
    ##                           gene_id         tx_name
    ##                       <character>     <character>
    ##   ENST00000545779 ENSG00000266818 ENST00000545779
    ##   ENST00000408051 ENSG00000275604 ENST00000408051
    ##   ENST00000318050 ENSG00000176695 ENST00000318050
    ##   -------
    ##   seqinfo: 319 sequences (1 circular) from GRCh38 genome

<div id="cb7" class="sourceCode">

``` r
colData(geuv19)[1:3,1:4]
```

</div>

    ## DataFrame with 3 rows and 4 columns
    ##         snp_19_1392636 snp_19_1393723 snp_19_1394530 snp_19_1396952
    ##              <integer>      <integer>      <integer>      <integer>
    ## NA06984              1              0              0              0
    ## NA06985              0              0              0              0
    ## NA06986              0              0              0              0

From a computational viewpoint, the simplest assessment of molecular
QTLs involves fitting a linear model to assess additive association of
minor allele dose with the molecular response. From an organizational
viewpoint, a simple approach is to start and end with a
RangedSummarizedExperiment instance. `bind_Zs` carries this out for the
situation in which genotype calls are managed in a natural way in
colData.

<div id="cb9" class="sourceCode">

``` r
m = maf(colData(geuv19))
ok = which(m > .3)
limg = geuv19[1:100,]
colData(limg) = colData(limg)[,ok]
limg = bind_Zs(limg, 
   colselector=function(se) colnames(colData(se)))
limg
```

</div>

    ## class: RangedSummarizedExperiment 
    ## dim: 100 91 
    ## metadata(1): snplocs
    ## assays(1): counts
    ## rownames(100): ENST00000545779 ENST00000408051 ... ENST00000313093
    ##   ENST00000543365
    ## rowData names(243): tx_id tx_biotype ... snp_19_58269173
    ##   snp_19_58274425
    ## colnames(91): NA06984 NA06985 ... NA12889 NA12890
    ## colData names(237): snp_19_1400679 snp_19_1400766 ... snp_19_58269173
    ##   snp_19_58274425

All the Z-scores for tests of association are in the rowRanges of
`limg`. The first 6 elements of the mcols of the rowRanges are related
to annotation. A heatmap of the Z-scores shows a band of transcripts
with all zero counts. Very light tiles correspond to negative Z-scores,
dark tiles to positive.

<div id="cb11" class="sourceCode">

``` r
zs = data.matrix(as.data.frame(mcols(rowRanges(limg))[,-c(1:6)]))
zs[is.na(zs)] = 0
heatmap(zs, Colv=NA, scale="none")
```

</div>

![](biocXqtl_files/figure-html/lklll-1.png)

With some helper functions, the basic data layout can be seen,
illustrating inappropriateness standard linear modeling assumptions
underlying interpretation of the Z-score.

<div id="cb12" class="sourceCode">

``` r
mo = function(x) as.numeric(assay(limg[x,]))
sn = function(x) as.numeric(colData(limg)[,x])
beeswarm::beeswarm(jitter(mo("ENST00000529442"))~sn("snp_19_19631444"))
```

</div>

![](biocXqtl_files/figure-html/lkda-1.png)

</div>

</div>
