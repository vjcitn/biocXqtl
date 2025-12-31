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

## Working with variants from VCF

Thanks to the AWS Open Data project, tabix-indexed VCFs for genotyped
lymphoblastoid cell lines from the 1000 genomes project are
[available](https://registry.opendata.aws/1000-genomes/). Data from the
Multi-Ancestry analysis of Gene Expression are
[available](https://github.com/mccoy-lab/MAGE) via Zenodo or Dropbox. In
our experience the Dropbox link is more performant.

We have acquired the “DESeq2” gene quantifications from the MAGE
archive, and a selection of `release` genotypes from the associated cell
lines.

<div class="section level3">

### Data on expression and genotype

The expression data:

<div id="cb1" class="sourceCode">

``` r
library(biocXqtl)
```

</div>

    ## Warning: package 'GenomicRanges' was built under R version 4.5.2

<div id="cb3" class="sourceCode">

``` r
data(mageSE_19)
```

</div>

Genotype data:

<div id="cb4" class="sourceCode">

``` r
dv = demo_vcf()
h = VariantAnnotation::scanVcfHeader(dv)
head(samples(h))
```

</div>

    ## [1] "HG00096" "HG00100" "HG00105" "HG00108" "HG00110" "HG00113"

<div id="cb6" class="sourceCode">

``` r
length(intersect(samples(h), colnames(mageSE_19)))
```

</div>

    ## [1] 731

</div>

<div class="section level3">

### Building an XqtlExperiment

We will focus on SNVs but the genotype information has other types of
variants.

Bind the minor allele counts to the expression data:

<div id="cb8" class="sourceCode">

``` r
mins  = minorAlleleCounts(demo_vcf(), GRanges("19:1-50000000"))
mxx = XqtlExperiment(mageSE_19, mins)
colData(mxx) = NULL
```

</div>

The last command above allows us to compute crude measures of
association. If we populate colData with covariate information, test
procedures in the package will incorporate it.

</div>

<div class="section level3">

### Creating association statistics with visualizations

The following step uses C++ modules to compute association tests for all
genotypes and all gene expression measures in mageSE\_19.

<div id="cb9" class="sourceCode">

``` r
sds = rowSds(assay(mxx), na.rm=TRUE)
qq = quantile(sds, .95)
ok = which(sds > qq)
system.time(zzz <- bind_Zs(mxx[ok,]))
```

</div>

    ## some variants have MAF > 0.5

    ##     user   system  elapsed 
    ## 1364.280   26.647  184.635

Here’s a helper function to visualize one association.

<div id="cb12" class="sourceCode">

``` r
onebox = function(xse, mfeat="ENSG00000174837", vnt="rs4897932", title) {
if (missing(title)) title=""
boxplot(split(as.numeric(assay(xse[mfeat,])),
   as.numeric(data.matrix(mcols(getCalls(xse)[vnt,])))),
   ylab=mfeat, xlab=vnt, main=title)
}
onebox(mxx[ok,], title="MAGE eQTL")
```

</div>

![](biocXqtl_files/figure-html/helper-1.png)

We can also visualize interactively:

<div id="cb13" class="sourceCode">

``` r
viz_stats(zzz, midchop=5)
```

</div>

<div id="htmlwidget-ac96cb3ee4656e2e9ec3"
class="plotly html-widget html-fill-item"
style="width:700px;height:432.632880098888px;">

</div>

</div>

<div class="section level3">

### Covariate adjustment

First we reanalyze with adjustment for continental group.

<div id="cb14" class="sourceCode">

``` r
data(mageSE_19)
sds = rowSds(assay(mageSE_19), na.rm=TRUE)
qq = quantile(sds, .9)
ok = which(sds > qq)
cd = colData(mageSE_19)
mins  = minorAlleleCounts(demo_vcf(), GRanges("19:1-3000000"))
mxx = XqtlExperiment(mageSE_19[ok,], mins)
print(mxx)
```

</div>

    ## class: XqtlExperiment 
    ## dim: 138 731 
    ## metadata(0):
    ## assays(1): logcounts
    ## rownames(138): ENSG00000141934 ENSG00000099812 ... ENSG00000196867
    ##   ENSG00000268107
    ## rowData names(6): gene_id gene_name ... symbol entrezid
    ## colnames(731): HG00096 HG00100 ... NA21129 NA21130
    ## colData names(13): SRA_accession internal_libraryID ...
    ##   RNAQubitTotalAmount_ng RIN
    ##   9594 genotype calls present.
    ##   use getCalls() to see them with addresses.

<div id="cb16" class="sourceCode">

``` r
colData(mxx)=NULL
rowRanges(mxx) = rowRanges(mxx)[,1:6]
mxx$continent = cd$continentalGroup
system.time(zzz <- bind_Zs(mxx))
```

</div>

    ## some variants have MAF > 0.5

    ##    user  system elapsed 
    ## 276.801   1.613  35.633

<div id="cb19" class="sourceCode">

``` r
viz_stats(zzz, midchop=5)
```

</div>

<div id="htmlwidget-e5c8c404fe174e4c81bd"
class="plotly html-widget html-fill-item"
style="width:700px;height:432.632880098888px;">

</div>

Then we add sex as well.

<div id="cb20" class="sourceCode">

``` r
data(mageSE_19)
sds = rowSds(assay(mageSE_19), na.rm=TRUE)
qq = quantile(sds, .9)
ok = which(sds > qq)
cd = colData(mageSE_19)
mins  = minorAlleleCounts(demo_vcf(), GRanges("19:1-3000000"))
mxx = XqtlExperiment(mageSE_19[ok,], mins)
mxx$continent = cd$continentalGroup
mxx$sex = cd$sex
system.time(zzz <- bind_Zs(mxx))
```

</div>

    ## some variants have MAF > 0.5

    ##     user   system  elapsed 
    ## 1421.558    3.766  180.295

<div id="cb23" class="sourceCode">

``` r
viz_stats(zzz, midchop=7)
```

</div>

<div id="htmlwidget-36aa3d2a04d42bbc2145"
class="plotly html-widget html-fill-item"
style="width:700px;height:432.632880098888px;">

</div>

</div>

</div>

</div>
