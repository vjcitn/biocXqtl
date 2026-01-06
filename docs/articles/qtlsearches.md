<div id="main" class="col-md-9" role="main">

# eQTL searching with biocXqtl

<div class="section level2">

## Introduction

[eQTLs](https://en.wikipedia.org/wiki/Expression_quantitative_trait_loci)
are genetic loci exhibiting association with gene expression.

Consortia have produced catalogs of eQTLs based on various cohorts and
experimental procedures.

-   [EBI eQTL catalog](https://www.ebi.ac.uk/eqtl/)
-   [eQTLGen](https://www.eqtlgen.org/); a related project with
    single-cell experiments is [sc-eqtlGen](https://www.eqtlgen.org/sc/)

biocXqtl facilitates flexible analyses for genetics of variation in
molecular phenotypes. In this vignette we review analytic results for
455 GEUVADIS samples.

</div>

<div class="section level2">

## tensorQTL example

<div id="cb1" class="sourceCode">

``` r
library(biocXqtl)
```

</div>

    ## Warning: package 'GenomicRanges' was built under R version 4.5.2

    ## Warning: replacing previous import 'BiocGenerics::type' by 'arrow::type' when
    ## loading 'biocXqtl'

<div id="cb4" class="sourceCode">

``` r
library(arrow)
```

</div>

    ## Warning: package 'arrow' was built under R version 4.5.2

<div id="cb6" class="sourceCode">

``` r
library(dplyr)
library(DT)
tens = example_tensorQTL()
tens |> arrange(pval_nominal) |> head(10) |> 
    as.data.frame() |> datatable()
```

</div>

<div id="htmlwidget-ac96cb3ee4656e2e9ec3"
class="datatables html-widget html-fill-item"
style="width:100%;height:auto;">

</div>

</div>

</div>
