<div id="main" class="col-md-9" role="main">

# produce a plotly display of statistics of xQTL association

<div class="ref-description section level2">

produce a plotly display of statistics of xQTL association

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
viz_stats(
  se,
  jitter_fac = 500,
  ptcolor = "blue",
  midchop = 2,
  xlabel = "SNP addr",
  ylabel = "xQTL association Z"
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   se:

    SummarizedExperiment post "bind\_Zs"

-   jitter\_fac:

    plot employs jittering of SNP address by this factor

-   ptcolor:

    point color

-   midchop:

    stats with absolute value less than this value are omitted from
    display

-   xlabel:

    axis label

-   ylabel:

    axis label

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
data(geuv19)
sds = rowSds(assay(geuv19), na.rm=TRUE)
qq = quantile(sds, .8)
ok = which(sds > qq)
lk = geuv19[ok,]
mafs = maf(colData(lk)) # only snps here
mins = apply(data.matrix(as.data.frame(colData(lk))), 2, min, na.rm=TRUE) # some -1 values
colData(lk) = colData(lk)[,which(mafs>.25 & mins > -1)]
lk <- bind_Zs(lk, colselector = function(se) colnames(colData(se)))
viz_stats(lk)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rowData': object 'lk' not found
```

</div>

</div>

</div>
