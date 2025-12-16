<div id="main" class="col-md-9" role="main">

# create a matrix of Z statistics for an additive effect of rare allele count in calls on molecular phenotype in molec

<div class="ref-description section level2">

create a matrix of Z statistics for an additive effect of rare allele
count in calls on molecular phenotype in molec

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
getzs(molec, calls)
```

</div>

</div>

<div class="section level2">

## Arguments

-   molec:

    numeric matrix of molecular phenotype values, columns are samples

-   calls:

    numeric matrix of rare allele counts, rows are samples, columns are
    loci

</div>

<div class="section level2">

## Value

A matrix with one row per

</div>

<div class="section level2">

## Note

For a SNP with MAF 0, NA is returned.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
data(geuv19)
mol = assay(geuv19)
calls = data.matrix(as.data.frame(colData(geuv19)))
m = maf(calls)
allz = getzs(mol[1:100,], calls[,m>.3])
summary(as.numeric(allz))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#> -4.4767 -0.5831  0.1399  0.1047  0.7790  4.3300    1988 
```

</div>

</div>

</div>
