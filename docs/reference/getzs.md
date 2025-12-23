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
getzs(
  molec,
  calls,
  statfun = function(x, y) RcppEigen::fastLmPure(X = x, y = y)
)
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

-   statfun:

    a function with arguments x, y, with the intention that x is a
    design matrix lacking a column of 1s and y is a response vector with
    nrow(x) elements.

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
sds = rowSds(mol)
mol = mol[which(sds>0), ] # drop constant features
calls = data.matrix(as.data.frame(colData(geuv19)))
csds = apply(calls,2,sd, na.rm=TRUE)
mins = apply(calls,2,min, na.rm=TRUE)  # some snps include -1 values
m = maf(calls)
allz = getzs(mol[1:100,], calls[,m>.3 & csds>0 & mins > -1])
summary(as.numeric(allz))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -4.4767 -0.5657  0.1605  0.1210  0.7876  4.3300 
  # 90 sec?
 if (requireNamespace("MASS")) {
 nbz = function(x,y) { 
 # error messages will be thrown from glm ... maybe condition to allow warning/error
   f = tryCatch(MASS::glm.nb(y~x, data=data.frame(y=y, x=x[,2])), 
            error=function(e) return(list(coefficients=NA, se=NA) ))  # getzs adds column of 1s
   #if (inherits(f, "try-error")) return(list(coefficients=NA, se=NA))
   dat = summary(f)$coefficients
   list(coefficients=dat[,1], se=dat[,2])
 }  
 # for now NA is noisily returned
 allz2 = suppressWarnings(getzs(mol[1:100,], calls[,m>.3 & csds>0 & mins > -1], statfun = nbz))
 # do NB results differ substantially from OLS?
 plot(allz2[1,], allz[1,], xlab="NB", ylab="OLS", main=sprintf("txQTL Zs for %s", rownames(mol)[1]))
 }
#> Loading required namespace: MASS

```

</div>

</div>

</div>
