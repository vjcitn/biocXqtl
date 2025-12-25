<div id="main" class="col-md-9" role="main">

# Fit multiple linear models with common design matrix

<div class="ref-description section level2">

Fit multiple linear models with common design matrix

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
fastLmMany_R(X, Y, intercept = TRUE)
```

</div>

</div>

<div class="section level2">

## Arguments

-   X:

    Design matrix (n x p)

-   Y:

    Response matrix (n x m), each column is a separate response

-   intercept:

    Logical, add intercept column to X?

</div>

<div class="section level2">

## Value

List with coefficients and standard errors matrices

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
set.seed(123)  # thanks claude
n <- 50
p <- 3
m <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- matrix(rnorm(n * m), n, m)
result <- fastLmMany_R(X, Y)
str(result)
#> List of 3
#>  $ coefficients: num [1:4, 1:5] 0.061 -0.12816 -0.17016 -0.0281 0.00798 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "(Intercept)" "X1" "X2" "X3"
#>   .. ..$ : chr [1:5] "Y1" "Y2" "Y3" "Y4" ...
#>  $ se          : num [1:4, 1:5] 0.139 0.145 0.15 0.138 0.139 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "(Intercept)" "X1" "X2" "X3"
#>   .. ..$ : chr [1:5] "Y1" "Y2" "Y3" "Y4" ...
#>  $ rank        : int 4
```

</div>

</div>

</div>
