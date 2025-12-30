<div id="main" class="col-md-9" role="main">

# internal function to combine new xQTL test results with existing XqtlExperiment instance

<div class="ref-description section level2">

internal function to combine new xQTL test results with existing
XqtlExperiment instance

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
add_scores(xse, mat, check.cols = TRUE)
```

</div>

</div>

<div class="section level2">

## Arguments

-   xse:

    XqtlExperiment instance

-   mat:

    a matrix of test results conformant with the XqtlExperiment:
    rownames agree,

-   check.cols:

    logical(1) will message if redundant variants are found when new
    scores are bound in columns will be added to the
    mcols(rowRanges(xse))

</div>

</div>
