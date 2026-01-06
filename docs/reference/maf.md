<div id="main" class="col-md-9" role="main">

# compute putative minor allele frequency for XqtlExperiment

<div class="ref-description section level2">

compute putative minor allele frequency for XqtlExperiment

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
maf(xse)
```

</div>

</div>

<div class="section level2">

## Arguments

-   xse:

    XqtlExperiment instance

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
example(XqtlExperiment)  # makes nn
#> Warning: no help found for 'XqtlExperiment'
summary(maf(nn))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'summary': error in evaluating the argument 'x' in selecting a method for function 'mcols': object 'nn' not found
```

</div>

</div>

</div>
