<div id="main" class="col-md-9" role="main">

# provide path to a demonstration VCF file

<div class="ref-description section level2">

provide path to a demonstration VCF file

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
demo_vcf()
```

</div>

</div>

<div class="section level2">

## Note

derived from output of \`aws s3 cp –no-sign-request
s3://1000genomes/release/20130502/ALL.chr19.phase3\_shapeit2\_mvncall\_integrated\_v5a.20130502.genotypes.vcf.gz
.\` 26 Dec 2025.

Variants were filtered by restricting to samples in \`mageSE\_19\` and
variants with at least 10 individuals in each genotype class (00, 01,
11).

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
v = demo_vcf()
VariantAnnotation::scanVcfHeader(v)
#> class: VCFHeader 
#> samples(731): HG00096 HG00100 ... NA21129 NA21130
#> meta(5): fileDate fileformat reference source contig
#> fixed(2): FILTER ALT
#> info(27): CIEND CIPOS ... EX_TARGET MULTI_ALLELIC
#> geno(1): GT
```

</div>

</div>

</div>
