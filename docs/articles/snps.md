<div id="main" class="col-md-9" role="main">

# Handling variants in biocXqtl

<div class="section level2">

## Introduction

Diverse strategies have been proposed for managing archives of DNA
variants in humans.

-   [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) is widely
    used and many toolsets are available for working with this format.  
-   [Plink](https://www.cog-genomics.org/plink/2.0/) has a format in
    wide use.
-   Hail’s MatrixTable and associated tooling offer [“easy ways to
    slice, dice, query, and
    summarize”](https://hail.is/docs/0.2/tutorials/01-genome-wide-association-study.html)
    DNA variants.

We’ll focus on VCF and Plink in this document.

</div>

<div class="section level2">

## VCF access and manipulations

We will work with bgzipped, tabix-indexed resources.

A demonstration VCF archive is available with biocXqtl. The associated
vcf.gz file was formed by filtering publicly available variants from the
1000 genomes project. See the Appendix for information on how the
demonstration VCF was produced. A key filtering step involved limiting
inclusion to those variants for which at least 10 individuals were
present in each of the classes: homozygous rare, heterozygous,
homozygous common.

Setup:

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
library(VariantAnnotation)
vpath = demo_vcf()
```

</div>

First we examine the header.

<div id="cb5" class="sourceCode">

``` r
h = scanVcfHeader(vpath)
h
```

</div>

    ## class: VCFHeader 
    ## samples(731): HG00096 HG00100 ... NA21129 NA21130
    ## meta(5): fileDate fileformat reference source contig
    ## fixed(2): FILTER ALT
    ## info(27): CIEND CIPOS ... EX_TARGET MULTI_ALLELIC
    ## geno(1): GT

The header report indicates that genotypes are available for 731
individuals.

Direct ingestion of VCF can consume considerable RAM. We can select
regions of interest and scan with inline filtering. See the
documentation on `ScanVcfParam` in VariantAnnotation.

<div id="cb7" class="sourceCode">

``` r
library(GenomicRanges)
p = GRanges("19:1-3000000")
ss = ScanVcfParam(which=p)
sc1 = scanVcf(vpath, param=ss)
```

</div>

Each range in the `ScanVcfParam` will produce a list element in the
result of `scanVcf`.

<div id="cb8" class="sourceCode">

``` r
names(sc1[[1]])
```

</div>

    ## [1] "rowRanges" "REF"       "ALT"       "QUAL"      "FILTER"    "INFO"     
    ## [7] "GENO"

<div id="cb10" class="sourceCode">

``` r
dim(sc1[[1]]$GENO[[1]])
```

</div>

    ## [1] 9594  731

<div id="cb12" class="sourceCode">

``` r
sc1[[1]]$rowRanges
```

</div>

    ## GRanges object with 9594 ranges and 0 metadata columns:
    ##               seqnames          ranges strand
    ##                  <Rle>       <IRanges>  <Rle>
    ##   rs567986644       19          160565      *
    ##   rs548347106       19   183097-183102      *
    ##   rs527871033       19   212785-212786      *
    ##   rs533243783       19          219620      *
    ##    rs62102979       19          226776      *
    ##           ...      ...             ...    ...
    ##    rs12982139       19         2999553      *
    ##   rs199521660       19         2999724      *
    ##   rs200315688       19 2999734-2999738      *
    ##    rs13345197       19         2999759      *
    ##    rs62125443       19         2999804      *
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

For purposes of statistical modeling, we have a helper function that
retrieves minor allele counts from a tabix-indexed VCF. This returns a
GRanges with counts in the mcols. The addresses are bound with the
counts.

<div id="cb14" class="sourceCode">

``` r
cnts = minorAlleleCounts(vpath, region=GRanges("19:1-3000000"))
cnts[,1:2]
```

</div>

    ## GRanges object with 9594 ranges and 2 metadata columns:
    ##               seqnames          ranges strand |   HG00096   HG00100
    ##                  <Rle>       <IRanges>  <Rle> | <numeric> <numeric>
    ##   rs567986644       19          160565      * |         0         1
    ##   rs548347106       19   183097-183102      * |         1         1
    ##   rs527871033       19   212785-212786      * |         2         1
    ##   rs533243783       19          219620      * |         1         1
    ##    rs62102979       19          226776      * |         1         1
    ##           ...      ...             ...    ... .       ...       ...
    ##    rs12982139       19         2999553      * |         1         0
    ##   rs199521660       19         2999724      * |         0         1
    ##   rs200315688       19 2999734-2999738      * |         1         0
    ##    rs13345197       19         2999759      * |         1         0
    ##    rs62125443       19         2999804      * |         0         2
    ##   -------
    ##   seqinfo: 1 sequence from GRCh38 genome; no seqlengths

</div>

<div class="section level2">

## PLINK

<div class="section level3">

### Example data

The data were extracted from the GWAS tutorial.

<div id="cb16" class="sourceCode">

``` r
zfi = system.file("plink", "chr19_1kglim.zip", package="biocXqtl")
tdir = tempdir()
unzip(zfi, exdir=tdir)
dir(tdir)
```

</div>

    ##  [1] "BiocStyle"                      "chr19_limsamples_maf10.bed"    
    ##  [3] "chr19_limsamples_maf10.bim"     "chr19_limsamples_maf10.fam"    
    ##  [5] "chr19_limsamples_maf10.log"     "file91142cab4512"              
    ##  [7] "file9114377a83a4"               "file91143dce457d"              
    ##  [9] "file9114462d36e0"               "file91144e8d4038"              
    ## [11] "make_bed.sh"                    "rmarkdown-str91142db56957.html"

</div>

<div class="section level3">

### Using genio

This package seems to be capable only of ‘full reads’.

<div id="cb18" class="sourceCode">

``` r
library(genio)
fi = grep("bed$", dir(tdir, full=TRUE), value=TRUE)
dem = read_plink(fi)
dim(dem$X)
```

</div>

    ## [1] 163640    104

</div>

<div class="section level3">

### Using python `bed-reader`

Here we can make targeted reads from the resource.

<div id="cb20" class="sourceCode">

``` r
library(reticulate)
```

</div>

    ## Warning: package 'reticulate' was built under R version 4.5.2

<div id="cb22" class="sourceCode">

``` r
py_require("bed-reader")
br = import("bed_reader")
res = br$open_bed(fi)
toget = tuple(0:4, 0:5)
limread = res$read(toget)
dim(limread)
```

</div>

    ## [1] 5 6

<div id="cb24" class="sourceCode">

``` r
head(res$bp_position)
```

</div>

    ## [1] 80840 81039 81806 90854 93234 95981

<div id="cb26" class="sourceCode">

``` r
head(res$chromosome)
```

</div>

    ## [1] "19" "19" "19" "19" "19" "19"

<div id="cb28" class="sourceCode">

``` r
head(res$sid)
```

</div>

    ## [1] "rs201816663" "rs879919317" "rs2432259"   "rs200254023" "rs151275984"
    ## [6] "rs8110113"

<div id="cb30" class="sourceCode">

``` r
head(res$iid)
```

</div>

    ## [1] "HG00096" "HG00097" "HG00099" "HG00100" "HG00101" "HG00171"

</div>

</div>

<div class="section level2">

## Appendix

<div class="section level3">

### VCF production

Code to produce the demonstration VCF is

    library(biocXqtl)
    library(VariantAnnotation)
    data(mageSE_19) # produced from MAGE "DESeq2" quantifications
    cc = colnames(mageSE_19)
    data(tiles38)
    library(GenomeInfoDb)
    seqlevelsStyle(tiles38) = "Ensembl"
    t19 = tiles38[which(seqnames(tiles38)=="19")]
    bigvcf = filterVcfByGenotypeCounts("~/MAGE/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", samples=cc, min_per_genotype=10, regions=range(t19))

The vcf.gz file was obtained following instructions at
`https://registry.opendata.aws/1000-genomes/`.

The code given above runs on a 32 GB macbook in a few minutes, with the
following report:

    Total samples in VCF: 2504
    Analyzing 731 samples
    Reading VCF...
    [W::hts_idx_load3] The index file is older than the data file: 
     .../MAGE/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
    Read 1816698 variants
    Counting genotypes...
    Variants passing filter: 156943 / 1816698
      Mean counts - HomRef: 690.8, Het: 25.4, HomAlt: 14.4

</div>

<div class="section level3">

### PLINK resource production

Folder `GWAS/project1/data/raw` in the tutorial download has the
relevant contents.

An extract from chr19 is produced using

    plink2 \
      --pfile all_hg38 vzs \
      --chr 19 \
      --maf 0.1 \
      --keep kp.txt \
      --make-pgen vzs \
      --allow-extra-chr \
      --out chr19_limsamples_pp1

`kp.txt` is in `inst/plink` and is a selection of samples with diverse
ancestries.

PLINK bed-format files are produced using

    plink2 \
      --pfile chr19_limsamples_pp1 vzs \
      --make-bed \
      --out chr19_limsamples_maf10

</div>

</div>

</div>
