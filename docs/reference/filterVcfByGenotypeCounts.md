<div id="main" class="col-md-9" role="main">

# Filter VCF to variants with minimum counts in each genotype group

<div class="ref-description section level2">

Filter VCF to variants with minimum counts in each genotype group

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
filterVcfByGenotypeCounts(
  vcf_file,
  samples = NULL,
  min_per_genotype = 1,
  genome = "hg38",
  regions = NULL,
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   vcf\_file:

    Path to tabix-indexed VCF file (.vcf.gz with .tbi index)

-   samples:

    Character vector of sample IDs to retain

-   min\_per\_genotype:

    Minimum number of individuals required in each genotype group
    (homozygous reference, heterozygous, homozygous alternate)

-   genome:

    Genome identifier (e.g., "hg19", "hg38")

-   regions:

    Obligatory GRanges object to restrict to specific genomic regions

-   verbose:

    Logical, print progress messages

</div>

<div class="section level2">

## Value

A filtered VCF object

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
filtvcf <- filterVcfByGenotypeCounts(
  vcf_file = demo_vcf(),
  samples = c("HG00096", "HG00100", "HG00105", "HG00108", "HG00110", "HG00113", 
     "HG00114", "HG00121", "HG00122", "HG00127", "HG00130", "HG00132", 
     "HG00141", "HG00142", "HG00146", "HG00148", "HG00151", "HG00160", 
     "HG00177", "HG00179"),
  min_per_genotype = 2,
  regions = GenomicRanges::GRanges("19:1-50000000"),
  genome = "hg38"
)
#> Total samples in VCF: 731
#> Analyzing 20 samples
#> Reading VCF...
#> Read 8465 variants
#> Counting genotypes...
#> Variants passing filter: 2922 / 8465
#>   Mean counts - HomRef: 10.8, Het: 6, HomAlt: 3.1
filtvcf
#> class: CollapsedVCF 
#> dim: 2922 20 
#> rowRanges(vcf):
#>   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
#> info(vcf):
#>   DataFrame with 27 columns: CIEND, CIPOS, CS, END, IMPRECISE, MC, MEINFO, M...
#> info(header(vcf)):
#>                  Number Type    Description                                    
#>    CIEND         2      Integer Confidence interval around END for imprecise...
#>    CIPOS         2      Integer Confidence interval around POS for imprecise...
#>    CS            1      String  Source call set.                               
#>    END           1      Integer End coordinate of this variant                 
#>    IMPRECISE     0      Flag    Imprecise structural variation                 
#>    MC            .      String  Merged calls.                                  
#>    MEINFO        4      String  Mobile element info of the form NAME,START,E...
#>    MEND          1      Integer Mitochondrial end coordinate of inserted seq...
#>    MLEN          1      Integer Estimated length of mitochondrial insert       
#>    MSTART        1      Integer Mitochondrial start coordinate of inserted s...
#>    SVLEN         .      Integer SV length. It is only calculated for structu...
#>    SVTYPE        1      String  Type of structural variant                     
#>    TSD           1      String  Precise Target Site Duplication for bases, i...
#>    AC            A      Integer Total number of alternate alleles in called ...
#>    AF            A      Float   Estimated allele frequency in the range (0,1)  
#>    NS            1      Integer Number of samples with data                    
#>    AN            1      Integer Total number of alleles in called genotypes    
#>    EAS_AF        A      Float   Allele frequency in the EAS populations calc...
#>    EUR_AF        A      Float   Allele frequency in the EUR populations calc...
#>    AFR_AF        A      Float   Allele frequency in the AFR populations calc...
#>    AMR_AF        A      Float   Allele frequency in the AMR populations calc...
#>    SAS_AF        A      Float   Allele frequency in the SAS populations calc...
#>    DP            1      Integer Total read depth; only low coverage data wer...
#>    AA            1      String  Ancestral Allele. Format: AA|REF|ALT|IndelTy...
#>    VT            .      String  indicates what type of variant the line repr...
#>    EX_TARGET     0      Flag    indicates whether a variant is within the ex...
#>    MULTI_ALLELIC 0      Flag    indicates whether a site is multi-allelic      
#> geno(vcf):
#>   List of length 1: GT
#> geno(header(vcf)):
#>       Number Type   Description
#>    GT 1      String Genotype   
```

</div>

</div>

</div>
