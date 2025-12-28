

#' Filter VCF to variants with minimum counts in each genotype group
#' @import VariantAnnotation
#' @importFrom Rsamtools TabixFile
#'
#' @param vcf_file Path to tabix-indexed VCF file (.vcf.gz with .tbi index)
#' @param samples Character vector of sample IDs to retain
#' @param min_per_genotype Minimum number of individuals required in each 
#'        genotype group (homozygous reference, heterozygous, homozygous alternate)
#' @param genome Genome identifier (e.g., "hg19", "hg38")
#' @param regions Obligatory GRanges object to restrict to specific genomic regions
#' @param verbose Logical, print progress messages
#'
#' @return A filtered VCF object
#'
#' @examples
#' filtvcf <- filterVcfByGenotypeCounts(
#'   vcf_file = demo_vcf(),
#'   samples = c("HG00096", "HG00100", "HG00105", "HG00108", "HG00110", "HG00113", 
#'      "HG00114", "HG00121", "HG00122", "HG00127", "HG00130", "HG00132", 
#'      "HG00141", "HG00142", "HG00146", "HG00148", "HG00151", "HG00160", 
#'      "HG00177", "HG00179"),
#'   min_per_genotype = 2,
#'   regions = GenomicRanges::GRanges("19:1-50000000"),
#'   genome = "hg38"
#' )
#' filtvcf
#' @export
filterVcfByGenotypeCounts <- function(
    vcf_file, 
    samples = NULL,
    min_per_genotype = 1,
    genome = "hg38",
    regions = NULL,
    verbose = TRUE
) {
    if (missing(regions)) stop("regions must be supplied as a GRanges instance.")
    
    # Validate inputs
    if (!file.exists(vcf_file)) {
        stop("VCF file not found: ", vcf_file)
    }
    
    tbi_file <- paste0(vcf_file, ".tbi")
    if (!file.exists(tbi_file)) {
        stop("Tabix index not found. Run: indexVcf('", vcf_file, "')")
    }
    
    # Create TabixFile
    tbx <- Rsamtools::TabixFile(vcf_file)
    
    # Check available samples
    hdr <- scanVcfHeader(tbx)
    available_samples <- samples(hdr)
    
    if (verbose) {
        message("Total samples in VCF: ", length(available_samples))
    }
    
    # Validate requested samples
    if (!is.null(samples)) {
        missing_samples <- setdiff(samples, available_samples)
        if (length(missing_samples) > 0) {
            stop("Samples not found in VCF: ", 
                 paste(missing_samples, collapse=", "))
        }
        selected_samples <- samples
    } else {
        selected_samples <- available_samples
    }
    
    if (verbose) {
        message("Analyzing ", length(selected_samples), " samples")
    }
    
    # Create ScanVcfParam
    param <- ScanVcfParam(
        samples = selected_samples,
        geno = "GT",
        which = regions
    )
    
    # Read VCF   # FIXME -- could blow up ... may want to introduce yield discipline; regions control may be good enough
    if (verbose) message("Reading VCF...")
    vcf <- readVcf(tbx, genome = genome, param = param)
    
    if (verbose) {
        message("Read ", nrow(vcf), " variants")
    }
    
    # Extract genotypes
    gt <- geno(vcf)$GT
    
    # Count genotypes for each variant
    if (verbose) message("Counting genotypes...")
    genotype_counts <- countGenotypes(gt)
    
    # Filter variants
    keep <- (genotype_counts$n_hom_ref >= min_per_genotype) &
            (genotype_counts$n_het >= min_per_genotype) &
            (genotype_counts$n_hom_alt >= min_per_genotype)
    
    if (verbose) {
        message("Variants passing filter: ", sum(keep), " / ", length(keep))
        message("  Mean counts - HomRef: ", round(mean(genotype_counts$n_hom_ref), 1),
                ", Het: ", round(mean(genotype_counts$n_het), 1),
                ", HomAlt: ", round(mean(genotype_counts$n_hom_alt), 1))
    }
    
    # Subset VCF
    vcf_filtered <- vcf[keep, ]
    
    return(vcf_filtered)
}

#' Count genotypes across samples for each variant
#'
#' @param gt matrix of genotype calls (variants x samples) with form "0|0" or "0/0",
#' as produced by VariantAnnotion::geno applied to a VCF object
#' @return DataFrame with counts for each genotype category
#' @export
countGenotypes <- function(gt) {
    
    # Initialize counts
    n_variants <- nrow(gt)
    n_hom_ref <- integer(n_variants)
    n_het <- integer(n_variants)
    n_hom_alt <- integer(n_variants)
    n_missing <- integer(n_variants)
    
    # Count for each variant
    for (i in seq_len(n_variants)) {
        gts <- gt[i, ]
        
        # Count each genotype type
        # 0/0 or 0|0 = homozygous reference
        n_hom_ref[i] <- sum(gts %in% c("0/0", "0|0"), na.rm = TRUE)
        
        # 0/1, 1/0, 0|1, 1|0 = heterozygous
        n_het[i] <- sum(gts %in% c("0/1", "1/0", "0|1", "1|0"), na.rm = TRUE)
        
        # 1/1 or 1|1 = homozygous alternate
        n_hom_alt[i] <- sum(gts %in% c("1/1", "1|1"), na.rm = TRUE)
        
        # Missing
        n_missing[i] <- sum(is.na(gts) | gts == "./." | gts == ".|.")
    }
    
    DataFrame(
        n_hom_ref = n_hom_ref,
        n_het = n_het,
        n_hom_alt = n_hom_alt,
        n_missing = n_missing,
        n_called = n_hom_ref + n_het + n_hom_alt
    )
}

#' More efficient vectorized version for large datasets
#' @param gt matrix of genotype calls (variants x samples) with form "0|0" or "0/0"
#' as produced by VariantAnnotion::geno applied to a VCF object
#' @return DataFrame with counts for each genotype category
#' @export
countGenotypesVectorized <- function(gt) {
    
    # Convert to character matrix if needed
    if (is(gt, "matrix")) {
        gt <- as.matrix(gt)
    }
    
    # Count genotypes using apply
    counts <- t(apply(gt, 1, function(row) {
        c(
            hom_ref = sum(row %in% c("0/0", "0|0"), na.rm = TRUE),
            het = sum(row %in% c("0/1", "1/0", "0|1", "1|0"), na.rm = TRUE),
            hom_alt = sum(row %in% c("1/1", "1|1"), na.rm = TRUE),
            missing = sum(is.na(row) | row %in% c("./.", ".|."))
        )
    }))
    
    DataFrame(
        n_hom_ref = counts[, "hom_ref"],
        n_het = counts[, "het"],
        n_hom_alt = counts[, "hom_alt"],
        n_missing = counts[, "missing"],
        n_called = counts[, "hom_ref"] + counts[, "het"] + counts[, "hom_alt"]
    )
}

## Example usage:  (Most code produced by Claude, 26 Dec 2025)

## Basic usage - filter to variants with at least 5 individuals in each group
# vcf_filtered <- filterVcfByGenotypeCounts(
#     vcf_file = "mydata.vcf.gz",
#     samples = c("sample1", "sample2", "sample3", "sample4", "sample5"),
#     min_per_genotype = 5,
#     genome = "hg38"
# )

## With specific genomic regions
# library(GenomicRanges)
# regions <- GRanges("chr1", IRanges(start = 1000000, end = 2000000))
# vcf_filtered <- filterVcfByGenotypeCounts(
#     vcf_file = "mydata.vcf.gz",
#     samples = NULL,  # use all samples
#     min_per_genotype = 3,
#     genome = "hg38",
#     regions = regions
# )

## Access the filtered data
# rowRanges(vcf_filtered)  # Genomic positions
# geno(vcf_filtered)$GT    # Genotypes
# ref(vcf_filtered)        # Reference alleles
# alt(vcf_filtered)        # Alternate alleles


#' Get genotype summary statistics
#' 
#' @param vcf A VCF object
#' @return DataFrame with genotype statistics per variant
getGenotypeSummary <- function(vcf) {
    gt <- geno(vcf)$GT
    counts <- countGenotypesVectorized(gt)
    
    # Add MAF (minor allele frequency)
    n_alleles <- 2 * counts$n_called
    n_alt_alleles <- 2 * counts$n_hom_alt + counts$n_het
    maf <- pmin(n_alt_alleles / n_alleles, 
                1 - n_alt_alleles / n_alleles)
    
    cbind(counts, 
          DataFrame(
              maf = maf,
              call_rate = counts$n_called / (counts$n_called + counts$n_missing)
          ))
}

#' Filter VCF with multiple criteria
#'
#' @param vcf VCF object
#' @param min_per_genotype Minimum samples per genotype group
#' @param min_maf Minimum minor allele frequency
#' @param min_call_rate Minimum call rate
#' @return Filtered VCF object
filterVcfMultipleCriteria <- function(
    vcf,
    min_per_genotype = NULL,
    min_maf = NULL,
    min_call_rate = NULL
) {
    
    stats <- getGenotypeSummary(vcf)
    keep <- rep(TRUE, nrow(vcf))
    
    if (!is.null(min_per_genotype)) {
        keep <- keep & 
                (stats$n_hom_ref >= min_per_genotype) &
                (stats$n_het >= min_per_genotype) &
                (stats$n_hom_alt >= min_per_genotype)
    }
    
    if (!is.null(min_maf)) {
        keep <- keep & (stats$maf >= min_maf)
    }
    
    if (!is.null(min_call_rate)) {
        keep <- keep & (stats$call_rate >= min_call_rate)
    }
    
    message("Keeping ", sum(keep), " / ", length(keep), " variants")
    
    vcf[keep, ]
}
