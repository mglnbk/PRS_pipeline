library(data.table)
library(magrittr)
library(argparse)

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("--input", required = TRUE)
parser$add_argument("--output", required = TRUE)
args <- parser$parse_args()


an_gsr <- fread(args$input)

######## QC of base data ########
# Heritability check
# LD score regression and SumHer
# default

# Genome Build
# GRCh37 hg19

# Standard GWAS QC
# AN
#[1] "CHROM"    "POS"      "ID"       "REF"      "ALT"      "BETA"
# [7] "SE"       "PVAL"     "NGT"      "IMPINFO"  "NEFFDIV2" "NCAS"
#[13] "NCON"     "DIRE"
# IMPINFO > .8 and MAF > .1
result <- an_gsr[IMPINFO > .8]

# Unknown SNP
unknown_snp <- result[!grep("^rs", ID)][, ID]
res <- result[!ID %in% unknown_snp]

# Mismatching SNP
# exclude duplication
duplicated_snps <- res[, .N, by = ID][N > 2, ID]
res <- res[!ID %in% duplicated_snps]

# Ambiguous SNPs removal
# Only A1 A,G, C, T; A2 T,C,G,A
res <- res[!((REF == "T" & ALT == "A") |
               (REF == "A" & ALT == "T") |
               (REF == "C" & ALT == "G") |
               (REF == "G" & ALT == "C"))]
data.table::fwrite(res, args$output, sep = "\t")
