library(data.table)
library(magrittr)
library(argparse)

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("--input", required = TRUE)
parser$add_argument("--output", required = TRUE)
parser$add_argument("--keys", nargs = "+", required = TRUE)
parser$add_argument("--values", nargs = "+", required = TRUE)

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
#[1] "CHROM"    "POS"      "ID"       "A1"      "A2"      "BETA"
# [7] "SE"       "PVAL"     "NGT"      "IMPINFO"  "NEFFDIV2" "NCAS"
#[13] "NCON"     "DIRE"
# IMPINFO > .8 and MAF > .1
setnames(
  an_gsr, old = args$keys,
  new = args$values,
  skip_absent = TRUE
)
result <- an_gsr[INFO > .8]

# Unknown SNP
unknown_snp <- result[!grep("^rs", SNP)][, SNP]
res <- result[!SNP %in% unknown_snp]

# Mismatching SNP
# exclude duplication
duplicated_snps <- res[, .N, by = SNP][N > 2, SNP]
res <- res[!SNP %in% duplicated_snps]

# Ambiguous SNPs removal
# Only A1 A,G, C, T; A2 T,C,G,A
res <- res[!((A1 == "T" & A2 == "A") |
               (A1 == "A" & A2 == "T") |
               (A1 == "C" & A2 == "G") |
               (A1 == "G" & A2 == "C"))]
data.table::fwrite(res, args$output, sep = "\t")
