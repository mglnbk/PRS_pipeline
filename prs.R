library(data.table)
library(magrittr)

an_gsr <- fread("./pgcAN2.2019-07.vcf.tsv")
ptsd_gsr <- fread("./pts_eur_freeze2_overall.results")

colnames(an_gsr)
colnames(ptsd_gsr)

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
nrow(res)

# Ambiguous SNPs removal
# Only A1 A,G, C, T; A2 T,C,G,A
res <- res[!((REF == "T" & ALT == "A") |
             (REF == "A" & ALT == "T") |
             (REF == "C" & ALT == "G") |
             (REF == "G" & ALT == "C"))]
nrow(res)
QC_GSR_output_file <- "AN_GSR_BASE.QC.gz"
data.table::fwrite(res, QC_GSR_output_file, sep = "\t")

######## QC of Target Data #########
# standard GWAS QC
# removing SNPs with low genotyping rate, low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing individuals with low genotyping rate

# input file should include .bed file accompanied with .bim and .fam files
input_file <- "/proj/tengfei/CHRX/autosome/ukb_imp_allchr_v3_40k_fmri"
output_file <- "./res.QC"

cmd1 <- paste("/proj/htzhu/UKB_GWAS/phase1and2/plink --bfile", input_file, 
             "--geno 0.01",   # allowable largest proportion of NA genes 
             "--maf 0.01",    # freq of minor allele
             "--hwe 1e-6",    # Hardy-Weinberg p-vlaue
             "--mind 0.01",   # Excludes individuals who have a high rate of genotype missingness
             "--write-snplist",
             "--make-just-fam",
             "--out", output_file)

system2(cmd1, stdout = TRUE, stderr = TRUE)

# remove highly correlated SNPs
# this will generate two files, prune.in. and prune.out.
# prune.in contains all of SNPs r2<0.25
input_file <- "res.QC"
output_file <- "res.QC"
cmd2 <- paste("/proj/htzhu/UKB_GWAS/phase1and2/plink --bfile", input_file, 
             "--extract res.QC.snplist",   # use QCed snps
             "--keep res.QC.fam", # use QCed fam files
             "--indep-pairwise 200 50 0.25", # slide windows to filter those snp r2>0.25
             "--out", output_file)


# Calculate heterozygosity rate
input_file <- "res.QC"
output_file <- "res.QC"
cmd2 <- paste("/proj/htzhu/UKB_GWAS/phase1and2/plink --bfile", input_file, 
             "--extract res.QC.prune.in",   # output prune.in
             "--keep res.QC.fam",
             "--het",
             "--out", output_file)

# Read in file
dat <- fread("res.QC.het")
# Get samples with F coefficient within 3 SD of the population mean
valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 
# print FID and IID for valid samples
fwrite(valid[,c("FID","IID")], "res.valid.sample", sep="\t") 


# magrittr allow us to do piping, which help to reduce the 
# amount of intermediate data types
library(data.table)
library(magrittr)
# Read in bim file 
bim <- fread("EUR.bim") %>%
    # Note: . represents the output from previous step
    # The syntax here means, setnames of the data read from
    # the bim file, and replace the original column names by 
    # the new names
    setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
    # And immediately change the alleles to upper cases
    .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]
# Read in summary statistic data (require data.table v1.12.0+)
height <- fread("Height.QC.gz") %>%
    # And immediately change the alleles to upper cases
    .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
# Read in QCed SNPs
qc <- fread("EUR.QC.snplist", header=F)
