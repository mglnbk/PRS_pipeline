library(data.table)
library(magrittr)
library(argparse)

# Set up argument parser
parser <- ArgumentParser(description='Process SNPs data.')

# Define arguments
parser$add_argument('--bim', required=TRUE, help='Path to the BIM file')
parser$add_argument('--target', required=TRUE, help='Path to the summary statistics file')
parser$add_argument('--snp', required=TRUE, help='Path to the QCed SNPs file')
parser$add_argument('--output_a1', required=TRUE, help='Path to the output A1 file')
parser$add_argument('--output_mismatch', required=TRUE, help='Path to the output mismatch file')

# Parse arguments
args <- parser$parse_args()
# Read in bim file
bim <- fread(args$bim) %>%
    # Set column names
    setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
    # Change alleles to upper cases
    .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]

# Read in summary statistic data
target <- fread(args$target) %>%
# Change alleles to upper cases
    .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]

# Read in QCed SNPs
qc <- fread(args$snp, header=FALSE)

# Merge summary statistic with target and filter out QCed SNPs
info <- merge(bim, target, by=c("SNP", "CHR", "BP")) %>%
    .[SNP %in% qc[, V1]]

# Function for calculating the complementary allele
complement <- function(x) {
    switch (x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}

# Get SNPs that have the same alleles across base and target
info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]

# Identify SNPs that are complementary between base and target
com.snps <- info[sapply(B.A1, complement) == A1 &
                 sapply(B.A2, complement) == A2, SNP]

# Update the bim file for complementary SNPs
bim[SNP %in% com.snps, c("B.A1", "B.A2") := list(sapply(B.A1, complement), sapply(B.A2, complement))]

# Identify SNPs that need recoding
recode.snps <- info[B.A1 == A2 & B.A2 == A1, SNP]
# Update the bim file for recoding SNPs
bim[SNP %in% recode.snps, c("B.A1", "B.A2") := list(B.A2, B.A1)]

# Identify SNPs that need recoding & complement
com.recode <- info[sapply(B.A1, complement) == A2 &
                   sapply(B.A2, complement) == A1, SNP]
# Update the bim file for recoding & complementary SNPs
bim[SNP %in% com.recode,
    c("B.A1", "B.A2") := list(sapply(B.A2, complement),
                              sapply(B.A1, complement))]

# Write the updated bim file
fwrite(bim[, c("SNP", "B.A1")],
       args$output_a1,
       col.names = FALSE,
       sep = "\t")

# Identify SNPs that have different allele in base and target
mismatch <- bim[!(SNP %in% info.match |
                    SNP %in% com.snps |
                    SNP %in% recode.snps |
                    SNP %in% com.recode), SNP]
write.table(mismatch,
            args$output_mismatch,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
