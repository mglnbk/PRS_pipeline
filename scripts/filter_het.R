library(data.table)
args <- commandArgs(trailingOnly = TRUE)

# 输入和输出文件路径
input_file <- args[1]
output_file <- args[2]

# Read in file
dat <- fread(input_file)
# Get samples with F coefficient within 3 SD of the population mean
valid <- dat[F <= mean(F) + 3 * sd(F) & F >= mean(F) - 3 * sd(F)]
# print FID and IID for valid samples
fwrite(valid[,c("FID","IID")], output_file, sep="\t") 
