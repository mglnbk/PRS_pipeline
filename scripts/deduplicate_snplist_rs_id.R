# Load necessary libraries
library(data.table)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# 读取输入文件
snplist <- fread(input_file, header = FALSE)

# 去除重复的 rs-ID
dedup_snplist <- unique(snplist)

# 写入去重后的 SNP 列表
fwrite(dedup_snplist, file = output_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
