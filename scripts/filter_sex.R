library(data.table)

args <- commandArgs(trailingOnly = TRUE)
valid_file <- args[1]
sexcheck_file <- args[2]
output_file <- args[3]

valid <- fread(valid_file)

sexcheck <- fread(sexcheck_file)
valid_sexcheck <- sexcheck[FID %in% valid$FID]

output_data <- valid_sexcheck[STATUS == "OK", .(FID, IID)]

fwrite(output_data, file = output_file, sep = "\t", quote = FALSE, col.names = FALSE)
