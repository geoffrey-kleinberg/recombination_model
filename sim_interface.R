source('simulator.R')

args <- commandArgs(trailingOnly = TRUE)

len <- as.numeric(args[1])
q_val <- as.numeric(args[2])
m <- as.numeric(args[3])
n <- as.numeric(args[4])
file_name <- args[5]

get_and_write_estimate(len, q_val, m, n, file_name)