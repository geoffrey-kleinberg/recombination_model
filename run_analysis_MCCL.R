#######################
# Author: G Rhodes    #
# Date: Mar. 20, 2023 #
#######################

# This script is wrapped code for analyzing data 
# with the Recombination Model

###############################
# SECTION 1                   #
# load functions from scripts #
###############################

# functions needed to estimate ancestor
source("ancestor_bitwise.R")
source("ancestor_general.R")
# functions needed to generate samples of descendants
source("descendant_sequences.R")
# functions needed for m=0 marginal estimation
source("onewise_marginal.R")
# functions needed for m=1 marginal estimation
source("pairwise_marginal.R")
# functions needed for m=2 marginal estimation
source("threewise_marginal.R")
# functions needed for m=3 marginal estimation
source("fourwise_marginal_corrected.R")
# functions needed for joint estimation
source("ancestor_estimation.R")

###############################
# SECTION 2                   #
# Data wrangling into binary  #
###############################

### SET DATA TO BE USED HERE ###
#SNP_data <- read.csv2(file = "INSERT PATH HERE", header = TRUE, sep = "")

## select everything except for the SNP rsIDs and their locations
# snps <- SNP_data[,-c(1,2)]
## select the SNP rsIDs and their locations
# snps_id_loc <- SNP_data[,c(1,2)]

# check that each row has at most 2 letters, not less than 1 letter, and no NA's
# unique_letter <- apply(snps, 1, unique)
# num_letters <- sapply(unique_letter, length)
# print("Num rows w/ invalid number of unique letters:")
# sum(num_letters >2 & num_letters < 1)
# print("Num rows w/ NA unique letters:")
# sum(is.na(num_letters))
# print("NOTE both these should be zero!")

## Outline for following code -->
# Step one: count the number of each allele for each SNP
# Step two: identify the major allele
# Step three: check if each allele for the SNP is that allele; if yes, code 0 and code 1 otherwise
# Step four: repeat steps one through three for each row of the dataset and store into new dataset

# Step one/two
## This is a self-written function that takes a SNP row and returns a string of the major allele in that row
### parameter s passed to this function is a SNP, a row of hapmap data
find_major_allele <- function(s){
  ## create a vector of the alleles in this row
  alleles_in_snp <- unique(s)
  ## initialize a vector where the length is the number of alleles in the SNP
  num_alleles <- rep(NA, length(alleles_in_snp))
  
  ## calculate the number of each allele in the SNP
  for (j in 1:length(num_alleles)){
    num_alleles[j] <- sum(s %in% alleles_in_snp[j])
  }
  
  ## pick the major allele based on counts
  major_allele <- alleles_in_snp[which.max(num_alleles)]
  return (major_allele)
}

# Step three
## This is a self written function that takes a SNP row and a string of the major allele in that row and returns the row converted
##  to a binary sequence where the major allele is 0 and the minor allele is 1.
### The parameter snp_row is a SNP, a row of hapmap_data
### The parameter allele_m is a string of one letter, ATCG, indicating the major allele for that SNP
binary_seq <- function(snp_row, allele_m){
  seq <- rep(NA, length = length(snp_row))
  for (i in 1:length(snp_row)){
    if (snp_row[i] == allele_m){
      seq[i] <- 0
    }else{
      seq[i] <- 1
    }
  }
  return (seq)
}

# Step four
## Write a for loop iterating over all SNPs in the dataset and calling both functions on each and saving the generated binary sequence
##  rows to a new dataset

# initialize new dataset as an empty matrix
# SNP_dat_binary <- matrix(NA, nrow = dim(snps)[1], ncol = dim(snps)[2])

# find the major allele for each row
# alleles_major <- apply(snps, 1, find_major_allele)

# loop over each SNP
# for (i in 1:dim(snps)[1]){
  # generate the binary sequence for that SNP
#   binary <- binary_seq(snps[i,], alleles_major[i])
  # attach the recoded row into the initialized dataset
#   SNP_dat_binary[i,] <- binary}

# use the rsIDs from original dataset as the row names
# rownames(SNP_dat_binary) <- snps_id_loc[,1]
# use the genotype IDs from the original dataset as the column names
# colnames(SNP_dat_binary) <- colnames(snps)

# take transpose so each column is a SNP
# SNP_data_final <- t(SNP_dat_binary)
######################

###############################
# SECTION 3                   #
# Run analysis                #
###############################

## SET q-VALUE HERE ##
# q_val <- 0.05

## SET SEQUENCE LENGTH HERE
# Len <- 2

## SET SAMPLE SIZE
# nn <- 10

# SET SNP_data_final
# SNP_data_final <- data.frame(c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1), c(0, 0, 0, 0, 1, 1, 1, 0, 0, 1))

# marginal estimates, m=0
# m0_marg_est <- estimates_m0(descendents = SNP_data_final, L = Len, n = nn)

# marginal estimates, m=1
# m1_marg_est <- estimates_m1(L = Len, q = q_val, n = nn, d = SNP_data_final)

# marginal estimates, m=2
# m2_marg_est <- estimates_m2(L = Len, q = q_val, n = nn, descend = SNP_data_final)

# marginal estimates, m=3
# m3_marg_est <- estimates_m3(L = Len, q = q_val, n = nn, d = SNP_data_final)

# joint estimates from m=1 (PAIRWISE)
# m1_joint_est <- ancestor_pair_estimation(L = Len, m=1, pairs_est = m1_marg_est, ones_est = m0_marg_est)

# joint estimates from m=2 (THREEWISE)
# m2_joint_est <- ancestor_three(L = Len, m=2, threes_est = m2_marg_est, pairs_est = m1_marg_est)

# joint estimates from m=3 (FOURWISE)
# m3_joint_est <- ancestor_four(L = Len, m=3, four_est = m3_marg_est, three_est = m2_marg_est)
