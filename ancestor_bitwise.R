# date created: June 29, 2021
# date edited: June 30, 2021
# this script reapproaches generating a simulated ancestral distribution for L=3 from HapMap using the bitwise form

###############
# load packages
library(base)
###############

###############
## write a function that calculates the decimal value based on general form (x_s * 2^(L-s))
### parameter x is the binary value at a site on the sequence, parameter e is the exponent (L-s)
bitwise <- function(x, e){
  decimal_val <- x * 2^e
  return(decimal_val)
}
###############

###############
## write a function to create a vector of the decimal forms of the 8 sequences
### parameter b is a a vector of sequences length L such that each item in the vector is of the form "x1, x2, x3, ..."
bitwise_vector <- function(b, l){
  # initialize an empty vector to hold the decimal forms
  sequence_bitwise <- vector(length = length(b))
  # iterate over the 8 sequences
  for (s in 1:length(b)){
    # for each sequence initialize empty vector to hold the three terms
    bit_wise <- vector(length = l)
    # split sequence string
    sequence_split <- unlist(strsplit(b[s], ", "))
    # print(sequence_split)
    # iterate over 3 terms in the sequence
    for (i in 1:l){
      # print(sequence_split[i])
      x_s <- as.numeric(sequence_split[i])
      bit <- bitwise(x_s, (l-i))
      bit_wise[i] <- bit
    }
    decimal <- sum(bit_wise)
    sequence_bitwise[s] <- decimal 
  }
  return(sequence_bitwise)}
###############


###############
# write a function that calculates the sequence proportions / simulated ancestor distribution
## parameter possible is a vector of the decimal values for the possible ancestor sequences
## parameter pool is a vector of the decimal values for the ancestor pool generated from HapMap
sim_ancestor <- function(possible, pool){
  #initialize vector to hold counts
  sequencecounts <- vector(length = length(possible))
  # initialize all counts at 0
  for (v in 1:length(possible)){
    sequencecounts[v] <- 0
  }
  # iterate over each decimal in the pool
  for (x in 1:length(pool)){
    # for each decimal in the pool, check if it matches each possible and increment counts
    for (s in 1:length(possible)){
      if (pool[x] == possible[s]){
        sequencecounts[s] <- sequencecounts[s] + 1
      }
    }
  }
  # print(sequencecounts)
  # initialize empty vector to hold the proportions
  sequenceprop <- vector(length = length(possible))
  # calculate each proportion as the number of obs in pool with that decimal over the total number of observations
  for (g in 1:length(possible)){
    sequenceprop[g] <- sequencecounts[g] / length(pool)
  }
  # return the vector of proportions
  return(sequenceprop)
}
###############

###############
# call bitwise_vector function on possible ancestor sequences and pool of sequences from HapMap
# possible_ancest_bit <- bitwise_vector(sequences, 3)
# pool_ancest_bit <- bitwise_vector(haps_length3, 3)
###############

###############
# calculate the ancestor distribution 
# ancestor_3 <- sim_ancestor(possible_ancest_bit, pool_ancest_bit)
###############
