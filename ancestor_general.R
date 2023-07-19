# date created: June 30, 2021
# date edited: July 19, 2021

# this script calculates the ancestor distribution based upon the HapMap data as a pool of ancestor sequences for a sequence of length L
## see ancestor_L3_bitwise.R for functions bitwise, bitwise_vector, sim_ancestor

###############
# SET L HERE
# L <- 10 #or whatever sequence length
###############

ancestor <- function(L, hapmap_data){
  ###############
  # select a subset of the data for this sequence length
  hapmap_L <- hapmap_data[1:L,]
  ###############
  
  ###############
  # add a row to new dataset to store full sequence as one value --> string of form "x1, x2, ...., xL"
  
  # initialize empty row
  haplotypeL <- rep(NA, length = dim(hapmap_L)[2])
  # combine new row with dataset of first three SNPs
  hapmap_L <- rbind(hapmap_L, haplotypeL)
  
  # iterate over COLUMNS
  for (i in 1:dim(hapmap_L)[2]){
    hap <- vector(length = L)
    # iterate over ROWS
    for (j in 1:(dim(hapmap_L)[1]-1)){
      hap[j] <- hapmap_L[j,i]
      # print(hap)
      # print(unr_bin_L3[j,i])
    }
    full_hap <- toString(hap)
    hapmap_L[(L+1),i] <- full_hap}
  
  # vector to store just the strings
  pool_ancestor <- hapmap_L[(L+1),]
  ###############
  
  ###############
  # create a vector of the possible ancestor sequences in bitwise form
  ## it can be shown that the decimal values for sequences of length L are 0, 1, ... , (2^L - 1)
  
  possible_ancestor <- vector(length = (2^L))
  for (j in 0:((2^L)-1)){
    possible_ancestor[j+1] <- j
  }
  ###############
  
  ###############
  # convert pool of ancestor sequences to bitwise form
  ## use function bitwise_vector from other script (this involves a function call to bitwise)
  
  pool_ancestor_decimal <- bitwise_vector(pool_ancestor, L)
  ###############
  
  ###############
  # generate simulated ancestor distribution
  ## use function sim_ancestor from other script
  
  ancestor_L <- sim_ancestor(possible_ancestor, pool_ancestor_decimal)
  ###############
  return(as.matrix(ancestor_L))
}

# ancestor_5 <- ancestor(5, hapmap_binary)
