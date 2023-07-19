# date created: June 30, 2021
# date edited: July 1, July 2

# this script generates descendant sequences from the simulated ancestor distribution


#######
## Pseudo-code -->
# Based on the probabilities stored in the simulated ancestral distribution, select an ancestor sequence
# Use a random function to decide whether recombination occurs between sites 1 and 2 (using probability q) --> model as a binomial
## If there is no recombination, keep the selected ancestral distribution and continue
## If there is recombination, keep the selected site 1 and then randomly select another ancestor distribution for sites 2-L
# Use a random function to decide whether recombination occurs between sites 1 and 2
## If there is no recombination, keep the selected ancestral distribution and continue
## If there is recombination, keep the selected site 1,2 and then randomly select another ancestor distribution for sites 3-L
# 
###Continue through recombination between sites L-1 and L
#######

library(sfsmisc)

#######
# set example L, q, n
# L = 5
# q = 0.1
# n = 100
#######

#######
# write a function to obtain the descendent sequences
descendent_sample <- function (L, q, n, seed, hapmap){
  ######
  # REMEMBER TO SET SEED
  set.seed(seed)
  ######
  
  #######
  # utilize function ancestor defined elsewhere to get the ancestor dist'n for the entered length
  ancestor_pi <- ancestor(L, hapmap)
  colnames(ancestor_pi) <- "proportion"
  #######
  
  #######
  # initialize matrix to hold descendant sequences
  descendant_seq <- matrix(NA, nrow = n, ncol = L)
  #######
  
  #######
  # select whether recombination occurs at each of the possible locations; 0 no recombination, 1 recombination
  ## we need n*L samples because we need one for each descendant and then one for each site on the sequence for each descendant 
  ## probability is probability of recombination
  recombo_indx <- rbinom(n = L*n , size = 1, prob = q )
  head(recombo_indx)
  #######
  
  #######
  # select ancestor sequences
  ## sample from range 1:2^L, indices of rows in ancestor_pi
  ## number we are sampling same logic as above
  ancest_indx <- sample((1:(2^L)), size = L*n , prob = ancestor_pi, replace = TRUE)
  head(ancest_indx)
  #######
  
  #######
  # use a for loop to pick the n descendants
  ## iterate over the rows
  for(i in 1:n){
    # start with first column for each row
    start <- ((i-1)*L )
    # as starting point, set each descendant to initial ancestor sequence
    descendant_seq[i,] <- drop(digitsBase(x = ancest_indx[start + 1] - 1, base = 2, ndigits = L))
    
    # iterate over locations of possible recombinations
    for (k in 2:L){
      descendant_seq[i, k:L] <- (1 - recombo_indx[start + k])*descendant_seq[i, k:L] + recombo_indx[start + k]*drop(digitsBase(ancest_indx[k+start] - 1, base = 2, ndigits = L)[k:L])
    }
  }
  return(descendant_seq)
  #######
}
#######

#######
# testing function
# descendent1 <- descendent_sample(L, q, n, seed = 7121, hapmap_binary)
#######
