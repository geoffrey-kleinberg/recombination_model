# date created: July 2, 2021

# this script calculates the onewise-marginal distributions as part of the hierarchical estimator

########
# for onewise margins, m=0, the MLE is the sample proportions
# the sample refers to the generated descendant sequences
########

#######
# set example L, q, n
# L = 5
# q = 0.1
# n = 100
# descendant sequences for these examples is descendent1
#######

########
estimates_m0 <- function(d, L, n, q){
  onewise_margins <- matrix(NA, nrow = 2, ncol = L)

  
  rownames(onewise_margins) <- paste("pi_",c(0,1), sep="")
  colnames(onewise_margins) <- paste("site_", (1:(L)), sep="")
  for (i in 1:(dim(d)[2])){
    onewise_margins[1,i] <- (n - sum(d[,i]))/n
    onewise_margins[2,i] <- sum(d[,i])/n
  }
  return(onewise_margins)
}
########

########
# calculate for set examples
# onewise <- estimates_m0(descendent_sample(L,q,n,seed=7121), L)
########