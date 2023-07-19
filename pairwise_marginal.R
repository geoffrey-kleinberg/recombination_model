# date created: July 6, 2021
## edited on August 30, 2021 to account for the case where the bounds on phi are equal/ approx equal

# this script defines a function that calculates the pairwise marginal distributions, holding fixed the onewise marginal distributions (see onewise_marginal.R for function estimates_m0)

# #######
# # set example L, q, n
# L = 5
# q = 0.1
# n = 100
# # descendant sequences for these examples
# d <- descendent_sample(L,q,n,seed=7121)
# # onewise pi for these examples is onewise
# #######

#######
estimates_m1 <- function(L, q, n, d){
  # fix lower order marginals
  pi_m0 <- estimates_m0(d, L, n) 

  # initialize matrix for pairwise estimates; this will be filled in and returned
  pairwise_margins <- matrix(NA, nrow = 4, ncol = (L-1))
  rownames(pairwise_margins) <- paste("pi_", c("0,0", "0,1" , "1,0" , "1,1") )
  colnames(pairwise_margins) <- paste( "site_" , seq(1:(L-1)))
  
  # iterate over the pairs
  for (j in 1:(L-1)){
    # select relevant onewise marginals, pi_s(0)
    pi_zeroes <- pi_m0[1, j:(j+1)]
    
    # print("site:")
    # print(j)
    
    # calculate upper and lower bounds
    l_bound <- max(0, (sum(pi_zeroes)-1))
    u_bound <- min(pi_zeroes)
    # print("Bounds on phi:")
    # print(l_bound)
    # print(u_bound)
    
    #check if there is space bw bounds; if not set estimate to be avg of the bounds
    if (abs(l_bound - u_bound) <= 10^(-10)){
      phi <- (l_bound + u_bound)/2
    }
    if (abs(l_bound - u_bound) > 10^(-10)){
      # calculate unconstrained MLE
      n_00 <- sum(d[,j]==0 & d[,j+1]==0)
      mle <- (1 / (1-q))*((n_00/n) - q*pi_zeroes[1]*pi_zeroes[2])
      # print("unconstrained mle:")
      # print(mle)
      
      # find constrained estimate
      phi <- ((l_bound <= mle) & (mle <= u_bound))*mle + l_bound*(mle < l_bound) + u_bound*(mle > u_bound)
      # print("constrained estimate:")
      # print(phi)
    }
    
    # calculate pi's from phi
    pi_00 <- phi
    pi_01 <- pi_zeroes[1] - phi
    pi_10 <- pi_zeroes[2] - phi
    pi_11 <- 1 - pi_zeroes[1] - pi_zeroes[2] + phi
    
    #attach to matrix
    pairwise_margins[,j] <- round(c(pi_00, pi_01, pi_10, pi_11), 10)
  }
  
  return (pairwise_margins)
}
####### 
# pairwise <- estimates_m1(L, q, n, d)