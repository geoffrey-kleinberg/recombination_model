# date created: July 6, 2021
## edited on August 30, 2021 to replace polyroot() with cubic()

# this script defines a function that calculates the threewise marginal ancestor distributions

# #######
# # set example L, q, n
# L = 5
# q = 0.1
# n = 100
# # descendant sequences for these examples
# descend <- descendent_sample(L, q, n, seed = 7121)
# # onewise pi for these examples is onewise
# #######

#######
# write functions to formulate the cubic equation

# find the coefficient of phi
a_s <- function(q){
  # order: a_00, a_01, a_10, a_11
  a <- c((1-q)^2, -(1-q)^2, -(1-q)^2 , (1-q)^2)
  return(a)
}

# find the constant b for phi(0)
b_s0 <- function(q, pi_one, pi_two){
  b_00 <- q*(1-q)*(pi_one[1,1]*pi_two[1,2] + pi_two[1,1]*pi_one[1,3]) + (q^2)*pi_one[1,1]*pi_one[1,2]*pi_one[1,3]
  b_01 <- ((1-q)^2)*pi_two[1,1] + q*(1-q)*(pi_one[1,1]*pi_two[2,2] + pi_two[1,1]*pi_one[2,3]) + (q^2)*pi_one[1,1]*pi_one[1,2]*pi_one[2,3]
  b_10 <- ((1-q)^2)*pi_two[1,2] + q*(1-q)*(pi_one[2,1]*pi_two[1,2] + pi_two[3,1]*pi_one[1,3]) + (q^2)*pi_one[2,1]*pi_one[1,2]*pi_one[1,3]
  b_11 <- ((1-q)^2)*(pi_one[1,2] - pi_two[1,1] - pi_two[1,2]) + q*(1-q)*(pi_one[2,1]*pi_two[2,2] + pi_two[3,1]*pi_one[2,3]) + (q^2)*pi_one[2,1]*pi_one[1,2]*pi_one[2,3]
  
  return(c(b_00, b_01, b_10, b_11))
}

# find the coefficient of cubic term
A_1 <- function(n, a){
  A1 <- sum(n)*prod(a)
  return(A1)
}

# find the coefficient of squared term
A_2 <- function(n, a, b){
  t1 <- prod(a[1], a[2], a[3], b[4])*sum(n[1], n[2], n[3])
  t2 <- prod(a[1], a[2], b[3], a[4])*sum(n[1], n[2], n[4])
  t3 <- prod(a[1], b[2], a[3], a[4])*sum(n[1], n[3], n[4])
  t4 <- prod(b[1], a[2], a[3], a[4])*sum(n[2], n[3], n[4])
  A2 <- sum(c(t1,t2,t3,t4))
  
  return(A2)
}

# find the coefficient of linear term
A_3 <- function(n, a, b){
  t1 <- prod(a[1], a[2], b[3], b[4])* sum(n[1], n[2])
  t2 <- prod(a[1], b[2], a[3], b[4])* sum(n[1], n[3])
  t3 <- prod(a[1], b[2], b[3], a[4])* sum(n[1], n[4])
  t4 <- prod(b[1], a[2], a[3], b[4])* sum(n[2], n[3])
  t5 <- prod(b[1], a[2], b[3], a[4])* sum(n[2], n[4])
  t6 <- prod(b[1], b[2], a[3], a[4])* sum(n[3], n[4])
  A3 <- sum(c(t1, t2, t3, t4, t5, t6))
  
  return(A3)
}

# find the constant in cubic expression
A_4 <- function(n, a, b){
  t1 <- prod(a[1], b[2], b[3], b[4])*n[1]
  t2 <- prod(b[1], a[2], b[3], b[4])*n[2]
  t3 <- prod(b[1], b[2], a[3], b[4])*n[3]
  t4 <- prod(b[1], b[2], b[3], a[4])*n[4]
  A4 <- sum(c(t1, t2, t3, t4))
  
  return(A4)
}

# function to use to pick the real roots
pick_real = function(aaa){
  ## sequence of complex numbers to real numbers
  
  idx  = (abs(Im(aaa))>1e-4);
  return(Re(aaa)[!idx])
}

# find roots of the cubic equation for phi(0)
roots_threewise_0 <- function(q, n, pi_one, pi_two, lower_b, upper_b){
  
  # find a's and b's
  a_0 <- a_s(q)
  # print(a_0)
  b_0 <- b_s0(q, pi_one, pi_two)
  # print(b_0)
  
  # find A's
  A1_0 <- A_1(n, a_0)
  A2_0 <- A_2(n, a_0, b_0)
  A3_0 <- A_3(n, a_0, b_0)
  A4_0 <- A_4(n, a_0, b_0)
  # print("phi0 cubic coefficients:")
  # print(A1_0)
  # print(A2_0)
  # print(A3_0)
  # print(A4_0)
  
  # calculate the roots
  ## USING OLD FUNCTION
  # all_roots <- polyroot(c(A4_0, A3_0, A2_0, A1_0))
  ## USING NEW FUNCTION
  all_roots <- cubic(A1_0, A2_0, A3_0, A4_0)[c(1:3)]
  #coef <- c(A1_0, A2_0, A3_0, A4_0) #change by JSUN 05-01-2023
  #all_roots <- cubic(coef)[c(1:3)]
  # print("phi0 all roots")
  # print(all_roots)
  
  # find the real roots 
  real_roots <- pick_real(all_roots)
  # print("phi0 real roots")
  # print(real_roots)
  
  # if there are 3 real roots, pick the middle one
  if (length(real_roots) == 1){
    candidate <- real_roots
  } else {
    candidate <- sort(real_roots)[2]
  }
  
  
  # make final choice of estimate based on bounds
  estimate <- candidate*(lower_b <= candidate & candidate <= upper_b) + upper_b*(candidate > upper_b) + lower_b*(candidate < lower_b)
  
  return(round(estimate, 10))
}

# this function finds the constant b for phi(1)
b_s1 <- function(q, pi_one, pi_two){
  b_00 <- q*(1-q)*(pi_one[1,1]*pi_two[3,2] + pi_two[2,1]*pi_one[1,3]) + (q^2)*pi_one[1,1]*pi_one[2,2]*pi_one[1,3]
  b_01 <- ((1-q)^2)*pi_two[2,1] + q*(1-q)*(pi_one[1,1]*pi_two[4,2] + pi_two[2,1]*pi_one[2,3]) + (q^2)*pi_one[1,1]*pi_one[2,2]*pi_one[2,3]
  b_10 <- ((1-q)^2)*pi_two[3,2] + q*(1-q)*(pi_one[2,1]*pi_two[3,2] + pi_two[4,1]*pi_one[1,3]) + (q^2)*pi_one[2,1]*pi_one[2,2]*pi_one[1,3]
  b_11 <- ((1-q)^2)*(pi_one[2,2] - pi_two[2,1] - pi_two[3,2]) + q*(1-q)*(pi_one[2,1]*pi_two[4,2] + pi_two[4,1]*pi_one[2,3]) + (q^2)*pi_one[2,1]*pi_one[2,2]*pi_one[2,3]
  
  return(c(b_00, b_01, b_10, b_11))
}

# find roots of the cubic equation for phi(1)
roots_threewise_1 <- function(q, n, pi_one, pi_two, lower_b, upper_b){
  
  # find a's and b's
  a_1 <- a_s(q)
  b_1 <- b_s1(q, pi_one, pi_two)
  # print(b_1)
  
  # find A's
  A1_1 <- A_1(n, a_1)
  A2_1 <- A_2(n, a_1, b_1)
  A3_1 <- A_3(n, a_1, b_1)
  A4_1 <- A_4(n, a_1, b_1)
  # print(A1_1)
  # print(A2_1)
  # print(A3_1)
  # print(A4_1)
  
  # calculate the roots
  ### USING OLD FUNCTION
  # all_roots <- polyroot(c(A4_1, A3_1, A2_1, A1_1))
  ### USING NEW FUNCTION
  all_roots <- cubic(A1_1, A2_1, A3_1, A4_1)[c(1:3)] 
  #coef <- c(A1_1, A2_1, A3_1, A4_1) # change by JSUN,05/01/23
  #all_roots <- cubic(coef)[c(1:3)]
  # print("phi1 all roots")
  # print(all_roots)
  
  # find the real roots 
  real_roots <- pick_real(all_roots)
  # print("phi1 real roots")
  # print(real_roots)
  
  # if there are 3 real roots, pick the middle one
  if (length(real_roots) == 1){
    candidate <- real_roots
  } else {
    candidate <- sort(real_roots)[2]
  }
  
  # make final choice of estimate based on bounds
  estimate <- candidate*(lower_b <= candidate & candidate <= upper_b) + upper_b*(candidate > upper_b) + lower_b*(candidate < lower_b)
  
  return(round(estimate, 10))
}
#######



#######
# write a function to calculate the threewise marginal estimates
estimates_m2 <- function(L, q, n, d){
  # fix onewise estimates
  pi_one <- estimates_m0(d, L,n)
  # fix pairwise estimates
  pi_two <- estimates_m1(L, q, n, d)
  
  # initialize matrix to put the threewise margin estimates into
  threewise_margins <- matrix(NA, nrow = 8, ncol = (L-2))
  
  ## all m+1-wise ancestral sequences
  an_seq <- matrix(c(0:(2^(2+1)-1)), ncol=1)
  an_seq <- t( apply(an_seq,1,digitsBase, base=2, ndigits=2+1) ) 
  an_seq <- matrix( as.character(an_seq), nrow=2^{2+1}, ncol=2+1)
  an_seq <- apply(an_seq, 1, paste, collapse="")
  
  rownames(threewise_margins) <- paste("pi_", an_seq ,sep="")
  colnames(threewise_margins) <- paste("site_",seq(1:(L-2)),sep="")
  
  #### ESTIMATE PHIS ####
  # iterate over groups of 3 sites
  for (j in 1:(L-2)){
    pi_ones <- pi_one[,j:(j+2)]
    pi_pairs <- pi_two[,j:(j+1)]
    
    ### first estimate phi(0) ###
    # calculate n(xs, 0, xs+2)
    n_000 <- sum(d[,j]==0 & d[,j+1]==0 & d[,j+2]==0)
    n_001 <- sum(d[,j]==0 & d[,j+1]==0 & d[,j+2]==1)
    n_100 <- sum(d[,j]==1 & d[,j+1]==0 & d[,j+2]==0)
    n_101 <- sum(d[,j]==1 & d[,j+1]==0 & d[,j+2]==1)
    n_0 <- c(n_000, n_001, n_100, n_101)
    # print(n_0)
    
    # find bounds for phi(0)
    lower_phi0 <- max(0, (pi_pairs[1,1] + pi_pairs[1,2] - pi_ones[1,2]))
    upper_phi0 <- min(pi_pairs[1,1], pi_pairs[1,2])
    # print(lower_phi0)
    # print(upper_phi0)
    
    # estimate phi(0)
    if (abs(lower_phi0 - upper_phi0) <= 10^(-10)){
      hat_phi0 <- (lower_phi0 + upper_phi0)/2
    }
    if (abs(lower_phi0 - upper_phi0) > 10^(-10)){
      hat_phi0 <- roots_threewise_0(q, n_0, pi_ones, pi_pairs, lower_phi0, upper_phi0)
    }
    
    # calculate other pi's based on estimated phi
    pi_000 <- hat_phi0
    pi_001 <- round(pi_pairs[1,1] - hat_phi0, 10)
    pi_100 <- round(pi_pairs[1,2] - hat_phi0, 10)
    pi_101 <- round(pi_ones[1,2] - pi_pairs[1,2] - pi_pairs[1,1] + hat_phi0, 10)
    
    ### now estimate phi(1) ###
    # calculate the n(xs, 1, xs+2)
    n_010 <- sum(d[,j]==0 & d[,j+1]==1 & d[,j+2]==0)
    n_011 <- sum(d[,j]==0 & d[,j+1]==1 & d[,j+2]==1)
    n_110 <- sum(d[,j]==1 & d[,j+1]==1 & d[,j+2]==0)
    n_111 <- sum(d[,j]==1 & d[,j+1]==1 & d[,j+2]==1)
    n_1 <- c(n_010, n_011, n_110, n_111)
    
    # find bounds for phi(1)
    lower_phi1 <- max(0, (pi_pairs[2,1] + pi_pairs[3,2] - pi_ones[2,2]))
    upper_phi1 <- min(pi_pairs[2,1], pi_pairs[3,2])
    # print(lower_phi1)
    # print(upper_phi1)
    
    #estimate phi(1)
    if (abs(lower_phi1 - upper_phi1) <= 10^(-10)){
      hat_phi1 <- (lower_phi1 + upper_phi1)/2
    }
    if (abs(lower_phi1 - upper_phi1) > 10^(-10)){
      hat_phi1 <- roots_threewise_1(q, n_1, pi_ones, pi_pairs, lower_phi1, upper_phi1)
    }
    
    # calculate other pi's
    pi_010 <- hat_phi1
    pi_011 <- round(pi_pairs[2,1] - hat_phi1, 10)
    pi_110 <- round(pi_pairs[3,2] - hat_phi1, 10)
    pi_111 <- round(pi_ones[2,2] - pi_pairs[3,2] - pi_pairs[2,1] + hat_phi1, 10)
    
    ### append estimates into matrix 
    threewise_margins[,j] <- c(pi_000,pi_001,pi_010,pi_011,pi_100,pi_101,pi_110,pi_111)
    
  }
  return(threewise_margins)
}
#######

#######
# find threewise margins for the set example:
# threewise <- estimates_m2(L, q, n, descend)

