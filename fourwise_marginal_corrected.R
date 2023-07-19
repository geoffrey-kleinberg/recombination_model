# date created: July 8, 2021
# date edited: June 30, 2022

# this script calculates the fourwise marginal distributions (CORRECTED VERSION)
## functions estimates_m0, estimates_m1, estimates_m2, and descendent_sample are defined in their own respective scripts

## the four free parameters are:
# phi(0,0) = pi_{s,...,s+3}(0,0,0,0) GROUP 1
# phi(0,1) = pi_{s,...,s+3}(0,0,1,0) GROUP 2
# phi(1,0) = pi_{s,...,s+3}(0,1,0,0) GROUP 3
# phi(1,1) = pi_{s,...,s+3}(0,1,1,0) GROUP 4

#### SUPPORTING FUNCTIONS ####
########
# define functions for coefficients of phi

# this function is for the coefficients a_{xs, xs+2}^s
## same function for all four phi's
a_s_four <- function(q){
  a_00 <- (1-q)^3
  a_01 <- -(1-q)^3
  a_10 <- -(1-q)^3
  a_11 <- (1-q)^3
  return(c(a_00, a_01, a_10, a_11))
}

# this function is for the constants b_{xs, xs+2}^s
## this function is for phi(0,0)
b_s_gr1 <- function(q, pi_one, pi_two, pi_three){
  b_00 <- (q*((1-q)^2)*sum(pi_one[1,1]*pi_three[1,2], pi_two[1,1]*pi_two[1,3], pi_three[1,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[1,2]*pi_two[1,3], pi_two[1,1]*pi_one[1,3]*pi_one[1,4], pi_one[1,1]*pi_two[1,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[1,2]*pi_one[1,3]*pi_one[1,4]))
  b_01 <- ((1-q)^3*pi_three[1,1] +
             q*((1-q)^2)*sum(pi_one[1,1]*pi_three[2,2], pi_two[1,1]*pi_two[2,3], pi_three[1,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[1,2]*pi_two[2,3], pi_two[1,1]*pi_one[1,3]*pi_one[2,4], pi_one[1,1]*pi_two[1,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[1,2]*pi_one[1,3]*pi_one[2,4]))
  b_10 <- ((1-q)^3*pi_three[1,2] +
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[1,2], pi_two[3,1]*pi_two[1,3], pi_three[5,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[1,2]*pi_two[1,3], pi_two[3,1]*pi_one[1,3]*pi_one[1,4], pi_one[2,1]*pi_two[1,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[1,2]*pi_one[1,3]*pi_one[1,4]))
  b_11 <- ((1-q)^3*(pi_two[1,2] - pi_three[1,1] - pi_three[1,2]) +
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[2,2], pi_two[3,1]*pi_two[2,3], pi_three[5,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[1,2]*pi_two[2,3], pi_two[3,1]*pi_one[1,3]*pi_one[2,4], pi_one[2,1]*pi_two[1,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[1,2]*pi_one[1,3]*pi_one[2,4]))
  return(c(b_00, b_01, b_10, b_11))
}

# this function is for the constants b_{xs, xs+2}^s
## this function is for phi(0,1)
b_s_gr2 <- function(q, pi_one, pi_two, pi_three){
  b_00 <- (q*((1-q)^2)*sum(pi_one[1,1]*pi_three[3,2], pi_two[1,1]*pi_two[3,3], pi_three[2,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[1,2]*pi_two[3,3], pi_two[1,1]*pi_one[2,3]*pi_one[1,4], pi_one[1,1]*pi_two[2,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[1,2]*pi_one[2,3]*pi_one[1,4]))

  b_01 <- ((1-q)^3*pi_three[2,1] +
             q*((1-q)^2)*sum(pi_one[1,1]*pi_three[4,2], pi_two[1,1]*pi_two[4,3], pi_three[2,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[1,2]*pi_two[4,3], pi_two[1,1]*pi_one[2,3]*pi_one[2,4], pi_one[1,1]*pi_two[2,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[1,2]*pi_one[2,3]*pi_one[2,4]))


  b_10 <- ((1-q)^3*pi_three[3,2] +
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[3,2], pi_two[3,1]*pi_two[3,3], pi_three[6,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[1,2]*pi_two[3,3], pi_two[3,1]*pi_one[2,3]*pi_one[1,4], pi_one[2,1]*pi_two[2,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[1,2]*pi_one[2,3]*pi_one[1,4]))


  b_11 <- ((1-q)^3*(pi_two[2,2] - pi_three[2,1] - pi_three[3,2]) +
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[4,2], pi_two[3,1]*pi_two[4,3], pi_three[6,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[1,2]*pi_two[4,3], pi_two[3,1]*pi_one[2,3]*pi_one[2,4], pi_one[2,1]*pi_two[2,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[1,2]*pi_one[2,3]*pi_one[2,4]))
  return(c(b_00, b_01, b_10, b_11))
}

# this function is for the constants b_{xs, xs+2}^s
## this function is for phi(1,0)
b_s_gr3 <- function(q, pi_one, pi_two, pi_three){
  b_00 <- (q*((1-q)^2)*sum(pi_one[1,1]*pi_three[5,2], pi_two[2,1]*pi_two[1,3], pi_three[3,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[2,2]*pi_two[1,3], pi_two[2,1]*pi_one[1,3]*pi_one[1,4], pi_one[1,1]*pi_two[3,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[2,2]*pi_one[1,3]*pi_one[1,4]))


  b_01 <- ((1-q)^3*pi_three[3,1] +
             q*((1-q)^2)*sum(pi_one[1,1]*pi_three[6,2], pi_two[2,1]*pi_two[2,3], pi_three[3,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[2,2]*pi_two[2,3], pi_two[2,1]*pi_one[1,3]*pi_one[2,4], pi_one[1,1]*pi_two[3,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[2,2]*pi_one[1,3]*pi_one[2,4]))

  b_10 <- ((1-q)^3*pi_three[5,2] +
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[5,2], pi_two[4,1]*pi_two[1,3], pi_three[7,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[2,2]*pi_two[1,3], pi_two[4,1]*pi_one[1,3]*pi_one[1,4], pi_one[2,1]*pi_two[3,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[2,2]*pi_one[1,3]*pi_one[1,4]))


  b_11 <- ((1-q)^3*(pi_two[3,2] - pi_three[3,1] - pi_three[5,2]) +
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[6,2], pi_two[4,1]*pi_two[2,3], pi_three[7,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[2,2]*pi_two[2,3], pi_two[4,1]*pi_one[1,3]*pi_one[2,4], pi_one[2,1]*pi_two[3,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[2,2]*pi_one[1,3]*pi_one[2,4]))
  return(c(b_00, b_01, b_10, b_11))
}

# this function is for the constants b_{xs, xs+2}^s
## this function is for phi(1,1)
b_s_gr4 <- function(q, pi_one, pi_two, pi_three){
  b_00 <- (q*((1-q)^2)*sum(pi_one[1,1]*pi_three[7,2], pi_two[2,1]*pi_two[3,3], pi_three[4,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[2,2]*pi_two[3,3], pi_two[2,1]*pi_one[2,3]*pi_one[1,4], pi_one[1,1]*pi_two[4,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[2,2]*pi_one[2,3]*pi_one[1,4]))

  b_01 <- ((1-q)^3*pi_three[4,1] +
             q*((1-q)^2)*sum(pi_one[1,1]*pi_three[8,2], pi_two[2,1]*pi_two[4,3], pi_three[4,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[2,2]*pi_two[4,3], pi_two[2,1]*pi_one[2,3]*pi_one[2,4], pi_one[1,1]*pi_two[4,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[2,2]*pi_one[2,3]*pi_one[2,4]))

  b_10 <- ((1-q)^3*pi_three[7,2] +
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[7,2], pi_two[4,1]*pi_two[3,3], pi_three[8,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[2,2]*pi_two[3,3], pi_two[4,1]*pi_one[2,3]*pi_one[1,4], pi_one[2,1]*pi_two[4,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[2,2]*pi_one[2,3]*pi_one[1,4]))


  b_11 <- ((1-q)^3*(pi_two[4,2] - pi_three[4,1] - pi_three[7,2]) +
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[8,2], pi_two[4,1]*pi_two[4,3], pi_three[8,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[2,2]*pi_two[4,3], pi_two[4,1]*pi_one[2,3]*pi_one[2,4], pi_one[2,1]*pi_two[4,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[2,2]*pi_one[2,3]*pi_one[2,4]))
  return(c(b_00, b_01, b_10, b_11))
}
########

########
# define functions for coefficients in cubic expression

## A1,2,3,4 for phi(0,0); group 1
A1_four <- function(n, a){
  A1 <- sum(n)*prod(a)
  return(A1)
}

A2_four <-function(n, a, b){
  t1 <- sum(n[1], n[2], n[3])*prod(a[1], a[2], a[3], b[4])
  t2 <- sum(n[1], n[2], n[4])*prod(a[1], a[2], b[3], a[4])
  t3 <- sum(n[1], n[3], n[4])*prod(a[1], b[2], a[3], a[4])
  t4 <- sum(n[2], n[3], n[4])*prod(b[1], a[2], a[3], a[4])
  return(sum(t1, t2, t3, t4))
}

A3_four <- function(n, a, b){
  t1 <- sum(n[1], n[2])*prod(a[1], a[2], b[3], b[4])
  t2 <- sum(n[1], n[3])*prod(a[1], b[2], a[3], b[4])
  t3 <- sum(n[1], n[4])*prod(a[1], b[2], b[3], a[4])
  t4 <- sum(n[2], n[3])*prod(b[1], a[2], a[3], b[4])
  t5 <- sum(n[2], n[4])*prod(b[1], a[2], b[3], a[4])
  t6 <- sum(n[3], n[4])*prod(b[1], b[2], a[3], a[4])
  return(sum(t1, t2, t3, t4, t5, t6))
}

A4_four <- function(n, a, b){
  t1 <- prod(n[1], a[1], b[2], b[3], b[4])
  t2 <- prod(n[2], b[1], a[2], b[3], b[4])
  t3 <- prod(n[3], b[1], b[2], a[3], b[4])
  t4 <- prod(n[4], b[1], b[2], b[3], a[4])
  return(sum(t1, t2, t3, t4))
}
########

########
# define functions for finding roots of cubic expression

## FUNCTION FOR PHI(0,0)=PI_{S,...,S+3}
roots_fourwise_gr1 <- function(q, n_vec, pi_one, pi_two, pi_three, lower_b, upper_b){
  # find a coefficients
  a <- a_s_four(q)
  # print("phi00 a's:")
  # print(a)

  # find b constants
  b_gr1 <- b_s_gr1(q, pi_one, pi_two, pi_three)
  # print("phi00 b's:")
  # print(b_gr1)

  # find cubic coefficients
  A1_gr1 <- A1_four(n_vec, a)
  A2_gr1 <- A2_four(n_vec, a, b_gr1)
  A3_gr1 <- A3_four(n_vec, a, b_gr1)
  A4_gr1 <- A4_four(n_vec, a, b_gr1)
  # print("phi00 A's:")
  # print(A1_gr1)
  # print(A2_gr1)
  # print(A3_gr1)
  # print(A4_gr1)
  # As_gr1 <- c(A1_four(n_vec, a), A2_four(n_vec, a, b_gr1), A3_four(n_vec, a, b_gr1),A4_four(n_vec, a, b_gr1))
  # print("phi00 A's:")
  # print(As_gr1)

  # find all roots of the cubic expression (real and complex)
  all_roots <- cubic(A1_gr1, A2_gr1, A3_gr1, A4_gr1)[c(1:3)] 
  #all_roots <- cubic(c(A1_gr1, A2_gr1, A3_gr1, A4_gr1))[c(1:3)] # change by JSUN 05-01-2023
  # all_roots <- polyroot(c(A4_gr1, A3_gr1, A2_gr1, A1_gr1))
  # print("phi00 all root:")
  # print(all_roots)

  # select real roots
  real_roots <- pick_real(all_roots)
  # print("phi00 real root:")
  # print(real_roots)

  # pick candidate estimate from real roots
  if (length(real_roots)==1){
    candidate <- real_roots
  }
  else{
    candidate <- sort(real_roots)[2]
  }
  # print("candidate:")
  # print(candidate)

  # select final estimate based on bounds
  estimate <- candidate*((candidate >= lower_b) & (candidate <= upper_b)) + lower_b*(candidate < lower_b) + upper_b*(candidate > upper_b)

  # return the final estimate of phi(0,0)
  return(round(estimate, 10))
}


## FUNCTION FOR PHI(0,1)##
roots_fourwise_gr2 <- function(q, n_vec, pi_one, pi_two, pi_three, lower_b, upper_b){
  # find a coefficients
  a <- a_s_four(q)
  # print("phi01 a's")
  # print(a)

  # find b constants
  b_gr2 <- b_s_gr2(q, pi_one, pi_two, pi_three)
  # print("phi01 b's")
  # print(b_gr2)

  # find cubic coefficients
  # As_gr2 <- c(A1_four(n_vec, a), A2_four(n_vec, a, b_gr2), A3_four(n_vec, a, b_gr2),A4_four(n_vec, a, b_gr2))

  # find cubic coefficients
  A1_gr2 <- A1_four(n_vec, a)
  A2_gr2 <- A2_four(n_vec, a, b_gr2)
  A3_gr2 <- A3_four(n_vec, a, b_gr2)
  A4_gr2 <- A4_four(n_vec, a, b_gr2)
  # print("phi01 A's:")
  # print(A1_gr2)
  # print(A2_gr2)
  # print(A3_gr2)
  # print(A4_gr2)

  # find all roots of the cubic expression (real and complex)
  all_roots <- cubic(A1_gr2, A2_gr2, A3_gr2, A4_gr2)[c(1:3)] 
  #all_roots <- cubic(c(A1_gr2, A2_gr2, A3_gr2, A4_gr2))[c(1:3)] # change by JSUN 05-03-2023
  # all_roots <- polyroot(c(A4_gr2, A3_gr2, A2_gr2, A1_gr2))
  # print("phi01 all roots:")
  # print(all_roots)

  # select real roots
  real_roots <- pick_real(all_roots)
  # print("phi01 real roots:")
  # print(real_roots)

  # pick candidate estimate from real roots
  if (length(real_roots)==1){
    candidate <- real_roots
  }
  else{
    candidate <- sort(real_roots)[2]
  }

  # select final estimate based on bounds
  estimate <- candidate*(candidate >= lower_b & candidate <= upper_b) + lower_b*(candidate < lower_b) + upper_b*(candidate > upper_b)

  # return the final estimate of phi(0,1)
  return(round(estimate, 10))
}

## FUNCTION FOR PHI(1,0)##
roots_fourwise_gr3 <- function(q, n_vec, pi_one, pi_two, pi_three, lower_b, upper_b){
  # find a coefficients
  a <- a_s_four(q)
  # print("phi10 a's:")
  # print(a)

  # find b constants
  b_gr3 <- b_s_gr3(q, pi_one, pi_two, pi_three)
  # print("phi10 b's")
  # print(b_gr3)

  # find cubic coefficients
  # As_gr3 <- c(A1_four(n_vec, a), A2_four(n_vec, a, b_gr3), A3_four(n_vec, a, b_gr3),A4_four(n_vec, a, b_gr3))
  # find cubic coefficients
  A1_gr3 <- A1_four(n_vec, a)
  A2_gr3 <- A2_four(n_vec, a, b_gr3)
  A3_gr3 <- A3_four(n_vec, a, b_gr3)
  A4_gr3 <- A4_four(n_vec, a, b_gr3)
  # print("phi10 A's:")
  # print(A1_gr3)
  # print(A2_gr3)
  # print(A3_gr3)
  # print(A4_gr3)

  # find all roots of the cubic expression (real and complex)
  all_roots <- cubic(A1_gr3, A2_gr3, A3_gr3, A4_gr3)[c(1:3)]
  # all_roots <- polyroot(c(A4_gr3, A3_gr3, A2_gr3, A1_gr3))
  # print("phi10 all roots:")
  # print(all_roots)

  # select real roots
  real_roots <- pick_real(all_roots)
  # print("phi10 real roots:")
  # print(real_roots)

  # pick candidate estimate from real roots
  if (length(real_roots)==1){
    candidate <- real_roots
  }
  else{
    candidate <- sort(real_roots)[2]
  }

  # select final estimate based on bounds
  estimate <- candidate*(candidate >= lower_b & candidate <= upper_b) + lower_b*(candidate < lower_b) + upper_b*(candidate > upper_b)

  # return the final estimate of phi(1,0)
  return(round(estimate, 10))
}

## FUNCTION FOR PHI(1,1)
roots_fourwise_gr4 <- function(q, n_vec, pi_one, pi_two, pi_three, lower_b, upper_b){
  # find a coefficients
  a <- a_s_four(q)
  # print("phi11 a's:")
  # print(a)

  # find b constants
  b_gr4 <- b_s_gr4(q, pi_one, pi_two, pi_three)
  # print("phi11 b's")
  # print(b_gr4)

  # find cubic coefficients
  # As_gr4 <- c(A1_four(n_vec, a), A2_four(n_vec, a, b_gr4), A3_four(n_vec, a, b_gr4),A4_four(n_vec, a, b_gr4))
  # find cubic coefficients
  A1_gr4 <- A1_four(n_vec, a)
  A2_gr4 <- A2_four(n_vec, a, b_gr4)
  A3_gr4 <- A3_four(n_vec, a, b_gr4)
  A4_gr4 <- A4_four(n_vec, a, b_gr4)
  # print("phi11 A's:")
  # print(A1_gr4)
  # print(A2_gr4)
  # print(A3_gr4)
  # print(A4_gr4)

  # find all roots of the cubic expression (real and complex)
  all_roots <- cubic(A1_gr4, A2_gr4, A3_gr4, A4_gr4)[c(1:3)]
  # all_roots <- polyroot(c(A4_gr4, A3_gr4, A2_gr4, A1_gr4))
  # print("phi11 all roots:")
  # print(all_roots)

  # select real roots
  real_roots <- pick_real(all_roots)
  # print("phi11 real roots:")
  # print(real_roots)

  # pick candidate estimate from real roots
  if (length(real_roots)==1){
    candidate <- real_roots
  }
  else{
    candidate <- sort(real_roots)[2]
  }

  # select final estimate based on bounds
  estimate <- candidate*(candidate >= lower_b & candidate <= upper_b) + lower_b*(candidate < lower_b) + upper_b*(candidate > upper_b)

  # return the final estimate of phi(1,1)
  return(round(estimate, 10))
}
########

#### MAIN FUNCTION ####
########
# define function that carries out the full fourwise estimation

estimates_m3 <- function(L, q, n, d){
  # fix onewise
  onewise <- estimates_m0(d, L, n)
  # fix pairwise
  pairwise <- estimates_m1(L, q, n , d)
  # fix threewise
  threewise <- estimates_m2(L, q, n, d)

  # initialize matrix to hold fourwise estimates
  fourwise_marginals <- matrix(NA, nrow = 16, ncol = L-3)

  ## all m+1-wise ancestral sequences
  an_seq <- matrix(c(0:(2^(3+1)-1)), ncol=1)
  an_seq <- t( apply(an_seq,1,digitsBase, base=2, ndigits=3+1) )
  an_seq <- matrix( as.character(an_seq), nrow=2^{3+1}, ncol=3+1)
  an_seq <- apply(an_seq, 1, paste, collapse="")

  rownames(fourwise_marginals) <- paste("pi_", an_seq ,sep="")
  colnames(fourwise_marginals) <- paste("site_",seq(1:(L-3)),sep="")

  # iterate over the fourwise combos
  for (s in 1:(L-3)){
    # find relevant onewise
    pi_one <- onewise[,s:(s+3)]
    # print(pi_one)

    # find relevant pairwise
    pi_two <- pairwise[,s:(s+2)]

    # find relevant threewise
    pi_three <- threewise[,s:(s+1)]
    # print(pi_three)

    ### FIND PHI(0,0)=PI_{S,...,S+3}(0,0,0,0) ###

    # calculate num observations for ea sequence
    n_0000 <- sum(d[,s]==0 & d[,s+1]==0 & d[,s+2]==0 & d[,s+3]==0) #pi(0,0,0,0)
    n_0001 <- sum(d[,s]==0 & d[,s+1]==0 & d[,s+2]==0 & d[,s+3]==1) #pi(0,0,0,1)
    n_1000 <- sum(d[,s]==1 & d[,s+1]==0 & d[,s+2]==0 & d[,s+3]==0) #pi(1,0,0,0)
    n_1001 <- sum(d[,s]==1 & d[,s+1]==0 & d[,s+2]==0 & d[,s+3]==1) #pi(1,0,0,1)
    n_gr1 <- c(n_0000, n_0001, n_1000, n_1001)
    # print("phi00 n's:")
    # print(n_gr1)

    # calculate bounds on phi(0,0)
    lower_gr1 <- max(0, (pi_three[1,2]+ pi_three[1,1] - pi_two[1,2]))
    upper_gr1 <- min(pi_three[1,2], pi_three[1,1])
    # print("bounds on phi00")
    # print(lower_gr1)
    # print(upper_gr1)

    # estimate phi(0,0)
    if (abs(lower_gr1 - upper_gr1) <= 10^(-10)){
      hat_phi_gr1 <- (lower_gr1 + upper_gr1)/2
    }
    if (abs(lower_gr1 - upper_gr1) > 10^(-10)){
      hat_phi_gr1 <- roots_fourwise_gr1(q, n_gr1, pi_one, pi_two, pi_three, lower_b = lower_gr1, upper_b = upper_gr1)
    }

    # estimate other pi's in group 1 based on constrained MLE for phi(0,0)
    pi_0000 <- hat_phi_gr1
    pi_0001 <- round((pi_three[1,1] - hat_phi_gr1), 10)
    pi_1000 <- round((pi_three[1,2] - hat_phi_gr1), 10)
    pi_1001 <- round((pi_two[1,2] - pi_three[1,1] - pi_three[1,2] + hat_phi_gr1), 10)


    ### FIND PHI(0,1)=PI_{S,...,S+3}(0,0,1,0) ###

    # calculate num observations for ea sequence
    n_0010 <- sum(d[,s]==0 & d[,s+1]==0 & d[,s+2]==1 & d[,s+3]==0) #pi(0,0,1,0)
    n_0011 <- sum(d[,s]==0 & d[,s+1]==0 & d[,s+2]==1 & d[,s+3]==1) #pi(0,0,1,1)
    n_1010 <- sum(d[,s]==1 & d[,s+1]==0 & d[,s+2]==1 & d[,s+3]==0) #pi(1,0,1,0)
    n_1011 <- sum(d[,s]==1 & d[,s+1]==0 & d[,s+2]==1 & d[,s+3]==1) #pi(1,0,1,1)
    n_gr2 <- c(n_0010, n_0011, n_1010, n_1011)
    # print("phi01 n's:")
    # print(n_gr2)

    # calculate bounds on phi(0,1)
    lower_gr2 <- max(0, (pi_three[3,2] + pi_three[2,1] - pi_two[2,2]))
    upper_gr2 <- min(pi_three[2,1], pi_three[3,2])
    # print("bounds on phi01")
    # print(lower_gr2)
    # print(upper_gr2)

    # estimate phi(0,1)
    if (abs(lower_gr2 - upper_gr2) <= 10^(-10)){
      hat_phi_gr2 <- (lower_gr2 + upper_gr2)/2
    }
    if (abs(lower_gr2 - upper_gr2) > 10^(-10)){
      hat_phi_gr2 <- roots_fourwise_gr2(q, n_gr2, pi_one, pi_two, pi_three, lower_gr2, upper_gr2)
    }

    # calculate estimates for group 2 based on constrained MLE of phi(0,1)
    pi_0010 <- hat_phi_gr2
    pi_0011 <- round((pi_three[2,1] - hat_phi_gr2), 10)
    pi_1010 <- round((pi_three[3,2] - hat_phi_gr2), 10)
    pi_1011 <- round((pi_two[2,2] - pi_three[2,1] - pi_three[3,2] + hat_phi_gr2), 10)

    ### FIND PHI(1,0)=PI_{S,...,S+3}(0,1,0,0) ###

    # calculate num observations for ea sequence
    n_0100 <- sum(d[,s]==0 & d[,s+1]==1 & d[,s+2]==0 & d[,s+3]==0) #pi(0,1,0,0)
    n_0101 <- sum(d[,s]==0 & d[,s+1]==1 & d[,s+2]==0 & d[,s+3]==1) #pi(0,1,0,1)
    n_1100 <- sum(d[,s]==1 & d[,s+1]==1 & d[,s+2]==0 & d[,s+3]==0) #pi(1,1,0,0)
    n_1101 <- sum(d[,s]==1 & d[,s+1]==1 & d[,s+2]==0 & d[,s+3]==1) #pi(1,1,0,1)
    n_gr3 <- c(n_0100, n_0101, n_1100, n_1101)
    # print("phi10 n's:")
    # print(n_gr3)

    # calculate bounds on phi(1,0)
    lower_gr3 <- max(0, (pi_three[5,2] + pi_three[3,1] - pi_two[3,2]))
    upper_gr3 <- min(pi_three[3,1], pi_three[5,2])
    # print("bounds on phi10")
    # print(lower_gr3)
    # print(upper_gr3)

    # estimate phi(1,0)
    if (abs(lower_gr3 - upper_gr3) <= 10^(-10)){
      hat_phi_gr3 <- (lower_gr3 + upper_gr3)/2
    }
    if (abs(lower_gr3 - upper_gr3) > 10^(-10)){
      hat_phi_gr3 <- roots_fourwise_gr3(q, n_gr3, pi_one, pi_two, pi_three, lower_gr3, upper_gr3)
    }

    # calculate estimates for group 3 based on constrained MLE for phi(1,0)
    pi_0100 <- hat_phi_gr3
    pi_0101 <- round((pi_three[3,1] - hat_phi_gr3), 10)
    pi_1100 <- round((pi_three[5,2] - hat_phi_gr3), 10)
    pi_1101 <- round((pi_two[3,2] - pi_three[3,1] - pi_three[5,2] + hat_phi_gr3), 10)

    ### FIND PHI(0,0)=PI_{S,...,S+3}(0,1,1,0) ###

    # calculate num observations for ea sequence
    n_0110 <- sum(d[,s]==0 & d[,s+1]==1 & d[,s+2]==1 & d[,s+3]==0) #pi(0,1,1,0)
    n_0111 <- sum(d[,s]==0 & d[,s+1]==1 & d[,s+2]==1 & d[,s+3]==1) #pi(0,1,1,1)
    n_1110 <- sum(d[,s]==1 & d[,s+1]==1 & d[,s+2]==1 & d[,s+3]==0) #pi(1,1,1,0)
    n_1111 <- sum(d[,s]==1 & d[,s+1]==1 & d[,s+2]==1 & d[,s+3]==1) #pi(1,1,1,1)
    n_gr4 <- c(n_0110, n_0111, n_1110, n_1111)
    # print("phi11 n's:")
    # print(n_gr4)

    # calculate bounds on phi(1,1)
    lower_gr4 <- max(0, (pi_three[7,2] + pi_three[4,1] - pi_two[4,2]))
    upper_gr4 <- min(pi_three[4,1], pi_three[7,2])
    # print("bounds on phi11")
    # print(lower_gr4)
    # print(upper_gr4)

    # estimate phi(1,1)
    if (abs(lower_gr4 - upper_gr4) <= 10^(-10)){
      hat_phi_gr4 <- (lower_gr4 + upper_gr4)/2
    }
    if (abs(lower_gr4 - upper_gr4) > 10^(-10)){
      hat_phi_gr4 <- roots_fourwise_gr4(q, n_gr4, pi_one, pi_two, pi_three, lower_gr4, upper_gr4)
    }

    # estimate others in group 4 based on constrained MLE phi(1,1)
    pi_0110 <- hat_phi_gr4
    pi_0111 <- round((pi_three[4,1] - hat_phi_gr4), 10)
    pi_1110 <- round((pi_three[7,2] - hat_phi_gr4), 10)
    pi_1111 <- round((pi_two[4,2] - pi_three[4,1] - pi_three[7,2] + hat_phi_gr4), 10)

    # append to fourwise matrix
    fourwise_marginals[,s] <- c(pi_0000, pi_0001, pi_0010, pi_0011, pi_0100, pi_0101, pi_0110, pi_0111, pi_1000, pi_1001, pi_1010, pi_1011, pi_1100, pi_1101, pi_1110, pi_1111)
  }
  return(fourwise_marginals)
}
#######
