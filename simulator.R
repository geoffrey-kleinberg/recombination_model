library('stringr')

# actual MCCL estimator
source("run_analysis_MCCL.R")

marg_ests <- c(estimates_m0, estimates_m1, estimates_m2, estimates_m3)

joint_ests <- c(ancestor_pair_estimation, ancestor_three_estimation, ancestor_four_estimation)

get_and_write_estimate <- function(len, q_val, m, n, file_name) {
  
  cat('', file=str_glue('rawResultsMCCLm{m}/{file_name}.txt'))
  
  for (run in 1:50) {
    
    name <- str_glue("simData/{file_name}/s{run}.txt")
    
    start_time <- Sys.time()
    
    m1_joint_est <- get_one_estimate(len, q_val, m, n, name)
    
    total_time <- Sys.time() - start_time
    
    result_file <- str_glue("rawResultsMCCLm{m}/{file_name}.txt")
    
    cat(paste(m1_joint_est, collapse='\t'),file=result_file, append=TRUE)
    cat('\t', file=result_file, append=TRUE)
    cat(total_time, file=result_file, append=TRUE)
    cat("\n",file=result_file, append=TRUE)
    
  }
  
}

get_one_estimate <- function(len, q_val, m, n, file_name) {
  SNP_data <- readLines(file_name)
  
  SNP_data <- str_split_fixed(SNP_data, '', len)
  SNP_data <- as.data.frame(SNP_data)
  SNP_data <- sapply(SNP_data, as.numeric)
  
  # previous marginal estimates, m-1
  m0_marg_est <- marg_ests[[m]](d = SNP_data, L = len, n = n, q = q_val)
  
  # marginal estimates, m
  m1_marg_est <- marg_ests[[m + 1]](L = len, q = q_val, n = n, d = SNP_data)
  
  # return the joint estimate
  joint_ests[[m]](L = len, m=m, m1_est = m1_marg_est, m_est = m0_marg_est)
}


# replace parameters to do one estimate
print(get_one_estimate(len = 3, q_val = 0.001, m = 1, n = 100, file_name = 'simData/l3q1pi1/s1.txt'))
