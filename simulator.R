library('stringr')

# actual MCCL estimator
source("run_analysis_MCCL.R")

lengths <- c(3, 4, 5)
qs <- c(0.001, 0.01, 0.05, 0.1)
nn <- 100

marg_ests <- c(estimates_m0, estimates_m1, estimates_m2, estimates_m3)

joint_ests <- c(ancestor_pair_estimation, ancestor_three_estimation, ancestor_four_estimation)

for (m in 1:3) {
  Len = m + 2
  while (Len <= 5) {
    for (j in 1:4) {
      for (k in 1:3) {
        
        print(Len)
      
        q_val <- qs[j]
        
        file_name <- str_glue('l{Len}q{j}pi{k}')
        cat('', file=str_glue('rawResultsMCCLm{m}/{file_name}.txt'))
        
        for (run in 1:50) {
          
          start_time <- Sys.time()
          SNP_data <- readLines(str_glue("simData/{file_name}/s{run}.txt"))
          
          SNP_data <- str_split_fixed(SNP_data, '', Len)
          SNP_data <- as.data.frame(SNP_data)
          SNP_data <- sapply(SNP_data, as.numeric)
          
          # previous marginal estimates, m-1
          m0_marg_est <- marg_ests[[m]](d = SNP_data, L = Len, n = nn, q = q_val)
          
          # marginal estimates, m
          m1_marg_est <- marg_ests[[m + 1]](L = Len, q = q_val, n = nn, d = SNP_data)
          
          m1_joint_est <- joint_ests[[m]](L = Len, m=m, m1_est = m1_marg_est, m_est = m0_marg_est)
          
          total_time <- Sys.time() - start_time
          
          result_file <- str_glue("rawResultsMCCLm{m}/{file_name}.txt")
          
          cat(paste(m1_joint_est, collapse='\t'),file=result_file, append=TRUE)
          cat('\t', file=result_file, append=TRUE)
          cat(total_time, file=result_file, append=TRUE)
          cat("\n",file=result_file, append=TRUE)
          
        }
        
      }
    }
    Len <- Len + 1
  }
}
