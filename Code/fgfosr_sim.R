library(here)
library(fastGFPCA)
library(mgcv)
library(ggplot2)
library(splines)
library(tidyverse)
library(refund)
library(lme4)
library(parallel)
library(R.utils)

source(here::here("fgfosr_funks.R"))

seed_start = 1134

set.seed(seed_start)

theta_0 = rnorm(10 + 4, sd=1)
theta_1 = rnorm(10 + 4, sd=1)

#set true fixed effects
f_0 = function(s) bs(s, knots=seq(0.1, 0.9, len=10), Boundary.knots=c(0,1), intercept=TRUE) %*% theta_0
f_1 = function(s) bs(s, knots=seq(0.1, 0.9, len=10), Boundary.knots=c(0,1), intercept=TRUE) %*% theta_1

N_list = c(200)
J_list = c(200)
#c(200,500,1000)
bin_count_list = c(20,30,40)
N_iter = 100     # number of iterations for each simulation scenario

params = expand.grid(seed_start = seed_start, J = J_list, N = N_list)

final_results_df = data.frame("seed_start" = integer(), "N" = integer(), "J"=integer(), "bin_count" = integer(), "bias" = double(), "coverage" = double(), "MSE" = double(), runtime_min = double())

for(i in 1:nrow(params)){
  ########### Make le data
  set.seed(params$seed_start[i])
  gen_data_list = lapply(1:N_iter, function(x) gen_data(params$N[i], params$J[i], x))
  
  # Run gfosr with varying  bin counts
  for(r in bin_count_list){
    print(paste0("N: ", params$N[i], " J: ",params$J[i], " bin count: ", r))
    
    ############ Process le data
    fgfosr_start_t = Sys.time() #Time that sucka
    # if run on windows remove cores (default is 1) or set to 1 to remove parallel running
    mod_list = lapply(gen_data_list, function (x) get_fgfosr_mod(x , r, cores))
    fgfosr_end_t = Sys.time() #The cake is done!
    fgfosr_time_diff =  as.double(difftime(fgfosr_end_t, fgfosr_start_t, units="mins"))
    
    # Process results to pred and true values in a dataframe
    results = results_data(mod_list, gen_data_list, "fgfosr")
    
    # pull bias, coverage and MSE
    comp_mets = bias_coverage_mse(results)
    
    # aggregate it to put in a final df
    bin_out_df = params[i,] %>% mutate(bin_count = r, bias = comp_mets[1], coverage = comp_mets[2], MSE = comp_mets[3], runtime_min = fgfosr_time_diff)
    
    #append to a final output dataset
    final_results_df= rbind(final_results_df, bin_out_df)
    
    #make unique name and save results to a df for every bin run
    save_df = paste0("N_", params$N[i], "J_",params$J[i], "bc_", r,".rds")
    saveRDS(results, save_df)
    
    time_df = paste0("N_", params$N[i], "J_",params$J[i], "bc", r,"_summary.rds")        
    saveRDS(bin_out_df, time_df)

    #keep memory clear
    rm(results)
    
  }
  
}
final_results_df %>% mutate("num_iters" = N_iter)

write.csv(final_results_df, "fgfosr_output.csv")


final_results_df = data.frame("seed_start" = integer(), "N" = integer(), "J"=integer(), "bias" = character(), "coverage" = character(), "MSE" = character(), runtime_min = character())


for(i in 1:nrow(params)){
  print(paste0("N: ", params$N[i], " J: ",params$J[i]))
  
  N = params$N[i]
  J = params$J[i]
  
  gen_data_list = lapply(1:N_iter, function(x) gen_data(N = N, J = J, x))
  
  pffr_start_t = Sys.time() #Time that sucka
    
  ## PFFR
  lapply(gen_data_list, function(x) get_pffr_mod(x, N = N, J = J))
    
  pffr_end_t = Sys.time() #The cake is done!
  
  
  pffr_time_diff =  as.double(difftime(pffr_end_t, pffr_start_t, units="mins"))
  
  # Process results to pred and true values in a dataframe
  results = results_data(pffr_mod_list, gen_data_list, "pffr")
  
  # pull bias, coverage and MSE
  comp_mets = bias_coverage_mse(results)
  
  # aggregate it to put in a final df
  bin_out_df = params[i,] %>% mutate(bias = as.character(comp_mets[1]), coverage = as.character(comp_mets[2]), MSE = as.character(comp_mets[3]), runtime_min = as.character(pffr_time_diff))
  
  #append to a final output dataset
  final_results_df= rbind(final_results_df, bin_out_df)
  
  #make unique name and save results to a df for every bin run
  save_df = paste0("N_", params$N[i], "J_",params$J[i],"_pffr.rds")
  saveRDS(results, save_df)
  
  time_df = paste0("N_", params$N[i], "J_",params$J[i], "bc_pffr","_summary.rds")        
  saveRDS(bin_out_df, time_df)
  
  #keep memory clear
  rm(results)
    
}
  

final_results_df = final_results_df %>% mutate("num_iters" = N_iter)

write.csv(final_results_df, "pffr_output.csv")

