

 gen_data = function(N, J, run_num){
  
  ######################################################
  # Input: N - number of subjects                      #
  #        J - number of observations for each subject #
  #                                                    #
  # Output: list  [1] simulated data  [2] true data    #
  ######################################################
  
  #number of eigenvectors
  K =4
  
  # functional domain, equally spaced grid on 0,1
  sind <- seq(0,1,len=J)
  
  # create dataframe of true values at each observation point
  true_val = data.frame(sind = sind, true_val = f_0(sind)) %>%
    mutate(X = 0) %>%
    rbind(., data.frame(sind = sind, true_val = f_1(sind), X = 1)) %>%
    mutate(run = run_num)
  
  # Simulate X values
  X_i <- rbinom(N, size=1, prob=0.5)
  
  ## random effects
  # eigenfunctions \phi_k(s) evaluated on the observed grid
  phi <- sqrt(2)*cbind(sin(2*pi*sind),cos(2*pi*sind),
                       sin(4*pi*sind),cos(4*pi*sind))
  # eigenvalues \lambda
  lambda <- 0.5^(0:(K-1))
  # subject-specific weights/coefficients \xi_ik
  xi <- matrix(rnorm(N*K),N,K);
  xi <- xi %*% diag(sqrt(lambda))
  b_i <- xi %*% t(phi); # of size N by J
  
  ## linear predictor \eta_i(s)
  eta_i <- t(vapply(1:N, function(x){
    f_0(sind) + f_1(sind)*X_i[x] + b_i[x,]
  }, numeric(J)))
  
  ## response Y_i(s) \sim Binom(p = expit(eta_i))
  Y_i <- matrix(rbinom(N*J, size=1, prob=plogis(eta_i)), 
                N, J, byrow=FALSE)
  
  ## create data matrix in long format for estimation 
  df <- data.frame(id = rep(1:N, each=J),
                   sind = rep(sind, N), 
                   Y = as.vector(t(Y_i)), 
                   X = rep(X_i, each=J),
                   run = run_num)
  
  out =list(df, true_val)
  
  return(out)
}


get_coefs = function(midp, data){
  
  ##################################################################
  # Input: bin_seg - list of upper, lower and midpoint of each bin #
  #        data    - simulated data                                #
  #                                                                #
  # Output: dataframe with bin estimates for each id and errors    #
  ##################################################################
  
  # function to fit glmm in a bin width and return the mean at the midpoint of the bin
  mod = glmer(data = data ,Y ~ X + (1|id), family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  
  warn = ifelse(length(mod@optinfo$conv$lme4$messages[[1]])==0, "None", mod@optinfo$conv$lme4$messages[[1]])
  
  return( data.frame("id" = rownames(coef(mod)$id), "eta_i"= coef(mod)$id[[1]] + ranef(mod)[[1]]$"(Intercept)", "s_m" = midp, "error" = warn, row.names = NULL) )
  
}



get_fgfosr_mod = function(df_list, bin_count){
  
  ################################################################################################################
  # Input: df        - simulated dataframe with [1] simulated data and [2] true values at each observation point #
  #        bin_count - number of bins to split the domain into                                                   #
  #                                                                                                              #
  # Output: list, [1] fit model  [2] domain (sind)  [3] run nubmer                                               #
  ################################################################################################################
  
  
  #isolate the data from the list of dataframes passed in
  df = df_list[[1]]
  
  # Step 1
  #bin count & size
  bin_size = (max(df$sind)-min(df$sind))/bin_count
  
  #get a list of bin boundaries
  bin_bound = seq(min(df$sind), max(df$sind), bin_size)
  mid_points = rowMeans(embed(bin_bound, 2))
  bin_groups = data.frame(lb = head(bin_bound, -1), ub = tail(bin_bound, -1), mp = mid_points)
  bin_seg_list = as.list(as.data.frame(t(bin_groups)))
  
  # Step 2
  #fit glmm data in each bin with get_coef function, pass in filtered data to reduce memory overhead when running in parallel
  fit_df = bind_rows(lapply(bin_seg_list, function(x) get_coefs(x[[3]], df[df$sind >= x[[1]] & df$sind <= x[[2]],])))

  
  # Step 3
  #get number of ids to make a wide matrix of the fit_df
  N = length(unique(fit_df$id))
  
  K = 14
  
  #get eigen functions
  fpca_latent = fpca.face(matrix(fit_df$eta_i, N, bin_count, byrow=FALSE), pve=0.95, argvals=mid_points, knots= K,lower=0)
  
  #if there is only 1 remaining eigenfunction, need to adjust to before sending to reeval_efunctions
  efun = if (is.null(ncol(fpca_latent$efunctions))) matrix(fpca_latent$efunctions) else fpca_latent$efunctions
  
  ## project latent eigen functions onto the original grid
  ef_ogrid = fastGFPCA:::reeval_efunctions(knots = K, argvals_bin = mid_points, argvals = unique(df$sind), efunctions =  efun, npc = ncol(efun)) %>%
    data.frame() %>%
    mutate(sind = unique(df$sind))
  
  df$id = as.factor(df$id)
  
  df_w_ef = left_join(df, ef_ogrid, by = "sind")
  
  ### Dynamically change column names based on number of eigenfunctions
  colnames(df_w_ef) = gsub(x =colnames(df_w_ef), pattern = "\\.", replacement = "phi1")  
  colnames(df_w_ef) = gsub(x =colnames(df_w_ef), pattern = "X(?=[0-9])", replacement = "phi", perl=T)
  
  phi_num = sum(startsWith(colnames(df_w_ef), "phi"))
  
  ### Dynamically change gam formula names based on number of eigenfunctions
  add_ef = dplyr::case_when(
    phi_num == 0 ~ "",
    phi_num != 0 ~ paste0(lapply(1:phi_num, function(x) paste0("s(id, by = phi", x, ", bs=\"re\")")), collapse=" + "))
  
  
  dynamic_form = formula(paste0("Y ~ s(sind, k=30, bs = 'cr') + 
                          s(sind, by=X, k = 30, bs = 'cr') + ", add_ef))
  
  
  #### fit gam with formula
  gam_mod = mgcv::bam(dynamic_form, method = "fREML", data = df_w_ef, discrete = T,  family=binomial, gc.level = 1)
  
  #tie simulated data run number to the model fit
  run = unique(df$run)
  
  cov_err_df = fit_df %>% filter(error != "None") %>% mutate(run = run)
  
  gam_mod_list = list(gam_mod, df_w_ef$sind, run, cov_err_df)
  
  return(gam_mod_list)
  
}


gfosr_get_pred = function (gam_mod) {
  
  ###############################################################################
  # Input: gam mod - takes in list with [1] gam mod, [2] domain, [3] run number #
  #                                                                             #
  # Output: dataframe with sind, predicteve values, se, and X indicator         #
  ###############################################################################
  
  sind_t = unique(gam_mod[[2]])
  
  form = toString(gam_mod[[1]]$formula)
  
  #vary base prediction dataframe based on presense of number of eigencvectors from the gam mod
  df_pred_f1=case_when(
    grepl("phi4", form, fixed = TRUE) ~ "data.frame(sind = sind_t, phi1 = 0, phi2 = 0, phi3 = 0, phi4 = 0,id = 1, X = 1)",
    grepl("phi3", form, fixed = TRUE) ~ "data.frame(sind = sind_t, phi1 = 0, phi2 = 0, phi3 = 0,id = 1, X = 1)",
    grepl("phi2", form, fixed = TRUE) ~ "data.frame(sind = sind_t, phi1 = 0, phi2 = 0,id = 1, X = 1)",
    grepl("phi1", form, fixed = TRUE) ~ "data.frame(sind = sind_t, phi1 = 0, id = 1, X = 1)",
    TRUE ~ "data.frame(sind = sind_t, id = 1 X = 1)",
  )
  
  # excutes dynamic dataframe creation from previous step
  df_pred_f1 = eval(parse(text = df_pred_f1))
  
  #predicts based on game mod and created df
  pred_fit_f1 = predict(gam_mod[[1]], newdata=df_pred_f1, se.fit=TRUE, type = "terms")
  
  #f*0 has a constraint of a constant, f_0(s) = \beta_0 + f*_0(s)
  b_0 = coef(gam_mod[[1]])[1]
  f0_star = data.frame(pred_fit_f1$fit)[[1]]
  f0 = b_0+f0_star
  
  f1 = data.frame(pred_fit_f1$fit)[[2]]
  
  pred_x0 = data.frame(sind = sind_t, pred_val = f0, X=0, se = data.frame(pred_fit_f1$se.fit)[[1]])
  pred_x1 = data.frame(sind = sind_t, pred_val = f1, X=1, se = data.frame(pred_fit_f1$se.fit)[[2]])
  
  pred_vals = rbind(pred_x0, pred_x1) %>%
    mutate(run = unique(gam_mod[[3]]))
  
  return(pred_vals)
}


results_data = function(mod_list, df_list, meth){
  
  #################################################################################################
  # Input: gam mod - takes in [1] fgfosr list [2] generated data list                             #
  #                                                                                               #
  # Output: dataframe with joined pred and true values with bias, CI, and if pred value is in CI  #
  #################################################################################################

  #Use get_pred to get predicted values
  if(meth == "fgfosr"){
    pred_dfs = lapply(mod_list, function(x) gfosr_get_pred(x))  
  }else{
    pred_dfs = lapply(mod_list, function(x) pffr_get_pred(x))  
  }
  
  
  #combine the predicted values for each run into 1 df and add a run identifier column
  all_pred_data = bind_rows(pred_dfs) %>%
    mutate(run = as.integer(run), X = as.factor(X))
  
  ############### Organize True data into 1 df and combine with predicted values by run
  true_val_list = do.call(rbind, lapply(df_list, function(x) x[[2]])) %>%
    mutate(run = as.integer(run), X = as.factor(X))
  
  
  # Calculate bias and CIs 
  all_data = full_join(all_pred_data, true_val_list, by = c("run", "sind", "X")) %>%
    mutate(bias = true_val - pred_val) %>%  #calculate bias
    mutate(CI_U = pred_val + 1.96*se, CI_L = pred_val - 1.96*se) %>%
    mutate(CI_check = ifelse(true_val <= CI_U & true_val >= CI_L, 1, 0))
  
  
  return(all_data)
}


bias_coverage_mse = function(df){
  
  #####################################
  # Input: Reuslts dataframe          #
  #                                   #
  # Output : bias, coverage and MSE   #
  #####################################
  
  bias_out = mean(df$bias)
  
  coverage_out = df %>% 
    group_by(run) %>%
    summarize(run_cov = mean(CI_check)) %>%
    summarise(cov = mean(run_cov)) %>%
    pull(cov)
  
  MSE = df %>%
    summarize( mean(true_val - pred_val)^2) %>%
    pull()
  
  out = c(bias_out, coverage_out, MSE)
  return(out)
}







############################ pffr functions


get_pffr_mod = function(df_list, N, J){
  
  ############################################
  # Input: simulated data, N, J              #
  #                                          #
  # Output : pffr model, domain (sind), run  #
  ############################################
  
  df = df_list[[1]]
  
  #pffr setup
  run = unique(df$run)
  sind = unique(df$sind)
  df_pffr <- data.frame(Y = I(matrix(df$Y, N, J, byrow=TRUE)), 
                        id = factor(1:N), 
                        X = df$X[!duplicated(df$id)])
  
  # Penalized flexible functional regression
  fit_pffr <- pffr(Y ~ X + s(id, bs="re"), 
                   # use mgcv::bam with fastREML smoothing parameter selection to estimate the model
                   algorithm="bam", method="fREML", discrete=TRUE,
                   # specify the bases used for estimating f_0(s), f_1(s) via bs.int and bs.yindex, respectively
                   bs.yindex=list(bs="cr", k=30), bs.int=list(bs="cr",k=30), data=df_pffr,
                   # specify outcome distribution and the functional domain (yind)
                   family="binomial", yind=sind)
  
  
  out = list(fit_pffr, sind, run)
  
  return(out)

}


pffr_get_pred = function(mod_list){
  
  #################################################################################
  # Input: pffr mod - takes in list with [1] pffr mod, [2] domain, [3] run number #
  #                                                                               #
  # Output: dataframe with sind, predicteve values, se, and X indicator           #
  #################################################################################
  
  sind = mod_list[[2]]
  fit_pffr = mod_list[[1]]
  run = mod_list[[3]]
  
  ## extract estiamted coefficients
  # this is slightly different than when using mgcv directly as refund::pffr has it's own naming convention for variables
  df_pred <- data.frame("sind.vec"=sind, X=1, id=1)
  
  # note that we need to explicity call the predict method from the mgcv package
  # (refund::pffr has it's own predict method that doesnt do exactly what we want it to)
  beta_hat <- mgcv::predict.gam(fit_pffr, newdata=df_pred, type='iterms', se.fit=TRUE)
  
  # add back constant removed from f_0(s) due to identifiability constraint 
  beta_hat$fit[,1] <- beta_hat$fit[,1] + attributes(beta_hat)$constant
  
  # put results in a data frame
  results_df <- data.frame("pred_val" = as.vector(beta_hat$fit[,1:2]), 
                           "se" = as.vector(beta_hat$se.fit[,1:2]), 
                           "sind" = rep(sind, 2), 
                           "X" = rep(c("0", "1"), each=length(sind))) %>%
                  mutate(run = run)
  
  return(results_df)
  
  
}


