# Wald confidence intervals
standard_wald = function(simulated_data){
  model_fit = lm(y ~ x, data = simulated_data)
  broom::tidy(model_fit, conf.int = TRUE) |> 
    filter(term == "x") |> 
    rename(est_se = std.error) |> 
    select(est_se, conf.low, conf.high) 
}



# Nonparametric bootstrap percentile intervals
boot_quantile = function(simulated_data, nboot, alpha = 0.05){
  boot_beta <- rep(NA, nboot)
  sample_size <- nrow(simulated_data)
  
  ## bootstrap ***********************
  for(b in 1:nboot){
    # get the sample for each bootstrap
    sample_index = sample(1:sample_size, size = sample_size, replace = TRUE)
    boot_sample = simulated_data[sample_index, ]
    
    # run linear regression
    model = lm(y ~ x, data = boot_sample)
    
    # get estimated beta
    boot_beta[b] = model$coefficients[["x"]]
  }
  
  ## output    ***********************
  data.frame(est_se = sd(boot_beta, na.rm = T),
             conf.low = quantile(boot_beta, alpha/2, na.rm = T),
             conf.high = quantile(boot_beta, 1 - alpha/2, na.rm = T)
  )
}





# Nonparametric bootstrap t intervals
boot_t = function(simulated_data, est_beta, est_se, nboot_t, nboot_se, alpha = 0.05){
  t_star <- rep(NA, nboot_t)
  sample_size <- nrow(simulated_data)
  
  ## bootstrap ***********************
  for(b in 1:nboot_t){
    # get the sample for each bootstrap
    sample_index = sample(1:sample_size, size = sample_size, replace = TRUE)
    boot_sample = simulated_data[sample_index, ]
    
    # run linear regression
    model = lm(y ~ x, data = boot_sample)
    
    # get estimated beta
    boot_beta = model$coefficients[["x"]]
    
    # get se_beta of each bootstrap using bootstrap
    # ! use boot_sample as input
    inner_boot <- boot_quantile(boot_sample, nboot_se, alpha = 0.05)
    boot_se <- inner_boot$est_se
    
    # get boot_t
    t_star[b] <- (boot_beta - est_beta)/boot_se
  }
  
  
  ## output    ***********************
  boot_t_low = quantile(t_star, alpha/2, na.rm = T)
  boot_t_up = quantile(t_star, 1 - alpha/2, na.rm = T)
  
  data.frame(est_se = est_se,
             conf.low = est_beta - boot_t_up*est_se,
             conf.high = est_beta - boot_t_low*est_se
  )
  
}


