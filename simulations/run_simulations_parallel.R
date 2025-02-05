pacman::p_load(tidyverse, tictoc, future, furrr)

source(here::here("source", "01_simulate_data.R"))
source(here::here("source", "02_apply_methods.R"))


###############################################################
## set simulation design elements
###############################################################

# Parameters for 12 scenarios
nsim = 475
n = c(10, 50, 500)
beta_true = c(0, 0.5, 2)
epsilon_distr = c("Normal", "logNormal")
CI_method = c("Wald", "boot quantile", "boot t")

params.all = expand.grid(
  n = n,
  beta_true = beta_true,
  epsilon_distr = epsilon_distr,
  CI_method = CI_method
)

params.all$scenario_id = 1:nrow(params.all)

###############################################################
## start simulation code
###############################################################

# Define simulation function for each scenario
simulate_scenario <- function(params, seeds, nsim, 
                              nboot , nboot_t , nboot_se ) {
  # Generate seeds for reproducibility
  results <- vector("list", nsim)
  
  for (i in seq_len(nsim)) {
    set.seed(seeds[i])
    
    # ************************************
    # Simulate data  *********************
    simdata <- get_simdata(
      n = params$n,
      beta_treat = params$beta_true,
      epsilon_distr = params$epsilon_distr
    )
    
    # **************************************
    # Apply method(s)  *********************
    model <- fit_model(simdata, params$epsilon_distr)
    est_beta <- model$coefficients[["x"]]
    est_se <- summary(model)$coefficients["x", "Std. Error"]
    
    
    ## Wald
    if(params$CI_method == "Wald"){
      tic()
      est_out <- standard_wald(simdata, params$epsilon_distr) 
      time_stamp = toc(quiet = TRUE)
    }
    
    ## boot quantile
    if(params$CI_method == "boot quantile"){
      tic()
      est_out <- boot_quantile(simdata, params$epsilon_distr, nboot)
      time_stamp = toc(quiet = TRUE)
    } 
    
    ## boot t
    if(params$CI_method == "boot t"){
      tic()
      est_out <- boot_t(simdata, params$epsilon_distr, est_beta, est_se, nboot_t, nboot_se)
      time_stamp = toc(quiet = TRUE)
    } 
    
    # record time 
    cal_CI_time = time_stamp$toc - time_stamp$tic
    
    
    
    # ******************************************
    # summarise estimates  *********************
    
    # Store results
    result <- data.frame(params, seed = seeds[i], cal_CI_time, est_beta, est_out) |> 
      mutate(error = est_beta - beta_true, 
             coverage = ifelse(beta_true >= conf.low & beta_true <= conf.high, 1, 0))
    
    results[[i]] = result
  }
  
  # Combine results of nsim for each scenario
  results_dt <- do.call(rbind, results)
  
  # Save results for the current scenario
  filename <- paste0("scenario_", params$scenario_id, ".Rdata")
  save(results_dt, file = here::here("data", filename))
  
  return(results_dt)
}




###############################################################
##  parallel calculation for each simulation scenario
###############################################################

# ****************************************************
#  set parallel cores         ************************
plan(multisession, workers = parallel::detectCores() - 1)

# set seeds
set.seed(1)
seeds = sample(1:10000, nsim) 

# for boot quantile, nboot = 500 
# for boot t interval, nboot_t = 500, and nboot_se = 200 for nested bootstrap
nboot = 500
nboot_t = 500
nboot_se = 200

# ****************************************************
# Run simulations for all scenarios in parallel ******

results_parallel <- future_map(
  seq_len(nrow(params.all)),
  function(scenario_id) {
    params <- params.all[scenario_id, ]
    simulate_scenario(params, seeds, nsim, nboot, nboot_t, nboot_se)
  }, .options = furrr_options(seed = TRUE)
)


####################
# save results
#
# Combine results from all scenarios if needed
final_results <- do.call(rbind, results_parallel)
rownames(final_results) <- NULL

# Save combined results
save(final_results, file = here::here("data", "all_scenarios.Rdata"))

