# load libraries
library(tidyverse)

# Two distribution
get_simdata = function(n, beta_treat, epsilon_distr){
  beta0 = 1
  x = rbinom(n, 1, prob = 0.5)
  
  if(epsilon_distr == "Normal"){
    epsilon = rnorm(n, 0, sd = 2)
  }
  
  if(epsilon_distr == "logNormal"){
    epsilon = rlnorm(n, 0, sdlog = 2)
  }
  
  y = beta0 + beta_treat * x + epsilon
  
  tibble(
    x = x,
    y = y
  )
}



