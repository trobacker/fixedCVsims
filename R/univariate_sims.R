## Simulations for univariate classical times series
## comparing rat of estimator corrections, MSE, bandwidth selection rules, and
## testing performance (Type I Error)
library(Matrix)
library(MASS)
library(distr)
library(beepr)

repo_path <- here::here() 
setwd(paste0(repo_path, "/R")) 
source("helper_functions.R")
source("kernels.R")
source("lugsail.R")
source("get_cv.R")
source("get_b.R")

seed_one <- 1234

#### Type 1 Error Rates  ####

## Scenario 1 ##
set.seed(seed_one)
nsim <- 1000    # Number of simulations for Type errors
alpha <- 0.05   # Significance level
big_T <- 200  # Time Series length
#phi = 0 # 0.9       # Autocorrelation parameter
d = 1        # X dimension (univariate Y for now)

# True parameter values for testing
null_means = rep(0, d)     #  add d+1 for intercept?

phi_vec <- c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99)

type1_all <- rep(NA, length(phi_vec))
b_mean_vec <- rep(NA, length(phi_vec))

phi_index = 0

### Store the get_b values to explore how b's vary across settings!!

# Select Kernel(s)
#the_kernel <- bartlett    # Define kernel function
the_kernel <- parzen
#the_kernel <- qs
#the_kernel <- th

# Parzen, Bartlett, QS, TH
#the_kernel_string <- 'Bartlett' # String for printing
the_kernel_string <- 'Parzen'
#the_kernel_string <- 'QS'
#the_kernel_string <- 'TH'


for(phi in phi_vec){
  phi_index <- phi_index + 1
  
  type1_vec <- rep(NA, nsim)   # Store type1 freq
  test_stats_vec <- rep(NA, nsim) # store stats
  
  b_vec <- rep(NA, nsim)
  b_sd_vec <- rep(NA, nsim)
  
  for(i in 1:nsim){
    # Simulate data
    data_homo <- AR1_AR_u(big_T = big_T, phi = phi, d = d, theta = rep(0,d))
    Y <- data_homo$Y
    X <- data_homo$X
    big_T <- nrow(X) # Sample Size
    
    # The means (+ intercept? Not now)
    # consider autocor parameter vs. ols param estimates, if centered then always phi?
    theta_hat = ols(y = Y, X = X)
    
    # LRV Estimate
    # Bandwidth decision rule
    #b = 0.1     # There exist better rules
    
    # Such as:
    # the_data = both X and Y?
    # if it's not bartlett, just not "bartlett" (lowercase!)
    # X * Y instead of just X/Y
    b <- get_b(the_data = X*Y, the_kernel = tolower(the_kernel_string)) # INPUT is error terms
    b_vec[i] <- b
    
    # Have to weight the autocovariances by kernel
    all_autocovariances <- sapply(0:(big_T-1),
                                  autocovariance_multivariate,
                                  Y = as.matrix(X * Y)) # ERRORS! # X * (Y - X %*% theta_hat))
    
    # Apply kernel to autocovariances (matrix may be more involved) 
    #weighted_autocov <- sapply(X = all_autocovariances, FUN = the_kernel)
    
    # Lugsail Kernel test
    # When r = 1, c = 0, lugsail == mother kernel
    lugsail_parameters <- list(r = 1, c = 0)
    #lugsail_parameters <- list(r = 0.5, c = 0.5)
    
    # !!!
    #lugsail_parameters <- get_lugsail_parameters(big_T = big_T, q = 1, method = "Zero")
    
    # Apply bandwidth and kernel rules to get LRV (Mother kernels)
    omega_hat_LRV <- LRV_mother_estimator(b = b, all_autocovariances = all_autocovariances, the_kernel = the_kernel)
    
    
    # Construct omega-hat (LRV)
    #omega_hat <- sum(weighted_autocov)
    
    # Data term for covariance matrix of model
    M <- solve(t(X) %*% X / big_T )  # Data Term, d by d
    
    # Covariance matrix for model
    #sigma_covariance <- solve(M) * sum(weighted_autocov) * solve(M)
    sigma_covariance <- M * omega_hat_LRV[1] * M
    
    # TEST STATISTICs
    # !!!
    # See `second_regressor.R` line 405+: F_stat
    # Generally, an F_stat
    # Testing THETA's, not rho/phi autocorr.
    # NO: t_stat <- sqrt(big_T)*(var_model$ar[1] - phi) / (sqrt(omega_hat))
    
    # Here, d = 4, p = 1 (t-stat?)
    
    diff <- (theta_hat - null_means) # add +1 for intercept? 
    
    # F stat (YAY!)
    F_stat <- t(diff) %*% solve(sigma_covariance) %*% diff * big_T / d
    
    # d = 1
    #t_stat <- diff/sqrt(sigma_covariance/big_T) # DIVIDE BY SIGMA_COVARIANCE
    test_stats_vec[i] <- F_stat
    # record omegahats, 
    
    
    
    # VAR - multiple Ys
    # Issue here, dim(omega_hat_LRV) = square p by p, 
    # so t(diff) %*% omega^-1 %*% diff doesn't work but that can't be right.
    
    # Square for chi-square.. Use F-stat
    #test_statistic <- t_stat ** 2 # Chi-square for p = 1?
    test_statistic <- F_stat
    
    # CRITICAL VALUE
    critval <- get_cv(new_b = b, d = d, alpha = alpha, 
                      the_kernel =  the_kernel_string,
                      lugsail = "Mother",
                      method = "analytical")
    
    if(i %% 100 == 0){cat("(n = ", i, ") Test Statistic: ", test_statistic, 
                          "   Crit Value: ", critval, "\n")}
    
    
    # Compare critical value and test statistic 
    # But for t, there's one for each component of X. 
    type1_error <- ifelse(test_statistic < critval, 0, 1) # 1 if reject H0
    type1_vec[i] <- type1_error # Store 0's, 1's

  }

  # Store type 1 error rate for each phi
  type1_all[phi_index] <- mean(type1_vec) 
  b_mean_vec[phi_index] <- mean(b_vec)
  b_sd_vec[phi_index] <- sd(b_vec)
  
}
beep(sound = "complete") 


## Some output for reference


## get_b warnings.. check in high/low correlation settings

# big_T = 200
# Kernel = bartlett using X*Y (errors) in get_b()
# d = 1
# phi_vec <- c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99)
# type1_all
# [1] 0.049 0.054 0.058 0.059 0.055 0.097 0.408 0.522 0.646 0.741 0.826
# b_mean_vec
# [1] 0.011380 0.017760 0.040975 0.095530 0.160010 0.275110 0.165565 0.121400
# [9] 0.080750 0.034785 0.018305



