## Simulations for univariate classical times series
## comparing rat of estimator corrections, MSE, bandwidth selection rules, and
## testing performance (Type I Error)
here::i_am("R/univariate_sims.R")
#proj_path <- here::here()
#setwd("./R")
library(Matrix)
library(MASS)
library(distr)

setwd("./R")
source("kernels.R")
source("lugsail.R")
source("get_cv.R")
#source("R/get_cv.R") #
#source("get_b.R") # Has an error right now, Set b = 0.1 for now!
#source("estimate_LRV.R")


#### Helper Functions ####

## Good old OLS
ols <- function(y, X) {
  # Add a column of ones to the predictor matrix for the intercept term
  #X <- cbind(1, X) # Add intercept

  # Calculate the OLS estimates using the formula (X'X)^(-1)X'y
  beta <- solve(t(X) %*% X) %*% t(X) %*% y

  # Return the estimated coefficients
  return(beta)
}

## Get OLS Errors/Residuals
ols_errors <- function(y, X) {
  # Add a column of ones to the predictor matrix for the intercept term
  X <- cbind(1, X)

  # Calculate the OLS estimates using the formula (X'X)^(-1)X'y
  beta <- solve(t(X) %*% X) %*% t(X) %*% y

  # Calculate the errors/residuals
  errors <- y - X %*% beta

  return(errors)
}

# AUTOCOVARIANCE #

#' Function to calculate the autocovariance of a multivariate time series
#'
#' @param Y matrix of the Y-values of the time-series
#' @param h The lag at which to calculate the autocovariance
#' @return Returns the value of the autocovariance at lag h
autocovariance_multivariate <- function(Y, h) {
  big_T <- nrow(Y)
  k <- ncol(Y)
  Y_mean <- colMeans(Y)

  # Initialize the autocovariance matrix
  gamma_h <- matrix(0, nrow = k, ncol = k)

  # Calculate the autocovariance for lag h
  for (i in 1:(big_T-h)) {
    gamma_h <- gamma_h + (Y[i, ] - Y_mean) %*% t(Y[i+h, ] - Y_mean)
  }

  #not for lag 0
  if(h==0){
    gamma_h <- gamma_h / big_T
  }
  else{
    # add transpose for h > 1
    gamma_h <- 2 * gamma_h / big_T  
  }
  

  return(gamma_h)
}


##### Simulating Data Functions #####

#' Generates AR(1) correlated time series X values, of [d]-dimension, with
#' autocorrelation parameter [phi] with Gaussian noise.
#'
#' @param big_T numeric, Length of the time series
#' @param phi numeric, Autocorrelation coefficient
#' @param d numeric, Dimension and number of columns of X of the time series to create
#' @return Return a matrix X that has dimension big_T by d
AR_X <- function(big_T, phi, d){
  # Make a random normal draw
  X <- matrix(0, nrow = big_T, ncol = d)

  for(i in 2:big_T){
    X[i, ] <- phi*X[(i-1), ] + mvrnorm(n = 1, rep(0, d), diag(d))
  }
  # Return matrix of dim big_T by d
  return(X)
}

#' AR(1) correlated error/noise generation with autocorrelation parameter [phi]
#'
#' @param big_T numeric, Length of time series
#' @param phi numeric, Autocorrelation coefficient
#' @return Returns a numeric vector of length big_T
AR_u <- function(big_T, phi){
  u <- rep(0, big_T)
  for(i in 2:big_T){
    u[i] <- phi*u[(i-1)] + rnorm(1)
  }
  return(u)
}

#' Generate Sinusoidal error/noise vector (not AR)
#'
#' @param big_T numeric, Length of time series
#' @param phi numeric, Autocorrelation coefficient
#' @return Returns a numeric vector of length big_T
sine_u <- function(big_T, phi) {
  u <- rep(0, big_T)
  for(i in 2:big_T){
    #u[i] <- rnorm(1, mean = 0, sd = 1 + 0.5 * sin(2*pi*i / 50))
    u[i] <- 1 + 0.5 * sin(2*pi*i / 50) + rnorm(1)
  }
  return(u)  # Example: sinusoidal variation
}

#' AR(1) Heteroskedastic correlated observations [theta]%*%[AR_X]
#' plus AR(1) correlated errors[AR_u]:
#' Y = X %*%theta + u
#'
#' @param big_T numeric, Length of time series
#' @param phi numeric, Autocorrelation coefficient
#' @param d numeric, Dimension and number of columns of X of the time series to create
#' @param theta numeric vector length d, default is 0's, modifies the weight to AR_X
#' @return Returns a list of response Y  and X
AR1_AR_u <- function(big_T, phi, d, theta = rep(0, d)){
  X <- AR_X(big_T, phi, d)
  u <- AR_u(big_T, phi)
  Y <- X%*%theta + u
  return(list(Y = c(Y), X = X))
}


#' AR(1) Heteroskedastic correlated observations [theta]%*%[AR_X]
#' plus AR(1) correlated errors[AR_u]:
#' Y = X %*%theta + u
#'
#' @param big_T numeric, Length of time series
#' @param phi numeric, Autocorrelation coefficient
#' @param d numeric, Dimension and number of columns of X of the time series to create
#' @param theta numeric vector length d, default is 0's, modifies the weight to AR_X
#' @return Returns a list of response Y  and X
AR1_SINE <- function(big_T, phi, d, theta = rep(0, d)){
  X <- AR_X(big_T, phi, d)
  u <- sine_u(big_T, phi)
  Y <- X%*%theta + u
  return(list(Y = c(Y), X = X))
}


#### Examples ####

## Test autocovariance function
# Generate some example data
set.seed(123)
big_T <- 100 # length of time series
k <- 3 # number of columns
Y <- matrix(rnorm(big_T * k), big_T, k)

# Calculate the autocovariance for lag h = 1
gamma_h <- autocovariance_multivariate(Y, h = 1)
print(gamma_h)


###
## Using correlated observation and errors (univariate response Y)
###


# Generate data
set.seed(1234)
phi = 0.9
d = 4
data_sine <- AR1_SINE(big_T = big_T, phi = phi, d = d)
gamma_h <- autocovariance_multivariate(Y = as.matrix(data_sine$Y), h = 1)
print(gamma_h)

# LRV Estimate (Omega-hat)
#  SV estimator utilizes a kernel function k() and a bandwidth parameter b
#  to estimate the LRV via weighted average of the sample auto-covariances Gamma-hat(h)
#  Omega-hat = (1/T sum x%*%x')^(-1) %*%  \sum_{h=-(T-1)}^{T-1} k(h/(bt))*Gamma-hat(h)   %*% (1/T sum x%*%x')^(-1)

# Note: when Y is VAR, (multiple Y's), this will need to be modified to sum over
# the autocovariances I believe.
Y <- data_sine$Y
X <- data_sine$X
big_T <- nrow(X) # Sample Size
M <- solve(t(X) %*% X / big_T) # Data Term, d by d
b = 0.1

# Have to weight the autocovariances by kernel
all_autocovariances <- sapply(0:(big_T-1),
                                  autocovariance_multivariate,
                                  Y = as.matrix(data_sine$Y))
# Plot autocovs
plot(all_autocovariances~ c(1:100), 
     main = "Plot of Autocovariance vs Time",
     xlab = "Time", 
     ylab = "Auto Cov", 
     col = "orange", 
     pch = 16)


# Implementing bandwidth rules/functions to apply with 'all_autocovariances'

LRV_mother_estimator <- function(b, all_autocovariances, the_kernel){
  big_T <- nrow(as.matrix(all_autocovariances))
  
  # Make the weights that correspond to the autocovariances
  # Each row corresponds to the weights for a specific weight.
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios.
  W <- sapply(b, function(b){
    if(b == 0){
      new_weights <- c(1, rep(0, c(big_T-1)))
    } else{
      M <- b*big_T
      new_weights <- sapply(0:(big_T -1)/(M), the_kernel)
    }
  })
  W <- t(W)
  rownames(W) = paste("b=", round(b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(nrow(as.matrix(all_autocovariances))-1), sep = "")
  
  
  # [#, ] the b value
  # [ , #]  a component of the estimated LRV
  #    Omega_11, Omega_12,  ..., Omega_1d; ...; Omega_d1, ...Omega_dd
  omega_hats <- W %*% all_autocovariances
  rownames(omega_hats) <- paste("b=", round(b, 3))
  
  return(omega_hats)
}

# Apply kernel to autocovariances
# kernels: bartlett, qs, parzen, th
# the_kernel <- bartlett
# the_kernel_string <- 'bartlett'
## Do LRV_estimator to use bandwidth
#weighted_autocov <- sapply(X = all_autocovariances, FUN = the_kernel)

# Plot weighted autocov
# lines(weighted_autocov~ c(1:100), 
#      col = "dodgerblue", 
#      pch = 16)

# plot(weighted_autocov~ c(1:100), 
#      main = "Plot of Weighted Autocovariance vs Time",
#      xlab = "Time", 
#      ylab = "Auto Cov", 
#      col = "dodgerblue", 
#      pch = 16)

# Lugsail Kernel test
# When r = 1, c = 0, lugsail == mother kernel

# lugsail_parameters <- list(r = 1, c = 0)
# lugsail_parameters <- list(r = 0.5, c = 0.5)
# lugsail_parameters <- get_lugsail_parameters(big_T = big_T, q = 1, method = "Zero")
# 
# weighted_autocov <- sapply(X = all_autocovariances, FUN = lugsail, lugsail_parameters = lugsail_parameters, 
#                            the_kernel = the_kernel)
# 
# # Plot autocov and weighted autocov
# plot(all_autocovariances ~ c(1:100), 
#      main = "Plot of Autocovariance vs Time \n (Zero Lugsail)",
#      xlab = "Time", 
#      ylab = "Auto Cov", 
#      col = "orange", 
#      pch = 16)
# lines(weighted_autocov ~ c(1:100),
#       col = "dodgerblue")


# Construct omega-hat (LRV)
# omega_hat <- sum(weighted_autocov)
# cat("Kernel:", the_kernel_string)
# cat("omega hat:", omega_hat)
# cat("sqrt(omega-hat)", sqrt(omega_hat))

# Covariance of model parameters
# sigma_covariance = solve(M) %*% kernel(h, b, big_T) * autocovariance(lag h) %*% solve(M)
# sigma_covariance <- solve(M) * sum(weighted_autocov) * solve(M)
# print(sigma_covariance)


#### Type 1 Error Example ####
set.seed(1234)
nsim <- 1000    # Number of simulations for Type errors
alpha <- 0.05   # Significance level
big_T <- 1000    # Time Series length
phi = 0 # 0.9       # Autocorrelation parameter
d = 1           # X dimension (univariate Y for now)

# True parameter values for testing
null_means = rep(0, d)     #  add d+1 for intercept?

type1_vec <- rep(NA, nsim)   # Store type1 freq
test_stats_vec <- rep(NA, nsim) # store stats


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
  b = 0.1     # There exist better rules
  
  # Have to weight the autocovariances by kernel
  all_autocovariances <- sapply(0:(big_T-1),
                              autocovariance_multivariate,
                              Y = as.matrix(X * Y)) # ERRORS! # X * (Y - X %*% theta_hat))
  
  # Select Kernel(s)
  the_kernel <- bartlett    # Define kernel function
  the_kernel_string <- 'bartlett' # String for printing
  
  # Apply kernel to autocovariances (matrix may be more involved) 
  #weighted_autocov <- sapply(X = all_autocovariances, FUN = the_kernel)
  
  # Lugsail Kernel test
  # When r = 1, c = 0, lugsail == mother kernel
  lugsail_parameters <- list(r = 1, c = 0)
  #lugsail_parameters <- list(r = 0.5, c = 0.5)
  #lugsail_parameters <- get_lugsail_parameters(big_T = big_T, q = 1, method = "Zero")
  
  # weighted_autocov <- sapply(X = all_autocovariances, FUN = lugsail,
  #                            lugsail_parameters = lugsail_parameters, 
  #                            the_kernel = the_kernel)
  
  # Apply bandwidth and kernel rules to get LRV (Mother kernels)
  omega_hat_LRV <- LRV_mother_estimator(b = 0.1, all_autocovariances = all_autocovariances, the_kernel = the_kernel)
  
  
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
  
  # d = 1
  t_stat <- diff/sqrt(sigma_covariance/big_T) # DIVIDE BY SIGMA_COVARIANCE
  test_stats_vec[i] <- t_stat
  # record omegahats, 
  
  
  
  # VAR - multiple Ys
  # Issue here, dim(omega_hat_LRV) = square p by p, 
  # so t(diff) %*% omega^-1 %*% diff doesn't work but that can't be right.
  
  # F Stat
  #F_stat = big_T * solve(omega_hat_LRV) * t(diff)  %*% diff / d
  
  
  # Square for chi-square.. Use F-stat
  test_statistic <- t_stat ** 2 # Chi-square for p = 1?
  #test_statistic <- F_stat
  
  # CRITICAL VALUE
  critval <- get_cv(new_b = b, d = d, alpha = alpha, 
                    the_kernel =  "Bartlett",
                    lugsail = "Mother",
                    method = "analytical")
  
  if(i %% 100 == 0){cat("(n = ", i, ") Test Statistic: ", test_statistic, 
                        "   Crit Value: ", critval, "\n")}
  
  
  # Compare critical value and test statistic 
  # But for t, there's one for each component of X. 
  type1_error <- ifelse(test_statistic < critval, 0, 1) # 1 if reject H0
  type1_vec[i] <- type1_error # Store 0's, 1's
  
  
  # Printing/debug
  # cat("Kernel:", the_kernel_string)
  # cat("omega hat:", omega_hat)
  # cat("sqrt(omega-hat)", sqrt(omega_hat))

}

# Type 1 error
mean(type1_vec)




# Plot data
plot(data_homo$Y ~ c(1:big_T), type = "l",
     main = paste0("AR1_HOMO Simulated Data, phi=", phi),
     xlab = "Time", ylab = "Response Y")

# # Plot data
# plot(data_sine$Y ~ c(1:big_T), type = "l",
#      main = paste0("AR1_HOMO Simulated Data, phi=", phi),
#      xlab = "Time", ylab = "Response Y")


# ## Easy version ## 
# #### Testing Type 1 Error ####
# 
# # Set parameters
# set.seed(123)
# n <- 30        # Sample size
# mu0 <- 0       # Null hypothesis mean
# mu1 <- 1       # Alternative hypothesis mean
# sigma <- 1     # Standard deviation
# 
# nsim <- 10000  # number of simulations
# type1_vec <- rep(NA, nsim) # Store type1 freq
# type2_vec <- rep(NA, nsim) # Store type2 freq
# 
# for(i in 1:nsim){
#   # Simulate data
#   data_null <- rnorm(n, mean = mu0, sd = sigma)
#   data_alt <- rnorm(n, mean = mu1, sd = sigma)
# 
#   # Perform t-test
#   t_test_null <- t.test(data_null, mu = mu0)
#   t_test_alt <- t.test(data_alt, mu = mu0)
# 
#   # Type I error
#   type1_error <- ifelse(t_test_null$p.value < alpha, 1, 0)
#   type1_vec[i] <- type1_error
# 
#   # Type II error
#   type2_error <- ifelse(t_test_alt$p.value >= alpha, 1, 0)
#   type2_vec[i] <- type2_error
# 
# }
# # Type 1 error
# type1 <- mean(type1_vec)
# 
# # Type 2 error
# type2 <- mean(type2_vec)
# 
# # Power
# power <- 1 - mean(type2_vec)
# 
# # Results
# cat("Type I Error:", type1, "\n")
# cat("Type II Error:", type2, "\n")
# cat("Power:", power, "\n")
# 
