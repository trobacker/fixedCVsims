## Functions for simulation use
here::i_am("R/helper_functions.R")
#proj_path <- here::here()
#setwd("./R")
library(Matrix)
library(MASS)
library(distr)

#setwd("./R")
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