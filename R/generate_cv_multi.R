# It can also be cleaned/sped up more overall.

# Load Functions ----------------------------------------------------------
source("kernels.R")
source("lugsail.R")
source("R.R")
library(Matrix)

# LRV Estimator -----------------------------------------------------------


# For mother kernels
LRV_mother_estimator <- function(new_b, all_autocovariances, the_kernel, d){
  big_T <- nrow(all_autocovariances)

  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the wieghts for a specific wieght.
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios.
  W <- sapply(new_b, function(b){
    M <- b*big_T
    if(M == 0){
      new_weights <- c(1, rep(0, big_T-1))
    } else{
      new_weights <- sapply(0:(big_T -1)/(M), the_kernel)
    }
  })
  W <- t(W)
  rownames(W) = paste("b=", round(new_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(big_T-1), sep = "")


  omega_hats <- W %*% all_autocovariances
  rownames(omega_hats) <- paste("b=", round(new_b, 3))
  return(omega_hats)
}



# For Lugsail Kernels
# Corrects the non-positive diagonal elements.

LRV_estimator <- function(new_b, all_autocovariances,
                          the_kernel, lugsail_parameters = list(r = 1, c= 0),
                          mother_omega, d){
  big_T <- nrow(all_autocovariances)

  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the weights for a specific weight.
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios.
  W <- sapply(new_b, function(b){
    M <- b*big_T
    if(M == 0){
      new_weights <- c(1, rep(0, big_T-1))
    } else{
      new_weights <- sapply(0:(big_T -1)/(M), lugsail, the_kernel = the_kernel,
                            lugsail_parameters = lugsail_parameters)
    }
  })

  W <- t(W)
  rownames(W) = paste("b=", round(new_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(big_T-1), sep = "")

  # [#, ] the new_b value
  # [ , #]  a component of the estimated LRV
  #    Omega_11, Omega_12,  ..., Omega_1d; ...; Omega_d1, ...Omega_dd
  omega_hats <- W %*% all_autocovariances
  rownames(omega_hats) <- paste("b=", round(new_b, 3))

  #  ---- Fix PSD Issues ----
  # diagonal indices, the columns to check for omega_hats
  check_index <- seq(1, d^2, by = (d+1))

  # Matrix with each row corresponding to a b, and each column corresponding
  # to a diagonal element
  # [#,  ]: new_b
  # [ , #]: diagonal element
  counts <- matrix(0, nrow = nrow(omega_hats), ncol = length(check_index))

  for(i in 1:length(check_index)){
    col_index <- check_index[i]
    index_fix <- omega_hats[, col_index] <= 0
    omega_hats[index_fix, col_index] <- mother_omega[index_fix,col_index]
  }

  return(omega_hats)
}






# Calculate F-Statistics --------------------------------------------------
# Calculate the F test statistic for dimensions 1, ..., d

# b: the bandwidth, proportion of estimated autocovs that have non-zero weight
# the_sim_data:
#     - Contains `big_T` simulated dependent random vectors of dimension (d).
#     - Should already be centered with hypothesized or estimated means.
# the_means:
#     - The estimated mean vector.
#     - Should be of length d.
#     - Need this because the_sim_data is already centered.
# all_autocovariances:
#     - Rows are correspond to estimated autocov (R) at lag [0, ..., (big_T -1)]
#     - Columns correspond to the vectorization of the estimated autocov matrix:
#          R11, R12, R13, ..., R1d, R21, R22, ..., R2d, ..., Rd1, ...Rdd
# kernel: the name of the kernel function
# lugsail_parameters:
#     - a named list, list(r = 1, c= 0) that contains the lugsail parameters
#     - default is non-lugsail

F_stats <- function(the_means, omega_hats, d = 1, big_T){
    null_means = rep(0, d)

    # Vector of F-statistics for each new_b value
    F_stat_by_b <- apply(omega_hats, 1, function(one_omega_hat){
      one_omega_hat <- matrix(one_omega_hat, nrow = sqrt(length(one_omega_hat)))
      one_omega_hat <- one_omega_hat[1:d, 1:d]
      one_omega_hat <- as.matrix(nearPD(one_omega_hat)$mat)
      inv_one_omega <- solve(one_omega_hat)
      num <- (the_means - null_means)[1: d]
      F_stat <-  big_T*(num%*% inv_one_omega  %*%num)/ d
      return(F_stat)
    })
  return(F_stat_by_b)
}


get_kernel_F_stats <- function(new_b, the_means, d, all_autocovariances,
                               the_kernel, lugsail_type, big_T, q = 1){
  # Need this for psd corrections.
  omega_mother <- LRV_mother_estimator(new_b, all_autocovariances, the_kernel)

  # Mother Kernel
  if(lugsail_type == "Mother"){
    simulated_F_stats <- t(F_stats(the_means, omega_mother, d, big_T))
  }

  # Zero Lugsail
  else if(lugsail_type == "Zero"){
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "Zero")
    omega_zero <- LRV_estimator(new_b, all_autocovariances,
                                the_kernel = the_kernel,
                                lugsail_parameters = lug_para,
                                mother_omega= omega_mother, d=d)
    simulated_F_stats <- t(F_stats(the_means, omega_zero, d, big_T))
  }

  # Over Lugsail
  else if (lugsail_type == "Over"){
    # Over Lugsail
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "Over")
    omega_over <- LRV_estimator(new_b, all_autocovariances,
                                the_kernel= the_kernel,
                                lugsail_parameters = lug_para,
                                mother_omega = omega_mother, d=d)
    simulated_F_stats <- t(F_stats(the_means, omega_over, d, big_T))
  }

  return(simulated_F_stats)
}



# Simulate the F-Statistics for this Setting ------------------------------


# Generates a null data set.  Calculates the autocovariance matrices.
# Calculates the F-statistic
# Need to repeat this simulation many times to then find what asymptotic CV is.


simulate_f_stat <- function(big_T = 1000, d = 1, the_kernel = bartlett, q=1, lugsail_type = "Mother", new_b = b){

  # Simulate the data
  sim_data <- matrix(rnorm(d*big_T), nrow = big_T, ncol = d)
  the_means <- colMeans(sim_data)
  sim_data <- apply(sim_data, 1, function(row) row - the_means)
  sim_data <- t(sim_data)

  # ------- AutoCovariance Matrices  -------
  # [#, ] the lag (0, ..., big_T-1)
  # [ , #]  a component of the vectorized autocov matrix.
  #    R11, R12,  ..., R1d; R21, R22, ..., R2d; ...; Rd1, ...Rdd
  max_lag <- ceiling(max(new_b)*big_T)
  all_autocovariances <- matrix(0, nrow = d^2, ncol = big_T)
  all_autocovariances[1:d^2, 1:(max_lag+1)] <-sapply(0:max_lag, R,the_sim_data = sim_data, big_T = big_T)
  all_autocovariances <- t(all_autocovariances)


  # ------- F-statistics for various b values settings (2D array) -------
  F_kernel <- get_kernel_F_stats(new_b, the_means, d,
                                 all_autocovariances,
                                 the_kernel, lugsail_type, big_T, q = q)

  return(c(F_kernel))
}

# Main --------------------------------------------------------------------

generate_cv_multi <- function(b, d = 2, alpha = 0.05,
                        the_kernel = bartlett, lugsail_type = "Mother",  q=1,
                        return_F_stats = F,
                        num_replicates = 50000, replicate_size = 1000){

  big_T <- replicate_size

  # Each row is a b values
  # Each column is a simulated data set
  F_stats <- replicate(num_replicates,
                            simulate_f_stat(big_T = big_T, d = d,
                                            the_kernel = the_kernel, q=q,  new_b = b))

  F_stats <- matrix(c(t(F_stats)), nrow = num_replicates, ncol = length(b))
  rownames(F_stats) <- paste("Sim", 1:num_replicates, sep = "")
  colnames(F_stats) <- paste("b = ", b, sep = "")

  critical_values <- sapply(alpha, function(the_alpha) {
    apply(F_stats, 2, quantile, probs = (1-the_alpha))
  })

  critical_values <- matrix(c(critical_values), nrow = length(b), ncol = length(alpha))
  rownames(critical_values) <- paste("b = ", b, sep = "")
  colnames(critical_values) <- paste("alpha = ", alpha, sep = "")


  if(return_F_stats){
    return(list(CV = critical_values, F_stats = F_stats))
  } else{
    return(critical_values)
  }
}


# Example
set.seed(26)
generate_cv_multi(b = c(0.1, 0.05), d = 2, the_kernel = bartlett,
                  return_F_stats = T, num_replicates = 10)



generate_cv_multi(b = c(0.1), d = 2, alpha = c(0.05, 0.01),
                  the_kernel = bartlett, return_F_stats = T, num_replicates = 10)





