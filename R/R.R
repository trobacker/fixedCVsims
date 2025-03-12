
R <- function(h, the_sim_data, big_T){
  index <- 1:(big_T -h)
  the_sim_data <- as.matrix(the_sim_data)

  # Already centered
  est <- lapply(index, function(i){
    est <- the_sim_data[i, ]%*%t(the_sim_data[(i+h), ])
  })

  # Sum together
  autocov_s <- Reduce('+', est)/big_T

  # Need the negative lags too.
  # Because of symmetry, we can do this.
  if(h!=0){
    autocov_s <- autocov_s + t(autocov_s)
  }

  return(autocov_s)
}

# For 1-dimentional data s
R_d1 <- function(h, one_sim_data){
  big_T <- length(one_sim_data)
  est_mean <- mean(one_sim_data)
  index <- 1:(big_T -h)
  est <- (one_sim_data[index] - est_mean)*(one_sim_data[(index+h)]-est_mean)

  autocov_s <- sum(est)/big_T

  # Multiply by 2 because we want the negative and positive autocovariance.
  if(h!=0){
    autocov_s <- 2*autocov_s
  }

  return(autocov_s)
}

# Calculate all the autocovariances for a 1-dimensional simulation
all_R <- function(one_sim_data){
  big_T <- length(one_sim_data)
  all_auto_cov <- sapply(0:(big_T-1), R_d1, one_sim_data= one_sim_data)
  return(all_auto_cov)
}
