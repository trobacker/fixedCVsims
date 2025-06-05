library(Matrix)

# For mother kernels
LRV_mother_estimator <- function(b, all_autocovariances, the_kernel,
                                 big_T= nrow(all_autocovariances),
                                 d = ncol(all_autocovariances)){

  # Make the weights that correspond to the autocovariances
  M <- b*big_T
  if(M == 0){
    new_weights <- c(1, rep(0, big_T-1))
  } else{
    new_weights <- sapply(0:(big_T -1)/(M), the_kernel)
  }
  W <- new_weights
  W <- t(as.matrix(W))
  colnames(W) = paste("Lag=", 0:(big_T-1), sep = "")


  # LRV estimate
  omega_hats <- W %*% all_autocovariances
  omega_hats <- matrix(omega_hats, nrow = sqrt(d))
  return(omega_hats)
}




# For Lugsail Kernels
# Corrects the non-positive diagonal elements.
LRV_lugsail_estimator <- function(b, all_autocovariances,
                          the_kernel, lugsail_parameters = list(r = 1, c= 0),
                          mother_omega, big_T= nrow(all_autocovariances),
                          d = ncol(all_autocovariances)){

  # Make the weights that correspond to the autocovariances
  M <- b*big_T
  if(M == 0){
    new_weights <- c(1, rep(0, big_T-1))
  } else{
    new_weights <- sapply(0:(big_T -1)/(M), lugsail, the_kernel = the_kernel,
                          lugsail_parameters = lugsail_parameters)
  }
  W <- new_weights
  W <- t(as.matrix(W))
  colnames(W) = paste("Lag=", 0:(big_T-1), sep = "")

  # LRV estimate
  omega_hats <- W %*% all_autocovariances
  omega_hats <- matrix(omega_hats, nrow = d)

  #  ---- Fix PSD Issues ----
  # diagonal indices, the columns to check for omega_hats
  check_index <- seq(1, d^2, by = (d+1))

  for(i in 1:d){
    if(omega_hats[i, i] <= 0){
      omega_hats[i, i] <- mother_omega[i,i]
    }
  }

  return(omega_hats)
}



# Master
LRV_estimator <- function(b, all_autocovariances,
                          the_kernel, lugsail_type,
                          big_T= nrow(all_autocovariances),
                          d = ncol(all_autocovariances)){

  # Start with making the mother estimator
  omega_mother <- LRV_mother_estimator(b, all_autocovariances, kernel_fct,
                                       big_T = big_T, d= ncol(all_autocovariances))
  # Mother Kernel
  if(lugsail_type == "Mother"){
    omega <- omega_mother
  }

  # Zero Lugsail
  else if(lugsail_type == "Zero"){
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "Zero")
    omega <- LRV_lugsail_estimator(b, all_autocovariances,
                                the_kernel = the_kernel,
                                lugsail_parameters = lug_para,
                                mother_omega= omega_mother, d=d)
  }

  # Over Lugsail
  else if (lugsail_type == "Over"){
    # Over Lugsail
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "Over")
    omega <- LRV_lugsail_estimator(b, all_autocovariances,
                                   the_kernel = the_kernel,
                                   lugsail_parameters = lug_para,
                                   mother_omega= omega_mother, d=d)
  }
  omega <- nearPD(omega)$mat # Check if computationally PD
  return(omega)
}
