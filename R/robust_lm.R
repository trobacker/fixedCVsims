library(Matrix)


# Description -------------------------------------------------------------

# Input:
# - fit: lm object that needs robust inference
# - the_kernel: the mother kernel
# - lugsail: the lugsail setting setting
# - method: how do you want to conduct the inference

# Returns:
# - robust_summary: a summary table of the lm object with robust inference
# - rhos_bs: the estimated rho's and b's for the t-statistics. Will be a matrix where
#            each row is a coefficient (and F-stat), and there is two columns, the rhos and bs selected
# - F_stat: this will be the F-statistic and corresponding p-value
# - model_details: will return the the_kernel, lugsail, and CV method.



# simulate an example data set --------------------------------------------


# Multivariate
set.seed(62)
d <- 5
big_T <- 200
model_rho <- 0.7
rho_matrix <- matrix(0, nrow = d, ncol = d)
diag(rho_matrix) <- model_rho
sim_data <- matrix(0, nrow = big_T, ncol = d)
sim_data[1, ] <- rnorm(d)
the_sd <- 1/4
for(i in 2:big_T){
  sim_data[i,] <- sim_data[i-1, ]%*%rho_matrix + rnorm(d, sd = the_sd)
}

disturbance <- rnorm(1)
for(i in 2:big_T){
  disturbance[i] <- disturbance[i-1]*model_rho + rnorm(1, sd = the_sd)
}

y <- apply(sim_data, 1, sum) + disturbance
the_data <- data.frame(y, sim_data)
colnames(the_data) <- c("y", paste("x", 1:d, sep = ""))

fit <- lm(y ~. , the_data)





# support functions -------------------------------------------------------
p_values <- function(test_stat, the_b = 0, the_d= 1,  the_kernel = "Bartlett",
                     lugsail = "Mother", method = "simulated"){

  # Get the appropriate CVs
  CV_01  <- get_cv(the_b, the_d, 0.01, the_kernel, lugsail, method)*the_d
  CV_025 <- get_cv(the_b, the_d, 0.025, the_kernel, lugsail, method)*the_d
  CV_05  <- get_cv(the_b, the_d, 0.05, the_kernel, lugsail, method)*the_d
  CV_10  <- get_cv(the_b, the_d, .10, the_kernel, lugsail, method)*the_d

  # Use t-statistic if d==1
  if(d == 1){
    CV_01 <- sqrt(CV_01)
    CV_025 <- sqrt(CV_025)
    CV_05 <- sqrt(CV_05)
    CV_10 <- sqrt(CV_10)
  }

  # Check the p_value code
  p_value <- ">=0.10"
  if(abs(test_stat)>CV_01){
    p_value <- "<0.01*"

  }else if(abs(test_stat)>CV_025){
    p_value <- "<0.025."

  }else if(abs(test_stat)>CV_05){
    p_value <- "<0.05."
  }else if(abs(test_stat)>CV_10){
    p_value <- "<0.10"
  }
  cv_table <- c(the_b, the_d, CV_10, CV_05, CV_025, CV_01)
  return(list(p_value = p_value, cv_table = cv_table))
}

# main  -------------------------------------------------------------------

robust_lm <- function(fit, the_kernel = "Bartlett", lugsail= "Mother",
                      method = "simulated", tau = 0.05*.15){


  # ------- Basic statistics needed from the LM object -------
  kernel_fct <- bartlett
  if(the_kernel == "QS"){
    kernel_fct <- qs
  } else if(the_kernel == "TH"){
    kernel_fct <- th
  } else if(the_kernel == "Parzen"){
    kernel_fct <- parzen
  }

  # ------- Basic statistics needed from the LM object -------
  X <- model.matrix(fit)
  coefs <- fit$coefficients
  residuals <- fit$residuals
  errors <- apply(X, 2, function(col) col*residuals) # not normal errors, moment conditions
  errors <- errors - apply(errors, 2, mean) # center immediately
  M <- t(X)%*%X/big_T
  big_T <- nrow(fit$model)

  # ------- AutoCovariance Matrices  -------
  # [#, ] the lag (0, ..., big_T-1)
  # [ , #]  a component of the vectorized autocov matrix.
  #    R11, R12,  ..., R1d; R21, R22, ..., R2d; ...; Rd1, ...Rdd
  all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = errors, big_T = big_T)
  all_autocovariances <- t(all_autocovariances)

  # ------- Summary Object -------
  # Also calculates b-values that will be moved somewhere else
  summary <- matrix(0, nrow = length(coefs), ncol = 4)
  for(i in 1:length(coefs)){

    # get the beta cov matrix
    the_b <- get_b(errors[,i], the_kernel = tolower(the_kernel),
                   lugsail=lugsail, tau = tau) # only consider the correlation levels for coef you want
    omega <- LRV_estimator(the_b, all_autocovariances, kernel_fct, lugsail,
                           big_T, d = length(coefs))
    beta_cov <- (solve(M)%*%omega%*%solve(M))

    # Standard error (se)
    one_beta_cov <- beta_cov[i,i]
    se <- sqrt(one_beta_cov/big_T)

    # Test statistic
    test_stat <- coefs[i]/se

    #return
    summary[i,] <- c(coefs[i], se, test_stat, the_b)
  }

  cv_table <- matrix(0, nrow = length(coefs), ncol = 6)
  coefs_p_values <- rep(NA, length(coefs))
  for(i in 1:length(coefs)){
    keep <- p_values(summary[i, 3], the_b = summary[i, 4],
                     the_d= 1,  the_kernel = the_kernel,
                     lugsail = lugsail, method = method)
    coefs_p_values[i] <- keep$p_value
    cv_table[i, ] <- keep$cv_table
  }
  summary <- data.frame(summary, coefs_p_values)
  rownames(summary) <- names(coefs)
  colnames(summary) <- c("Estimate", "Std. Error", "t value","b", "P(>|t|)")

  all_rhos <- rep(NA, 4)
  for(i in 1:ncol(errors)){
    all_rhos[i] <- stats::acf(errors[,i], plot = F)$acf[2]
  }
  summary$b <- NULL # delete the b column from the summary table

  # ------- F-test -------
  the_b <- get_b(errors, the_kernel = tolower(the_kernel),
                 lugsail = lugsail, tau = tau)
  omega <- LRV_estimator(the_b, all_autocovariances, kernel_fct, lugsail,
                         big_T, d = length(coefs))
  omega  <- as.matrix(omega) # Check if computationally PD
  beta_cov <- (solve(M)%*%omega%*%solve(M))[-1, -1]
  F_stat <- big_T * coefs[-1] %*% solve(beta_cov) %*% coefs[-1]/(length(coefs)-1)
  keep <- p_values(F_stat, the_b = the_b, the_d =c(length(coefs)-1),
                   the_kernel = the_kernel,
                   lugsail = lugsail, method = method)
  F_stat <- data.frame(F_stat, keep$p_value)
  names(F_stat) <- c("F Statistic", "P-Value")

  # ------- CV-Table -------
  cv_table <- rbind(cv_table, keep$cv_table)
  cv_table <- data.frame(rho = c(all_rhos, mean(all_rhos)),
                                 cv_table)
  rownames(cv_table) <- c(names(coefs), "F-test")
  colnames(cv_table) <- c("rho","b", "dim", ".10", ".05", ".025", ".01")

  # ------- Return Values -------
  return_me <- list("Summary_Table" = summary,
                    "F_test" = F_stat,
                    "CV_table" = cv_table)
  return(return_me)
}



robust_lm(fit)
robust_lm(fit, lugsail = "Zero")
robust_lm(fit, the_kernel = "QS")
robust_lm(fit, the_kernel = "QS", tau = 0.0117 )
robust_lm(fit, the_kernel = "QS", lugsail = "Zero")
robust_lm(fit, the_kernel = "QS", tau = alpha*.5)


robust_lm(fit, the_kernel = "QS")
robust_lm(fit, the_kernel = "QS", method = "analytical")

robust_lm(fit, tau = -sqrt(alpha)/ (big_T * log(rho)))


robust_lm(fit, the_kernel = "Bartlett", lugsail = "Mother", method = "analytical")

robust_lm(fit, the_kernel = "Bartlett", lugsail = "Mother", method = "simulated")

