# mother kernel g_q
# zero lugsail g_q = 0
# over lugsail is the negative of mother
g_q <- list("bartlett" = 1, "parzen" = 6, "th" = pi^2/4, "qs" = 1.421223)

# Suns: Mother bandwidth rule
mother_b_rule <- function(big_T, alpha, d, w_q, g_q, q = 1, tau = 0.15*alpha){
  cv <- qchisq((1-alpha), d)
  num <- dchisq(cv, d)*cv*g_q*w_q
  den <- tau
  opt_b <- (num/den)^(1/q)/big_T
  return(opt_b)
}

# Zero lugsail bandwidth rule
zero_b_rule <- function(rho, big_T, alpha = 0.05, d =1, tau = -alpha^(1/(2*d))/(big_T*log(rho))){
  cv <- qchisq((1-alpha), d)
  num <-  tau*(1+ rho)
  den <- dchisq(cv, d)*cv*2*rho^2
  opt_b <- log(num/den)/(big_T*log(rho))
  opt_b <- max(0, opt_b)
  if(is.na(opt_b)){opt_b <- 0 }

  return(opt_b)
}

# Over bandwidth rule
over_b_rule <- function(rho, big_T, alpha, d, w_q, g_q, q=1, tau = -alpha^(1/(2*d))/(big_T*log(rho))){
  try_b <- (1:(big_T/2))/big_T
  cv <- qchisq((1-alpha), d)

  # Type 1 error
  distortion <- dchisq(cv, d)*cv*2*rho^2*rho^(try_b*big_T)/(1 + rho)+ (try_b*big_T)^(-q)*dchisq(cv, d)*cv*g_q*w_q
  type_1 <- alpha + distortion
  opt_b <- try_b[which(abs(distortion) <= tau)]
  opt_b <- min(opt_b)
  if(opt_b == Inf){
    warning("No bandwidth meets criteria, it is recommended to increase tolerance level.")
    opt_b <- 0
  }
  return(opt_b)
}


# Andrews rule (just as a reference, will delete later)
mse_rule <- function(big_T, rho = .7, q = 1, d = 1){
  w_1 <- (2*rho/(1-rho^2))
  the_b  <- 1.1447*(w_1/big_T)^(2/3)
  if(is.na(the_b)){the_b <- 0}
  return(the_b)
}


# Main ---------------------------------------------------------

get_b <- function(the_data, alpha = 0.05, the_kernel ="bartlett", lugsail="Mother", tau = NA){

  # dimensions
  big_T <- nrow(the_data)
  d <- ncol(the_data)

  # calculate the rho
  # average the AR(1) coefficient for all dimensions
  all_rhos <- rep(0, d)
  for(i in 1:d){
    all_rhos[i] <- stats::acf(the_data[,i], plot = F)$acf[2]
  }
  rho <- mean(all_rhos)

  # tau, the neighborhood
  if(is.na(tau)){
    #tau <- -alpha^(1/(2*d))/(big_T*log(rho))
    tau <- alpha*.15
  }

  # kernel statistic information
  q <- 1
  if(the_kernel != "bartlett"){q <- 2}
  if(q == 1){
    w_q <- 2*rho/(1-rho^2)
  } else {w_q <- 2*rho/(1-rho)^2 }
  g_q <- g_q[[the_kernel]]

  # g_q based on lugsail type
  if(lugsail == "Zero"){
    g_q <- 0
  } else if (lugsail == "Over"){
    g_q <- - g_q
  }

  zero_b <- zero_b_rule(rho, big_T, alpha = 0.05, d =d, tau = tau)
  mother_b <- mother_b_rule(big_T, alpha, d, w_q, g_q, q, tau = tau)
  over_b <- over_b_rule(rho, big_T, alpha = 0.05, d =d, w_q, -g_q, tau = tau)
  mse_b <- mse_rule(big_T, rho, q, d)
  all_b <- c(mse_b, mother_b, zero_b, over_b)
  names(all_b) <- c("MSE", "Sun", "Zero", "Over")

  return(data.frame(all_b))
}



set.seed(62)
d <- 1
big_T <- 200
rho_matrix <- matrix(0, nrow = d, ncol = d)
diag(rho_matrix) <- 0.7
sim_data <- matrix(0, nrow = big_T, ncol = d)
sim_data[1, ] <- rnorm(d)
for(i in 2:big_T){
  sim_data[i,] <- sim_data[i-1, ]%*%rho_matrix + rnorm(d)
}


get_b(sim_data)



