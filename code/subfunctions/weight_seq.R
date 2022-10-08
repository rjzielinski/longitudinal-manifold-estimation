weight_seq <- function(x_obs, mu, sigma, epsilon = 0.001, max_iter = 1000) {
  
  # x_obs is the dataset of interest
  # There are n observations, with each observation being a D-dimensional point
  # x_obs is a n x D matrix
  # mu is a vector of the knots in a mixture density estimation
  # sigma is a bandwidth of this density estimation
  # epsilon is a predetermined tolerance of the Euclidean distance between thetas in two consecutive steps
  # max_iter is a predetermined upper bound of the number of steps in this iteration
  
  n_D <- dim(x_obs)
  n <- n_D[1]
  D <- n_D[2]
  N <- dim(mu)[1]
  
  A <- matrix(NA, nrow = n, ncol = N)
  for (j in 1:N) {
    A_prepare <- function(x) {
      return(ker(x, mu[j, ], sigma))
    }
    A[, j] <- apply(x_obs, 1, A_prepare) # A[i, j] is \psi_sigma (x_i - mu_j)
  }
  
  # The initial guess for weights
  theta_old <- rep(1 / N, N)
  
  # absolute value of the difference between theta_new and theta_old
  abs_diff <- 10 * epsilon
  
  count <- 0
  
  # the initial guess of the Lagrangian multiplier
  lambda_hat_old <- c(n, rep(-1, D))
  
  while ((abs_diff > epsilon) & (count <= max_iter)) {
    
    W <- t(t(A) * theta_old)
    W <- W / apply(W, 1, sum) # use rowsums as substitute for apply?
    w <- apply(W, 2, sum) # colsums?
    
    flambda <- function(lambda) {
      # function for computing Lagrangian multipliers
      denom_temp <- apply(
        t(t(cbind(rep(1, dim(mu)[1]), mu)) * lambda),
        1,
        sum
      ) # again use rowsums
      
      num_temp <- mu * w
      
      f1 <- sum(w / denom_temp)
      f2 <- apply(
        num_temp * (as.vector(1 / denom_temp)), # is as.vector necessary?
        2,
        sum
      ) # use colsums
      f <- dist_euclidean(f1, 1) + dist_euclidean(f2, apply(x_obs, 2, mean)) # use colmeans
      return(f)
    }
    
    lambda_hat <- nlm(flambda, lambda_hat_old, iterlim = 1000)$estimate
    
    theta_new <- w / apply(
      t(t(cbind(rep(1, dim(mu)[1]), mu)) * lambda_hat),
      1,
      sum
    ) # rowsums
    abs_diff <- dist_euclidean(theta_new, theta_old)
    if (is.na(abs_diff)) {
      abs_diff <- 0
      theta_new <- theta_old
    }
    
    # pmax and pmin guarantee that theta_j's are in [0, 1]
    theta_old <- pmax(theta_new, 0)
    theta_old <- pmin(theta_old, 1) 
    count <- count + 1
    lambda_hat_old <- lambda_hat
  }
  
  theta_hat <- pmax(theta_new, 0)
  theta_hat <- pmin(theta_hat, 1)
  return(theta_hat)
}