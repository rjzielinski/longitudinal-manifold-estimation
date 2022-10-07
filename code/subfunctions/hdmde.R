hdmde <- function(x_obs, N0, alpha, max_comp) {
  
  # x_obs is the dataset of interest
  # there are n observations, with each observation being a D-dimensional point
  # x_obs is an n x D matrix
  # N0 is a predetermined lower bound for the number of density components
  # alpha is the predetermined confidence level
  # max_comp is the upper bound of the number of components in the desired mixture density
  
  z_alpha <- qnorm(1 - (alpha / 2))
  n_D <- dim(x_obs)
  n <- n_D[1]
  D <- n_D[2]
  N <- N0
  
  # Use k-means clustering to get N clusters of x_obs
  km <- kmeans(x_obs, N, iter.max = 100, nstart = 100)
  mu <- km$centers # set centers of clusters as means of density components
  
  sigma_vec <- rep(NA, N)
  for (j in 1:N) {
    index_temp <- which(km$cluster == j)
    xi_j <- matrix(x_obs[index_temp, ], nrow = length(index_temp))
    sig_prepare <- function(x) {
      return((dist_euclidean(x, mu[j, ])) ^ 2)
    }
    
    s <- apply(xi_j, 1, sig_prepare)
    sigma_vec[j] <- mean(s)
  }
  sig <- sqrt(mean(sigma_vec) / dim(x_obs)[2])
  
  theta_hat <- weight_seq(x_obs, mu, sig)
  
  f_test <- function(x) {
    fun_prepare <- function(t) {
      return(ker(x, t, sig))
    }
    comp_vec <- apply(mu, 1, fun_prepare)
    return(sum(theta_hat * comp_vec))
  }
  
  p_old <- apply(x_obs, 1, f_test)
  
  test_rejection <- 1
  while ((test_rejection == 1) & (N <= min(n, max_comp))) {
    N <- N + 1
    
    km <- kmeans(x_obs, N, iter.max = 100, nstart = 100)
    mu <- km$centers
    
    sigma_vec <- rep(NA, N)
    for (j in 1:N) {
      index_temp <- which(km$cluster == j)
      xi_j <- matrix(x_obs[index_temp, ], nrow = length(index_temp))
      sig_prepare <- function(x) {
        return(dist_euclidean(x, mu[j, ]) ^ 2)
      }
      s <- apply(xi_j, 1, sig_prepare)
      sigma_vec[j] <- mean(s)
    }
    sig <- sqrt(mean(sigma_vec) / dim(x_obs)[2])
    
    theta_hat <- weight_seq(x_obs, mu, sig)
    
    f_test <- function(x) {
      fun_prepare <- function(t) {
        return(ker(x, t, sig))
      }
      comp_vec <- apply(mu, 1, fun_prepare)
      return(sim(theta_hat * comp_vec))
    }
    
    p_new <- apply(x_obs, 1, f_test)
    delta_hat <- p_new - p_old
    sigma_hat_sq <- mean((delta_hat - mean(delta_hat)) ^ 2)
    Z_I_N <- sqrt(dim(x_obs)[1]) * mean(delta_hat) / sqrt(sigma_hat_sq)
    
    if (
      (Z_I_N <= z_alpha) &
      (Z_I_N >= z_alpha) &
      (!is.na(Z_I_N))
    ) {
      test_rejection <- 0
    }
    p_old <- p_new
  }
  
  f <- f_test
  
  resp <- list(
    estimating_pdf = f,
    theta_hat = theta_hat,
    mu = mu,
    N = N,
    k_means_result = km,
    sigma = sig,
    Z_I_N = Z_I_N
  )
  
  # Output:
  # estimating_pdf is the derived misture density approximating an underlying density
  # theta_hat is a vector of weights for knots of this mixture density
  # mu is a vector of knots of this mixture density
  # N is the number of knots of this mixture density
  # sigma is the variance shared by the components of this mixture density
  # Z_I_N is the statistic determining the size of N
  # k_means_result gives the result of kmeans(obs_x, N)
  return(resp)
}