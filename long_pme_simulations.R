library(tidyverse)
library(plotly)
library(profvis)

source("code/lpme.R")
source("code/pme.R")

sim_D2d1_case1 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise, time_noise) { 
  I <- 1000
  t <- rnorm(I, mean = 0, sd = 1)
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 2)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, vertical_noise, horizontal_noise) {
    return(
      c(
        tau + (horizontal_multiplier * sin(time_val)) + horizontal_noise,
        sin(tau + pi / 2) + (vertical_multiplier * sin(time_val)) + vertical_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier, 
      vertical_noise, 
      horizontal_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  data.points <- X + cbind(e1, e2)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

sim_D2d1_case2 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise, time_noise) { 
  I <- 1000
  t <- runif(I, min = -3 * pi, max = 3 * pi) 
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 2)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, vertical_noise, horizontal_noise) {
    return(
      c(
        tau + (horizontal_multiplier * sin(time_val)) + horizontal_noise,
        sin(tau) + (vertical_multiplier * sin(time_val)) + vertical_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier, 
      vertical_noise, 
      horizontal_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  data.points <- X + cbind(e1, e2)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

sim_D2d1_case3 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise, time_noise) { 
  I <- 1000
  t <- runif(I, min = 0, max = 1.5 * pi)
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 2)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, vertical_noise, horizontal_noise) {
    return(
      c(
        cos(tau) + (horizontal_multiplier * sin(time_val)) + horizontal_noise,
        sin(tau) + (vertical_multiplier * sin(time_val)) + vertical_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier, 
      vertical_noise, 
      horizontal_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  data.points <- X + cbind(e1, e2)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

sim_D2d1_case4 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise, time_noise) { 
  I <- 1000
  t <- runif(
    I, 
    min = 0 + (time_val * (pi / 4)), 
    max = (1.5 * pi) + (time_val * (pi / 4))
  )
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 2)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, vertical_noise, horizontal_noise) {
    return(
      c(
        cos(tau) + (horizontal_multiplier * sin(time_val)) + horizontal_noise,
        sin(tau) + (vertical_multiplier * sin(time_val)) + vertical_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier, 
      vertical_noise, 
      horizontal_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  data.points <- X + cbind(e1, e2)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

sim_D3d1_case1 <- function(time_val, vertical_multiplier, horizontal_multiplier, depth_multiplier, noise, time_noise) { 
  I <- 1000
  t <- runif(I, min = -1, max = 1)
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  depth_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  e3 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 3)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, depth_multiplier, vertical_noise, horizontal_noise, depth_noise) {
    return(
      c(
        tau + (horizontal_multiplier * sin(time_val)) + horizontal_noise,
        (tau ^ 2) + (vertical_multiplier * sin(time_val)) + vertical_noise,
        (tau ^ 3) + (depth_multiplier * sin(time_val)) + depth_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier,
      depth_multiplier,
      vertical_noise, 
      horizontal_noise,
      depth_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 3, byrow = TRUE)
  data.points <- X + cbind(e1, e2, e3)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

sim_D3d1_case2 <- function(time_val, vertical_multiplier, horizontal_multiplier, depth_multiplier, noise, time_noise) { 
  I <- 1000
  t <- runif(I, min = 0, max = 3 * pi)
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  depth_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  e3 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 3)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, depth_multiplier, vertical_noise, horizontal_noise, depth_noise) {
    return(
      c(
        tau + (horizontal_multiplier * sin(time_val)) + horizontal_noise,
        cos(tau) + (vertical_multiplier * sin(time_val)) + vertical_noise,
        sin(tau) + (depth_multiplier * sin(time_val)) + depth_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier,
      depth_multiplier,
      vertical_noise, 
      horizontal_noise,
      depth_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 3, byrow = TRUE)
  data.points <- X + cbind(e1, e2, e3)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

sim_D3d2_case1 <- function(time_val, vertical_multiplier, horizontal_multiplier, depth_multiplier, noise, time_noise) { 
  I <- 1000
  t1 <- runif(I, min = -1, max = 1)
  t2 <- runif(I, min = -1, max = 1)
  t <- cbind(t1, t2)
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  depth_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  e3 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 3)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, depth_multiplier, vertical_noise, horizontal_noise, depth_noise) {
    return(
      c(
        tau[1] + (horizontal_multiplier * sin(time_val)) + horizontal_noise,
        tau[2] + (vertical_multiplier * sin(time_val)) + vertical_noise,
        (norm_euclidean(tau) ^ 2) + (depth_multiplier * sin(time_val)) + depth_noise
      )
    )
  }
  
  X <- apply(
    t, 
    1,
    manifold,
    time_val = time_val, 
    vertical_multiplier = vertical_multiplier, 
    horizontal_multiplier = horizontal_multiplier,
    depth_multiplier = depth_multiplier,
    vertical_noise = vertical_noise,
    horizontal_noise = horizontal_noise,
    depth_noise = depth_noise
  ) %>% 
    unlist() %>% 
    matrix(ncol = 3, byrow = TRUE)
  data.points <- X + cbind(e1, e2, e3)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

### Simulation Case 1

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals, 
  sim_D2d1_case1, 
  vertical_multiplier = 1,
  horizontal_multiplier = 1, 
  noise = 0.15, 
  time_noise = 1
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) & 
    (sim_pred[, dim_idx] < idx_max)
}
    
r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

### Simulation Case 2

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals, 
  sim_D2d1_case2, 
  vertical_multiplier = 1, 
  horizontal_multiplier = 1, 
  noise = 0.15, 
  time_noise = 0.5
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) & 
    (sim_pred[, dim_idx] < idx_max)
}
    
r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

### Simulation Case 3

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals, 
  sim_D2d1_case3, 
  vertical_multiplier = 1, 
  horizontal_multiplier = 1, 
  noise = 0.15, 
  time_noise = 0.5
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) & 
    (sim_pred[, dim_idx] < idx_max)
}
    
r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

## Simulation Case 4

time_vals <- seq(0, 1, 0.5)

set.seed(100)
df_list <- lapply(
  time_vals, 
  sim_D2d1_case4, 
  vertical_multiplier = 1, 
  horizontal_multiplier = 1, 
  noise = 0.15, 
  time_noise = 0.5
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) & 
    (sim_pred[, dim_idx] < idx_max)
}
    
r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

### D3, d1 Simulation Case 1

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals, 
  sim_D3d1_case1, 
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  depth_multiplier = 1,
  noise = 0.15, 
  time_noise = 1
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) & 
    (sim_pred[, dim_idx] < idx_max)
}
    
r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

### D3, d1 Simulation Case 2

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals, 
  sim_D3d1_case2,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  depth_multiplier = 1,
  noise = 0.15, 
  time_noise = 1
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(data_points) - 1)
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) & 
    (sim_pred[, dim_idx] < idx_max)
}
    
r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(data_points) - 1)
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y", "z")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~z,
  frame = ~time,
  opacity = 0.5,
  type = "scatter3d",
  mode = "markers"
)

### D3, d2 Simulation Case 1

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals, 
  sim_D3d2_case1, 
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  depth_multiplier = 1,
  noise = 0.15, 
  time_noise = 1
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 2)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 1)

r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>% 
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) & 
    (sim_pred[, dim_idx] < idx_max)
}
    
r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2:3][r_inrange,]))
r_max <- max(unlist(grid_mat[, 2:3][r_inrange,]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.2
)
r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>% 
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
  progress(i, nrow(sim_pred))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r1", "r2", "x1", "x2", "x3")

range_x1 <- max(sim_pred_full_df$x1) - min(sim_pred_full_df$x1)
range_x2 <- max(sim_pred_full_df$x2) - min(sim_pred_full_df$x2)
range_x3 <- max(sim_pred_full_df$x3) - min(sim_pred_full_df$x3)

plot_ly(
  sim_pred_full_df,
  x = ~x1,
  y = ~x2,
  z = ~x3,
  frame = ~time,
  opacity = 0.5,
  type = "scatter3d",
  mode = "markers"
) %>% 
  layout(
    scene = list(
      xaxis = list(
        range = c(
          min(sim_pred_full_df$x1) - (0.1 * range_x1),
          max(sim_pred_full_df$x1) + (0.1 * range_x1)
        )
      ),
      yaxis = list(
        range = c(
          min(sim_pred_full_df$x2) - (0.1 * range_x2),
          max(sim_pred_full_df$x2) + (0.1 * range_x2)
        )
      ),
      zaxis = list(
        range = c(
          min(sim_pred_full_df$x3) - (0.1 * range_x3),
          max(sim_pred_full_df$x3) + (0.1 * range_x3)
        )
      )
    )
  )
