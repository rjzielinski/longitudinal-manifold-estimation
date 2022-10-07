lpme <- function(df, d, tuning.para.seq = exp(-15:5), alpha = 0.05, max.comp = 100, epsilon = 0.05, max.iter = 100, print.MSDs = TRUE, print_plots = TRUE, SSD_ratio_threshold = 100) {
  # df is an N x (D + 1) matrix, with the first column corresponding
  # to the time point at which each observation was collected
  # this matrix should include the observations from all time points
  
  source("pme.R")
  require(plotly)
  require(svMisc)
  # source("Principal_Manifold_Estimation.R")
  
  time_points <- df[, 1] %>% 
    unique()
  
  pme_results <- list()
  funcs <- list()
  x_test <- list()
  r <- list()
  
  for (idx in 1:length(time_points)) {
    df_temp <- df[df[, 1] == time_points[idx], ]
    pme_results[[idx]] <- pme(
      x.obs = df_temp[, -1],
      d = d
    )
    funcs[[idx]] <- pme_results[[idx]]$embedding.map
    r_test <- seq(
      from = -10,
      to = 10,
      length.out = dim(pme_results[[idx]]$knots)[1]
    )
    r_list <- lapply(numeric(d), function(x) r_test)
    r_mat <- as.matrix(expand.grid(r_list))
    r_length <- dim(r_mat)[1]
    x_temp <- apply(
      r_mat,
      1,
      funcs[[idx]]
    ) %>% 
      matrix(nrow = dim(r_mat)[1], byrow = TRUE)
    # x_temp <- map(r_mat, ~ funcs[[idx]](.x)) %>% 
    #   unlist() %>% 
    #   matrix(nrow = r_length, byrow = TRUE)
    idx_inrange <- matrix(nrow = dim(x_temp)[1], ncol = dim(x_temp)[2])
    for (dim_idx in 1:dim(x_temp)[2]) {
      idx_range <- max(df_temp[, dim_idx + 1]) - min(df_temp[, dim_idx + 1])
      idx_min <- min(df_temp[, dim_idx + 1]) - (0.2 * idx_range)
      idx_max <- max(df_temp[, dim_idx + 1]) + (0.2 * idx_range)
      idx_inrange[, dim_idx] <- (x_temp[, dim_idx] > idx_min) & (x_temp[, dim_idx] < idx_max)
    }
    
    r_inrange <- rowSums(idx_inrange) == dim(x_temp)[2]
    if (sum(r_inrange) == 0) {
      r_inrange <- rowSums(idx_inrange) > 0
    }
    r_min <- min(r_mat[r_inrange, ])
    r_min <- ifelse(is.na(r_min), -10, r_min)
    r_max <- max(r_mat[r_inrange, ])
    r_max <- ifelse(is.na(r_max), 10, r_max)
    r_test <- seq(
      r_min,
      r_max,
      length.out = max(10, round((dim(pme_results[[idx]]$knots)[1]) ^ (1 / d)))
    )
    r_list <- lapply(numeric(d), function(x) r_test)
    r_mat <- as.matrix(expand.grid(r_list))
    r_length <- dim(r_mat)[1]
    x_test[[idx]] <- apply(
      r_mat,
      1,
      funcs[[idx]]
    ) %>% 
      matrix(nrow = dim(r_mat)[1], byrow = TRUE)
    
    # x_test[[idx]] <- map(r_test, ~ funcs[[idx]](.x)) %>% 
    #   unlist() %>% 
    #   matrix(nrow = r_length, byrow = TRUE)
    if (norm_euclidean(x_test[[idx]][nrow(x_test[[idx]]), ]) - norm_euclidean(x_test[[idx]][1, ]) < 0) {
      r_mat <- r_mat[nrow(r_mat):1, ]
    }
    r[[idx]] <- cbind(time_points[idx], r_mat)
  }
  
  r <- reduce(r, rbind)
  r_df <- data.frame(r)
  rnames <- paste0("r", 1:(dim(r)[2] - 1))
  names(r_df) <- c("time", rnames)
  x_test <- reduce(x_test, rbind) 
  x_test_df <- data.frame(x_test)
  xnames <- paste0("x", 1:(dim(x_test)[2]))
  names(x_test_df) <- xnames
  
  D_new <- dim(x_test)[2]
  d_new <- dim(r)[2]
  n_new <- dim(x_test)[1]
  gamma <- 4 - d_new

  sig_new <- 0.001
  theta_hat_new <- weight.seq(x_test, x_test, sig_new)
  centers_new <- x_test
  sigma_new <- 0.001
  W_new <- diag(theta_hat_new)
  X_new <- x_test
  I_new <- length(theta_hat_new)

  dissimilarity_matrix_new <- as.matrix(dist(X_new))
  # isomap_initial <- isomap(
  #   dissimilarity_matrix_new,
  #   ndim = d_new,
  #   k = 10
  # )
  # t_initial <- isomap_initial$points
  t_initial <- r
  
  MSE_seq_new <- vector()
  SOL_new <- list()
  TNEW_new <- list()

  for (tuning.ind in 1:length(tuning.para.seq)) {
    print(
      paste0(
        "The tuning parameter is gamma[",
        as.character(tuning.ind),
        "] = ",
        as.character(tuning.para.seq[tuning.ind]),
        "."
      )
    )

    w_new <- tuning.para.seq[tuning.ind]
    t_new <- t_initial
    T_new <- cbind(rep(1, I_new), t_new)

    E_new <- matrix(NA, ncol = I_new, nrow = I_new)
    for (j in 1:I_new) {
      E_prepare <- function(t) {
        eta.kernel(t - t_new[j, ], gamma)
      }
      E_new[, j] <- apply(t_new, 1, E_prepare)
    }

    M1_new <- cbind(
      2 * E_new %*% W_new %*% E_new + 2 * w_new * E_new,
      2 * E_new %*% W_new %*% T_new,
      T_new
    )
    M2_new <- cbind(
      2 * t(T_new) %*% W_new %*% E_new,
      2 * t(T_new) %*% W_new %*% T_new,
      matrix(0, ncol = d_new + 1, nrow = d_new + 1)
    )
    M3_new <- cbind(
      t(T_new),
      matrix(0, ncol = d_new + 1, nrow = d_new + 1),
      matrix(0, ncol = d_new + 1, nrow = d_new + 1)
    )
    M_new <- rbind(
      M1_new,
      M2_new,
      M3_new
    )

    b_new <- rbind(
      2 * E_new %*% W_new %*% X_new,
      2 * t(T_new) %*% W_new %*% X_new,
      matrix(0, nrow = d_new + 1, ncol = D_new)
    )
    sol_new <- ginv(M_new) %*% b_new

    eta.func <- function(t) {
      eta.func.prepare <- function(tau) {
        return(eta.kernel(t - tau, gamma))
      }
      return(
        matrix(
          apply(
            t_new,
            1,
            eta.func.prepare
          ),
          ncol = 1
        )
      )
    }

    f_new <- function(t) {
      return(
        as.vector(
          t(sol_new[1:I_new, ]) %*%
            eta.func(t) + t(sol_new[(I_new + 1):(I_new + d_new + 1),]) %*%
            matrix(c(1, t), ncol = 1)
        )
      )
    }

    f0_new <- f_new

    X_initial_guess <- cbind(X_new, t_new)
    projection.index.f0 <- function(x.init) {
      projection(
        x.init[1:D_new],
        f0_new,
        x.init[(D_new + 1):(D_new + d_new)]
      )
    }

    t_new <- matrix(
      t(apply(X_initial_guess, 1, projection.index.f0)),
      nrow = I_new
    )

    SSD.prepare <- function(x.prin, f) {
      return(
        dist.euclidean(
          x.prin[1:D_new],
          f(x.prin[(D_new + 1):(D_new + d_new)])
        ) ^ 2
      )
    }

    X_projection_index <- cbind(X_new, t_new)
    SSD.prepare.again <- function(x.init) {
      return(SSD.prepare(x.init, f_new))
    }
    SSD_new <- sum(as.vector(apply(X_projection_index, 1, SSD.prepare.again)))
    
    if (print_plots == TRUE) {
      time_vals <- seq(
        min(time_points),
        max(time_points),
        0.1
      )
      
      r_vals <- seq(
        from = -10,
        to = 10,
        by = 1
      )
      r_list <- lapply(numeric(d_new - 1), function(x) r_vals)
      r_mat <- as.matrix(expand.grid(r_list))
      r_length <- dim(r_mat)[1]
      pred_grid <- expand_grid(r_vals, r_mat) %>% 
        as.matrix()
      
      f_pred <- matrix(nrow = nrow(pred_grid), ncol = ncol(X_new))
      for (i in 1:nrow(pred_grid)) {
        f_pred[i, ] <- f_new(unlist(as.vector(pred_grid[i, ])))
        progress(i, max.value = nrow(pred_grid), console = FALSE)
      }
      
      idx_inrange <- matrix(nrow = dim(f_pred)[1], ncol = dim(f_pred)[2])
      for (dim_idx in 1:dim(f_pred)[2]) {
        idx_range <- max(X_new[, dim_idx]) - min(X_new[, dim_idx])
        idx_min <- min(X_new[, dim_idx]) - (0.2 * idx_range)
        idx_max <- max(X_new[, dim_idx]) + (0.2 * idx_range)
        idx_inrange[, dim_idx] <- (f_pred[, dim_idx] > idx_min) & 
          (f_pred[, dim_idx] < idx_max)
      }
      
      r_inrange <- rowSums(idx_inrange) == dim(f_pred)[2]
      if (sum(r_inrange) == 0) {
        r_inrange <- rowSums(idx_inrange) > 0
      }
      r_min <- min(pred_grid[r_inrange, -1])
      r_min <- ifelse(is.na(r_min), -10, r_min)
      r_max <- max(pred_grid[r_inrange, -1])
      r_max <- ifelse(is.na(r_max), 10, r_max)
      r_vals <- seq(
        r_min,
        r_max,
        by = (r_max - r_min) / 40
      )
      r_list <- lapply(numeric(d_new - 1), function(x) r_vals)
      r_mat <- as.matrix(expand.grid(r_list))
      pred_grid <- expand_grid(time_vals, r_mat) %>% 
        as.matrix()
      f_pred <- matrix(nrow = nrow(pred_grid), ncol = ncol(X_new))
      for (i in 1:nrow(pred_grid)) {
        f_pred[i, ] <- f_new(unlist(as.vector(pred_grid[i, ])))
        progress(i, max.value = nrow(pred_grid), console = FALSE)
      }
      
      f_pred_full <- cbind(pred_grid, f_pred)
      
      if (D_new == 2) {
        plt <- plot_ly(
          x = f_pred_full[, d_new + 1],
          y = f_pred_full[, d_new + 2],
          z = f_pred_full[, 1],
          color = f_pred_full[, 1],
          opacity = 0.5,
          type = "scatter3d",
          mode = "lines"
        )
        print(plt)
      } else if (D_new == 3) {
        plt <- plot_ly(
          x = f_pred_full[, d_new + 1],
          y = f_pred_full[, d_new + 2],
          z = f_pred_full[, d_new + 3],
          frame = f_pred_full[, 1],
          type = "scatter3d",
          mode = "markers",
          opacity = 0.5
        )
        print(plt)
      }  
    }
    

    count <- 1
    SSD_ratio <- 10 * epsilon

    while ((SSD_ratio > epsilon) & (SSD_ratio <= SSD_ratio_threshold) & (count <= (max.iter - 1))) {
      SSD_old <- SSD_new
      f0 <- f_new

      T_new <- cbind(rep(1, I_new), t_new)
      E_new <- matrix(NA, ncol = I_new, nrow = I_new)
      for (j in 1:I_new) {
        E.prepare <- function(t) {
          eta.kernel(t - t_new[j, ], gamma)
        }
        E_new[, j] <- apply(t_new, 1, E.prepare)
      }

      M1_new <- cbind(
        2 * E_new %*% W_new %*% E_new + 2 * w_new * E_new,
        2 * E_new %*% W_new %*% T_new,
        T_new
      )
      M2_new <- cbind(
        2 * t(T_new) %*% W_new %*% E_new,
        2 * t(T_new) %*% W_new %*% T_new,
        matrix(0, ncol = d_new + 1, nrow = d_new + 1)
      )
      M3_new <- cbind(
        t(T_new),
        matrix(0, ncol = d_new + 1, nrow = d_new + 1),
        matrix(0, ncol = d_new + 1, nrow = d_new + 1)
      )
      M_new <- rbind(
        M1_new,
        M2_new,
        M3_new
      )

      b_new <- rbind(
        2 * E_new %*% W_new %*% X_new,
        2 * t(T_new) %*% W_new %*% X_new,
        matrix(0, nrow = d_new + 1, ncol = D_new)
      )
      sol_new <- ginv(M_new) %*% b_new

      eta.func <- function(t) {
        eta.func.prepare <- function(tau) {
          return(eta.kernel(t - tau, gamma))
        }
        return(matrix(apply(t_new, 1, eta.func.prepare), ncol = 1))
      }

      f_new <- function(t) {
        return(
          as.vector(
            t(sol_new[1:I_new, ]) %*%
              eta.func(t) + t(sol_new[(I_new + 1):(I_new + d_new + 1), ]) %*%
              matrix(c(1, t), ncol = 1)
          )
        )
      }

      t_old <- t_new

      X_initial_guess <- cbind(X_new, t_new)
      projection.index.f0 <- function(x.init) {
        projection(x.init[1:D_new], f0, x.init[(D_new + 1):(D_new + d_new)])
      }
      t_new <- matrix(
        t(apply(X_initial_guess, 1, projection.index.f0)),
        nrow = I_new
      )

      X_projection_index <- cbind(X_new, t_new)
      SSD.prepare.again <- function(x.init) {
        return(SSD.prepare(x.init, f_new))
      }
      
      if (print_plots == TRUE) {
        time_vals <- seq(
          min(time_points),
          max(time_points),
          0.1
        )
      
        r_vals <- seq(
          from = -10,
          to = 10,
          by = 1
        )
        r_list <- lapply(numeric(d_new - 1), function(x) r_vals)
        r_mat <- as.matrix(expand.grid(r_list))
        r_length <- dim(r_mat)[1]
        pred_grid <- expand_grid(r_vals, r_mat) %>% 
          as.matrix()
        
        f_pred <- matrix(nrow = nrow(pred_grid), ncol = ncol(X_new))
        for (i in 1:nrow(pred_grid)) {
          f_pred[i, ] <- f_new(unlist(as.vector(pred_grid[i, ])))
          progress(i, max.value = nrow(pred_grid), console = FALSE)
        }
        
        idx_inrange <- matrix(nrow = dim(f_pred)[1], ncol = dim(f_pred)[2])
        for (dim_idx in 1:dim(f_pred)[2]) {
          idx_range <- max(X_new[, dim_idx]) - min(X_new[, dim_idx])
          idx_min <- min(X_new[, dim_idx]) - (0.2 * idx_range)
          idx_max <- max(X_new[, dim_idx]) + (0.2 * idx_range)
          idx_inrange[, dim_idx] <- (f_pred[, dim_idx] > idx_min) & 
            (f_pred[, dim_idx] < idx_max)
        }
        
        r_inrange <- rowSums(idx_inrange) == dim(f_pred)[2]
        if (sum(r_inrange) == 0) {
          r_inrange <- rowSums(idx_inrange) > 0
        }
        r_min <- min(pred_grid[r_inrange, -1])
        r_min <- ifelse(is.na(r_min), -10, r_min)
        r_max <- max(pred_grid[r_inrange, -1])
        r_max <- ifelse(is.na(r_max), 10, r_max)
        r_vals <- seq(
          r_min,
          r_max,
          by = (r_max - r_min) / 40
        )
        r_list <- lapply(numeric(d_new - 1), function(x) r_vals)
        r_mat <- as.matrix(expand.grid(r_list))
        pred_grid <- expand_grid(time_vals, r_mat) %>% 
          as.matrix()
        f_pred <- matrix(nrow = nrow(pred_grid), ncol = ncol(X_new))
        for (i in 1:nrow(pred_grid)) {
          f_pred[i, ] <- f_new(unlist(as.vector(pred_grid[i, ])))
          progress(i, max.value = nrow(pred_grid), console = FALSE)
        }
        
        f_pred_full <- cbind(pred_grid, f_pred)
        
        if (D_new == 2) {
          plt <- plot_ly(
            x = f_pred_full[, d_new + 1],
            y = f_pred_full[, d_new + 2],
            z = f_pred_full[, 1],
            color = f_pred_full[, 1],
            type = "scatter3d",
            mode = "lines",
            opacity = 0.5
          )
          print(plt)
        } else if (D_new == 3) {
          plt <- plot_ly(
            x = f_pred_full[, d_new + 1],
            y = f_pred_full[, d_new + 2],
            z = f_pred_full[, d_new + 3],
            frame = f_pred_full[, 1],
            type = "scatter3d",
            mode = "markers",
            opacity = 0.5
          )
          print(plt)
        }
      }
      
      
      SSD_new <- sum(as.vector(apply(X_projection_index, 1, SSD.prepare.again)))
      SSD_ratio <- abs(SSD_new - SSD_old) / SSD_old
      count <- count + 1

      print(
        paste0(
          "SSD = ",
          as.character(round(SSD_new)),
          ", SSD.ratio is ",
          as.character(round(SSD_ratio, 4)),
          ", and this is the ",
          as.character(count),
          "th step of iteration."
        )
      )
    }

    data_initial <- matrix(0, nrow = 1, ncol = D_new + d_new)
    for (i in 1:I_new) {
      length.temp <- 1
      X.i <- matrix(X_new[i, ], nrow = 1)
      t_temp <- matrix(rep(t_new[i, 1], length.temp))
      for (j in 1:d_new) {
        t_temp <- cbind(
          t_temp,
          rep(t_new[i, j], length.temp)
        )
      }
      t_temp <- matrix(t_temp[, -1], nrow = length.temp)
      data_initial <- rbind(data_initial, cbind(X.i, t_temp))
    }

    data_initial <- data_initial[-1, ]
    proj.para.prepare <- function(data.init) {
      return(
        projection(
          data.init[1:D_new],
          f_new,
          data.init[(D_new + 1):(D_new + d_new)]
        )
      )
    }
    proj_para <- matrix(
      t(apply(data_initial, 1, proj.para.prepare)),
      ncol = d_new
    )
    proj_points <- t(apply(proj_para, 1, f_new))
    diff_data_fit <- apply(
      data_initial[, 1:D_new] - proj_points,
      1,
      norm.euclidean
    )
    MSE_new <- mean(diff_data_fit ^ 2)
    MSE_seq_new[tuning.ind] <- MSE_new
    print(
      paste(
        "When gamma = ",
        as.character(w_new),
        ", MSD = ",
        as.character(MSE_new),
        "."
      )
    )
    SOL_new[[tuning.ind]] <- sol_new
    TNEW_new[[tuning.ind]] <- t_new

    if (tuning.ind >= 4) {
      if (
        (MSE_seq_new[tuning.ind] > MSE_seq_new[tuning.ind - 1]) &
        (MSE_seq_new[tuning.ind - 1] > MSE_seq_new[tuning.ind - 2]) &
        (MSE_seq_new[tuning.ind - 2] > MSE_seq_new[tuning.ind - 3])
      ) {
        break
      }
    }
  }

  optimal_ind <- min(which(MSE_seq_new == min(MSE_seq_new)))
  sol_opt <- SOL_new[[optimal_ind]]
  t_new_opt <- TNEW_new[[optimal_ind]]
  eta.func <- function(t) {
    eta.func.prepare <- function(tau) {
      return(eta.kernel(t - tau, gamma))
    }
    return(
      matrix(
        apply(t_new_opt, 1, eta.func.prepare),
        ncol = 1
      )
    )
  }
  f.optimal <- function(t) {
    return(
      as.vector(
        t(sol_opt[1:I_new, ]) %*%
          eta.func(t) + t(sol_opt[(I_new + 1):(I_new + d_new + 1), ]) %*%
          matrix(c(1, t), ncol = 1)
      )
    )
  }

  if (print.MSDs == TRUE) {
    plot(
      log(tuning.para.seq[1:length(MSE_seq_new)]),
      MSE_seq_new,
      xlab = "Log Gamma",
      ylab = "MSD",
      type = "l"
    )
    lines(
      log(tuning.para.seq[1:length(MSE_seq_new)]),
      MSE_seq_new,
      type = "p",
      pch = 20,
      col = "orange",
      cex = 2
    )
    abline(
      v = log(tuning.para.seq[optimal_ind]),
      lwd = 1.5,
      col = "darkgreen",
      lty = 2
    )
    print(
      paste(
        "The optimal tuning parameter is ",
        as.character(tuning.para.seq[optimal_ind]),
        ", and the MSD of the optimal fit is ",
        as.character(MSE_seq_new[optimal_ind], ".")
      )
    )
  }

  resp <- list(
    embedding_map = f.optimal,
    MSD = MSE_seq_new,
    knots = centers_new,
    weights_of_knots = theta_hat_new,
    coe_kernel = sol_opt[1:I_new, ],
    coe_poly = sol_opt[(I_new + 1):(I_new + d_new + 1), ],
    SOL = SOL_new,
    TNEW = TNEW_new,
    T_parameter = sol_opt,
    Coef = t_new_opt
  )
  return(resp)
}
  