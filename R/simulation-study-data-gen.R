library(cmdstanr)
library(dplyr)
library(stringr)
rdirichlet <- function(n, alpha) {
  d <- length(alpha)
  rates <- rep(1, d)
  gammas <- replicate(n = n,
                      rgamma(d, shape = alpha,
                             rate = rates)
                      )
  draws <- sweep(gammas, MARGIN = 2, STATS = colSums(gammas), FUN = "/")
  return(t(draws))
}
to_bin <- function(N, str, pad) {
  if (N == 1) {
    str <- paste0(str, "1")
    str <- str_pad(str, pad, side = "right", pad = "0")
    return(str)
  }
  if (N == 0) {
    str <- paste0(str, "0")
    str <- str_pad(str, pad, side = "right", pad = "0")
    return(str)
  }
  div <- N %/% 2
  rem <- N %% 2
  str <- paste0(str, rem)
  to_bin(div, str, pad)
}
softmax <- function(x) {
  x_app <- c(x,0)
  return(exp(x_app) / sum(exp(x_app)))
}
inv_softmax <- function(x) {
  l_x <- log(x) -
    log(x[length(x)])
  return(head(l_x, -1))
}


gen_param_set <- function(N_p_r, N_r, N_z, N_a, N_x, Y_const_a = TRUE,
                          restricted_to_zero = NULL, null_hypo = FALSE) {
  N_U <- 2^N_z
  names_C <- sapply(0:(N_U - 1), \(x) to_bin(x, "", N_z))
  if (!is.null(restricted_to_zero)) {
    names_C <- setdiff(names_C, restricted_to_zero)
  }
  N_U <- length(names_C)
  S_cond_U_Z <- sapply(names_C, \(x) as.numeric(strsplit(x, "")[[1]])) |> t()
  if (N_U == 4) {
    theta_glob <- c(0.91, 0.05, 0.005, 0.035)
  } else if (N_U == 8) {
    theta_glob <- c(0.91, 0.05, 0.001, 0.001,  0.001, 0.001, 0.001)
    theta_glob <- c(theta_glob, 1 - sum(theta_glob))
    p_1 <- theta_glob[which(S_cond_U_Z[, 1] == 1)] |> sum()
    p_2 <- theta_glob[which(S_cond_U_Z[, 2] == 1)] |> sum()
    p_3 <- theta_glob[which(S_cond_U_Z[, 3] == 1)] |> sum()
  } else {
    theta_glob <- rdirichlet(1, rep(3, N_U))[1, ]
  }
  thetas <- rdirichlet(N_r, theta_glob*100)
  thetas_raw <- matrix(0,N_r, N_U - 1)
  for (r in 1:N_r) {
    thetas_raw[r,] <- inv_softmax(thetas[r,])
  }
  beta_X <- matrix(rnorm((N_x - 1) * (N_U - 1)), N_x - 1, N_U - 1)
  thetas_w_X <- array(0, c(N_r, N_x, N_U))
  for (r in 1:N_r) {
    for (x in 1:N_x) {
      if (x == 1) {
        thetas_w_X[r, x, ] <- softmax(thetas_raw[r, ])
      } else {
        thetas_w_X[r, x, ] <- softmax(thetas_raw[r, ] + beta_X[x-1,])
      }
    }
  }

  a_vec_dummy <- rdirichlet(N_U, rep(2, N_a)) 
  a_vec_raw <- matrix(0, N_U, N_a - 1)
  beta_a <- matrix(rnorm((N_x - 1) * (N_a - 1)), N_x - 1, N_a - 1)
  for (u in 1:N_U)
    a_vec_raw[u,] <- inv_softmax(a_vec_dummy[u,])
  a_vec <- array(0, c(N_U, N_x, N_a))
  
  for (u in 1:N_U) {
    for (x in 1:N_x) {
      if (x == 1) {
        a_vec[u, x, ] <- softmax(a_vec_raw[u, ])
      } else {
        a_vec[u, x, ] <- softmax(a_vec_raw[u, ] + beta_a[x - 1, ])
      }
    }
  }


  beta_y <- log(1.1^(1:(N_x - 1)))
  N <- N_p_r * N_r
  slopes_C <- array(NA_real_, c(N_a, N_U, N_z))
  for (C_i in seq_along(names_C)) {
    names_i <- S_cond_U_Z[C_i, ]
    sel_vec <- which(names_i == 1)
    n_pos <- sum(names_i)
    slopes_C[1, C_i, sel_vec] <- rnorm(n_pos)
    if (names_C[C_i] == "11") {
      slopes_C[1, C_i, sel_vec] <- c(
        log(0.3 / 0.7),
        log(0.3 / 0.7 * 0.4)
      )
      if (null_hypo) {
        slopes_C[1, C_i, sel_vec] <- c(
          log(0.22 / 0.78),
          log(0.22 / 0.78)
        )
      }
    } else if (names_C[C_i] == "10") {
      slopes_C[1, C_i, sel_vec] <- c(
        log(0.15 / 0.85)
      )
      if (null_hypo) {
        slopes_C[1, C_i, sel_vec] <- c(
          log(0.25 / 0.75)
        )
      }
    } else if (names_C[C_i] == "01") {
      slopes_C[1, C_i, sel_vec] <- c(
        log(0.2 / 0.8)
      )
    } else if (names_C[C_i] == "111") {
      o_0 <- 0.3 / 0.7
      x <- o_0 * 0.4
      slopes_C[1, C_i, sel_vec] <- c(
        log(o_0),
        log(o_0),
        log(x)
      )
      if (null_hypo) {
        slopes_C[1, C_i, sel_vec] <- c(
          log(0.22 / 0.78),
          log(0.22 / 0.78),
          log(0.22 / 0.78)
        )
      }
    } else if (names_C[C_i] == "101") {
      slopes_C[1, C_i, sel_vec] <- c(
        log(0.2 / 0.8),
        log(0.10 / 0.9)
      )
    } else if (names_C[C_i] == "110") {
      slopes_C[1, C_i, sel_vec] <- c(
        log(0.3 / 0.7),
        log(0.15 / 0.85)
      )
    } else if (names_C[C_i] == "011") {
      slopes_C[1, C_i, sel_vec] <- c(
        log(0.25 / 0.75),
        log(0.08 / 0.92)
      )
    } else if (names_C[C_i] == "010") {
      slopes_C[1, C_i, sel_vec] <- c(log(0.25 / 0.75))
    } else if (names_C[C_i] == "100") {
      slopes_C[1, C_i, sel_vec] <- c(log(0.10 / 0.9))
      if (null_hypo) {
        slopes_C[1, C_i, sel_vec] <- c(
          log(0.25 / 0.8)
        )
      }
    } else if (names_C[C_i] == "001") {
      slopes_C[1, C_i, sel_vec] <- c(log(0.25 / 0.75))
    }
    if (!Y_const_a) {
      for (a in 2:N_a) {
        slopes_C[a, C_i, sel_vec] <- rep(0,n_pos)
        if (names_C[C_i] == "11") {
          slopes_C[a, C_i, sel_vec] <- c(log(0.925^(a - 1)), log(0.825^(a - 1)))
          if (null_hypo) {
            slopes_C[a, C_i, sel_vec] <- c(log(0.925^(a - 1)), log(0.925^(a - 1)))
          }
        }
        if (names_C[C_i] == "10") {
          slopes_C[a, C_i, sel_vec] <- log(0.95^(a - 1))
        }
        if (names_C[C_i] == "01") {
          slopes_C[a, C_i, sel_vec] <- 0
        }
        if (names_C[C_i] == "111") {
          slopes_C[a, C_i, sel_vec] <- c(log(0.925^(a - 1)),
                                         log(0.925^(a - 1)),
                                         log(0.925^(a - 1)))
        }
      }
    } else {
      for (a in 2:N_a) {
        slopes_C[a, C_i, sel_vec] <- 0 
      }
    }
  }
  return(
    list(
      thetas_w_X = thetas_w_X,
      a_vec = a_vec,
      slopes_C = slopes_C,
      beta_y = beta_y,
      theta_glob = theta_glob,
      thetas_raw = thetas_raw,
      beta_X = beta_X,
      N_U = N_U,
      names_C = names_C,
      S_cond_U_Z = S_cond_U_Z
    )
  )
}

gen_waning_set <- function(N_p_r = 5e3, N_r = 4, N_z = 3, N_a = 3, N_x = 3, Y_const_a = TRUE,
                          restricted_to_zero = c("100","010","101","011"), null_hypo = FALSE) {
  ## z_1: Placebo x placebo
  ## z_2: placebo x vax
  ## z_3: vax x placebo
  ## 0: not eligible for second time period
  ## 1: elibible for second time period
  ## VE_1 
  N_U <- 2^3
  names_C <- sapply(0:(N_U - 1), \(x) to_bin(x, "", N_z))
  if (!is.null(restricted_to_zero)) {
    names_C <- setdiff(names_C, restricted_to_zero)
  }
  N_U <- length(names_C)
  S_cond_U_Z <- sapply(names_C, \(x) as.numeric(strsplit(x, "")[[1]])) |> t()
  theta_glob <- c(0.035 / 2, 0.005 / 2, 0.05 / 2)
  theta_glob <- c(theta_glob, 1 - sum(theta_glob))
  p_1 <- theta_glob[which(S_cond_U_Z[, 1] == 0)] |> sum()
  p_2 <- theta_glob[which(S_cond_U_Z[, 2] == 0)] |> sum()
  p_3 <- theta_glob[which(S_cond_U_Z[, 3] == 0)] |> sum()
  thetas <- rdirichlet(N_r, theta_glob*100)
  thetas_raw <- matrix(0,N_r, N_U - 1)
  for (r in 1:N_r) {
    thetas_raw[r,] <- inv_softmax(thetas[r,])
  }
  beta_X <- matrix(rnorm((N_x - 1) * (N_U - 1)), N_x - 1, N_U - 1)
  thetas_w_X <- array(0, c(N_r, N_x, N_U))
  for (r in 1:N_r) {
    for (x in 1:N_x) {
      if (x == 1) {
        thetas_w_X[r, x, ] <- softmax(thetas_raw[r, ])
      } else {
        thetas_w_X[r, x, ] <- softmax(thetas_raw[r, ] + beta_X[x-1,])
      }
    }
  }

  a_vec_dummy <- rdirichlet(N_U, rep(2, N_a)) 
  a_vec_raw <- matrix(0, N_U, N_a - 1)
  beta_a <- matrix(rnorm((N_x - 1) * (N_a - 1)), N_x - 1, N_a - 1)
  for (u in 1:N_U)
    a_vec_raw[u,] <- inv_softmax(a_vec_dummy[u,])
  a_vec <- array(0, c(N_U, N_x, N_a))
  
  for (u in 1:N_U) {
    for (x in 1:N_x) {
      if (x == 1) {
        a_vec[u, x, ] <- softmax(a_vec_raw[u, ])
      } else {
        a_vec[u, x, ] <- softmax(a_vec_raw[u, ] + beta_a[x - 1, ])
      }
    }
  }


  beta_y <- log(1.1^(1:(N_x - 1)))
  N <- N_p_r * N_r
  slopes_C <- array(NA_real_, c(N_a, N_U, N_z))
  for (C_i in seq_along(names_C)) {
    names_i <- S_cond_U_Z[C_i, ]
    sel_vec <- which(names_i == 1)
    n_pos <- sum(names_i)
    slopes_C[1, C_i, sel_vec] <- rnorm(n_pos)
    o_0 <- p_1 * 0.5 / (1 - p_1 * 0.5)
    x <- p_3 * 0.5 / (1 - p_3 * 0.5)
    x2 <- x * 1.7
    if (names_C[C_i] == "111") {
      slopes_C[1, C_i, sel_vec] <- c(
        log(o_0),
        log(x),
        log(x2)
      )
      if (null_hypo) {
        slopes_C[1, C_i, sel_vec] <- c(
          log(o_0),
          log(x),
          log(x)
        )
      }
    } else if (names_C[C_i] == "110") {
      slopes_C[1, C_i, sel_vec] <- c(
        log(p_1 / (1 - p_1)),
        log(p_3 / (1 - p_3))
      )
    } else if (names_C[C_i] == "001") {
      slopes_C[1, C_i, sel_vec] <- log(p_1 / (1 - p_1))
    }
    if (!Y_const_a) {
      for (a in 2:N_a) {
        slopes_C[a, C_i, sel_vec] <- rep(0,n_pos)
        if (names_C[C_i] == "111") {
          slopes_C[a, C_i, sel_vec] <- c(log(0.925^(a - 1)),
                                         log(0.925^(a - 1)),
                                         log(0.925^(a - 1)))
        }
      }
    } else {
      for (a in 2:N_a) {
        slopes_C[a, C_i, sel_vec] <- 0 
      }
    }
  }
  return(
    list(
      thetas_w_X = thetas_w_X,
      a_vec = a_vec,
      slopes_C = slopes_C,
      beta_y = beta_y,
      theta_glob = theta_glob,
      thetas_raw = thetas_raw,
      beta_X = beta_X,
      N_U = N_U,
      names_C = names_C,
      S_cond_U_Z = S_cond_U_Z,
      restricted = restricted_to_zero
    )
  )
}

gen_data <- function(N_p_r, N_r, N_z, N_a, N_x, Y_const_a = TRUE,
                     restricted_to_zero = NULL,
                     sn_Y = 0.95, sp_Y = 0.6, sn_S = 0.9, sp_S = 0.95,
                     par_list = NULL, wane = FALSE) {

  N_U <- par_list$N_U
  S_cond_U_Z <- par_list$S_cond_U_Z
  names_C <- par_list$names_C
  thetas_w_X <- par_list$thetas_w_X
  a_vec <- par_list$a_vec
  slopes_C <- par_list$slopes_C
  beta_y <- par_list$beta_y

  N <- N_p_r * N_r
  X <- sample(N_x, N, replace = TRUE)
  idx <- as.vector(sapply(1:N_r,rep,times = N_p_r))
  Z <- sample(1:N_z,N_p_r * N_r, replace = TRUE)

  N_pos_beta_y <- sum(S_cond_U_Z)

  df_X <- data.frame(X = factor(X, levels = 1:N_x))
  X_mat <- model.matrix(~X, df_X)[, -1]
  S_P_0 <- mapply(function(i,x) sample(names_C, 1, TRUE, thetas_w_X[i, x, ]),
                  idx, X) 
  S <- sapply(S_P_0, \(x) as.numeric(strsplit(x,"")[[1]])) |> t()
  U <- 1 + S %*% 2^(0:(N_z-1))
  if (length(unique(U)) < max(U)) {
    U <- as.integer(as.factor(U))
  }
  Y <- matrix(NA_integer_, N, N_z)
  P <- array(NA_integer_, c(N, N_U, N_z))
  A <- mapply(
    function(u, x) sample(1:N_a, 1, replace = TRUE, a_vec[u, x, ]),
    U, X
  )
  for (C_i in seq_along(names_C)) {
    C <- names_C[C_i]
    sel_vec <- which(S_cond_U_Z[C_i, ] == 1)
    N_C <- sum(S_P_0 == C)
    df_A <- data.frame(A = factor(A[S_P_0 == C], levels = 1:N_a))
    A_mat <- model.matrix(~ A, df_A)
    X_mat_i <- X_mat[S_P_0 == C,]
    N_pos <- length(sel_vec)
    if (C_i > 1) {
      c_sum <- sum(S_cond_U_Z[1:(C_i - 1), ])
      if (N_pos == 1) {
        P[S_P_0 == C, C_i, sel_vec] <- plogis(A_mat %*% slopes_C[,C_i, sel_vec]
                                              + X_mat_i %*% beta_y)
        Y[S_P_0 == C, sel_vec] <- rbinom(N_C, 1, P[S_P_0 == C, C_i, sel_vec])
      } else {
        P_dummy <- matrix(NA_real_, N_C, N_pos)
        for (j in 1:N_pos) {
          P_dummy[, j] <- plogis(A_mat %*% slopes_C[, C_i, sel_vec[j]]
                                 + X_mat_i %*% beta_y)
          Y[S_P_0 == C, sel_vec[j]] <- rbinom(N_C, 1, P_dummy[,j])
        }
        P[S_P_0 == C, C_i, sel_vec] <- P_dummy
      }
    }
  }
  S_o <- sapply(1:N, \(x) S[x, Z[x]])
  Y_o <- sapply(1:N, \(x) Y[x, Z[x]])
  Y_o[is.na(Y_o)] <- 0
  S_tilde <- sapply(S_o,
                    \(x)
                    ifelse(x == 1,
                           rbinom(1, 1, sn_S),
                           rbinom(1, 1, (1 - sp_S))))
  Y_tilde <- sapply(Y_o,
                    \(x)
                    ifelse(x == 1,
                           rbinom(1, 1, sn_Y),
                           rbinom(1, 1, (1 - sp_Y))))
  p_x_given_r <- prop.table(table(idx, X), margin = 1)
  p_r_given_x <- prop.table(table(X, idx), margin = 1)
  p_x_given_a <- prop.table(table(A, X), margin = 1)
  prop_r_x <- prop.table(table(idx, X))
  prop_x <- prop.table(table(X))

  # Note to future Rob: The indices for beta_y don't correspond
  # to the indices for p_y_1, need to fix this for future research
  # Possible fix is to compute these qtys up in the loop over C_i
  # instead of down here
  
  ve <- matrix(NA_real_, N_a, N_x)
  p_y <- array(NA_real_, c(N_a, N_x, sum(S_cond_U_Z)))
  ve_a <- matrix(NA_real_, N_a, N_z - 1)
  ve_S <- matrix(NA_real_, N_r, N_z - 1)
  marg_ve_S_by_x <- matrix(NA_real_, N_x, N_z - 1)
  marg_ve <- matrix(NA_real_, N_x, N_z - 1)
  p_S <- matrix(NA_real_, N_x, N_U)
  p_y_1 <- array(NA_real_, c(N_a, sum(S_cond_U_Z)))
  for (a in 1:N_a) {
    id_i <- 1
    for (z in 1:N_z) {
      for (u in 1:N_U) {
        if (S_cond_U_Z[u, z] == 1) {
          p_y_1[a, id_i] <- plogis(slopes_C[1, u, z] + ifelse(a > 1, slopes_C[a, u, z], 0))
          id_i <- id_i + 1
        }
      }
    }
  }

  id_s <- as.integer(!wane)
  for (z in 1:(N_z - 1)) {
    for (r in 1:N_r) {
      ve_S[r,z] = 1 -
        sum(rowSums(thetas_w_X[r, , which(S_cond_U_Z[, z+1] == id_s)]) * p_x_given_r[r]) /
          sum(rowSums(thetas_w_X[r, , which(S_cond_U_Z[, 1] == id_s)]) * p_x_given_r[r])
    }
  }
  z_idx <- cumsum(colSums(S_cond_U_Z))
  marg_marg_ve <- rep(NA_real_, N_z - 1)
  marg_marg_diff_ve <- rep(NA_real_, N_z - 1)
  for (z in 1:(N_z - 1)) {
    denom_x <- num_x <- rep(NA_real_, N_x)
    for (x in 1:N_x) {
      marg_ve_S_by_x[x,z] <- 1 -
        sum(rowSums(thetas_w_X[, x, which(S_cond_U_Z[, z+1] == id_s)]) * p_r_given_x[x]) /
          sum(rowSums(thetas_w_X[, x, which(S_cond_U_Z[, 1] == id_s)]) * p_r_given_x[x])
      num_x[x] <- sum(p_y_1[, z_idx[z + 1]] * a_vec[N_U, x, ])
      denom_x[x] <- sum(p_y_1[, z_idx[1]] * a_vec[N_U, x, ])
      marg_ve[x,z] <- 1 - num_x[x] / denom_x[x] 
    }
    marg_marg_ve[z] <- 1 - sum(num_x * prop_x) / sum(denom_x * prop_x)
    marg_marg_diff_ve[z] <- sum(denom_x * prop_x) - sum(num_x * prop_x)

    for (a in 1:N_a)
      ve_a[a,z] = 1 - sum(p_y_1[a, z_idx[z + 1]] * p_x_given_a[a, ]) /
        sum(p_y_1[a, z_idx[1]] * p_x_given_a[a, ])
  }
  marg_ve_S <- rep(NA_real_, N_z-1)
  for (z in 1:(N_z-1))
    marg_ve_S[z] = 1 -
      sum(apply(thetas_w_X[, , which(S_cond_U_Z[,z+1] == id_s)],c(1,2),sum) * prop_r_x) / 
      sum(apply(thetas_w_X[, , which(S_cond_U_Z[,1] == id_s)],c(1,2),sum) * prop_r_x)
  for (x in 1:N_x)
    for (u in 1:N_U)
      p_S[x,u] <- sum(thetas_w_X[,x,u] * p_r_given_x[x])

  A_obs_mat <- matrix(
    c(
      c(0.5, 0.5, rep(0, 5)),
      c(0.25, 0.5, 0.25, rep(0, 4)),
      c(0, 0.25, 0.5, 0.25, rep(0, 3)),
      c(0, 0, 0.25, 0.5, 0.25, rep(0, 2)),
      c(0, 0, 0, 0.25, 0.5, 0.25, 0),
      c(0, 0, 0, 0, 0.25, 0.5, 0.25),
      c(0, 0, 0, 0, 0, 0.5, 0.5)
    ),
    7, 7,
    byrow = T
  )
  if (N_a == 3) {
    A_obs_mat <- matrix(
      c(
        c(0.75, 0.25, 0),
        c(0.125, 0.75, 0.125),
        c(0, 0.25, 0.75)
      ),
      3, 3,
      byrow = T
    )
  }
  A_obs <- sapply(
    A,
    \(x) sample(1:N_a, 1, prob = A_obs_mat[x, ])
  )

  return(list(S = S_o, Z = Z, Y = Y_o, A = A, A_obs = A_obs, S_P_0 = S_P_0,
              S_tilde = S_tilde,
              Y_tilde = Y_tilde,
              names_C = names_C,
              ## theta_glob = theta_glob,
              ## thetas_raw = thetas_raw,
              thetas_w_X = thetas_w_X,
              ## beta_X = beta_X,
              X = X,
              beta_y = beta_y,
              idx = idx,
              p_y_1 = p_y_1,
              S_cond_U_Z = S_cond_U_Z,
              a_vec = a_vec,
              slopes_C = slopes_C,
              sn_S = sn_S,
              sn_Y = sn_Y,
              sp_S = sp_S,
              sp_Y = sp_Y,
              ve = ve,
              ve_a = ve_a,
              ve_S = ve_S,
              marg_ve_S_by_x = marg_ve_S_by_x,
              p_S = p_S,
              marg_ve_S = marg_ve_S,
              marg_ve = marg_ve,
              marg_marg_ve = marg_marg_ve,
              p_x_given_r = p_x_given_r,
              p_r_given_x = p_r_given_x,
              p_x_given_a = p_x_given_a,
              prop_r_x = prop_r_x,
              prop_x = prop_x
              )
         )
}

proc_data <- function(S, Y, A, Z, idx,
                      S_cond_U_Z, X, metric) {
  N_a <- length(unique(A))
  N_z <- length(unique(Z))
  N_x <- length(unique(X))
  N_U <- nrow(S_cond_U_Z)
  N_Y_1 <- sum(S_cond_U_Z)
  N_r <- length(unique(idx))
  theta_raw_metric <- matrix(NA_real_,N_r,N_U - 1)
  beta_raw_metric <- matrix(NA_real_, N_x - 1, N_U - 1)
  p_a_raw_metric <- matrix(NA_real_, N_U, N_a - 1)
  beta_a_metric <- matrix(NA_real_, N_x - 1, N_a - 1)
  p_y_1_tilde_raw_metric <- matrix(NA_real_, N_a, N_Y_1)
  i <- 1
  for (r in 1:N_r) {
    for (u in 1:(N_U - 1)) {
      theta_raw_metric[r, u] <- metric[i]
      i <- i + 1
    }
  }
  for (x in 1:(N_x - 1)) {
    for (u in 1:(N_U - 1)) {
      beta_raw_metric[x,u] <- metric[i]
      i <- i + 1
    }
  }
  for (u in 1:N_U)
    for (a in 1:(N_a - 1)) {
      p_a_raw_metric[u,a] <- metric[i]
      i <- i + 1
    }

  beta_a_metric <- matrix(NA_real_, N_x - 1, N_a - 1)
  p_y_1_tilde_raw_metric <- matrix(NA_real_, N_a, N_Y_1)
  for (x in 1:(N_x - 1))
    for (a in 1:(N_a - 1)) {
      beta_a_metric[x,a] <- metric[i]
      i <- i + 1
    }
  sn_S_raw_metric <- metric[i]; i <- i + 1
  sp_S_raw_metric <- metric[i]; i <- i + 1
  sp_Y_raw_metric <- metric[i]; i <- i + 1
  for (a in 1:N_a)
    for (y in 1:N_Y_1) {
      p_y_1_tilde_raw_metric[a,y] <- metric[i]
      i <- i + 1
    }
  comb_data <- data.frame(S = S, Y = Y, A = A, Z = Z, idx = idx, X = X)
  tab <- table(
    comb_data$S,
    comb_data$Y,
    comb_data$A,
    comb_data$idx,
    comb_data$Z,
    comb_data$X
  )
  dimnames(tab) <- list(
    S = 0:1,
    Y = 0:1,
    A = 1:N_a,
    idx = seq_len(N_r),
    Z = 1:N_z,
    X = 1:N_x
  )
  new_array <- array(NA_integer_, dim = c(N_r, N_z, N_a, N_x, 4))
  for (z in 1:N_z) {
    for (r in seq_len(N_r)) {
      for (a in 1:N_a) {
        for (x in 
1:N_x)
        new_array[r, z, a, x, ] <- as.vector(tab[, , a, r, z, x])
      }
    }
  }

  S_0_cond_U_Z <- 1 - S_cond_U_Z
  n_Y_p_U <- rowSums(S_cond_U_Z)
  n_Y_p_z <- colSums(S_cond_U_Z)
  ordered_n_Y <- order(n_Y_p_U)
  N_p_S <- table(n_Y_p_U)
  N_p_S_mat <- data.frame(t(N_p_S)) %>%
    mutate(
      n_Y_p_U = as.numeric(as.character(n_Y_p_U))
    )
  N_p_S_mat <- as.matrix(N_p_S_mat[,2:3])
  z_idx <- apply(S_cond_U_Z, 1,
                 \(x) which(x == 1))
  z_idxs <- list()
  for (count_i in seq_along(N_p_S)) {
    count <- N_p_S_mat[count_i,1]
    sel_vec <- which(n_Y_p_U == count)
    N_count <- length(sel_vec)
    raw_idx <- z_idx[sel_vec] |>
      unlist()
    idx_mat <- matrix(raw_idx, byrow=T, N_count, count)
    z_idxs[[count_i]] <- idx_mat
  }
  idx_mat <- matrix(0,N_U,N_z)
  for (i in 2:N_U) {
    idx_sel <- which(S_cond_U_Z[i,] == 1)
    idx_mat[i, seq_len(length(idx_sel))] <- idx_sel
  }
  return(
    list(
      xtab = new_array,
      S1_cond_U_Z = S_cond_U_Z,
      S0_cond_U_Z = S_0_cond_U_Z,
      N_U = N_U,
      N_r = N_r,
      N_a = N_a,
      N_z = N_z,
      N_x = N_x,
      n_Y_p_z = n_Y_p_z,
      n_Y0_p_z = colSums(S_0_cond_U_Z),
      n_Y_p_U = n_Y_p_U,
      N_S = length(N_p_S),
      N_p_S = N_p_S_mat,
      theta_raw_metric = theta_raw_metric,
      beta_raw_metric = beta_raw_metric,
      p_a_raw_metric = p_a_raw_metric,
      beta_a_metric = beta_a_metric,
      sn_S_raw_metric = sn_S_raw_metric,
      sp_S_raw_metric = sp_S_raw_metric,
      sp_Y_raw_metric = sp_Y_raw_metric,
      p_y_1_tilde_raw_metric = p_y_1_tilde_raw_metric,
      sn_S_shape_1 = 1,
      sn_S_shape_2 = 1,
      sp_S_shape_1 = 1,
      sp_S_shape_2 = 1,
      sn_Y_shape_1 = 3,
      sn_Y_shape_2 = 2,
      sp_Y_shape_1 = 1,
      sp_Y_shape_2 = 1,
      prior_means = c(1, 0.5, rep(0, N_U - 3)),
      prior_sds = rep(1, N_U - 1),
      N_post = rowSums(S_cond_U_Z),
      N_gt_1 = sum(rowSums(S_cond_U_Z) > 1)
    )
  )
}

gen_obs_ve <- function(dat) {
  S_P_0 <- dat$S_P_0
  A <- dat$A
  Y <- dat$Y[S_P_0 == "111" & A == 1]
  Z <- dat$Z[S_P_0 == "111" & A == 1]
  p_1 <- mean(Y[Z == 1])
  p_2 <- mean(Y[Z == 2])
  p_3 <- mean(Y[Z == 3])
  ve_1 <- 1 - mean(Y[Z == 2]) / mean(Y[Z == 1])
  ve_2 <- 1 - mean(Y[Z == 3]) / mean(Y[Z == 1])
  return(list(ve_1 = ve_1, ve_2 = ve_2,
              p_1 = p_1, p_2 = p_2, p_3 = p_3))
}

set.seed(3)

par_list <- gen_param_set(N_p_r = 1e3, N_r = 8,
                          N_z = 3, N_a = 7,
                          N_x = 3,
                          Y_const_a = FALSE)
for (N_p_r_i in c(5, 10, 15)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 8, N_z = 3, N_a = 7, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.9, sn_S = 0.8, sn_Y = 0.99, par_list = par_list)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }

  saveRDS(dat_list, file = paste0(N_p_r_i, "k_design_Z_3_sparse_sp_Y_0_pt_9.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_3_sparse_noisy_sp_Y_0_pt_9.RDS"))
}

set.seed(5)

par_list <- gen_param_set(N_p_r = 1e3, N_r = 4,
                          N_z = 2, N_a = 3,
                          N_x = 3,
                          Y_const_a = FALSE)
for (N_p_r_i in c(1, 5, 10, 20)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 4, N_z = 2, N_a = 3, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.9, sn_S = 0.8, sn_Y = 0.99, par_list = par_list)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }

  saveRDS(dat_list, file = paste0(N_p_r_i,"k_design_Z_2_sparse_sp_Y_0_pt_9.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_2_sparse_noisy_sp_Y_0_pt_9.RDS"))
}

set.seed(5)

par_list <- gen_param_set(N_p_r = 1e3, N_r = 4,
                          N_z = 2, N_a = 3,
                          N_x = 3,
                          Y_const_a = FALSE)
for (N_p_r_i in c(1, 5, 10, 20)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 4, N_z = 2, N_a = 3, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.99, sn_S = 0.8, sn_Y = 0.8, par_list = par_list)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }

  saveRDS(dat_list, file = paste0(N_p_r_i,"k_design_Z_2_transmission_sparse.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_2_transmission_sparse_noisy.RDS"))
}

set.seed(5)

par_list <- gen_param_set(N_p_r = 1e3, N_r = 4,
                          N_z = 2, N_a = 3,
                          N_x = 3,
                          Y_const_a = FALSE,
                          null_hypo = TRUE)
for (N_p_r_i in c(1, 5, 10, 20)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 4, N_z = 2, N_a = 3, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.9, sn_S = 0.8, sn_Y = 0.99, par_list = par_list)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }

  saveRDS(dat_list, file = paste0(N_p_r_i, "k_design_Z_2_sparse_null_sp_Y_0_pt_9.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_2_sparse_noisy_null_sp_Y_0_pt_9.RDS"))
}

set.seed(3)

par_list <- gen_param_set(N_p_r = 1e3, N_r = 8,
                          N_z = 3, N_a = 7,
                          N_x = 3,
                          Y_const_a = FALSE,
                          null_hypo = TRUE)
for (N_p_r_i in c(5, 10, 15)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 8, N_z = 3, N_a = 7, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.9, sn_S = 0.8, sn_Y = 0.99, par_list = par_list)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }
  saveRDS(dat_list, file = paste0(N_p_r_i, "k_design_Z_3_sparse_null_sp_Y_0_pt_9.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_3_sparse_noisy_null_sp_Y_0_pt_9.RDS"))
}

set.seed(5)

par_list <- gen_param_set(N_p_r = 1e3, N_r = 4,
                          N_z = 2, N_a = 3,
                          N_x = 3,
                          Y_const_a = FALSE,
                          null_hypo = TRUE)
for (N_p_r_i in c(1, 5, 10, 20)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 4, N_z = 2, N_a = 3, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.99, sn_S = 0.8, sn_Y = 0.8, par_list = par_list)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }

  saveRDS(dat_list, file = paste0(N_p_r_i,"k_design_Z_2_transmission_sparse_null.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_2_transmission_sparse_noisy_null.RDS"))
}

set.seed(1)

par_list <- gen_waning_set(null_hypo = FALSE)
for (N_p_r_i in c(1, 5, 10, 20)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 4, N_z = 3, N_a = 3, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.99, sn_S = 0.8, sn_Y = 0.8, par_list = par_list,wane=TRUE)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }

  saveRDS(dat_list, file = paste0(N_p_r_i,"k_design_Z_2_wane_sparse.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_2_wane_sparse_noisy.RDS"))
}

set.seed(1)

par_list <- gen_waning_set(null_hypo = TRUE)
for (N_p_r_i in c(1, 5, 10, 20)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 4, N_z = 3, N_a = 3, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.99, sn_S = 0.8, sn_Y = 0.8, par_list = par_list,wane=TRUE)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }

  saveRDS(dat_list, file = paste0(N_p_r_i,"k_design_Z_2_wane_sparse_null.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_2_wane_sparse_noisy_null.RDS"))
}

set.seed(5)

par_list <- gen_param_set(N_p_r = 1e3, N_r = 4,
                          N_z = 2, N_a = 3,
                          N_x = 3,
                          Y_const_a = FALSE,
                          null_hypo = TRUE)
for (N_p_r_i in c(1, 5, 10, 20)) {
  dat_list <- list()
  dat_list_noisy <- list()
  N_p_r <- N_p_r_i * 1e3
  for (i in 1:100) {
    print(i)
    dat_1 <- gen_data(N_p_r = N_p_r, N_r = 4, N_z = 2, N_a = 3, N_x = 3, Y_const_a = FALSE, sp_S = 0.99, sp_Y = 0.99, sn_S = 0.8, sn_Y = 0.8, par_list = par_list)
    dat_proc <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    dat_proc_noise <- proc_data(dat_1$S_tilde, dat_1$Y_tilde, dat_1$A_obs, dat_1$Z, dat_1$idx,
      dat_1$S_cond_U_Z, dat_1$X,
      metric = rep(1, 217)
    )
    sel_nms <- c(
      "p_x_given_r",
      "p_r_given_x",
      "p_x_given_a",
      "prop_r_x",
      "prop_x"
    )
    dat_list[[i]] <- list(
      stan_data = c(dat_proc, dat_1[sel_nms]),
      pars = dat_1
    )
    dat_list_noisy[[i]] <- list(
      stan_data = c(dat_proc_noise, dat_1[sel_nms]),
      pars = dat_1
    )
  }

  saveRDS(dat_list, file = paste0(N_p_r_i,"k_design_Z_2_transmission_sparse_null.RDS"))
  saveRDS(dat_list_noisy, file = paste0(N_p_r_i, "k_design_Z_2_transmission_sparse_noisy_null.RDS"))
}


tt <- matrix(rnorm(1e5 * 3) * 2, 1e5, 3)
tt <- sweep(tt, 2, c(-2.5, -2.5, -2.5), FUN = "+")
tt_prior <- apply(tt, 1, softmax) |> t()
mean(1 - rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 3] == 0)]) /
  rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 0)]))
quantile(1 - rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 3] == 0)]) /
  rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 0)]))
hist(1 - rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 3] == 0)]) /
  rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 0)]))
apply(tt_prior, 2, quantile, c(0.01, 0.5, 0.99))
hist(tt_prior[, 1])
hist(tt_prior[, 2])
hist(tt_prior[, 3])
hist(tt_prior[, 4])
hist(rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 0)]))
hist(1 - rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 2] == 0)]) /
     rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 0)]))
hist(rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 3] == 0)]))
mean(1 - rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 2] == 1)]) /
  rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 1)]))
mean(1 - rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 3] == 1)]) /
  rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 1)]))
hist(1 - rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 2] == 1)]) /
  rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 1)]))
quantile(1 - rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 3] == 0)]) /
  rowSums(tt_prior[, which(dat_1$S_cond_U_Z[, 1] == 0)]))


tt_alt <- rdirichlet(10000, c(1,1,1,10))
hist(tt_alt[, 4])

tt_raw <- apply(tt_alt, 1, inv_softmax) |>
t()
hist(tt_raw[, 3], breaks = 100)
apply(tt_raw, 2, \(x) c(mean(x),sd(x)))

p_y_1 <- plogis(rnorm(10000) - 1)
p_y_2 <- plogis(rnorm(10000) - 1)
p_y_3 <- plogis(rnorm(10000) - 1)
hist(p_y_2)

mod_uk_wane <- cmdstan_model("covariate-meas-error-wane.stan")
mod_uk <- cmdstan_model("covariate-meas-error-a.stan")
mod_uk_o <- cmdstan_model("covariate-meas-error-a-equal-sn-sp.stan")
mod_uk_alt <- cmdstan_model("covariate-meas-error-a-trans-param.stan")
mod_uk_prior <- cmdstan_model("covariate-meas-error-equal-data-gen.stan")

gen_data_prior <- mod_uk_prior$sample(fixed_param = TRUE, iter_sampling = 100,
                                      data = dat_list[[1]]$stan_data)
tt <- gen_data_prior$draws(format = "draws_matrix")

dat <- readRDS("15k_design_Z_3_sparse.RDS")
fit <- mod_uk_wane$sample(data = dat_list[[1]]$stan_data,
                       parallel_chains = 4,
                       iter_warmup = 2000,
                       iter_sampling = 2000, 
                       init = 1,
                       refresh = 1e3,
                       adapt_delta = 0.9,
                       max_treedepth = 12)
for (i in 1:100) {
  tt <- mod_uk_o$optimize(data = dat[[2]]$stan_data, refresh=0)
  print(
    tt$mle()["sn_S"]
  )
  print(
    tt$lp()
  )
}

prop.table(table(dat[[92]]$pars$S, dat[[92]]$pars$S_tilde), 1)
prop.table(table(dat[[92]]$pars$Y, dat[[92]]$pars$Y_tilde), 1)
prop.table(table(dat[[2]]$pars$S, dat[[2]]$pars$S_tilde), 1)
prop.table(table(dat[[2]]$pars$Y, dat[[2]]$pars$Y_tilde), 1)
table(dat[[2]]$pars$S, dat[[2]]$pars$S_tilde)
table(dat[[2]]$pars$A, dat[[2]]$pars$S_tilde)
samps <- fit$draws(format = "draws_matrix")
fit$diagnostic_summary()
sum_fit <- fit$summary()
max_rhat <- max(sum_fit$rhat,na.rm=T)
min_ess_bulk <- min(sum_fit$ess_bulk,na.rm=T)
min_ess_bulk_i <- which.min(sum_fit$ess_bulk)
par_max_rhat <- sum_fit$variable[which.max(sum_fit$rhat)]
hist(samps[, "sn_S"])
plot(samps[,"sn_S"],samps[,"sp_Y"])
points(0.8, 0.99, col = "red", pch = 19)
plot(samps[, "sn_Y"], samps[, "lp__"])
plot(samps[, "theta[3,2,5]"], samps[, "p_a[3,1,1]"])
hist(samps[, "p_a[6,3,5]"], breaks =  100)
plot(samps[, "beta_raw[1,3]"], samps[, "theta_raw[3,1]"])
hist(samps[, "marg_ve_s[1]"])
abline(v = dat_list[[1]]$pars$marg_marg_ve[1], col = "red")
hist(samps[, "marg_ve_a_x[1]"])
abline(v = dat_list[[1]]$pars$marg_marg_ve[1], col="red")
hist(samps[,"marg_ve_a_x[2]"])
abline(v = dat_list[[1]]$pars$marg_marg_ve[2], col = "red")
hist(samps[, "marg_ve_a_x[2]"] - samps[, "marg_ve_a_x[1]"])
abline(v = dat_list[[1]]$pars$marg_marg_ve[2] - dat_list[[1]]$pars$marg_marg_ve[1], col = "red")
head(sort(sum_fit$ess_bulk))

theta <- rbeta(10000,0.7, 1)
