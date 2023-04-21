library(posterior)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
print(args)
fits_dir <- args[1]
data_fl <- args[2]

fixed_data <- readRDS(data_fl)

save_nm <- strsplit(fits_dir,"/")[[1]][1]
pattern <- regexpr("[0-9]+", save_nm)
parse_size <- regmatches(save_nm, pattern)
print(paste0("size_parse: ",parse_size))

read_fits <- function(fit_dir) {
  fls <- list.files(path = fit_dir, full.names = TRUE)
  ret_list <- list()
  for (fl_i in seq_along(fls)) {
    fl <- fls[fl_i]
    fl_idx <- strsplit(fl, split = "_")[[1]] |> tail(1)
    fl_idx <- strsplit(fl_idx, split = ".RDS")[[1]][1] |>
      as.numeric()
    fit <- readRDS(fl)
    ret_list[[fl_idx]] <- fit
  }
  return(ret_list)
}

cover_quantiles <- function(draws, true, p) {
  lower <- (1 - p) / 2
  upper <- 1 - lower
  quants <- apply(draws, 2, quantile, c(lower, upper))
  n_pars <- ncol(draws)
  cover <- rep(NA_integer_, n_pars)
  length <- quants[2,] - quants[1,]
  for (idx_par in 1:n_pars) {
    cover_i <- 0
    if (true[idx_par] >= quants[1, idx_par] &&
      true[idx_par] <= quants[2, idx_par]) {
      cover_i <- 1
    }
    cover[idx_par] <- cover_i
  }
  return(
	list(
	cover = cover,
        length = length
	)
  )
}

mse_pars <- function(draws, true) {
  post_means <- colMeans(draws)
  post_vars <- apply(draws, 2, var)
  bias <- (post_means - true)
  return(list(
    bias = bias,
    vars = post_vars,
    mse = bias^2 + post_vars,
    means = post_means
  ))
}

names_suffix <- function(N) {
  return(paste0("[", 1:N, "]"))
}

parse_fits <- function(ensemble_fits, fixed_data,
                       thin_step = 1,parse_size) {
  print("Parsing fits")
  N_z <- fixed_data[[1]]$stan_data$N_z
  nms_suffix <- names_suffix(N_z-1)

  ve_i_nms <- paste0("marg_ve_S",nms_suffix)
  ve_p_nms <- paste0("marg_ve_a_x",nms_suffix)
  S <- length(ensemble_fits)
  L_tot <- ensemble_fits[[1]]$samps %>% nrow()
  L_tot <- L_tot / 4
  print("assumes 4 chains")
  treedepths <- divergences <- rep(NA_integer_, S)
  max_rhats <- min_ess_bulk <- min_ess_tail <- rep(NA_real_, S)
  fits_par_nms <- ensemble_fits[[1]]$samps %>%
    select(
      -`.chain`,
      -`.iteration`,
      -`.draw`
    ) %>%
    names()
  all_nms <- c(ve_i_nms, 
               ve_p_nms)
  if (N_z > 2) {
    nms_diff_ve_i <- paste0("diff_ve_i_",3:N_z,"_v_",2)
    nms_diff_ve_p <- paste0("diff_ve_p_",3:N_z,"_v_",2)
    all_nms <- c(all_nms, nms_diff_ve_i, nms_diff_ve_p)
  }
  
  n_vars <- length(all_nms)
  mse_mat <- matrix(NA_real_, S, n_vars)
  metadata <- list()

  colnames(mse_mat) <- all_nms
  post_means_mat <- bias_mat <- cover_mat_95 <- cover_mat_80 <- par_mat <- mse_mat
  for (idx_fit in seq_along(ensemble_fits)) {
    print(paste0("Parsing fit: ", idx_fit))
    fit <- ensemble_fits[[idx_fit]]
    fit_data <- fixed_data[[idx_fit]]$pars
      if (!is.null(fit)) {
        samps_comp <- fit$samps 
        sum_fit <- fit$sum_fit
        if (thin_step > 1) {
          samps_comp <- samps_comp %>%
            subset_draws(
              iteration = seq(1, L_tot, by = thin_step)
            )  
        }
        samps <- samps_comp %>%
          select(
            -`.chain`,
            -`.iteration`,
            -`.draw`
          ) %>%
          as.data.frame()

        true <- c(fit_data$marg_ve_S, fit_data$marg_marg_ve)
        names(true) <- c(ve_i_nms, ve_p_nms)
        true_it <- true
        if (N_z > 2) {
          diff_ve_i <- true[ve_i_nms]
          diff_ve_i <- diff_ve_i[2:(N_z-1)] - diff_ve_i[1]
          names(diff_ve_i) <- nms_diff_ve_i

          diff_ve_p <- true[ve_p_nms]
          diff_ve_p <- diff_ve_p[2:(N_z-1)] - diff_ve_p[1]
          names(diff_ve_p) <- nms_diff_ve_p

          true_it <- c(true, diff_ve_i, diff_ve_p)

          draws_diff_ve_i <- samps[,ve_i_nms[2:(N_z-1)]] - samps[,ve_i_nms[1]] |>
            as.matrix()
          colnames(draws_diff_ve_i) <- names(diff_ve_i)

          draws_diff_ve_p <- samps[, ve_p_nms[2:(N_z - 1)]] - samps[, ve_p_nms[1]]|>
          as.matrix()
          colnames(draws_diff_ve_p) <- names(diff_ve_p)

          samps <- cbind(samps, draws_diff_ve_i, draws_diff_ve_p)
        }

        stopifnot(all(names(samps) == names(true_it)))
        stopifnot(all(names(samps) == names(mse_mat)))
        zs <- mse_pars(samps, true_it)
        cover_80 <- cover_quantiles(samps, true_it, 0.8)
        cover_95 <- cover_quantiles(samps, true_it, 0.95)
        bias_mat[idx_fit, ] <- zs$bias
        mse_mat[idx_fit, ] <- zs$mse
        post_means_mat[idx_fit, ] <- zs$means
        cover_mat_95[idx_fit, ] <- cover_95$cover
        cover_mat_80[idx_fit, ] <- cover_80$cover
        diags_s <- fit$diags
        div <- sum(diags_s$divergent__ == 1)
        treedepth <- sum(diags_s$treedepth__ == 14)
        par_mat[idx_fit, ] <- true_it
        treedepths[idx_fit] <- treedepth
        divergences[idx_fit] <- div
        max_rhats[idx_fit] <- sum_fit$rhat %>% max
        min_bulk <- sum_fit$ess_bulk %>% min()
        min_tail <- sum_fit$ess_tail %>% min(na.rm=T)
        samps_denom <- nrow(fit$samps)
        min_ess_bulk[idx_fit] <- min_bulk / samps_denom
        min_ess_tail[idx_fit] <- min_tail / samps_denom
      }
  }
  metadata$S <- nrow(samps)
  return(
    list(
      cover_95 = cover_mat_95,
      cover_80 = cover_mat_80,
      mse = mse_mat,
      true_pars = par_mat,
      bias = bias_mat,
      treedepths = treedepths,
      divergences = divergences,
      max_rhats = max_rhats,
      min_ess_tail = min_ess_tail,
      min_ess_bulk = min_ess_bulk,
      metadata = metadata,
      parse_size = parse_size,
      means = post_means_mat
    )
  )
}

fits <- read_fits(fits_dir)
fixed_data <- readRDS(data_fl)
parsed <- parse_fits(fits, fixed_data, thin_step = 1, parse_size)
saveRDS(parsed, file = paste0("parsed_fits_",save_nm,".RDS"))
