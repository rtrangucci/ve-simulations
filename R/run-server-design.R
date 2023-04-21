args <- commandArgs(trailingOnly = TRUE)
data_fl <- args[1]
data_idx <- as.numeric(args[2])
mod_nm <-  args[3]
prep <- args[4]
seed <- data_idx
set.seed(seed)
print(args)
data_fit <- readRDS(data_fl)[[data_idx]]
inf_model <- readRDS(mod_nm)

options(CMDSTANR_NO_VER_CHECK = TRUE)
print("loading cmdstanr")
library(cmdstanr, quietly = T)
print("loaded cmdstanr")
library(dplyr)
library(posterior)

print("loading fit_mod")
fit_mod <- function(model, data) {
  fit <- model$sample(
    data = data$stan_data,
    iter_warmup = 6000,
    iter_sampling = 6000,
    refresh = 3000,
    chains = 4, parallel_chains = 4,
    adapt_delta = 0.90,
    init = 0.1,
    max_treedepth = 12
  )
  pars <- posterior::as_draws_df(
    fit$draws(
                variables =   c("marg_ve_S","marg_ve_a_x")
    )
  )
  samp_pars <- fit$sampler_diagnostics() %>%
    as_draws_df()
  return(
    list(
      samps = pars,
      diags = samp_pars,
      metadata = fit$metadata(),
      sum_fit = fit$summary()
    )
  )
}

print("loaded fit_mod")

print("fitting-model")
fit <- fit_mod(inf_model, data_fit)

saveRDS(fit, paste0("/scratch/stats_dept_root/stats_dept1/trangucc/",prep,"_",data_idx,".RDS"))
