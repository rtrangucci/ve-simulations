library(cmdstanr)
library(dplyr)
library(stringr)
library(xtable)
library(ggplot2)
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
add_flair <- function(df, noise, study_type, z_dim, null, n) {
  df$noise <- noise
  df$study_type <- study_type
  df$z_dim <- z_dim
  df$null <- null
  df$n_p_site <- n * 1e3
  df$n_tot <- n * 1e3 * ifelse(z_dim == 3, 8, 4)
  df$idx <- 1:nrow(df)
  return(df)
}

read_files <- function(fl_dir) {
  fls <- list.files(fl_dir, full.names = TRUE)
  dfs_bias_z2 <- dfs_bias_sq_z2 <- dfs_cover_95_z2 <- dfs_post_pos_z2 <- list()
  dfs_bias_z3 <- dfs_bias_sq_z3 <- dfs_cover_95_z3 <- dfs_post_pos_z3 <- list()
  for (fl_i in seq_along(fls)) {
    fl <- fls[fl_i]
    noise <- ifelse(grepl("nois",fl), "A-w-noise", "A-wo-noise")
    study_type <- ifelse(grepl("trans", fl), "ve-t",
      ifelse(grepl("wane", fl), "ve-w", "ve-p")
    )
    null <- ifelse(grepl("null", fl), "null", "alt")
    z_dim <- ifelse(grepl("z3", fl), 3, 2)
    n <- readr::parse_number(fl)
    dfs <- readRDS(fl)
    print(fl)
    if (z_dim == 2) {
      dfs_bias_z2[[fl_i]] <- dfs[["bias"]] |>
      as.data.frame() |>
      add_flair(noise, study_type, z_dim, null, n)
      dfs_bias_sq_z2[[fl_i]] <- dfs[["bias"]]^2 |>
      as.data.frame() |>
      add_flair(noise, study_type, z_dim, null, n)
      dfs_cover_95_z2[[fl_i]] <- dfs[["cover_95"]] |>
      as.data.frame() |>
      add_flair(noise, study_type, z_dim, null, n)
      dfs_post_pos_z2[[fl_i]] <- data.frame(test_stat = dfs[["joint_hypo"]]) |>
      add_flair(noise, study_type, z_dim, null, n)
    } else {
      dfs_bias_z3[[fl_i]] <- dfs[["bias"]] |>
      as.data.frame() |>
      add_flair(noise, study_type, z_dim, null, n)
      dfs_bias_sq_z3[[fl_i]] <- dfs[["bias"]]^2 |>
      as.data.frame() |>
      add_flair(noise, study_type, z_dim, null, n)
      dfs_cover_95_z3[[fl_i]] <- dfs[["cover_95"]] |>
      as.data.frame() |>
      add_flair(noise, study_type, z_dim, null, n)
      dfs_post_pos_z3[[fl_i]] <- dfs[["joint_hypo"]] |>
      as.data.frame() |>
      add_flair(noise, study_type, z_dim, null, n)
    }
  }
  return(list(
    post_pos_z2 = do.call(rbind, dfs_post_pos_z2),
    bias_z2 = do.call(rbind, dfs_bias_z2),
    bias_sq_z2 = do.call(rbind, dfs_bias_sq_z2),
    cover_95_z2 = do.call(rbind, dfs_cover_95_z2),
    post_pos_z3 = do.call(rbind, dfs_post_pos_z3),
    bias_z3 = do.call(rbind, dfs_bias_z3),
    bias_sq_z3 = do.call(rbind, dfs_bias_sq_z3),
    cover_95_z3 = do.call(rbind, dfs_cover_95_z3)
  ))
}


dir <- "fits_new_sp"
res <- read_files(dir)
res$post_pos_z2 |>
  group_by(n_tot, noise, null) |>
  summarise(m = mean(test_stat > 0.9))

res$post_pos_z3 |>
group_by(n_tot, noise, null) |>
summarise(m = mean(V2 > 0.85))

res$post_pos_z2 |>
  filter(study_type == "ve-t") |>
  group_by(n_tot, noise, null) |>
  summarise(m = mean(test_stat > 0.925))  |>
  data.frame()

res$post_pos_z2 |>
filter(study_type == "ve-p", null == "alt") |>
group_by(n_tot, noise) |>
summarise(m = mean(`marg_ve_a_x[1]` > 0.9))  |>
  data.frame() |>
  tidyr::pivot_wider(values_from = "m",
                     names_from = "n_tot") |>
  xtable::xtable() |> print(include.rownames = FALSE)


res$post_pos_z3 |>
filter(study_type == "ve-p", null == "alt") |>
group_by(n_tot, noise) |>
summarise(m = mean(V2 > 0.9))  |>
  data.frame() |>
  tidyr::pivot_wider(values_from = "m",
                     names_from = "n_tot") |>
  xtable::xtable()|> print(include.rownames = FALSE)

res$post_pos_z2 |>
filter(study_type == "ve-p", null == "alt") |>
group_by(n_tot, noise) |>
  summarise(m = mean(test_stat > 0.9),
            sd = sd(test_stat > 0.9)/sqrt(n()),
            upper_ci = m + 1.96*sd,
            lower_ci = m - 1.96*sd)  |>
data.frame() |>
tidyr::pivot_wider(
  values_from = "m",
  names_from = "n_tot"
) |>
xtable::xtable()|> print(include.rownames = FALSE)

g <- res$post_pos_z2 |>
  filter(study_type == "ve-p", null == "alt") |>
  group_by(n_tot, noise) |>
  summarise(m = mean(test_stat > 0.9),
            sd = sd(test_stat > 0.9)/sqrt(n()),
            upper_ci = m + 1.96*sd,
            lower_ci = m - 1.96*sd)  |>
  left_join(data.frame(noise = c("A-w-noise","A-wo-noise"),
                       label = c("A'","A"))) |>
  ggplot(aes(x = n_tot, y = m, group = label, colour = label)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),width=10,size=1) +
  geom_hline(aes(yintercept = 0.8),linewidth=1,lty=2) +
  ylim(c(-0.01,1)) +
  theme_bw(base_size = 18) +
  xlab("Sample size") +
  ylab("Prop. trials rejecting the null") +
  scale_colour_manual(values = c("A'"="red","A"="blue"),name = "Cov. measurement") 
  
ggsave(filename = "Power.pdf",plot = g, width = 10, height = 10)

g <- res$post_pos_z2 |>
  filter(study_type == "ve-p", null == "alt") |>
  group_by(n_tot, noise) |>
  summarise(m = mean(test_stat > 0.9),
            sd = sd(test_stat > 0.9)/sqrt(n()),
            upper_ci = m + 1.96*sd,
            lower_ci = m - 1.96*sd)  |>
  left_join(data.frame(noise = c("A-w-noise","A-wo-noise"),
                       label = c("A'","A"))) |>
  ggplot(aes(x = n_tot, y = m, group = label, colour = label)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),width=10,size=1) +
  geom_hline(aes(yintercept = 0.8),linewidth=1,lty=2) +
  geom_vline(aes(xintercept = 4e4),linewidth=1,lty=1) +
  ylim(c(-0.01,1)) +
  theme_bw(base_size = 18) +
  xlab("Sample size") +
  ylab("Prop. trials rejecting the null") +
  scale_colour_manual(values = c("A'"="red","A"="blue"),name = "Cov. measurement") 

ggsave(filename = "Power-alt.pdf",plot = g, width = 10, height = 10)
  

g2 <- res$post_pos_z2 |>
  filter(study_type == "ve-p", null == "null") |>
  group_by(n_tot, noise) |>
  summarise(m = mean(test_stat > 0.9),
            sd = sd(test_stat > 0.9)/sqrt(n()),
            upper_ci = m + 1.96*sd,
            lower_ci = m - 1.96*sd)  |>
  left_join(data.frame(noise = c("A-w-noise","A-wo-noise"),
                       label = c("A'","A"))) |>
  ggplot(aes(x = n_tot, y = m, group = label, colour = label)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),width=10,size=1) +
  geom_hline(aes(yintercept = 0.05),linewidth=1,lty=2) +
  ylim(c(-0.01,1)) +
  theme_bw(base_size = 18) +
  xlab("Sample size") +
  ylab("Prop. trials rejecting the null") +
  scale_colour_manual(values = c("A'"="red","A"="blue"),name = "Cov. measurement") 
  
ggsave(filename = "Size.pdf",plot = g2, width = 10, height = 10)

res$post_pos_z3 |>
filter(study_type == "ve-p", null == "null") |>
group_by(n_tot, noise) |>
summarise(m = mean(V2 > 0.9))  |>
data.frame() |>
tidyr::pivot_wider(
  values_from = "m",
  names_from = "n_tot"
) |>
xtable::xtable()|> print(include.rownames = FALSE)

res$post_pos_z2 |>
filter(study_type == "ve-p", null == "null") |>
group_by(n_tot, noise) |>
summarise(m = mean(test_stat > 0.9))  |>
data.frame() |>
tidyr::pivot_wider(
  values_from = "m",
  names_from = "n_tot"
) |>
xtable::xtable()|> print(include.rownames = FALSE)

res$post_pos_z3 |>
filter(study_type == "ve-p", null == "null") |>
group_by(n_tot, noise) |>
summarise(m = mean(diff_ve_p_3_v_2 > 0.8))  |>
data.frame() |>
tidyr::pivot_wider(
  values_from = "m",
  names_from = "n_tot"
) |>
xtable::xtable()|> print(include.rownames = FALSE)

res$post_pos_z2 |>
filter(study_type == "ve-t", null == "alt") |>
group_by(n_tot, noise) |>
summarise(m = mean(test_stat > 0.925))  |>
data.frame()|>
tidyr::pivot_wider(
  values_from = "m",
  names_from = "n_tot"
) |>
xtable::xtable()|> print(include.rownames = FALSE)

res$post_pos_z2 |>
filter(study_type == "ve-t", null == "null") |>
group_by(n_tot, noise) |>
summarise(m = mean(test_stat > 0.925))  |>
data.frame()|>
tidyr::pivot_wider(
  values_from = "m",
  names_from = "n_tot"
) |>
xtable::xtable()|> print(include.rownames = FALSE)

res$post_pos_z2 |>
filter(study_type == "ve-t", null == "alt") |>
group_by(n_tot, noise) |>
  summarise(m1 = mean(`marg_ve_a_x[1]` > 0.925),
            m2 = mean(`marg_ve_S[1]` > 0.925))  |>
data.frame()|>
tidyr::pivot_wider(
  values_from = c("m1","m2"),
  names_from = "n_tot"
) |>
xtable::xtable()|> print(include.rownames = FALSE)
