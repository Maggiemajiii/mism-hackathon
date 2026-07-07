## =========================================================================
## fit_spatial_sir(): full pipeline
##   - builds transmission/recovery datasets from before/after-drug tracks
##   - fits MLE (bbmle::mle2), prints convergence code
##   - fits Bayesian posterior (mcmc::metrop), auto-tuning acceptance rate
##   - returns both sets of estimates
##   - produces 5 plots:
##       1. df_before, simulated (MLE params)     vs observed I(t)
##       2. df_before, simulated (Bayesian params) vs observed I(t)
##       3. df_after,  simulated (MLE params)     vs observed I(t)
##       4. df_after,  simulated (Bayesian params) vs observed I(t)
##       5. parameter estimates: MLE (dark grey) vs Bayesian (dark orange),
##          with a vertical dashed line at 0, one row per parameter
## =========================================================================

library(dplyr)
library(bbmle)
library(mcmc)
library(ggplot2)

fit_spatial_sir <- function(df_before, df_after, delta = 1,
                             n_mcmc = 10000, n_sims = 30, burn = 1000,
                             tune_target = c(0.2, 0.4), tune_max_iter = 15,
                             mult_start = 9) {

  ## ---- 0. Combine data ----
  df_before$drug <- 0
  df_after$drug  <- 1
  tracks <- bind_rows(df_before, df_after) %>% arrange(drug, cell, frame)

  ## ---- 1. Build transmission blocks (distance matrices, precomputed once) ----
  build_transmission_blocks <- function(tracks, delta) {
    blocks <- list(); k <- 1
    for (d in unique(tracks$drug)) {
      sub <- tracks %>% filter(drug == d)
      frames <- sort(unique(sub$frame))
      for (t in frames) {
        if (!(t + delta) %in% frames) next
        now <- sub %>% filter(frame == t)
        nxt <- sub %>% filter(frame == t + delta) %>% select(cell, state_next = state)
        susc <- now %>% filter(state == 1) %>% left_join(nxt, by = "cell")
        inf  <- now %>% filter(state == 2)
        if (nrow(susc) == 0 || nrow(inf) == 0) next
        D <- sqrt(outer(susc$x, inf$x, "-")^2 + outer(susc$y, inf$y, "-")^2)
        blocks[[k]] <- list(D = D, became_infected = as.numeric(susc$state_next == 2), drug = d)
        k <- k + 1
      }
    }
    blocks
  }

  build_recovery_data <- function(tracks, delta) {
    out <- list(); k <- 1
    for (d in unique(tracks$drug)) {
      sub <- tracks %>% filter(drug == d)
      frames <- sort(unique(sub$frame))
      for (t in frames) {
        if (!(t + delta) %in% frames) next
        now <- sub %>% filter(frame == t, state == 2)
        nxt <- sub %>% filter(frame == t + delta) %>% select(cell, state_next = state)
        if (nrow(now) == 0) next
        now <- now %>% left_join(nxt, by = "cell")
        out[[k]] <- data.frame(became_removed = as.numeric(now$state_next == 3), drug = d)
        k <- k + 1
      }
    }
    bind_rows(out)
  }

  trans_blocks <- build_transmission_blocks(tracks, delta)
  recov_df     <- build_recovery_data(tracks, delta)

  cat("Number of (frame, condition) transmission blocks:", length(trans_blocks), "\n")
  cat("Number of recovery observations:", nrow(recov_df), "\n")

  ## ---- 2. Negative log-likelihood ----
  nll <- function(log_beta1, log_phi, betaD, log_gamma, gammaD) {
    beta1 <- exp(log_beta1); phi <- exp(log_phi); gamma <- exp(log_gamma)

    ll_trans <- 0
    for (blk in trans_blocks) {
      kappa_sum <- rowSums(exp(-blk$D / phi))
      lambda    <- beta1 * exp(betaD * blk$drug) * kappa_sum
      p_inf     <- 1 - exp(-lambda * delta)
      p_inf     <- pmin(pmax(p_inf, 1e-10), 1 - 1e-10)
      ll_trans  <- ll_trans + sum(dbinom(blk$became_infected, 1, p_inf, log = TRUE))
    }

    p_rec  <- 1 - exp(-gamma * exp(gammaD * recov_df$drug) * delta)
    p_rec  <- pmin(pmax(p_rec, 1e-10), 1 - 1e-10)
    ll_rec <- sum(dbinom(recov_df$became_removed, 1, p_rec, log = TRUE))

    -(ll_trans + ll_rec)
  }

  ## ---- 3. MLE fit ----
  start_vals <- list(log_beta1 = log(0.01), log_phi = log(3), betaD = 0,
                      log_gamma = log(0.02), gammaD = 0)

  fit <- mle2(nll, start = start_vals, method = "BFGS", control = list(maxit = 500))

  cat("\nMLE convergence code (0 = converged):\n")
  print(fit@details$convergence)

  mle_est <- coef(fit)
  mle_se  <- sqrt(diag(vcov(fit)))

  ## ---- 4. Bayesian fit: Metropolis-Hastings via mcmc::metrop ----
  log_posterior <- function(par) {
    log_beta1 <- par[1]; log_phi <- par[2]; betaD <- par[3]
    log_gamma <- par[4]; gammaD <- par[5]

    ll <- -nll(log_beta1, log_phi, betaD, log_gamma, gammaD)

    lp <- dnorm(log_beta1, log(0.01), 3, log = TRUE) +
          dnorm(log_phi,   log(3),    3, log = TRUE) +
          dnorm(betaD,     0,         5, log = TRUE) +
          dnorm(log_gamma, log(0.02), 3, log = TRUE) +
          dnorm(gammaD,    0,         5, log = TRUE)

    ll + lp
  }

  init  <- coef(fit)
  Sigma <- vcov(fit)

  # Tune proposal scale toward target acceptance rate by CONTINUING the
  # chain (not restarting), matching the confirmed-working manual approach.
  # Tune proposal scale toward target acceptance rate by CONTINUING the
  # chain (not restarting). Starts from mult_start (default 8, since that's
  # been confirmed to work well across your datasets) rather than always
  # searching from scratch -- the loop still runs in case a given dataset
  # needs slight adjustment, but usually should exit on the first check.
  mult <- mult_start
  out <- metrop(log_posterior, initial = init, nbatch = 1000, scale = Sigma * mult)
  for (i in 1:tune_max_iter) {
    acc <- out$accept
    if (acc >= tune_target[1] && acc <= tune_target[2]) break
    mult <- mult * ifelse(acc > tune_target[2], 1.5, 0.7)
    out <- metrop(out, scale = Sigma * mult)
  }
  cat("\nTuned MCMC acceptance rate:", round(out$accept, 3),
      "(scale multiplier:", round(mult, 3), ")\n")

  # Production run: continue from the tuned chain (do NOT restart from init;
  # continuing is what's been confirmed to work).
  out_long <- metrop(out, nbatch = n_mcmc, scale = Sigma * mult)
  cat("Final long-run acceptance rate:", round(out_long$accept, 3), "\n")

  post <- as.data.frame(out_long$batch)
  colnames(post) <- names(init)
  post <- post[-(1:min(burn, nrow(post) - 1)), ]

  bayes_mean <- colMeans(post)
  bayes_ci   <- apply(post, 2, quantile, probs = c(0.025, 0.975))

  ## ---- 5. Simulator for posterior-predictive-style fit plots ----
  simulate_epidemic <- function(tracks_template, params, delta) {
    frames   <- sort(unique(tracks_template$frame))
    drug_val <- unique(tracks_template$drug)
    state <- tracks_template %>% filter(frame == frames[1]) %>% select(cell, x, y, state)

    sim_states <- list(state)
    for (t in frames[-length(frames)]) {
      coords_now <- tracks_template %>% filter(frame == t) %>% select(cell, x, y)
      state <- state %>% select(-x, -y) %>% left_join(coords_now, by = "cell")

      susc <- state %>% filter(state == 1)
      inf  <- state %>% filter(state == 2)
      if (nrow(susc) > 0 && nrow(inf) > 0) {
        D <- sqrt(outer(susc$x, inf$x, "-")^2 + outer(susc$y, inf$y, "-")^2)
        kappa_sum <- rowSums(exp(-D / params["phi"]))
        lambda    <- params["beta1"] * exp(params["betaD"] * drug_val) * kappa_sum
        p_inf     <- 1 - exp(-lambda * delta)
        newly_infected <- susc$cell[runif(nrow(susc)) < p_inf]
        state$state[state$cell %in% newly_infected] <- 2
      }

      infd <- state %>% filter(state == 2)
      if (nrow(infd) > 0) {
        p_rec <- 1 - exp(-params["gamma"] * exp(params["gammaD"] * drug_val) * delta)
        newly_removed <- infd$cell[runif(nrow(infd)) < p_rec]
        state$state[state$cell %in% newly_removed] <- 3
      }
      sim_states[[length(sim_states) + 1]] <- state
    }
    for (i in seq_along(sim_states)) sim_states[[i]]$frame <- frames[i]
    bind_rows(sim_states)
  }

  ## Fixed-parameter simulation set: used for MLE plots, where there is only
  ## a single point estimate (no posterior to sample from). Variability across
  ## the n_sims curves here reflects stochastic simulation noise only.
  sim_curve_set_fixed <- function(tracks_template, params, n_sims, delta) {
    results <- vector("list", n_sims)
    for (s in 1:n_sims) {
      sim <- simulate_epidemic(tracks_template, params, delta)
      curve <- sim %>% group_by(frame) %>% summarise(I = sum(state == 2), .groups = "drop")
      curve$sim_id <- s
      results[[s]] <- curve
    }
    bind_rows(results)
  }

  ## Posterior-sampling simulation set: used for Bayesian plots. A fresh
  ## parameter draw is taken from the posterior for EACH simulation, so the
  ## resulting band reflects both parameter uncertainty and stochastic
  ## simulation noise -- this is the version that matches what a proper
  ## posterior-predictive check should look like.
  sim_curve_set_posterior <- function(tracks_template, post, n_sims, delta) {
    results <- vector("list", n_sims)
    for (s in 1:n_sims) {
      draw <- post[sample(nrow(post), 1), ]
      params <- c(beta1 = exp(draw$log_beta1), phi = exp(draw$log_phi),
                  betaD = draw$betaD, gamma = exp(draw$log_gamma), gammaD = draw$gammaD)
      sim <- simulate_epidemic(tracks_template, params, delta)
      curve <- sim %>% group_by(frame) %>% summarise(I = sum(state == 2), .groups = "drop")
      curve$sim_id <- s
      results[[s]] <- curve
    }
    bind_rows(results)
  }

  make_fit_plot <- function(sims, tracks_template, title) {
    obs  <- tracks_template %>% group_by(frame) %>% summarise(I = sum(state == 2), .groups = "drop")
    ggplot() +
      geom_line(data = sims, aes(frame, I, group = sim_id), alpha = 0.15, color = "steelblue") +
      geom_line(data = obs, aes(frame, I), color = "black", linewidth = 1) +
      labs(title = title, x = "frame", y = "I(t)") +
      theme_minimal()
  }

  params_mle <- c(beta1 = unname(exp(mle_est["log_beta1"])), phi = unname(exp(mle_est["log_phi"])),
                   betaD = unname(mle_est["betaD"]), gamma = unname(exp(mle_est["log_gamma"])),
                   gammaD = unname(mle_est["gammaD"]))

  sims_before_mle   <- sim_curve_set_fixed(df_before, params_mle, n_sims, delta)
  sims_after_mle    <- sim_curve_set_fixed(df_after,  params_mle, n_sims, delta)
  sims_before_bayes <- sim_curve_set_posterior(df_before, post, n_sims, delta)
  sims_after_bayes  <- sim_curve_set_posterior(df_after,  post, n_sims, delta)

  plot1 <- make_fit_plot(sims_before_mle,   df_before, "Before drug \u2014 MLE fit")
  plot2 <- make_fit_plot(sims_before_bayes, df_before, "Before drug \u2014 Bayesian fit")
  plot3 <- make_fit_plot(sims_after_mle,    df_after,  "After drug \u2014 MLE fit")
  plot4 <- make_fit_plot(sims_after_bayes,  df_after,  "After drug \u2014 Bayesian fit")

  ## ---- 6. Parameter comparison plot ----
  param_names <- c("log_beta1", "log_phi", "betaD", "log_gamma", "gammaD")

  param_df <- data.frame(
    param    = rep(param_names, 2),
    method   = rep(c("MLE", "Bayesian"), each = 5),
    estimate = c(unname(mle_est[param_names]), unname(bayes_mean[param_names])),
    lower    = c(unname(mle_est[param_names] - 1.96 * mle_se[param_names]),
                 unname(bayes_ci[1, param_names])),
    upper    = c(unname(mle_est[param_names] + 1.96 * mle_se[param_names]),
                 unname(bayes_ci[2, param_names]))
  )
  param_df$param <- factor(param_df$param, levels = param_names)

  plot5 <- ggplot(param_df, aes(x = param, y = estimate, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), width = 0.15, linewidth = 1) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    scale_color_manual(values = c(MLE = "darkgrey", Bayesian = "darkorange")) +
    coord_flip() +
    labs(title = "Parameter estimates: MLE vs Bayesian (mean + 95% interval)",
         x = NULL, y = "Estimate (fitted scale)") +
    theme_minimal()

  ## ---- 7. Return everything ----
  list(
    mle_fit          = fit,
    mle_estimates    = mle_est,
    mle_se           = mle_se,
    bayes_posterior  = post,
    bayes_mean       = bayes_mean,
    bayes_ci         = bayes_ci,
    trans_blocks     = trans_blocks,
    recov_df         = recov_df,
    plots = list(
      before_mle        = plot1,
      before_bayes      = plot2,
      after_mle         = plot3,
      after_bayes       = plot4,
      param_comparison  = plot5
    )
  )
}

## =========================================================================
## Usage
## =========================================================================
df_before <- read.csv("trajectory_no_drug_0.csv")
df_after  <- read.csv("trajectory_with_drug_0.csv")

results <- fit_spatial_sir(df_before, df_after)

results$mle_estimates
results$bayes_mean
results$bayes_ci

results$plots$before_mle
results$plots$before_bayes
results$plots$after_mle
results$plots$after_bayes
results$plots$param_comparison
