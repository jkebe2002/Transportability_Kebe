# Estimators for the interpolation simulation
#
# dr_balance() returns, for ONE simulated dataset, a flat named list of fits
# (one per method x t). Names come from fit_labels(t_grid) and are interleaved:
#   OT_t_<t1>, OT_interpol_<t1>, OT_t_<t2>, OT_interpol_<t2>, ...
#
#   OT_t        : weights = t * 1 + (1 - t) * (OT weights to target)
#   OT_interpol : weights from solving  t * OT(source_w, source) + (1 - t) * OT(source_w, target)
#
# t = 1 -> no transport (all weights 1); t = 0 -> pure OT weights.

###=============================================================================
### Packages
###=============================================================================

library(eventglm)
library(survival)       # Surv() in the pseudo-observation formula
library(SuperLearner)
library(causalOT)

###=============================================================================
### Helpers
###=============================================================================

# diagnostics stored for every fit (order matters; used by the simulation script)
diag_names <- c("ot_target", "ot_source", "ot_source_w", "tipton",
                "minw", "maxw", "medw")

fit_labels <- function(t_grid) {
  as.vector(rbind(paste0("OT_t_", t_grid), paste0("OT_interpol_", t_grid)))
}

failed_fit <- function(n_cuts) {
  list(estimate = rep(NA_real_, n_cuts),
       variance = rep(NA_real_, n_cuts),
       wfit     = NULL,
       diag     = setNames(rep(NA_real_, length(diag_names)), diag_names))
}

# Binned Tipton generalizability (B) index: Bhattacharyya coefficient between the
# weighted trial and the target distributions of the logit of the sampling
# propensity. 1 = perfect overlap. logit_all is for everyone (S = 1 and S = 0).
tipton_index <- function(logit_all, S, w, nbins = 20) {
  if (anyNA(logit_all)) return(NA_real_)
  brks <- unique(quantile(logit_all, seq(0, 1, length.out = nbins + 1)))
  if (length(brks) < 3) return(NA_real_)
  brks[1] <- -Inf; brks[length(brks)] <- Inf
  p_t <- as.numeric(prop.table(table(cut(logit_all[S == 0], brks))))
  p_s <- as.numeric(tapply(w, cut(logit_all[S == 1], brks), sum))
  p_s[is.na(p_s)] <- 0
  p_s <- p_s / sum(p_s)
  sum(sqrt(p_s * p_t))
}

# OT distance between two weighted point clouds. Weights are normalized to
# probability vectors. CHECK ?ot_distance: confirm that this call returns the
# W2 quantity you want (cost exponent p = 2) and not only an entropic penalty term.
ot_dist <- function(x1, x2, a, b) {
  tryCatch(
    ot_distance(x1, x2, a = a / sum(a), b = b / sum(b), penalty = 1),
    error = function(e) NA_real_
  )
}

# Solve an OTProblem with the settings used throughout the project
solve_otp <- function(otp) {
  otp$setup_arguments(lambda = 4, delta = NULL, grid.length = 7L,
                      cost.function = NULL, p = 2, cost.online = "auto",
                      debias = TRUE, diameter = NULL,
                      ot_niter = 1000L, ot_tol = 1e-3)
  otp$solve(
    niter = 1000L,
    tol = 1e-5,
    optimizer = "torch",
    torch_optim = torch::optim_lbfgs,
    torch_scheduler = torch::lr_reduce_on_plateau,
    torch_args = list(line_search_fn = "strong_wolfe"),
    osqp_args = NULL,
    quick.balance.function = TRUE
  )
  invisible(otp)
}

# Weights (mean 1 over the n1 source rows) from the interpolated OT problem
#   tt * OT(source_w, source) + (1 - tt) * OT(source_w, target)
# Fresh Measure objects are built on every call so nothing is warm-started
# from a previous solve (Measure/OTProblem are reference objects).
fit_ot_weights <- function(Xmat1, Xmat0, tt) {
  n1 <- nrow(Xmat1)
  if (tt >= 1) return(rep(1, n1))
  
  m1_w <- causalOT::Measure(x = Xmat1, adapt = "weights")   # weights are learned
  m1   <- causalOT::Measure(x = Xmat1, adapt = "none")
  m0   <- causalOT::Measure(x = Xmat0, adapt = "none")
  
  otp_target <- causalOT::OTProblem(m1_w, m0)
  otp <- if (tt <= 0) otp_target else
    tt * causalOT::OTProblem(m1_w, m1) + (1 - tt) * otp_target
  
  solve_otp(otp)
  n1 * as.numeric(m1_w$weights)
}

###=============================================================================
### Doubly-robust transportability estimator over the interpolation grid
###=============================================================================

dr_balance <- function(S, X, Z1, time1, event1, cuts,
                       t_grid = seq(0, 1, by = 0.25),
                       sl.lib = c("SL.mean", "SL.glm")) {
  
  X  <- as.data.frame(X)
  n  <- nrow(X)
  n1 <- sum(S == 1)
  X1 <- X[S == 1, , drop = FALSE]
  X0 <- X[S == 0, , drop = FALSE]
  Xmat1 <- as.matrix(X1)
  Xmat0 <- as.matrix(X0)
  prZ <- mean(Z1)
  
  ## ---- Computed ONCE per replicate (none of this depends on the weights) -----
  
  # pseudo-observations for survival at each cut
  ps_dat <- data.frame(t_ = time1, e_ = event1)
  psi1.mat <- vapply(cuts, function(tt) {
    eventglm::pseudo_infjack(
      formula = Surv(t_, e_) ~ 1,
      time = tt,
      cause = 1,
      type = "survival",
      data = ps_dat
    )
  }, numeric(n1))
  
  # outcome-model predictions under Z = 0 / Z = 1 for everyone, per cut
  mu_list <- lapply(seq_along(cuts), function(i) {
    out <- SuperLearner(Y = psi1.mat[, i], X = data.frame(X1, Z = Z1),
                        SL.library = sl.lib, family = gaussian())
    cbind(predict(out, newdata = data.frame(X, Z = 0))$pred,
          predict(out, newdata = data.frame(X, Z = 1))$pred)
  })
  
  # logit of the sampling propensity for the Tipton index
  logit_all <- tryCatch(
    suppressWarnings(predict(glm(S ~ ., data = data.frame(S = S, X),
                                 family = binomial), type = "link")),
    error = function(e) rep(NA_real_, n)
  )
  
  # full-length treatment indicator (0 outside the trial)
  Zfull <- numeric(n); Zfull[S == 1] <- Z1
  
  ## ---- Weighted DR estimate for one set of source weights -------------------
  est_effect <- function(w) {
    samp_full <- numeric(n); samp_full[S == 1] <- w
    w0 <- (1 - Zfull) * S * samp_full / (1 - prZ)
    w1 <- Zfull * S * samp_full / prZ
    
    res <- lapply(seq_along(cuts), function(i) {
      psi <- numeric(n); psi[S == 1] <- psi1.mat[, i]
      mu_tmp <- mu_list[[i]]
      mu <- cbind(Zfull * mu_tmp[, 2] + (1 - Zfull) * mu_tmp[, 1], mu_tmp)
      
      aug_est <- sum((w1 - w0) * (psi - mu[, 1])) / n1 +
        mean(mu[S == 0, 3] - mu[S == 0, 2])
      eic <- ((w1 - w0) * (psi - mu[, 1])) / mean(S == 1) +
        ((S == 0) * (mu[, 3] - mu[, 2] - aug_est)) / mean(S == 0)
      
      c(estimate = aug_est, variance = var(eic) / n)
    })
    
    list(estimate = vapply(res, `[[`, numeric(1), "estimate"),
         variance = vapply(res, `[[`, numeric(1), "variance"))
  }
  
  ## ---- Pure OT weights (t = 0 end of both methods) --------------------------
  wts_OT <- fit_ot_weights(Xmat1, Xmat0, 0)
  
  ## ---- Estimate + diagnostics for one weight vector -------------------------
  make_fit <- function(w) {
    r <- est_effect(w)
    list(
      estimate = r$estimate,
      variance = r$variance,
      wfit = list(weights = w, ESS = sum(w)^2 / sum(w^2)),
      diag = c(
        ot_target   = ot_dist(Xmat1, Xmat0, w, rep(1, nrow(Xmat0))),  # vs target
        ot_source   = ot_dist(Xmat1, Xmat1, w, rep(1, n1)),           # vs unweighted source
        ot_source_w = ot_dist(Xmat1, Xmat1, w, wts_OT),               # vs OT-weighted source
        tipton      = tipton_index(logit_all, S, w),
        minw = min(w), maxw = max(w), medw = median(w)
      )
    )
  }
  
  safe_make <- function(expr_fun) {
    tryCatch(expr_fun(), error = function(e) { print(e); failed_fit(length(cuts)) })
  }
  
  ## ---- Loop over the interpolation grid -------------------------------------
  fits <- setNames(vector("list", 2 * length(t_grid)), fit_labels(t_grid))
  
  for (j in seq_along(t_grid)) {
    tt <- t_grid[j]
    
    # OT_t: linear interpolation between unit weights and OT weights
    fits[[2 * j - 1]] <- safe_make(function() make_fit(tt * 1 + (1 - tt) * wts_OT))
    
    # OT_interpol: interpolate the OT problems, then solve
    fits[[2 * j]] <- safe_make(function() {
      w_int <- if (tt <= 0) wts_OT else fit_ot_weights(Xmat1, Xmat0, tt)
      make_fit(w_int)
    })
  }
  
  fits
}