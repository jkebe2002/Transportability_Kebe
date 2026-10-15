###=============================================================================
### Data generating process for the LEADER transportability simulation
###=============================================================================

weibull_shape <- 1.5
lambda0 <- 0.2

selection_lp <- function(x1, x2, scenario, u = 0) {
  switch(scenario,
         baseline = 0.1 + 1.4 * x1 - 1.0 * x2,
         omitted_selection = 0.1 + 1.4 * x1 - 1.0 * x2 + 1.0 * u,
         omitted_outcome = 0.1 + 1.4 * x1 - 1.0 * x2,
         omitted_both = 0.1 + 1.4 * x1 - 1.0 * x2 + 1.0 * u,
         both_nonlinear = 0.15 + 1.2 * x1 - 0.9 * x2 - 0.35 * x1^2 - 0.4 * x1 * x2,
         positivity_practical = -0.4 + 2.0 * x1 - 2.0 * x2,
         positivity_structural = -0.4 + 0.5 * x1 - 0.5 * x2,
         outcome_extrapolation = -0.1 + 0.9 * x1 - 0.7 * x2
  )
}

outcome_lp <- function(x1, x2, A, scenario, u = 0) {
  base <- 0.3 * x1 - 0.4 * x2 - 0.5 * A
  switch(scenario,
         baseline = base,
         omitted_selection = base,
         omitted_outcome = base + 0.6 * u + 0.8 * u * A,
         omitted_both = base + 0.6 * u + 0.8 * u * A,
         both_nonlinear = base + 0.4 * x1^2 - 0.3 * x1 * x2 + 0.8 * x1 * A + 0.6 * x1^2 * A,
         positivity_practical = base,
         positivity_structural = base,
         outcome_extrapolation = base + 0.6 * pmax(x1 - 0.5, 0)
  )
}

true_surv <- function(t, A, x1, x2, scenario, u = 0) {
  lp <- outcome_lp(x1, x2, A, scenario, u)
  exp(-lambda0 * exp(lp) * t^weibull_shape)
}

generate_data <- function(n_total = 1000, cuts, scenario = "baseline") {

  x1 <- rnorm(n_total, mean = -0.5, sd = 1)
  x2 <- rbinom(n_total, 1, 0.5)

  u <- rnorm(n_total)

  if (scenario == "outcome_extrapolation") {
    shift <- rnorm(n_total, mean = 1.0, sd = 1.2)
    pi_trial <- plogis(selection_lp(x1, x2, scenario))
    S <- rbinom(n_total, 1, pi_trial)
    x1 <- ifelse(S == 0, shift, x1)
  } else {
    pi_trial <- plogis(selection_lp(x1, x2, scenario, u))
    S <- rbinom(n_total, 1, pi_trial)
    if (scenario == "positivity_structural") S[x1 > -0.5] <- 0
  }

  A <- rbinom(n_total, 1, 0.5)

  v <- runif(n_total)
  lp_t <- outcome_lp(x1, x2, A, scenario, u)
  T_event <- (-log(v) / (lambda0 * exp(lp_t)))^(1 / weibull_shape)
  C_time <- rexp(n_total, rate = exp(-3))

  tgt <- which(S == 0)
  truth <- do.call(rbind, lapply(cuts, function(t) {
    s0 <- mean(true_surv(t, A = 0, x1 = x1[tgt], x2 = x2[tgt], scenario, u[tgt]))
    s1 <- mean(true_surv(t, A = 1, x1 = x1[tgt], x2 = x2[tgt], scenario, u[tgt]))
    data.frame(time = t, surv0 = s0, surv1 = s1, diff = s1 - s0)
  }))

  time <- pmin(T_event, C_time)
  status <- as.numeric(T_event <= C_time)

  df <- data.frame(id = seq_len(n_total), S, A, x1, x2, time, status)
  df$A[df$S == 0] <- NA
  df$time[df$S == 0] <- NA
  df$status[df$S == 0] <- NA

  list(df = df, truth = truth)

}
