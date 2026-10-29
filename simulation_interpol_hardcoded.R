# interpolation simulation

# what are you doing at each run?

# sample both distributions

# calculate the OT weights 

# interpolate
#   what are you doing at each t?

#   create interpolated weights t*1 + (1-t) * OT weights

#   create OT_interpolated weights: t*otp1+(1-t)*otp0

#   Compare the weight distributions, the ESS, the W2 distance from the target, the source, and the OT-weighted source, Tipton index, outcome and variance 





###=============================================================================
### Run the simulation, aggregate, and plot
###=============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

out_dir <- file.path("~", "Documents/GitHub", "Transportability_Kebe")
out_dir <- normalizePath(out_dir, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(out_dir, "01_dgp.R"))
source(file.path(out_dir, "simulation_interpol_effects.R"))

scenarios <- c("baseline",
               "omitted_selection",
               "omitted_outcome",
               "omitted_both",
               "both_nonlinear",
               "positivity_practical",
               "positivity_structural") #,
#"outcome_extrapolation")
n_sim <- 500; n_total <- 1000; cuts <- 1:6
method_labels <- c("OT_t", "OT_interpol")
set.seed(420)

failed_fit <- function() list(estimate = rep(NA_real_, length(cuts)),
                              variance = rep(NA_real_, length(cuts)),
                              wfit = NULL)
safe_fit <- function(...) tryCatch(dr_balance_OT(...), error = function(e) { print(e); failed_fit() })

summarize_method <- function(label, scenario_name, est, vr, ess, truths, weights) {
  failure_rate <- colMeans(is.na(est))
  
  outlier <- abs(est) > 10 | vr > 100
  outlier[is.na(outlier)] <- FALSE
  est[outlier] <- NA; vr[outlier] <- NA
  
  avg_truth <- colMeans(truths, na.rm = TRUE)
  lo <- est - 1.96 * sqrt(vr); hi <- est + 1.96 * sqrt(vr)
  truth_mat <- matrix(avg_truth, nrow(est), ncol(est), byrow = TRUE)
  
  list(
    summary = data.frame(
      scenario = scenario_name, time = cuts, method = label,
      truth = avg_truth,
      bias  = colMeans(est, na.rm = TRUE) - avg_truth,
      mc_var = apply(est, 2, var, na.rm = TRUE),
      rmse  = sqrt(colMeans((est - truths)^2, na.rm = TRUE)),
      ESS   = mean(ess, na.rm = TRUE),
      coverage = colMeans(truth_mat >= lo & truth_mat <= hi, na.rm = TRUE),
      outlier_rate = colMeans(outlier),
      failure_rate = failure_rate),
    weights = if (is.null(weights)) NULL else
      data.frame(scenario = scenario_name, method = label, w = weights / mean(weights))
  )
}

run_it <- function(scenario_name) {
  K <- length(method_labels)
  est <- vr <- array(NA_real_, c(n_sim, length(cuts), K))
  ess <- matrix(NA_real_, n_sim, K)
  truths <- matrix(NA_real_, n_sim, length(cuts))
  wts <- vector("list", K)
  
  for (i in seq_len(n_sim)) {
    print(i)
    sim_data <- generate_data(n_total = n_total, scenario = scenario_name, cuts = cuts)
    df <- sim_data$df; S <- df$S
    # Wmat1 <- as.matrix(df[df$S==1, c("x1","x2")])
    # Wmat0 <- as.matrix(df[df$S==0,,c("x1","x2")])
    # 
    # # calculate the OT weights
    # m1 <- causalOT::Measure(x = Wmat1, adapt = "none")
    # m1_w <- causalOT::Measure(x = Wmat1, adapt = "weights")
    # m0 <- causalOT::Measure(x = Wmat0, adapt = "none")
    # otp_01 <- OTProblem(m1_test, m0)
    # otp_01$setup_arguments(lambda = 4, delta = NULL, grid.length = 7L,
    #                          cost.function = NULL, p = 2, cost.online = "auto",
    #                          debias = TRUE, diameter = NULL,
    #                          ot_niter = 1000L, ot_tol = 1e-3)
    # otp_01$solve(
    #   niter = 1000L, # maximum number of iterations
    #   tol = 1e-5, # tolerance for convergence
    #   optimizer = "torch", # which optimizer to use "torch" or "frank-wolfe"
    #   torch_optim = torch::optim_lbfgs, # torch optimizer to use if required
    #   torch_scheduler = torch::lr_reduce_on_plateau, # torch scheduler to use if required
    #   torch_args = list(line_search_fn = "strong_wolfe"), # args passed to the torch functions,
    #   osqp_args = NULL, #arguments passed to the osqp solver used for "frank-wolfe" and balance functions
    #   quick.balance.function = TRUE # if balance functions are also present, should an approximate value of the hyperparameter "delta" be found first
    # )
    # 
    # wts_OT <- as.numeric(m1_w$weights)
    # 

    args <- list(S = S, X = W, Z1 = df$A[S == 1], time1 = df$time[S == 1],
                 event1 = df$status[S == 1], cuts = cuts)
    
    fits <- c(
      list(do.call(safe_fit, c(args, list(method = "optweight", tols = 0.05)))),
      lapply(ot_lambdas, function(l)
        do.call(safe_fit, c(args, list(method = "OT", tols = 0.05, lambda = l))))
    )
    
    for (k in seq_len(K)) {
      f <- fits[[k]]
      est[i, , k] <- f$estimate
      vr[i, , k]  <- f$variance
      if (!is.null(f$wfit)) {
        ess[i, k] <- f$wfit$ESS
        if (is.null(wts[[k]])) {
          w <- f$wfit$weights
          wts[[k]] <- if (length(w) == length(S)) w[S == 1] else w
        }
      }
    }
    truths[i, ] <- sim_data$truth$diff
  }
  
  res <- lapply(seq_len(K), function(k)
    summarize_method(method_labels[k], scenario_name,
                     est[, , k], vr[, , k], ess[, k], truths, wts[[k]]))
  
  list(summary = do.call(rbind, lapply(res, `[[`, "summary")),
       weights = do.call(rbind, lapply(res, `[[`, "weights")))
}


results_by_scenario <- lapply(scenarios, run_it)
final_results <- do.call(rbind, lapply(results_by_scenario, `[[`, "summary"))
weights_long  <- do.call(rbind, lapply(results_by_scenario, `[[`, "weights"))

scen_labels <- c("Baseline", "Practical Positivity Violation", "Structural Positivity Violation")
final_results$scenario <- factor(final_results$scenario, levels = scenarios, labels = scen_labels)
weights_long$scenario  <- factor(weights_long$scenario,  levels = scenarios, labels = scen_labels)


library(ggplot2)



#####################################
# DO NOT CONSIDER EDITING PLOTS YET #
#####################################
p_bias <- ggplot(final_results,
                 aes(x = time, y = bias, color = method, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ scenario) +
  scale_x_continuous(breaks = cuts) +
  labs(title = "Bias by Follow-up Time",
       x = "Follow-up Time (t)",
       y = "Bias",
       color = "Method") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

p_bias
p_rmse<- ggplot(final_results,
                aes(x = time, y = rmse, color = method, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ scenario) +
  scale_x_continuous(breaks = cuts) +
  labs(title = "Bias by Follow-up Time",
       x = "Follow-up Time (t)",
       y = "RMSE",
       color = "Method") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

p_rmse

p_var<- ggplot(final_results,
               aes(x = time, y = mc_var, color = method, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ scenario) +
  scale_x_continuous(breaks = cuts) +
  labs(title = "Bias by Follow-up Time",
       x = "Follow-up Time (t)",
       y = "Var",
       color = "Method") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

p_var


dir.create(file.path(out_dir, "Figures"), showWarnings = FALSE)
ggsave(file.path(out_dir, "Figures/bias_by_method.png"), p_bias,
       width = 12, height = 6, dpi = 300)

results_long <- final_results %>%
  select(scenario, time, bias, coverage) %>%
  pivot_longer(cols = c("bias", "coverage"),
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric = case_when(metric == "bias" ~ "Bias",
                            metric == "coverage" ~ "Coverage Probability"))

write_csv(results_long, file.path(out_dir, "Results/simulation_table.csv"))

ref_lines <- data.frame(metric = c("Bias", "Coverage Probability"),
                        hline = c(0, 0.95))

faceted_plot <- ggplot(results_long, aes(x = time, y = value, color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.2) +
  facet_wrap(~ metric, scales = "free_y") +
  geom_hline(data = ref_lines,
             aes(yintercept = hline),
             linetype = "dashed",
             color = "black") +
  labs(title = "Estimator Performance by Scenario",
       x = "Follow-up Time (t)",
       y = "Value",
       color = "Scenario") +
  scale_x_continuous(breaks = cuts) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 21, face = "bold"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 19),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 19, face = "bold")) +
  guides(color = guide_legend(nrow = 2, override.aes = list(linewidth = 1.6, size = 4)))

ggsave(file.path(out_dir, "Figures/Figure1.png"), faceted_plot, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "Figures/Figure1.tiff"), faceted_plot, width = 14, height = 8, dpi = 300,  compression = "lzw")

###=============================================================================
### Weight-density panel (Reviewer 2 major comment 4)
###=============================================================================

weights_long$scenario <- factor(weights_long$scenario,
                                levels = scenarios,
                                labels = c("Baseline",
                                           "Incorrect Selection Model",
                                           "Incorrect Outcome Model",
                                           "Exchangeability Violation",
                                           "Nonlinear Misspecification",
                                           "Practical Positivity Violation",
                                           "Structural Positivity Violation",
                                           "Outcome Extrapolation"))

weight_plot <- ggplot(subset(weights_long, w > 0),
                      aes(x = log10(w))) +
  geom_density(fill = "grey40", 
               color = "grey20", 
               alpha = 0.5) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "black") +
  facet_wrap(~ scenario, 
             scales = "free_y",
             ncol = 4) +
  labs(x = expression(log[10] * "(trial weight)"),
       y = "Density") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(size = 10, face = "bold"))

ggsave(file.path(out_dir, "Figures/simulation_weight_distribution.png"), 
       weight_plot, width = 12, height = 5, dpi = 300)



################
# Copilot Plot #
################

library(ggplot2)
library(patchwork)

# Bias plot
p_bias <- ggplot(
  final_results_no_outlier_filter,
  aes(x = time, y = bias, color = method, group = method)
) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~scenario) +
  labs(
    title = "Bias by Scenario",
    x = "Time",
    y = "Bias",
    color = "Method"
  ) +
  theme_minimal()

# RMSE plot
p_rmse <- ggplot(
  final_results_no_outlier_filter,
  aes(x = time, y = rmse, color = method, group = method)
) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~scenario) +
  labs(
    title = "RMSE by Scenario",
    x = "Time",
    y = "RMSE",
    color = "Method"
  ) +
  theme_minimal()

# Display together
p_bias / p_rmse
