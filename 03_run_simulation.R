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
source(file.path(out_dir, "02_estimators.R"))

scenarios <- c("baseline",
               #"omitted_selection",
               #"omitted_outcome",
               #"omitted_both",
               #"both_nonlinear",
               "positivity_practical",
               "positivity_structural") #,
               #"outcome_extrapolation")

n_sim <- 500
n_total <- 1000
cuts <- c(1, 2, 3, 4, 5, 6)

set.seed(420)

summarize_data <- function(scenario_name, method, cuts = cuts, fit, truths = truths, lambda = NULL) {
  estimates[i, ] <- fit$estimate    
  variances[i, ] <- fit$variance
  if (method == "EB") {
    if (!is.null(fit$wfit)) {
      weights_list[[i]] <- fit$wfit$weights[S == 1]
    }}
  
  if (method != "EB" & !is.null(fit_ot$wfit)) {
    weights_list[[i]] <- fit$wfit$weights
  }
  failure_rate <- colMeans(is.na(estimates))
  
  outlier <- abs(estimates) > 10 | variances > 100
  estimates[outlier] <- NA
  variances[outlier] <- NA
  
  avg_truth <- colMeans(truths, na.rm = TRUE)
  bias <- colMeans(estimates, na.rm = TRUE) - avg_truth
  var <- colMeans(estimates, na.rm = TRUE) - avg_truth
  rmse <- sqrt(colMeans((estimates - truths)^2, na.rm = TRUE))
  
  lo <- estimates - 1.96 * sqrt(variances)
  hi <- estimates + 1.96 * sqrt(variances)
  truth_mat <- matrix(rep(avg_truth, n_sim), byrow = TRUE, nrow = n_sim)
  coverage <- colMeans(truth_mat >= lo & truth_mat <= hi, na.rm = TRUE)
  outlier_rate <- colMeans(outlier)
  
  rep_idx <- which(!vapply(weights_list, is.null, logical(1)))[1]
  
  weight_df <- 
    if (length(rep_idx) && !is.na(rep_idx)) {
      data.frame(
        scenario = scenario_name,
        method = method,
        w = weights_list[[rep_idx]]
      )
    } else {
      NULL
    }
  
  return(list(summary = data.frame(
    scenario = scenario_name,
    time = cuts,
    method = method,
    truth = avg_truth,
    bias = bias,
    rmse = rmse,
    coverage = coverage,
    outlier_rate = outlier_rate,
    failure_rate = failure_rate
    #outliers = outlier_eb
  ), weights = weights_df))
  
  
}


run_it <- function(scenario_name) {
  # MEASURE RUNTIME
  # approx OT package 
  # clue solve::LASP
  estimates_eb <- matrix(NA_real_, n_sim, length(cuts))
  variances_eb <- matrix(NA_real_, n_sim, length(cuts))
  truths <- matrix(NA_real_, n_sim, length(cuts))
  weights_list_eb <- vector("list", n_sim)
  
  estimates_ot <- matrix(NA_real_, n_sim, length(cuts))
  variances_ot <- matrix(NA_real_, n_sim, length(cuts))
  weights_list_ot <- vector("list", n_sim)
  
  for (i in seq_len(n_sim)) {
    print(i)
    sim_data <- generate_data(n_total = n_total, scenario = scenario_name,
                              cuts = cuts)
    df <- sim_data$df
    S <- df$S
    W <- data.frame(x1 = df$x1, x2 = df$x2)
    Z1 <- df$A[S == 1]
    t1 <- df$time[S == 1]
    e1 <- df$status[S == 1]
    
    fit_eb <- tryCatch(
      dr_balance(
        S = S,
        X = W,
        Z1 = Z1,
        time1 = t1,
        event1 = e1,
        cuts = cuts,
        method = "optweight",
        tols = 0.05
      ),
      error = function(e) {
        print(e)
        list(
          estimate = rep(NA, length(cuts)),
          variance = rep(NA, length(cuts)),
          wfit = NULL
        )
      }
    )
    fit_ot_0 <- tryCatch(
      dr_balance(S = S, X = W, Z1 = Z1, time1 = t1, event1 = e1,
                 cuts = cuts, method = "OT", tols = 0.05, lambda = 0),
      error = function(e) {
        print(e)
        list(estimate = rep(NA, length(cuts)),
                               variance = rep(NA, length(cuts)),
                               wfit = NULL)
      }
    )
    
    fit_ot_1 <- tryCatch(
      dr_balance(S = S, X = W, Z1 = Z1, time1 = t1, event1 = e1,
                 cuts = cuts, method = "OT", tols = 0.05, lambda = 1),
      error = function(e) {
        print(e)
        list(estimate = rep(NA, length(cuts)),
             variance = rep(NA, length(cuts)),
             wfit = NULL)
      }
    )
    
    fit_ot_10 <- tryCatch(
      dr_balance(S = S, X = W, Z1 = Z1, time1 = t1, event1 = e1,
                 cuts = cuts, method = "OT", tols = 0.05, lambda = 10),
      error = function(e) {
        print(e)
        list(estimate = rep(NA, length(cuts)),
             variance = rep(NA, length(cuts)),
             wfit = NULL)
      }
    )
    
    fit_ot_100 <- tryCatch(
      dr_balance(S = S, X = W, Z1 = Z1, time1 = t1, event1 = e1,
                 cuts = cuts, method = "OT", tols = 0.05, lambda = 100),
      error = function(e) {
        print(e)
        list(estimate = rep(NA, length(cuts)),
             variance = rep(NA, length(cuts)),
             wfit = NULL)
      }
    )
    
    fit_ot_1000 <- tryCatch(
      dr_balance(S = S, X = W, Z1 = Z1, time1 = t1, event1 = e1,
                 cuts = cuts, method = "OT", tols = 0.05, lambda = 1000),
      error = function(e) {
        print(e)
        list(estimate = rep(NA, length(cuts)),
             variance = rep(NA, length(cuts)),
             wfit = NULL)
      }
    )
    
    truths[i, ] <- sim_data$truth$diff
    
    
  
  

  # Flag extreme outliers (numerical blow-ups) as NA
  

  
 
  
  
  
  bias_ot <- colMeans(estimates_ot, na.rm = TRUE) - avg_truth
  rmse_ot <- sqrt(colMeans((estimates_ot - truths)^2, na.rm = TRUE))
  
  lo_ot <- estimates_ot - 1.96 * sqrt(variances_ot)
  hi_ot <- estimates_ot + 1.96 * sqrt(variances_ot)
  coverage_ot <- colMeans(truth_mat >= lo_ot & truth_mat <= hi_ot, na.rm = TRUE)
  
  outlier_rate_ot <- colMeans(outlier_ot)
  

  summary_df <- rbind(
    # (scenario_name, method, cuts = cuts, fit, truths = truths) 
    summarize_data(scenario_name = scenario_name,
                   method = "EB",
                   fit = fit_eb
                   ),
    
    summarize_data(scenario_name = scenario_name,
                   method = "OT",
                   fit = fit_ot_0, lambda = 0
    ),
    summarize_data(scenario_name = scenario_name,
                   method = "OT",
                   fit = fit_ot_1, lambda = 1
    ),
    summarize_data(scenario_name = scenario_name,
                   method = "OT",
                   fit = fit_ot_10, lambda = 10
    ),
    summarize_data(scenario_name = scenario_name,
                   method = "OT",
                   fit = fit_ot_100, lambda = 100
    ),
    summarize_data(scenario_name = scenario_name,
                   method = "OT",
                   fit = fit_ot_1000, lambda = 1000
    )
  )
  

  list(summary = summary_df, weights = weight_df)
  
}

results_by_scenario <- lapply(scenarios, run_it)
final_results_no_outlier_filter <- do.call(rbind, lapply(results_by_scenario, `[[`, "summary"))
weights_long <- do.call(rbind, lapply(results_by_scenario, `[[`, "weights"))

final_results$scenario <- factor(final_results$scenario,
                                 levels = scenarios,
                                 labels = c("Baseline",
                                            "Practical Positivity Violation",
                                            "Structural Positivity Violation"))

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