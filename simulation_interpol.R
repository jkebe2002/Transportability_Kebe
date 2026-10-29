# interpolation simulation
#
# At each run (replicate):
#   1. sample trial (S = 1) and target (S = 0) data from the DGP
#   2. dr_balance() (in simulation_interpol_effects.R) calculates the OT weights,
#      then for every t in t_grid builds
#        - OT_t        : t * 1 + (1 - t) * OT weights
#        - OT_interpol : weights from t * OT(source_w, source) + (1 - t) * OT(source_w, target)
#   3. for each method x t we record: estimate, variance, ESS, W2 distance to the
#      target / source / OT-weighted source, Tipton index, min/max/median weight

###=============================================================================
### Setup
###=============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

out_dir <- file.path("~", "Documents/GitHub", "Transportability_Kebe")
out_dir <- normalizePath(out_dir, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), showWarnings = FALSE)
dir.create(file.path(out_dir, "Results"), showWarnings = FALSE)

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
t_grid <- seq(0, 1, by = 0.25)

# Set TRUE for a fast smoke test before launching the full run
quick_test <- FALSE
if (quick_test) { n_sim <- 5; scenarios <- scenarios[1]; t_grid <- c(0, 0.5, 1) }

method_labels   <- fit_labels(t_grid)                     # e.g. "OT_t_0", "OT_interpol_0", ...
method_family   <- rep(c("OT_t", "OT_interpol"), times = length(t_grid))
method_t        <- rep(t_grid, each = 2)
K <- length(method_labels)

set.seed(420)

# One failed result per method if the whole replicate errors
safe_fit <- function(...) {
  tryCatch(dr_balance(...), error = function(e) {
    print(e)
    setNames(rep(list(failed_fit(length(cuts))), K), method_labels)
  })
}

###=============================================================================
### Summaries
###=============================================================================

summarize_method <- function(label, family, tt, scenario_name,
                             est, vr, ess, dg, truths, weights,
                             filter_outliers = TRUE) {
  failure_rate <- colMeans(is.na(est))
  
  outlier <- abs(est) > 10 | vr > 100
  outlier[is.na(outlier)] <- FALSE
  if (filter_outliers) { est[outlier] <- NA; vr[outlier] <- NA }
  
  avg_truth <- colMeans(truths, na.rm = TRUE)
  lo <- est - 1.96 * sqrt(vr); hi <- est + 1.96 * sqrt(vr)
  truth_mat <- matrix(avg_truth, nrow(est), ncol(est), byrow = TRUE)
  
  # replicate-averaged diagnostics (one value per method; recycled across times)
  dg_mean <- as.data.frame(as.list(colMeans(dg, na.rm = TRUE)))
  names(dg_mean) <- paste0("mean_", diag_names)
  
  list(
    summary = data.frame(
      scenario = scenario_name, time = cuts, method = label,
      family = family, t = tt,
      truth = avg_truth,
      bias  = colMeans(est, na.rm = TRUE) - avg_truth,
      mc_var = apply(est, 2, var, na.rm = TRUE),
      rmse  = sqrt(colMeans((est - truths)^2, na.rm = TRUE)),
      ESS   = mean(ess, na.rm = TRUE),
      coverage = colMeans(truth_mat >= lo & truth_mat <= hi, na.rm = TRUE),
      outlier_rate = colMeans(outlier),
      failure_rate = failure_rate,
      dg_mean),
    weights = if (is.null(weights)) NULL else
      data.frame(scenario = scenario_name, method = label, family = family, t = tt,
                 w = weights / mean(weights))
  )
}

###=============================================================================
### Run the simulation
###=============================================================================

run_it <- function(scenario_name) {
  est <- vr <- array(NA_real_, c(n_sim, length(cuts), K))
  ess    <- matrix(NA_real_, n_sim, K)
  dg     <- array(NA_real_, c(n_sim, K, length(diag_names)),
                  dimnames = list(NULL, NULL, diag_names))
  truths <- matrix(NA_real_, n_sim, length(cuts))
  wts <- vector("list", K)
  
  for (i in seq_len(n_sim)) {
    cat(scenario_name, "- replicate", i, "of", n_sim, "\n")
    sim_data <- generate_data(n_total = n_total, scenario = scenario_name, cuts = cuts)
    df <- sim_data$df; S <- df$S
    W <- df[, c("x1", "x2")]
    
    fits <- safe_fit(S = S, X = W, Z1 = df$A[S == 1],
                     time1 = df$time[S == 1], event1 = df$status[S == 1],
                     cuts = cuts, t_grid = t_grid)
    
    for (k in seq_len(K)) {
      f <- fits[[k]]
      est[i, , k] <- f$estimate
      vr[i, , k]  <- f$variance
      dg[i, k, ]  <- f$diag[diag_names]
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
  
  summarize_all <- function(filter_outliers) {
    lapply(seq_len(K), function(k)
      summarize_method(method_labels[k], method_family[k], method_t[k], scenario_name,
                       est[, , k], vr[, , k], ess[, k], dg[, k, ], truths, wts[[k]],
                       filter_outliers = filter_outliers))
  }
  res         <- summarize_all(TRUE)
  res_nofilt  <- summarize_all(FALSE)
  
  list(summary = do.call(rbind, lapply(res, `[[`, "summary")),
       summary_no_outlier_filter = do.call(rbind, lapply(res_nofilt, `[[`, "summary")),
       weights = do.call(rbind, lapply(res, `[[`, "weights")))
}

results_by_scenario <- lapply(scenarios, run_it)
final_results <- do.call(rbind, lapply(results_by_scenario, `[[`, "summary"))
final_results_no_outlier_filter <-
  do.call(rbind, lapply(results_by_scenario, `[[`, "summary_no_outlier_filter"))
weights_long  <- do.call(rbind, lapply(results_by_scenario, `[[`, "weights"))

# One set of scenario labels, applied once to every results table
scen_labels_all <- c(baseline              = "Baseline",
                     omitted_selection     = "Incorrect Selection Model",
                     omitted_outcome       = "Incorrect Outcome Model",
                     omitted_both          = "Exchangeability Violation",
                     both_nonlinear        = "Nonlinear Misspecification",
                     positivity_practical  = "Practical Positivity Violation",
                     positivity_structural = "Structural Positivity Violation",
                     outcome_extrapolation = "Outcome Extrapolation")
relabel_scenario <- function(d) {
  d$scenario <- factor(d$scenario, levels = scenarios,
                       labels = unname(scen_labels_all[scenarios]))
  d
}
final_results                   <- relabel_scenario(final_results)
final_results_no_outlier_filter <- relabel_scenario(final_results_no_outlier_filter)
weights_long                    <- relabel_scenario(weights_long)

# Full tables (all methods, all diagnostics) so nothing is lost to the plot reshaping below
write_csv(final_results, file.path(out_dir, "Results/simulation_full_summary.csv"))
write_csv(final_results_no_outlier_filter,
          file.path(out_dir, "Results/simulation_full_summary_no_outlier_filter.csv"))


#####################################
# DO NOT CONSIDER EDITING PLOTS YET #
#####################################
# (Plot code below is unchanged apart from removing the duplicate scenario relabel
#  and the stale scen_labels. NOTE: with the t-grid, `method` now has 2 x length(t_grid)
#  levels, and results_long / weights_long pool across methods, so these plots will
#  need rework once you get to them.)

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