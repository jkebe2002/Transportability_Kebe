# Cimulate trial data
# get summary statistics, outcome distribution for trial data
# simulate RWD
# get covariance structure from RWD
# reconstruct it, make reflection go halfway to reflection
# Use HIV data
#source("~/Documents/GitHub/Transportability_Kebe/HIVdata.R")

library(WeightIt)
library(optweight)
library(causalOT)
library(expm)
library(Matrix)
library('actuar') # For invgamma
library('copula') # For copula functions
library('dplyr')  # Piping function




vals <- c(
  18.84666876,   4.24102954,   0.01877518,  -0.123479512, -0.216549515,
  4.24102954,  76.11203450,  -2.43149491,   0.111415303,  0.067904772,
  0.01877518,  -2.43149491,  14.90460416,  -0.041154356,  0.674449853,
  -0.123479512,  0.111415303, -0.041154356,   0.075378151,  0.002052095,
  -0.216549515,  0.067904772,  0.674449853,   0.002052095,  1.534629076
)

# Covariance structure inspired by HIV dataset
#Here we make C which is the RWD covariance matrix for the 4 covariates and the outcome
C <- matrix(vals, nrow = 5, ncol = 5, byrow = TRUE,
              dimnames = list(paste0("X", 1:5), paste0("X", 1:5)))

R5 <- matrix(1, 5, 5, dimnames = list(vars5, vars5))

# change correlations for CT covariance
R5[upper.tri(R5)] <- c(
  -0.10,                 # BMI ~ age
  0.25,  0.20,          # DIAB ~ BMI, age
  -0.05, -0.15,  0.10,   # hard_drugs ~ BMI, age, DIAB
  0.10, -0.05,  0.15,  0.20   # diff_lvload ~ BMI, age, DIAB, hard_drugs
)
R5[lower.tri(R5)] <- t(R5)[lower.tri(R5)]

# SDs similar in scale to original
sds <- c(BMI = 4.34, age_0 = 8.74, DIAB_0 = 3.87,
         hard_drugs_01 = 0.268, diff_lvload = 1.23)

# Trial true covariance struct.
C_trial <- diag(sds) %*% R5 %*% diag(sds)






# code from Kevin EB paper
esteq_pate <- function(S, X, Y, Z, weights, base_weights, tau) {
  
  eq1 <- S * Z * weights * X -
    (1 - S) * base_weights * X
  
  eq2 <- S * (1 - Z) * weights * X -
    (1 - S) * base_weights * X
  
  eq3 <- S * weights *
    (Z * (Y - tau) - (1 - Z) * Y)
  
  eq <- c(eq1, eq2, eq3)
  
  return(eq)
  
}

# correlation matrix from covariance matrix
cor_from_cov <- function(Sigma) {
  D_inv <- diag(1 / sqrt(diag(Sigma)))
  R <- D_inv %*% Sigma %*% D_inv
  return(R)
}

#############################################################
#Next are wang-style data generating functions using copulas#
#############################################################

# Function to extract distribution and parameter arguments
extract_margin_distr_params = function(...) {
  args_list = list(...)
  n = length(args_list)
  
  margins_vector = character(length = n)
  params_list = vector("list", length = n)
  
  for (i in 1:n) {
    margins_vector[i] = args_list[[i]][[1]]
    params_list[[i]] = args_list[[i]][[2]]
  }
  
  return(list(margins_vector, params_list))
}
# Function to create a multivariate distribution
create_multivariate_dist = function(input, ...) {
  Margins = list(...)
  if (ncol(input) != length(Margins)) {
    stop('input cor matrix should have same number of args as ..., corresponding to marginal distributions and params')
  }
  
  distr_param = do.call(extract_margin_distr_params, Margins)
  
  marg = distr_param[[1]]
  params = distr_param[[2]]
  
  normal = normalCopula(param = P2p(input), dim = ncol(input), dispstr = "un")
  multivariate_dist = mvdc(copula = normal, margins = marg, paramMargins = params)
  return(multivariate_dist)
}

# Function to generate samples from a multivariate distribution
generate_sample = function(n, multivariate_distr, binomial_vectors = NULL, ...) {
  samples = as.data.frame(rMvdc(n, multivariate_distr))
  names(samples) = paste0('var', 1:ncol(samples))
  props = list(...)
  props = props[[1]]
  if (!is.null(binomial_vectors)) {
    if (length(binomial_vectors) != length(props)) {
      stop('need mean/sd for all binomial_vectors, ... additional arguments should be the same length as binomial_vectors')
    }
    for (i in 1:length(binomial_vectors)) {
      elements = props[[i]]
      samples[, binomial_vectors[i]] = ifelse(pnorm(samples[, binomial_vectors[i]], elements[[2]][['mean']], elements[[2]][['sd']]) <= elements[[2]][['mean']], 1, 0)
    }
  }
  return(samples)
}
#########################################
# geodesic transporting V to U at time t#
#########################################
geodesic_step <- function(V,U,t){
  #identity 
  E <-diag(nrow(V))
  Uhalf <- sqrtm(U)
  #map from V to U
  T <- Uhalf %*% sqrtm(solve((Uhalf %*% V %*% Uhalf))) %*% Uhalf
  
  Wt <- ((1-t) * E + t * T) %*% V %*% ((1-t) * E + t * T)
  return(Wt)
}

R_rwd <- cor_from_cov(C)

###############################################
# outcome model with noise taken from the data#
###############################################
outcome_model <- function(dat, trt) {
  -2.8 * trt - 0.1 * dat$X1 + 0.9 * dat$X2 + 0.5 * dat$X3 + 0.4 * dat$X4 + dat$eps
}

###############################################
# functions to generate data                  #
###############################################

trial_data_gen <- function(means_trial, C_trial, n_1, binomial_columns){
  dist_trial <- create_multivariate_dist(
    cor_from_cov(C_trial[1:4, 1:4]),
    list("norm", list(mean = means_trial[1], sd = sqrt(C_trial[1,1]))),
    list("norm", list(mean = means_trial[2], sd = sqrt(C_trial[2,2]))),
    list("norm", list(mean = means_trial[3], sd = sqrt(C_trial[3,3]))),
    list("norm", list(mean = means_trial[4], sd = sqrt(C_trial[4,4])))
  )
  #generate synthetic sample
  synth_samp_trial <- generate_sample(n_1, dist_trial, binomial_vectors = binomial_columns, #X3, X4 binary
                                      props = list(
                                        list("norm", list(mean = means_trial[3], sd = sqrt(C_trial[3,3]))),
                                        list("norm", list(mean = means_trial[4], sd = sqrt(C_trial[4,4])))
                                      ))
  
  # treatment assignment on the trial sample
  
  names(synth_samp_trial)[1:4] <- c("X1","X2","X3","X4")
  
  # treatment is random
  synth_samp_trial$trt <- rbinom(nrow(synth_samp_trial), 1, 0.5)
  
  sigma <- 1   # residual CD; adjust as needed
  
  # one noise draw per unit, shared across potential outcomes
  # normally distributed error for outcome model
  synth_samp_trial$eps <- rnorm(nrow(synth_samp_trial), 0, sigma)
  
  # potential outcomes
  synth_samp_trial$out_trt  <- outcome_model(synth_samp_trial, trt = 1)
  synth_samp_trial$out_ctrl <- outcome_model(synth_samp_trial, trt = 0)
  
  # observed outcome in the trial
  synth_samp_trial$out_real <- ifelse(synth_samp_trial$trt == 1,
                                      synth_samp_trial$out_trt,
                                      synth_samp_trial$out_ctrl)
  return(synth_samp_trial)
}


rwd_data_gen <- function(means_rws, C_rwd, n_0, binomial_columns) {
  
  dist_rwd <- create_multivariate_dist(
    cor_from_cov(C_rwd[1:4, 1:4]),
    list("norm", list(mean = means_rwd[1], sd = sqrt(C_rwd[1,1]))),
    list("norm", list(mean = means_rwd[2], sd = sqrt(C_rwd[2,2]))),
    list("norm", list(mean = means_rwd[3], sd = sqrt(C_rwd[3,3]))),
    list("norm", list(mean = means_rwd[4], sd = sqrt(C_rwd[4,4])))
  )
  
  #generate synthetic sample
  synth_samp_rwd <- generate_sample(n_0, dist_rwd, binomial_vectors = binomial_columns,
                                    props = list(
                                      list("norm", list(mean = means_rwd[3], sd = sqrt(C_rwd[3,3]))),
                                      list("norm", list(mean = means_rwd[4], sd = sqrt(C_rwd[4,4])))
                                    ))
  names(synth_samp_rwd)[1:4] <- c("X1","X2","X3","X4")
  
  sigma <- 1
  synth_samp_rwd$eps   <- rnorm(nrow(synth_samp_rwd),   0, sigma)
  
  #potential outcomes for RWD
  synth_samp_rwd$out_trt  <- outcome_model(synth_samp_rwd, trt = 1)
  synth_samp_rwd$out_ctrl <- outcome_model(synth_samp_rwd, trt = 0)
  
  return(synth_samp_rwd)
}




############################
# Reflect covariance matrix#
############################


flip_matrix <- function(C, method = "Higham") {
  
  R <- -C
  diag(R) <- diag(C)
  
  if (method == "Frobenius") {
    
    ev <- eigen(R, symmetric = TRUE)
    
    R_psd <- ev$vectors %*%
      diag(pmax(ev$values, 0)) %*%
      t(ev$vectors)
    
    R_psd <- (R_psd + t(R_psd))/2
    
    dimnames(R_psd) <- dimnames(C)
    
    return(R_psd)
  }
  
  # Higham keeps variances from original covariance matrix
  if (method == "Higham") {
    
    np <- nearPD(
      R,
      corr = FALSE,
      keepDiag = TRUE,
      base.matrix = TRUE
    )
    
    R_psd_diag <- np$mat
    
    dimnames(R_psd_diag) <- dimnames(C)
    
    return(R_psd_diag)
  }
  
  stop("Unknown method")
}

# generate synthetic trial data (OT)

synth_data_gen <- function(C = C_rwd,
                           means_trial, 
                           sd_trial, 
                           means_outcome, # one value for trt, one for control
                           sd_outcome, #one value for trt, one value for ctrl
                           binomial_columns = c(3,4), 
                           method,
                           t = .5) {
  C_refl <- flip_matrix(C, method)
  C_synth <- geodesic_step(C, C_refl, t)
  
  dist_trial_pseudo <- create_multivariate_dist(
    cor_from_cov(C_synth[1:4,1:4]),
    # now we are using the covariance matrix at time t between C and its reflection
    list("norm", list(mean = means_trial[1], sd = sd_trial[1])),
    list("norm", list(mean = means_trial[2], sd = sd_trial[2])),
    list("norm", list(mean = means_trial[3], sd = sd_trial[3])),
    list("norm", list(mean = means_trial[4], sd = sd_trial[4]))
  )
  #generate synthetic sample
  pseudo_samp_trial <- generate_sample(n_1, dist_trial_pseudo, binomial_vectors = binomial_columns,
                                      props = list(
                                        list("norm", list(mean = means_trial[3], sd = sd_trial[3])),
                                        list("norm", list(mean = means_trial[4], sd = sd_trial[4]))
                                      ))
  names(pseudo_samp_trial)[1:4] <- c("X1","X2","X3","X4")
  
  # treatment assignment on the trial sample
  pseudo_samp_trial$trt <- rbinom(nrow(pseudo_samp_trial), 1, 0.5)
  
  idx0 <- pseudo_samp_trial$trt == 0
  idx1 <- pseudo_samp_trial$trt == 1
  pseudo_samp_trial$out_real <- NA_real_
  
  #the pseudo data's outcomes come from normal distributions around the observed summary stats
  pseudo_samp_trial$out_real[idx0] <- rnorm(
    sum(idx0),
    mean = means_outcome[1],
    sd = sd_outcome[1]
  )
  
  pseudo_samp_trial$out_real[idx1] <- rnorm(
    sum(idx1),
    mean = means_outcome[2],
    sd = sd_outcome[2]
  )
  
  return(pseudo_samp_trial)
}


#Kevin Code:
# Estimation with individual-level data. Y and Z are arbitrary for the target sample, 
# but must be included. Set them to 0.
estimate_pate <- function(S, X, Y, Z,
                          weights,
                          base_weights = NULL, ...) {

  n_0 <- sum(1 - S)
  n_1 <- sum(S)
  n <- n_1 + n_0
  
  if (is.null(base_weights))
    base_weights <- rep(n_1/n_0, times = length(S))
  
  if (length(base_weights) != length(S))
    stop("base_weights must have the same length as S")
  
  tau <- sum(S*(weights*(2*Z - 1)*Y)/sum(S*Z*weights))
  
    m <- ncol(X)
    U <- matrix(0, ncol = 2*m, nrow = 2*m)
    v <- rep(0, times = 2*m + 1)
    meat <- matrix(0, ncol = 2*m + 1, nrow = 2*m + 1)
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - S[i]*Z[i] * weights[i] * (X[i,] %*% t(X[i,]))
      U[(m + 1):(2*m),(m + 1):(2*m)] <- U[(m + 1):(2*m),(m + 1):(2*m)] - S[i]*(1 - Z[i]) * weights[i] * X[i,] %*% t(X[i,])
      
      v[1:m] <- v[1:m] - S[i]*Z[i] * weights[i] * (Y[i] - tau) * X[i,]
      v[(m + 1):(2*m)] <- v[(m + 1):(2*m)] + S[i]*(1 - Z[i]) * weights[i] * Y[i] * X[i,]
      v[2*m + 1] <- v[2*m + 1] - S[i]*weights[i]*Z[i]
      s <- esteq_pate(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], weights = weights[i], 
                      base_weights = base_weights[i], tau = tau)
      meat <- meat + s %*% t(s)
      
    }
    
    invbread <- matrix(0, nrow = 2*m + 1, ncol = 2*m + 1)
    invbread[1:(2*m),1:(2*m)] <- U
    invbread[2*m + 1, ] <- v
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandwich <- NA
      variance <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[2*m + 1, 2*m + 1]
      
    }
    
  
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}

# Apply Kevin's code to transport pseudo trial to target at one iteration
estimate_pate_synth <- function(rwd_dat, synth_dat) {
  
  # Keep only variables needed downstream
  rwd_dat <- rwd_dat[, c("X1","X2","X3","X4")]
  synth_dat <- synth_dat[, c("X1","X2","X3","X4","trt","out_real")]
  
  # Kevin's original convention:
  # target sample has arbitrary Y and Z values
  rwd_dat$out_real <- 0
  rwd_dat$trt <- 0
  
  data_combined <- rbind(rwd_dat, synth_dat)
  
  # S = 1 for trial/pseudo-trial
  data_combined$S <- c(
    rep(0, nrow(rwd_dat)),
    rep(1, nrow(synth_dat))
  )
  
  S <- data_combined$S
  n <- nrow(data_combined)
  
  Xmat <- as.matrix(
    data_combined[, c("X1","X2","X3","X4")]
  )
  
  #################################################
  # Calibration / optweight
  #################################################
  
  opt_wt <- weightit(
    S ~ X1 + X2 + X3 + X4,
    data = data_combined,
    focal = 0,
    method = "optweight",
    estimand = "ATT",
    tols = .01
  )
  
  weights_cal <- rep(1, n)
  
  weights_cal[S == 1] <-
    opt_wt$weights[S == 1]
  
  #################################################
  # CausalOT
  #################################################
  
  weights_ot_unscaled <- calc_weight(
    x = Xmat,
    z = S,
    method = "COT",
    estimand = "ATC",
    options = list(lambda = 1)
  )
  
  weights_ot <- rep(1, n)
  
  weights_ot[S == 1] <-
    as.numeric(weights_ot_unscaled@w1) * sum(S == 1)
  
  stopifnot(length(weights_cal) == n)
  stopifnot(length(weights_ot) == n)
  
  #################################################
  # PATE estimation
  #################################################
  
  pate_cal <- estimate_pate(
    S = S,
    X = Xmat,
    Y = data_combined$out_real,
    Z = data_combined$trt,
    weights = weights_cal
  )
  
  # Note that here I am naively plugging the OT weights into the PATE estimation from Kevin's framework...
  pate_ot <- estimate_pate(
    S = S,
    X = Xmat,
    Y = data_combined$out_real,
    Z = data_combined$trt,
    weights = weights_ot
  )
  
  list(
    cal = pate_cal,
    ot  = pate_ot
  )
}



# what do you do for each simulation:
# generate RWD data
# generate trial data 
# generate pseudo data 
# transport RWD to pseudo data 
# transport RWD to trial data
# compare rwd outcomes and both transports

binomial_columns <- c(3,4)

#different covariance structures and mean structures for trial, RWD
means_trial <- c(22, 40, .5, .2)
means_rwd <- c(26, 50, .2, .6)
n_iter <- 500


results <- data.frame(
  iter = 1:n_iter,
  true_pate = NA,
  
  est_frob_cal = NA,
  est_frob_ot  = NA,
  
  est_hig_cal  = NA,
  est_hig_ot   = NA,
  
  se_frob_cal = NA,
  se_frob_ot  = NA,
  
  se_hig_cal = NA,
  se_hig_ot  = NA,
  
  cover_frob_cal = NA,
  cover_frob_ot  = NA,
  
  cover_hig_cal = NA,
  cover_hig_ot  = NA
)


###############
# Simulation. #
###############

for (i in 1:n_iter)
{
  n_1 = 500
  n_0 = 500
  n=n_1+n_0
  
  # generate samples using the 'true' means and covariances
  trial_dat <- trial_data_gen(means_trial, C_trial, n_1, binomial_columns)
  rwd_dat <- rwd_data_gen(means_rwd, C, n_0, binomial_columns)
  
  #sample RWD covariance for iteration
  C_sim <- cov(rwd_dat[, c("X1","X2","X3","X4")])
  
  # trial sample summary statistics
  means_trial_sim <- c(
    mean(trial_dat$X1),
    mean(trial_dat$X2),
    mean(trial_dat$X3),
    mean(trial_dat$X4)
  )
  sd_trial_sim <- c(
    sd(trial_dat$X1),
    sd(trial_dat$X2),
    sd(trial_dat$X3),
    sd(trial_dat$X4)
  )
  
  means_outcome <- c(
    mean(trial_dat$out_real[trial_dat$trt == 0]),
    mean(trial_dat$out_real[trial_dat$trt == 1])
  )
  sd_outcome <- c(
    sd(trial_dat$out_real[trial_dat$trt == 0]),
    sd(trial_dat$out_real[trial_dat$trt == 1])
  )
  
  # create pseudo trial data using both matrix flipping methods
  synth_dat_frob <- synth_data_gen(C = C_sim,
                                          means_trial = means_trial_sim, 
                                          sd_trial = sd_trial_sim, 
                                          means_outcome, # one value for trt, one for control
                                          sd_outcome, #one value probably makes most sense here 
                                          binomial_columns = c(3,4), 
                                          method = "Frobenius",
                                          t = .5) 
  
  synth_dat_hig <- synth_data_gen(C = C_sim,
                              means_trial = means_trial_sim, 
                              sd_trial = sd_trial_sim, 
                              means_outcome, # one value for trt, one for control
                              sd_outcome, #one value probably makes most sense here 
                              binomial_columns = c(3,4), 
                              method = "Higham",
                              t = .5) 
  
  # craft calibration weights
  # craft causalOT weights
  #run estimate PATE on both
  #save bias, RMSE, coverage
  
  true_pate <- mean(
    rwd_dat$out_trt - rwd_dat$out_ctrl
  )
  
  pate_frob <- estimate_pate_synth(rwd_dat, synth_dat_frob)
  pate_hig  <- estimate_pate_synth(rwd_dat, synth_dat_hig)
  
  # Calibration
  est_frob_cal <- pate_frob[[1]]$estimate
  se_frob_cal  <- sqrt(pate_frob[[1]]$variance)
  
  est_hig_cal <- pate_hig[[1]]$estimate
  se_hig_cal  <- sqrt(pate_hig[[1]]$variance)
  
  # OT
  est_frob_ot <- pate_frob[[2]]$estimate
  se_frob_ot  <- sqrt(pate_frob[[2]]$variance)
  
  est_hig_ot <- pate_hig[[2]]$estimate
  se_hig_ot  <- sqrt(pate_hig[[2]]$variance)
  
  results$true_pate[i] <- true_pate
  
  results$est_frob_cal[i] <- est_frob_cal
  results$est_frob_ot[i]  <- est_frob_ot
  
  results$est_hig_cal[i] <- est_hig_cal
  results$est_hig_ot[i]  <- est_hig_ot
  
  results$se_frob_cal[i] <- se_frob_cal
  results$se_frob_ot[i]  <- se_frob_ot
  
  results$se_hig_cal[i] <- se_hig_cal
  results$se_hig_ot[i]  <- se_hig_ot
  
  results$cover_frob_cal[i] <- ifelse(
    is.na(se_frob_cal),
    NA,
    as.numeric(
      results$true_pate[i] >= est_frob_cal - 1.96 * se_frob_cal &
        results$true_pate[i] <= est_frob_cal + 1.96 * se_frob_cal
    )
  )
  results$cover_frob_ot[i] <- ifelse(
    is.na(se_frob_ot),
    NA,
    as.numeric(
      results$true_pate[i] >= est_frob_ot - 1.96 * se_frob_ot &
        results$true_pate[i] <= est_frob_ot + 1.96 * se_frob_ot
    )
  )
  results$cover_hig_cal[i] <- ifelse(
    is.na(se_hig_cal),
    NA,
    as.numeric(
      results$true_pate[i] >= est_hig_cal - 1.96 * se_hig_cal &
        results$true_pate[i] <= est_hig_cal + 1.96 * se_hig_cal
    )
  )
  results$cover_hig_ot[i] <- ifelse(
    is.na(se_hig_ot),
    NA,
    as.numeric(
      results$true_pate[i] >= est_hig_ot - 1.96 * se_hig_ot &
        results$true_pate[i] <= est_hig_ot + 1.96 * se_hig_ot
    )
  )
  
}

summary_results <- data.frame(
  Method = c(
    "Frobenius-Cal",
    "Frobenius-OT",
    "Higham-Cal",
    "Higham-OT"
  ),
  
  Bias = c(
    mean(results$est_frob_cal - results$true_pate, na.rm=TRUE),
    mean(results$est_frob_ot  - results$true_pate, na.rm=TRUE),
    mean(results$est_hig_cal  - results$true_pate, na.rm=TRUE),
    mean(results$est_hig_ot   - results$true_pate, na.rm=TRUE)
  ),
  
  RMSE = c(
    sqrt(mean((results$est_frob_cal - results$true_pate)^2, na.rm=TRUE)),
    sqrt(mean((results$est_frob_ot  - results$true_pate)^2, na.rm=TRUE)),
    sqrt(mean((results$est_hig_cal  - results$true_pate)^2, na.rm=TRUE)),
    sqrt(mean((results$est_hig_ot   - results$true_pate)^2, na.rm=TRUE))
  ),
  
  Coverage = c(
    mean(results$cover_frob_cal, na.rm=TRUE),
    mean(results$cover_frob_ot,  na.rm=TRUE),
    mean(results$cover_hig_cal,  na.rm=TRUE),
    mean(results$cover_hig_ot,   na.rm=TRUE)
  )
)


# Bias
ggplot(summary_results,
       aes(x = Method, y = Bias, fill = Method)) +
  geom_col(width = .7) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(title = "Bias", y = "Bias")

# RMSE
ggplot(summary_results,
       aes(x = Method, y = RMSE, fill = Method)) +
  geom_col(width = .7) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(title = "RMSE", y = "RMSE")

# Coverage
ggplot(summary_results,
       aes(x = Method, y = Coverage, fill = Method)) +
  geom_col(width = .7) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(title = "Coverage", y = "Coverage")




