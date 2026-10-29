# estimtors for interpolation simulation

###=============================================================================
### Doubly-robust transportability estimator
###=============================================================================

library(eventglm)
library(WeightIt)
library(SuperLearner)
library(causalOT)

dr_balance_OT <- function(S, X, Z1, time1, event1, cuts, t_step, method = "OT_t",
                       sl.lib = c("SL.mean", "SL.glm"),
                       tols = 0.05) {

  
  X1 <- subset(X, S == 1)
  X0 <- subset(X, S == 0)
  

  Xmat1 <- as.matrix(X1)
  Xmat0 <- as.matrix(X0)
  
  # calculate the OT weights
  m1 <- causalOT::Measure(x = Xmat1, adapt = "none")
  m1_w <- causalOT::Measure(x = Xmat1, adapt = "weights")
  m0 <- causalOT::Measure(x = Xmat0, adapt = "none")
  otp01 <- OTProblem(m1_w, m0)
  otp01$setup_arguments(lambda = 4, delta = NULL, grid.length = 7L,
                         cost.function = NULL, p = 2, cost.online = "auto",
                         debias = TRUE, diameter = NULL,
                         ot_niter = 1000L, ot_tol = 1e-3)
  otp01$solve(
    niter = 1000L, # maximum number of iterations
    tol = 1e-5, # tolerance for convergence
    optimizer = "torch", # which optimizer to use "torch" or "frank-wolfe"
    torch_optim = torch::optim_lbfgs, # torch optimizer to use if required
    torch_scheduler = torch::lr_reduce_on_plateau, # torch scheduler to use if required
    torch_args = list(line_search_fn = "strong_wolfe"), # args passed to the torch functions,
    osqp_args = NULL, #arguments passed to the osqp solver used for "frank-wolfe" and balance functions
    quick.balance.function = TRUE # if balance functions are also present, should an approximate value of the hyperparameter "delta" be found first
  )
  
  otp11 <- OTProblem(m1_w, m1)
  
  n <- nrow(X)
  n1 <- sum(S == 1)
  
  #weights scaled 
  wts_OT <- n1*as.numeric(m1_w$weights)
  
  
  ps_dat <- data.frame(t_ = time1, e_ = event1)
  psi1.mat <- vapply(cuts, function(tt) {
    eventglm::pseudo_infjack(
      formula = Surv(t_, e_) ~ 1,
      time = tt,
      cause = 1,
      type = "survival",
      data = ps_dat
    )
  }, numeric(length(time1)))
  
  print("line 31")
  fmla <- as.formula(paste("s ~", paste(colnames(X), collapse = " + ")))
  d <- data.frame(s = S, X)
  
  #Interpolate by t

  #   what are you doing at each t?
  
  #   create interpolated weights t*1 + (1-t) * OT weights
  
  #   create OT_interpolated weights: t*otp1+(1-t)*otp0
  
  #   Compare the weight distributions, the ESS, the W2 distance from the target, the source, and the OT-weighted source, Tipton index, outcome and variance 
  prZ <- mean(Z1)
  est_effect <- function(weights) {
    samp_full <- numeric(n)
    samp_full[S == 1] <- weights
    
    pseudo.out <- lapply(seq_len(ncol(psi1.mat)), function(i) {
      
      psi1 <- psi1.mat[, i]
      out <- SuperLearner(Y = psi1, X = data.frame(X1, Z = Z1),
                          SL.library = sl.lib, family = gaussian())
      
      mu_tmp <- cbind(predict(out, newdata = data.frame(X, Z = 0))$pred,
                      predict(out, newdata = data.frame(X, Z = 1))$pred)
      
      psi <- numeric(n); psi[S == 1] <- psi1
      Z <- numeric(n); Z[S == 1] <- Z1
      
      w0 <- (1 - Z) * S * samp_full / (1 - prZ)
      w1 <- Z * S * samp_full / prZ    
      
      mu <- cbind(Z * mu_tmp[, 2] + (1 - Z) * mu_tmp[, 1], mu_tmp)
      
      aug_est <- sum((w1 - w0) * (psi - mu[, 1])) / n1 +
        mean(mu[S == 0, 3] - mu[S == 0, 2])
      eic <- ((w1 - w0) * (psi - mu[, 1])) / mean(S == 1) +
        ((S == 0) * (mu[, 3] - mu[, 2] - aug_est)) / mean(S == 0)
      
      list(estimate = aug_est, variance = var(eic) / n)
    })
    
  }
  for (t in seq(0,1, t_step)) {
    wts_OT_t <- t*1 + (1-t)* wts_OT
    
    otp_interpol <- t*otp11 + (1-t)*otp01
    otp_interpol$setup_arguments(lambda = 4, delta = NULL, grid.length = 7L,
                           cost.function = NULL, p = 2, cost.online = "auto",
                           debias = TRUE, diameter = NULL,
                           ot_niter = 1000L, ot_tol = 1e-3)
    otp_interpol$solve(
      niter = 1000L, # maximum number of iterations
      tol = 1e-5, # tolerance for convergence
      optimizer = "torch", # which optimizer to use "torch" or "frank-wolfe"
      torch_optim = torch::optim_lbfgs, # torch optimizer to use if required
      torch_scheduler = torch::lr_reduce_on_plateau, # torch scheduler to use if required
      torch_args = list(line_search_fn = "strong_wolfe"), # args passed to the torch functions,
      osqp_args = NULL, #arguments passed to the osqp solver used for "frank-wolfe" and balance functions
      quick.balance.function = TRUE # if balance functions are also present, should an approximate value of the hyperparameter "delta" be found first
    )
    
    wts_interpol <- n1*as.numeric(m1_test$weights)
    
    est_results <- function(method) {
    if (method == "OT_t") {
      #create new optimal transport weights (Dunipace GitHub)
      #d <- data.frame(s = S, X = X)
      res <- est_effect(wts_OT_t)
      estimate <- res$estimate
      variance <- res$variance
      ess <- sum(weights)^2 / sum(weights^2)
      #W2 distance from the target, the source, and the OT-weighted source
      
      ot_target <- ot_distance(m1_test$x, m0$x, 
                               a = weights, b = m0$init_weights, penalty = 1)
      ot_source <- ot_distance(m1_test$x, m1$x, 
                            a = weights, b = m1$init_weights, penalty = 1)
      ot_source_w <- ot_distance(m1_test$x, m1_test$x, 
                                 a = weights, b = wts_OT, penalty = 1)

      
      return(
        list(
          estimate = estimate,
          variance = variance,
          wfit = list(weights = wts_OT_t, 
                      ESS = ess
          )
        , method = method, t=t, minw = min(wts_OT_t), maxw = max(wts_OT_t), meanw = mean(wts_OT_t), medianw = median(wts_OT_t))
      )  }
    
    if (method == "OT_interpol") {
      #create new optimal transport weights (Dunipace GitHub)
      #d <- data.frame(s = S, X = X)
      res <- est_effect(wts_interpol)
      estimate <- res$estimate
      variance <- res$variance
      ess <- sum(weights)^2 / sum(weights^2)
      
      ot_target <- ot_distance(m1_test$x, m0$x, 
                               a = weights, b = m0$init_weights, penalty = 1)
      ot_source <- ot_distance(m1_test$x, m1$x, 
                               a = weights, b = m1$init_weights, penalty = 1)
      ot_source_w <- ot_distance(m1_test$x, m1_test$x, 
                                 a = weights, b = wts_OT, penalty = 1)
      
      return(
        list(
          estimate = estimate,
          variance = variance,
          wfit = list(weights = wts_OT_t, 
                      ESS = ess
          ), method = method, t=t, minw = min(wts_interpol), maxw = max(wts_interpol), meanw = mean(wts_interpol), medianw = median(wts_interpol),
          ot_target = ot_target,
          ot_source = ot_source,
          ot_source_w = ot_source_w
        )
      )  }
    }
    } #this bracket ends the for loop
}









