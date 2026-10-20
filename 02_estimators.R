###=============================================================================
### Doubly-robust transportability estimator
###=============================================================================

library(eventglm)
library(WeightIt)
library(SuperLearner)
library(causalOT)

dr_balance <- function(S, X, Z1, time1, event1, cuts,
                       method = c("optweight", "ebal", "super", "OT"),
                       sl.lib = c("SL.mean", "SL.glm"),
                       tols = 0.05, lambda = NULL) {
  print(method)
  
  method <- match.arg(method)
  
  print(method)
  
  X1 <- subset(X, S == 1)
  n <- nrow(X)
  n1 <- sum(S == 1)
  
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
  if (method == "optweight") {
    print("line 35")
    wfit <- switch(method,
                   super = weightit(fmla, 
                                    method = "super", 
                                    data = d,
                                    estimand = "ATC",
                                    SL.library = sl.lib, 
                                    SL.method = "method.balance",
                                    criterion = "ks.max"),
                   optweight = weightit(fmla, 
                                        method = "optweight", 
                                        data = d,
                                        estimand = "ATC", 
                                        tols = tols),
                   ebal = weightit(fmla, 
                                   method = "ebal", 
                                   data = d, 
                                   estimand = "ATC")
    )
    print("line 54")
    samp <- wfit$weights
    prZ <- mean(Z1)
    
    pseudo.out <- lapply(seq_len(ncol(psi1.mat)), function(i) {
      
      psi1 <- psi1.mat[, i]
      out <- SuperLearner(Y = psi1, X = data.frame(X1, Z = Z1),
                           SL.library = sl.lib, family = gaussian())
      
      mu_tmp <- cbind(predict(out, newdata = data.frame(X, Z = 0))$pred,
                      predict(out, newdata = data.frame(X, Z = 1))$pred)
      
      psi <- numeric(n); psi[S == 1] <- psi1
      Z <- numeric(n); Z[S == 1] <- Z1
      
      w0 <- (1 - Z) * S * samp / (1 - prZ)
      w1 <- Z * S * samp / prZ
      
      mu <- cbind(Z * mu_tmp[, 2] + (1 - Z) * mu_tmp[, 1], mu_tmp)
      
      aug_est <- sum((w1 - w0) * (psi - mu[, 1])) / n1 +
        mean(mu[S == 0, 3] - mu[S == 0, 2])
      eic <- ((w1 - w0) * (psi - mu[, 1])) / mean(S == 1) +
        ((S == 0) * (mu[, 3] - mu[, 2] - aug_est)) / mean(S == 0)
      
      list(estimate = aug_est, variance = var(eic) / n)
      
   
      
    })
    print("line 85")
    
    estimate <- vapply(pseudo.out, `[[`, numeric(1), "estimate")
    variance <- vapply(pseudo.out, `[[`, numeric(1), "variance")
    names(estimate) <- names(variance) <- cuts
    print("line 90")
    
    return(
      list(
        estimate = estimate,
        variance = variance,
        wfit = wfit
      )
    )  }
  if (method == "OT") {
    #create new optimal transport weights (Dunipace GitHub)
    #d <- data.frame(s = S, X = X)
    Xmat <- as.matrix(X)
    weights <- calc_weight(x = Xmat, z = S, 
                           method = "COT",
                           estimand = "ATC", 
                           options = list(lambda = lambda)
    )
  #   m_source <- causalOT::Measure(x = as.matrix(X[S==1, ]),adapt = "weights" )
  #   m_target <- causalOT::Measure(x = as.matrix(X[S==0, ]))
  #   
  #   otp <- causalOT::OTProblem(m_source, m_target)
  #   
  #   #set args
  #   otp$setup_arguments(
  #     lambda = 1, # penalty values of the optimal transport (OT) distances to try
  #     delta = 0, # constraint values to try for balancing functinos
  #     grid.length = 7L, # number of values of lambda and delta to try
  #     # if none are provided
  #     cost.function = NULL, # the ground cost to use between covariates
  #     # default is the Euclidean distance
  #     p = 2, # power to raise the cost by
  #     cost.online = "auto", #Should cost be calculated "online" or "tensorized" (stored in memory). "auto" will try to decide for you
  #     debias = TRUE, # use Sinkhorn divergences (debias = TRUE), i.e. debiased Sinkhorn distances, 
  #     # or use the Sinkhorn distances (debias = FALSE)
  #     diameter = NULL, # the diameter of the covariate space if known
  #     ot_niter = 1000L, # the number of iterations to run when solving OT distances
  #     ot_tol = 0.001 # the tolerance for convergance of OT distances
  #   )
  #   #solve otp prob
  #   otp$solve(
  #     niter = 1000L, # maximum number of iterations
  #     tol = 1e-5, # tolerance for convergence
  #     optimizer = "torch", # which optimizer to use "torch" or "frank-wolfe"
  #     torch_optim = torch::optim_lbfgs, # torch optimizer to use if required
  #     torch_scheduler = torch::lr_reduce_on_plateau, # torch scheduler to use if required
  #     torch_args = list(line_search_fn = "strong_wolfe"), # args passed to the torch functions,
  #     osqp_args = NULL, #arguments passed to the osqp solver used for "frank-wolfe" and balance functions
  #     quick.balance.function = TRUE # if balance functions are also present, should an approximate value of the hyperparameter "delta" be found first
  #   )
  #   otp$choose_hyperparameters(
  #     n_boot_lambda = 100L, #Number of bootstrap iterations to choose lambda
  #     n_boot_delta = 1000L, #Number of bootstrap iterations to choose delta
  #     lambda_bootstrap = Inf # penalty parameter to use for OT distances
  #   )
  #   
  # print(otp$selected_lambda)
    
 
    #sum to group size 
    samp_full <- numeric(n)
    samp_full[S == 1] <- as.numeric(weights@w1) * sum(S == 1)
    

    prZ <- mean(Z1)
    
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
    
    estimate <- vapply(pseudo.out, `[[`, numeric(1), "estimate")
    variance <- vapply(pseudo.out, `[[`, numeric(1), "variance")
    names(estimate) <- names(variance) <- cuts
    
    return(
      list(
        estimate = estimate,
        variance = variance,
        wfit = list(weights = as.numeric(weights@w1))
      )
    )  }
}