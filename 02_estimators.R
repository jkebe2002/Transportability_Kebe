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
      
      list(estimate = aug_est, variance = var(eic) / n )
      
   
      
    })
    print("line 85")
    
    estimate <- vapply(pseudo.out, `[[`, numeric(1), "estimate")
    variance <- vapply(pseudo.out, `[[`, numeric(1), "variance")
    names(estimate) <- names(variance) <- cuts
    print("line 90")
    w_S1 <- samp[S == 1] 
    ess  <- sum(w_S1)^2 / sum(w_S1^2)
    kish_ess <- function(w) sum(w)^2 / sum(w^2)
    return(
      list(
        estimate = estimate,
        variance = variance,
        wfit = list(weights = samp, ESS = kish_ess(samp[S == 1])),
        ESS = ess
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
    
    w1 <- as.numeric(weights@w1)
    ess <- sum(w1)^2 / sum(w1^2)
    
    return(
      list(
        estimate = estimate,
        variance = variance,
        wfit = list(weights = w1, 
                    ESS = ess
                      )
      )
    )  }
}