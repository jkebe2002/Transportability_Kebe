## Necessary libraries
library('actuar') # For invgamma
library('copula') # For copula functions
library('dplyr')  # Piping function
library(causalOT)
library(expm)
source("optimal_transport_fns.R")

#results df
results <- data.frame(
  run = integer(),
  t = numeric(),
  ESS_opt = numeric(),
  ESS_cot = numeric(),
  bures_t_0 = numeric(),
  bures_t_1 = numeric(),
  bures_t_mid = numeric(),
  Bindex_t = numeric()
)

for (i in 1:100)
{
  #Generate two Multivariate distribuitions
    #case: correlations have opposite sign
    #case: correlations are weaker
  
  #generate populations with low overlap:
  n_1 <- 1000
  n_0 <- 1000
  df1 = data.frame(var1 = rnorm(n_1,0,1),var2 = rbinom(n_1,size = 1, prob = .3),var3 = rbinom(n_1,size = 1, prob = .2), var4 = rnorm(n_1,1,5), var5 = rnorm(n_1,10,10))
  S1 = cov(df1)
  df0 = data.frame(var1 = rnorm(n_0,.5,1),var2 = rbinom(n_0,size = 1, prob = .7),var3 = rbinom(n_0,size = 1, prob = .6), var4 = rnorm(n_0,5,5), var5 = rnorm(n_0,15,10))
  S0 = cov(df0)
  
  binomial_columns <- c(2,3)
  
  
  means_0 <- colMeans(df0)
  means_1 <- colMeans(df1)
  #sds0 <- apply(df0, 2, sd)
  #sds1 <- apply(df1, 2, sd)
  
  #create an actual population that is even mixture of df0,df1
  n_mid <- (n_1+n_0) / 2
  
  df_mid <- rbind(
    
    df0[sample(nrow(df0), n_mid / 2, replace = TRUE), ],
    
    df1[sample(nrow(df1), n_mid / 2, replace = TRUE), ]
    
  )
  Smid <- cov(df_mid)
  
  
  for (i in seq(0,1,.01)) {
    #create synthetic population using correlation structures along the geodesic:
    
    #construct interpolated covariance
    St <- geodesic_step(S1,S0,t)
    #interpolated means
    means_t <- (1-t)*means_1 + t*means_0
    Rt <- cor_from_cov(St)
    dist_t <- create_multivariate_dist(
      Rt,
      list("norm", list(mean = means_t[1], sd = St[1,1])),
      list("norm", list(mean = means_t[2], sd = St[2,2])),
      list("norm", list(mean = means_t[3], sd = St[3,3])),
      list("norm", list(mean = means_t[4], sd = St[4,4])),
      list("norm", list(mean = means_t[5], sd = St[5,5]))
    )
    
    n_t <- round((1-t)*n_1 + t*n_0)
    
    #generate synthetic sample
    synth_samp_t <- generate_sample(n_t, dist_t, binomial_vectors = binomial_columns,
                                    props = list(
                                      list("norm", list(mean = means_t[2], sd = St[2,2])),
                                      list("norm", list(mean = means_t[3], sd = St[3,3]))
                                    ))
    
    #transport to df0
    df_0t <- combine_pops(synth_samp_t, df0)
    
    # call weightit on group var
    tryCatch({
      opt_w <- weightit(group ~ var1 + var2 + var3 + var4 + var5, data = df_0t, focal = 0,  
                        method = "optweight", estimand = "ATT", tols = 0.1)
     ESS_opt <- ESS(get.w(opt_wt))
    }, warning = function(w) {
    ESS_opt <- 0
    }, error = function(e) {
      ESS_opt <- 0
    })
    
    #call COT weight on group var
    z <- df_0t$group
    x <- as.matrix(df_0t[, c("var1", "var2", "var3", "var4", "var5")])
    weights <- calc_weight(
      x = x,
      z = z,
      method = "COT",
      estimand = "ATE",
      options = list(
        lambda.bootstrap = Inf,
        nboot = 1000L
      )
    )
    ESS_cot <- ESS(weights)
    
    #bures distance
    bures_t_0 <-   bures_dist(St,  S0)
    bures_t_1 <-   bures_dist(St,  S1)
    bures_t_mid <- bures_dist(St,Smid)
    
    #generalizability index (Tipton)
    fit <- glm(
      group ~ var1 + var2 + var3 + var4 + var5,
      data = df_0t,
      family = binomial()
    )
  
    # linear predictors (logits)
    logits <- predict(fit, type = "link")
    
    # separate by group
    logits_0 <- logits[df_0t$group == 0]
    logits_1 <- logits[df_0t$group == 1]
    
    Bindex_t <- Bindex(logits_0, logits_1)
    
    
    
    
  }

  
}

# B index function
Bindex <- function(dat1B,dat2B) {
  # bandwidth
  h = function(x){
    n = length(x)
    return((4*sqrt(var(x))^5/(3*n))^(1/5))
  }
  #kernel density
  kg = function(x,data){
    hb = h(data) #bin width
    k = r = length(x)
    for(i in 1:k) r[i] = mean(dnorm((x[i]-data)/hb))/hb
    return(r )
  }
  ##B index calculation
  return( as.numeric(integrate(function(x)
    sqrt(kg(x,dat1B)*kg(x,dat2B)),-Inf,Inf)$value))
}

# Correlation matrix from covariance matrix
cor_from_cov <- function(Sigma) {
  D_inv <- diag(1 / sqrt(diag(Sigma)))
  R <- D_inv %*% Sigma %*% D_inv
  return(R)
}
# geodesic transporting V to U at time t
geodesic_step <- function(V,U,t){
  #identity 
  E <-diag(nrow(V))
  Uhalf <- sqrtm(U)
  #map from V to U
  T <- Uhalf %*% sqrtm(solve((Uhalf %*% V %*% Uhalf))) %*% Uhalf
  
  Wt <- ((1-t) * E + t * T) %*% V %*% ((1-t) * E + t * T)
  return(Wt)
}

#compute bures distance metric
bures_dist <- function(A,B) {
  trace_A <- sum(diag(A))
  trace_B <- sum(diag(B))
  C <- sqrtm(A) %*% B %*% sqrtm(A)
  d <- sqrt(trace_A+trace_B-2*sqrt(sum(diag(C))))
  return(d)
}


# Function to combine sample and population data
combine_pops = function(pop, samp) {
  names(pop) = names(samp)
  pop[['group']] = 0
  samp[['group']] = 1
  data = rbind(pop, samp)
  return(data)
}
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
# save this idea for later 
gen_data <- function(n_0, n_1, prob, sig2, rho = 0, scenario = c("baseline", "positivity", "sample", "outcome", "missing", "sparse")) {
  
  # error variance
  R <- matrix(rho, nrow = 2, ncol = 2)
  diag(R) <- 1
  V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
  Sig <- V %*% R %*% V
  
  # treatment assignment
  z1 <- rbinom(n_1, 1, prob)
  
  if (scenario == "baseline") {
    
    # effect coefficients
    beta <- c(10, -3, -1, 1, 3, 2, 1, -2)
    alpha <- c(5, 3, -1, 1, -3, 2, 1, 2)
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 2)
    x01 <- rbinom(n_0, size = 1, prob = 0.6)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, size = 1, prob = 0.5)
    # added variables to induce dependence
    x04 <- x00+x02+rnorm(n_1, mean=0, sd = .5)
    x05 <- x01*x00
    x06 <- x03*x01
    
    # covariate values
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.4)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, size = 1, prob = 0.5)
    # added variables to induce dependence
    x14 <- x10+x12+rnorm(n_1, mean=0, sd = .5)
    x15 <- x11*x10
    x16 <- x13*x11
    
    theta <- c(1, -1, 0.6, 0, 0.5)
    X0 <- cbind(int0, x00, x01, x02, x03,x04,x05,x06)
    X1 <- cbind(int1, x10, x11, x12, x13,x14,x15,x16)
    
    # Outcome Mean
    mu_10 <- X1 %*% beta
    mu_11 <- X1 %*% beta + X1 %*% alpha
    
    SATE <- mean(X0 %*% alpha)
    PATE <- as.numeric(t(theta) %*% alpha)
    
  } else if (scenario == "positivity") {
    
    # effect coefficients
    beta <- c(10, -3, -1, 1, 3, 2, 1, -2)
    alpha <- c(5, 3, -1, 1, -3, 2, 1, 2)
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = sqrt(2))
    x01 <- rbinom(n_0, size = 1, prob = 0.8)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, size = 1, prob = 0.5)
    # added variables to induce dependence
    x04 <- x00+x02
    x05 <- x01*x00
    x06 <- x03*x01
    
    # covariate values
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = sqrt(2))
    x11 <- rbinom(n_1, size = 1, prob = 0.2)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, size = 1, prob = 0.5)
    # added variables to induce dependence
    x14 <- x10+x12
    x15 <- x11*x10
    x16 <- x13*x11
    
    theta <- c(1, -1, 0.8, 0, 0.5)
    X0 <- cbind(int0, x00, x01, x02, x03,x04,x05,x06)
    X1 <- cbind(int1, x10, x11, x12, x13,x14,x15,x16)
    
    # Outcome Mean
    mu_10 <- X1 %*% beta
    mu_11 <- X1 %*% beta + X1 %*% alpha
    
    SATE <- mean(X0 %*% alpha)
    PATE <- as.numeric(t(theta) %*% alpha)
    
  } else if (scenario == "sample") {
    
    # effect coefficients
    beta <- c(10, -3, -1, 1, 3)
    alpha <- c(5, 3, -1, 1, -3)
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 2)
    x01 <- rbinom(n_0, size = 1, prob = 0.6)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, 1, prob = 0.5)
  
    
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.4)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, 1, prob = 0.5)
    
    u00 <- exp(-0.25*x00 + 0.25*x02)
    u02 <- (0.5*x00 - 0.5*x02)^2
    u10 <- exp(-0.25*x10 + 0.25*x12)
    u12 <- (0.5*x10 - 0.5*x12)^2
    
    u0 <- scale(c(u00, u10))
    u2 <- scale(c(u02, u12))
    
    theta <- c(1, -1, 0.6, 0, 0.5)
    X0 <- cbind(int0, x00 = u0[1:n_0], x01, x02 = u2[1:n_0], x03)
    X0_tmp <- cbind(int0, x00, x01, x02, x03)
    X1 <- cbind(int1, x10 = u0[(n_0+1):(n_0+n_1)], x11, x12 = u2[(n_0+1):(n_0+n_1)], x13)
    X1_tmp <- cbind(int1, x10, x11, x12, x13)
    
    mu_10 <- X1_tmp%*%beta
    mu_11 <- X1_tmp%*%beta + X1_tmp%*%alpha
    
    SATE <- mean(X0_tmp %*% alpha)
    PATE <- NA
    
  } else if (scenario == "outcome") {
    
    # effect coefficients
    beta <- c(10, -3, -1, 1, 3)
    alpha <- c(5, 3, -1, 1, -3)
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 2)
    x01 <- rbinom(n_0, size = 1, prob = 0.6)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, size = 1, prob = 0.5)
    
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.4)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, size = 1, prob = 0.5)
    
    u00 <- exp(-0.25*x00 + 0.25*x02)
    u02 <- (0.5*x00 - 0.5*x02)^2
    u10 <- exp(-0.25*x10 + 0.25*x12)
    u12 <- (0.5*x10 - 0.5*x12)^2
    
    theta <- c(1, -1, 0.6, 0, 0.5)
    X0_tmp <- cbind(int0, u00, x01, u02, x03)
    X0 <- cbind(int0, x00, x01, x02, x03)
    X1_tmp <- cbind(int1, u10, x11, u12, x13)
    X1 <- cbind(int1, x10, x11, x12, x13)
    
    mu_10 <- X1_tmp %*% beta
    mu_11 <- X1_tmp %*% beta + X1_tmp %*% alpha
    
    SATE <- mean(X0_tmp %*% alpha)
    PATE <- NA
    
  } else if (scenario == "missing") {
    
    # effect coefficients
    beta <- c(10, -3, -1, 1, 3)
    alpha <- c(5, 3, -1, 1, -3)
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 2)
    x01 <- rbinom(n_0, size = 1, prob = 0.6)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, size = 1, prob = 0.5)
    
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.4)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, size = 1, prob = 0.5)
    
    theta <- c(1, -1, 0.6, 0, 0.5)
    X0_tmp <- cbind(int0, x00, x01, x02, x03)
    X0 <- cbind(int0, x00, x03)
    X1_tmp <- cbind(int1, x10, x11, x12, x13)
    X1 <- cbind(int1, x10, x13)
    
    mu_10 <- X1_tmp %*% beta
    mu_11 <- X1_tmp %*% beta + X1_tmp %*% alpha
    
    SATE <- mean(X0_tmp %*% alpha)
    PATE <- t(theta)%*%alpha
    
  } else if (scenario == "sparse") {
    
    # effect coefficients
    beta <- c(10, -3, -1, 1, 3, rep(0, times = 4))
    alpha <- c(5, 3, -1, 1, -3, rep(0, times = 4))
    
    # covariate values
    int0 <- rep(1, times = n_0)
    x00 <- rnorm(n_0, mean = -1, sd = 2)
    x01 <- rbinom(n_0, size = 1, prob = 0.6)
    x02 <- rnorm(n_0, mean = 0, sd = 1)
    x03 <- rbinom(n_0, 1, prob = 0.5)
    
    int1 <- rep(1, times = n_1)
    x10 <- rnorm(n_1, mean = 1, sd = 2)
    x11 <- rbinom(n_1, size = 1, prob = 0.4)
    x12 <- rnorm(n_1, mean = 0, sd = 1)
    x13 <- rbinom(n_1, 1, prob = 0.5)
    
    x04 <- rnorm(n_0, mean = 1, sd = 2)
    x05 <- rbinom(n_0, size = 1, prob = 0.4)
    x06 <- rnorm(n_0, mean = 0, sd = 1)
    x07 <- rbinom(n_0, size = 1, prob = 0.5)
    
    x14 <- rnorm(n_1, mean = -1, sd = 2)
    x15 <- rbinom(n_1, size = 1, prob = 0.6)
    x16 <- rnorm(n_1, mean = 0, sd = 1)
    x17 <- rbinom(n_1, size = 1, prob = 0.5)
    
    theta <- c(1, -1, 0.6, 0, 0.5, 1, 0.4, 0, 0.5)
    X0 <- cbind(int0, x00, x01, x02, x03, x04, x05, x06, x07)
    X1 <- cbind(int1, x10, x11, x12, x13, x14, x15, x16, x17)
    
    mu_10 <- X1%*%beta
    mu_11 <- X1%*%beta + X1%*%alpha
    
    SATE <- mean(X0 %*% alpha)
    PATE <- as.numeric(t(theta) %*% alpha)
    
  }
  
  # potential outcomes
  eval <- eigen(Sig, symmetric = TRUE)
  y1_init <- matrix(stats::rnorm(n_1*2, 0, 1), nrow = n_1, ncol = 2) # iid potential outcomes
  y1_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y1_init)) # SVD
  y1_pot <- y1_tmp + cbind(mu_10, mu_11) # include causal effect
  
  # observed outcome
  y1 <- z1*y1_pot[,2] + (1 - z1)*y1_pot[,1]
  
  out <- list(X0 = X0, X1 = X1, y1 = y1, z1 = z1, SATE = SATE, PATE = PATE)
  return(out)
  
}