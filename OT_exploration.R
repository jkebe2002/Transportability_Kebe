#generate populations with low overlap:
library(torch)
library(causalOT)
torch::install_torch()
#source("~/Documents/GitHub/Transportability_Kebe/cotOOP.R")
n_1 <- 1000
n_0 <- 1000
df1 = data.frame(var1 = rnorm(n_1,0,1),var2 = rbinom(n_1,size = 1, prob = .3),var3 = rbinom(n_1,size = 1, prob = .2), var4 = rnorm(n_1,1,5), var5 = rnorm(n_1,10,10), S = 1)
S1 = cov(df1)
df0 = data.frame(var1 = rnorm(n_0,.5,1),var2 = rbinom(n_0,size = 1, prob = .7),var3 = rbinom(n_0,size = 1, prob = .6), var4 = rnorm(n_0,5,5), var5 = rnorm(n_0,15,10), S = 0)
S0 = cov(df0)


df_all <- rbind(df1,df0)

x1 <- as.matrix(df1)
x0 <- as.matrix(df0)
x_all <- as.matrix(df_all)

m1 <- causalOT::Measure(x = x1[, -ncol(x1)], adapt = "none")
m1_test <- causalOT::Measure(x = x1[, -ncol(x1)], adapt = "weights")

m0 <- causalOT::Measure(x = x0[, -ncol(x0)], adapt = "none")
m_all <- causalOT::Measure(x = x_all[,-ncol(x0)], adapt = "weights")

otp0 <- OTProblem(m_all, m0)
otp1 <- OTProblem(m_all, m1)

m1 <- causalOT::Measure(x = x1[, -ncol(x1)], adapt = "none")
m1_test <- causalOT::Measure(x = x1[, -ncol(x1)], adapt = "weights")
otp_test <- OTProblem(m1_test, m1)
otp_test$setup_arguments(lambda = 100, delta = NULL, grid.length = 7L,
                      cost.function = NULL, p = 2, cost.online = "auto",
                      debias = TRUE, diameter = NULL,
                      ot_niter = 5000L, ot_tol = 1e-5)
otp_test$solve(
  niter = 1000L, # maximum number of iterations
  tol = 1e-5, # tolerance for convergence
  optimizer = "torch", # which optimizer to use "torch" or "frank-wolfe"
  torch_optim = torch::optim_lbfgs, # torch optimizer to use if required
  torch_scheduler = torch::lr_reduce_on_plateau, # torch scheduler to use if required
  torch_args = list(line_search_fn = "strong_wolfe"), # args passed to the torch functions,
  osqp_args = NULL, #arguments passed to the osqp solver used for "frank-wolfe" and balance functions
  quick.balance.function = TRUE # if balance functions are also present, should an approximate value of the hyperparameter "delta" be found first
)

t = 1
otp_t <- t * otp1 + (1-t)*otp0
otp_t$setup_arguments(lambda = 0, delta = NULL, grid.length = 7L,
                      cost.function = NULL, p = 2, cost.online = "auto",
                      debias = TRUE, diameter = NULL,
                      ot_niter = 10000L, ot_tol = 1e-7)
otp_t$solve(
  niter = 1000L, # maximum number of iterations
  tol = 1e-5, # tolerance for convergence
  optimizer = "torch", # which optimizer to use "torch" or "frank-wolfe"
  torch_optim = torch::optim_lbfgs, # torch optimizer to use if required
  torch_scheduler = torch::lr_reduce_on_plateau, # torch scheduler to use if required
  torch_args = list(line_search_fn = "strong_wolfe"), # args passed to the torch functions,
  osqp_args = NULL, #arguments passed to the osqp solver used for "frank-wolfe" and balance functions
  quick.balance.function = TRUE # if balance functions are also present, should an approximate value of the hyperparameter "delta" be found first
)

x_mat <- as.matrix(m_all$x$to(device = "cpu"))   # the $to() is harmless on CPU
colnames(x_mat) <- setdiff(colnames(df_all), "S")   # var1 ... var5
w <- as.numeric(m_all$weights$detach()$to(device = "cpu"))

new_df <- data.frame(x_mat, S = df_all$S, weights = w)



otp_t$loss

otp_t$


