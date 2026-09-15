# Install and load the package
library(expm)


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
  d <- sqrtm(trace_A+trace_B-2*sqrt(sum(diag(C))))
  return(d)
}

