# Reference computation in R for scripts/ch3_phase_type/ex_3_2_2_hyperexponential_lung.jl
# It uses the lung data of the R package survival and prints the numbers that the Julia script checks with @assert.
# Run with:  Rscript reference_checks_R/ex_3_2_2_hyperexponential_lung_check.R
library(survival)
y <- lung$time; d <- lung$status==2
nll <- function(p) { a <- plogis(p[1]); l1 <- exp(p[2]); l2 <- exp(p[3])
  S <- a*exp(-l1*y) + (1-a)*exp(-l2*y); f <- a*l1*exp(-l1*y) + (1-a)*l2*exp(-l2*y)
  -sum(ifelse(d, log(f), log(S))) }
best <- NULL
set.seed(1)
for (k in 1:200) {
  st <- c(rnorm(1), log(runif(1,1e-4,1e-2)), log(runif(1,1e-4,1e-2)))
  r <- try(optim(st, nll, method="Nelder-Mead", control=list(maxit=5000, reltol=1e-14)), silent=TRUE)
  if (!inherits(r,"try-error") && (is.null(best) || r$value < best$value)) best <- r
}
a <- plogis(best$par[1]); l <- exp(best$par[2:3])
cat("best loglik:", -best$value, " alpha:", a, " lambdas:", l, "\n")
cat("exponential loglik:", sum(d)*log(sum(d)/sum(y)) - sum(d), "\n")
