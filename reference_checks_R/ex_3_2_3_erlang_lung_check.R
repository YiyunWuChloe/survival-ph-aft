# Reference computation in R for scripts/ch3_phase_type/ex_3_2_3_erlang_lung.jl
# It uses the lung data of the R package survival and prints the numbers that the Julia script checks with @assert.
# Run with:  Rscript reference_checks_R/ex_3_2_3_erlang_lung_check.R
library(survival)
y <- lung$time; d <- lung$status==2
ll <- function(lam, m) {
  f <- dgamma(y, shape=m, rate=lam)               # Erlang(m, lam) density
  S <- pgamma(y, shape=m, rate=lam, lower.tail=FALSE)
  sum(ifelse(d, log(f), log(S)))
}
for (m in 1:4) {
  o <- optimize(function(l) ll(l, m), interval=c(1e-6, 1), maximum=TRUE, tol=1e-14)
  cat("m =", m, " lambda_hat =", o$maximum, " loglik =", o$objective, " mean m/lam =", m/o$maximum, " median =", qgamma(0.5, m, o$maximum), "\n")
  # multi-start check on a grid
  g <- 10^seq(-5, -1, length.out=4001); v <- sapply(g, ll, m=m); cat("   grid max:", g[which.max(v)], max(v), "\n")
}
