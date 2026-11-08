# Reference computation in R for scripts/ch4_ph_aft/ex_4_4_ph_aft_lung.jl
# It uses the lung data of the R package survival and prints the numbers that the Julia script checks with @assert.
# Run with:  Rscript reference_checks_R/ex_4_4_ph_aft_lung_check.R
library(survival)
d0 <- lung; y <- d0$time; d <- d0$status==2
X <- cbind(age62 = d0$age-62, female = as.numeric(d0$sex==2))
dat <- data.frame(time=y, status=d0$status, X)
sr <- function(dist) { f <- survreg(Surv(time, status==2) ~ age62 + female, data=dat, dist=dist)
  se <- sqrt(diag(vcov(f))); cat(dist, ": loglik =", f$loglik[2], " age =", coef(f)["age62"], "(", se["age62"], ") female =", coef(f)["female"], "(", se["female"], ") scale =", f$scale, "\n") }
sr("exponential"); sr("weibull")
Sf2 <- function(t, l1, mu1, mu2) { a <- l1 + mu1; b <- mu2
  if (abs(a-b) < 1e-9*max(a,b)) exp(-a*t)*(1 + l1*t) else exp(-a*t) + l1*(exp(-b*t) - exp(-a*t))/(a-b) }
ff2 <- function(t, l1, mu1, mu2) { a <- l1 + mu1; b <- mu2
  if (abs(a-b) < 1e-9*max(a,b)) exp(-a*t)*(mu1 + l1*mu2*t) else mu1*exp(-a*t) + l1*mu2*(exp(-b*t) - exp(-a*t))/(a-b) }
Sfun <- list(exp = function(u, r) exp(-r[1]*u), e2 = function(u, r) pgamma(u, 2, r[1], lower.tail=FALSE),
             e3 = function(u, r) pgamma(u, 3, r[1], lower.tail=FALSE), cx = function(u, r) sapply(u, Sf2, l1=r[1], mu1=r[2], mu2=r[3]))
ffun <- list(exp = function(u, r) r[1]*exp(-r[1]*u), e2 = function(u, r) dgamma(u, 2, r[1]),
             e3 = function(u, r) dgamma(u, 3, r[1]), cx = function(u, r) sapply(u, ff2, l1=r[1], mu1=r[2], mu2=r[3]))
nr <- c(exp=1, e2=1, e3=1, cx=3)
nll <- function(p, key) { beta <- p[1:2]; r <- exp(p[-(1:2)]); eta <- as.vector(X %*% beta); u <- y*exp(-eta)
  fi <- exp(-eta)*ffun[[key]](u, r); Si <- Sfun[[key]](u, r)
  -sum(ifelse(d, log(pmax(fi,1e-300)), log(pmax(Si,1e-300)))) }
set.seed(5)
for (key in c("exp","e2","e3","cx")) {
  best <- NULL
  for (k in 1:40) {
    st <- c(rnorm(2,0,.2), log(10^runif(nr[key],-3.3,-1.8)))
    r <- try(optim(st, nll, key=key, method="Nelder-Mead", control=list(maxit=40000, reltol=1e-15)), silent=TRUE)
    if (!inherits(r,"try-error") && (is.null(best) || r$value < best$value)) best <- r
  }
  H <- optimHess(best$par, nll, key=key); V <- try(solve(H)); se <- sqrt(diag(V))
  cat(key, ": loglik =", -best$value, " age =", best$par[1], "(", se[1], ") female =", best$par[2], "(", se[2], ") rates =", signif(exp(best$par[-(1:2)]),6), " eigen(H) min =", min(eigen(H)$values), "\n")
}
