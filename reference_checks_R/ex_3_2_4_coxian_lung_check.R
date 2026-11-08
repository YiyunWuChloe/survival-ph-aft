# Reference computation in R for scripts/ch3_phase_type/ex_3_2_4_coxian_lung.jl
# It uses the lung data of the R package survival and prints the numbers that the Julia script checks with @assert.
# Run with:  Rscript reference_checks_R/ex_3_2_4_coxian_lung_check.R
library(survival)
y <- lung$time; d <- lung$status==2
# Coxian m = 2: phase 1 -> phase 2 at rate l1, absorption from phase 1 at rate mu1, absorption from phase 2 at rate mu2
Sf <- function(t, l1, mu1, mu2) { a <- l1 + mu1; b <- mu2
  if (abs(a-b) < 1e-9*max(a,b)) exp(-a*t)*(1 + l1*t) else exp(-a*t) + l1*(exp(-b*t) - exp(-a*t))/(a-b) }
ff <- function(t, l1, mu1, mu2) { a <- l1 + mu1; b <- mu2
  if (abs(a-b) < 1e-9*max(a,b)) exp(-a*t)*(mu1 + l1*mu2*t) else mu1*exp(-a*t) + l1*mu2*(exp(-b*t) - exp(-a*t))/(a-b) }
# check f = -S'
tt <- c(50,200,500); h <- 1e-5
chk <- (Sf(tt-h,.004,.001,.006)-Sf(tt+h,.004,.001,.006))/(2*h); print(rbind(chk, ff(tt,.004,.001,.006)))
nll <- function(p) { l1 <- exp(p[1]); mu1 <- exp(p[2]); mu2 <- exp(p[3])
  f <- ff(y,l1,mu1,mu2); S <- Sf(y,l1,mu1,mu2)
  -sum(ifelse(d, log(pmax(f,1e-300)), log(pmax(S,1e-300)))) }
set.seed(2); res <- list()
for (k in 1:300) {
  st <- log(10^runif(3,-4,-1.5))
  r <- try(optim(st, nll, method="Nelder-Mead", control=list(maxit=20000, reltol=1e-15)), silent=TRUE)
  if (!inherits(r,"try-error")) res[[length(res)+1]] <- c(r$value, exp(r$par))
}
M <- do.call(rbind,res); M <- M[order(M[,1]),]
print(head(M,8), digits=7)
cat("number of starts within 1e-3 of the best:", sum(M[,1] < M[1,1]+1e-3), "of", nrow(M), "\n")

l1 <- .004; mu1 <- 0; mu2 <- .004; tt <- c(50,200,500)
print(rbind(num=(Sf(tt-h,l1,mu1,mu2)-Sf(tt+h,l1,mu1,mu2))/(2*h), eq=ff(tt,l1,mu1,mu2), erlang=dgamma(tt,2,.004)))
