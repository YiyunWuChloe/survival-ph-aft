# Reference computation in R for scripts/ch3_phase_type/ex_3_2_1_exponential_lung.jl
# It uses the lung data of the R package survival and prints the numbers that the Julia script checks with @assert.
# Run with:  Rscript reference_checks_R/ex_3_2_1_exponential_lung_check.R
library(survival)
d <- lung
f <- survreg(Surv(time, status==2) ~ 1, data=d, dist="exponential")
mu <- exp(coef(f)); cat("survreg: mean =", mu, " lambda =", 1/mu, " loglik =", f$loglik[1], " se(log mu) =", sqrt(vcov(f)[1,1]), "\n")
dd <- sum(d$status==2); T <- sum(d$time)
cat("closed form: lambda =", dd/T, " mean =", T/dd, " median =", log(2)*T/dd, " se =", (dd/T)/sqrt(dd), "\n")
lam <- dd/T; cat("loglik closed form:", dd*log(lam) - lam*T, "\n")
km <- survfit(Surv(time, status==2) ~ 1, data=d)
for (t in c(100,200,300,400,500,600,700,800)) { i <- max(which(km$time<=t)); cat("t =",t," KM =",km$surv[i]," exp =",exp(-lam*t),"\n") }
print(summary(km)$table[c("median")])
