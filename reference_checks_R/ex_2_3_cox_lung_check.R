# Reference computation in R for scripts/ch2_cox/ex_2_3_cox_lung.jl
# It uses the lung data of the R package survival and prints the numbers that the Julia script checks with @assert.
# Run with:  Rscript reference_checks_R/ex_2_3_cox_lung_check.R
library(survival)
d <- lung
d$age62 <- d$age - 62
d$female <- as.numeric(d$sex == 2)
cat("n =", nrow(d), " events =", sum(d$status==2), " mean age =", mean(d$age), " n women =", sum(d$female), "\n")
cat("distinct event times:", length(unique(d$time[d$status==2])), " events at tied times:", sum(duplicated(d$time[d$status==2])), "\n")
fit <- coxph(Surv(time, status==2) ~ age62 + female, data=d, ties="breslow")
print(summary(fit), digits=6)
cat("loglik:", fit$loglik, "\n")
bh <- basehaz(fit, centered=FALSE)
for (t in c(100,200,300,400,500,700,1000)) cat("H0(",t,") =", bh$hazard[max(which(bh$time<=t))], "\n")
nd <- data.frame(age62=c(0,0), female=c(0,1))
sf <- survfit(fit, newdata=nd)
for (t in c(100,200,300,500,700)) { i<-max(which(sf$time<=t)); cat("S(",t,"): men",sf$surv[i,1]," women",sf$surv[i,2],"\n") }
# median survival from model
cat("model medians:", summary(sf)$table[, "median"], "\n")
