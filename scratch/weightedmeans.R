a <- 12
b <- 5
tot <- a+b

a.surv <- 0.9
b.surv <- 0.5

weighted.mean.1 <- a.surv*(a/tot) + b.surv*b/tot
weighted.mean.2 <- sum(c(rep(a.surv,a),rep(b.surv,b)))/tot
                       