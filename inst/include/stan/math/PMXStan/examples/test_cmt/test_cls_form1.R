#R
require(rstan)
load("test_cls_form1.RData")
str(dat)
fit = stan("test_cls_form1.stan", data=dat, chains=1, iter=400)

