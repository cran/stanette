require(rstan)
source("ode1.data.R")
fit = stan("ode1.stan", data=.GlobalEnv, chains=1, iter=400)

require(wmisc)
p = as.vector(air(1))
plot(p, conc); abline(0,1, col="red")
