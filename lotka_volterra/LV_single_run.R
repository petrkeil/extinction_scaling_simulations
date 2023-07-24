# setwd("~/Dropbox/Projects/095_OldField2022/noodling/petr/")

rm(list=ls())

require(deSolve)
require(viridis)

system("R CMD SHLIB LVmod_metacom_log.c") # compile C code (needs Rtools)
source("LV_misc_functions.R")

#-------------------------------------------------------------------------------

n0mu = 0.01 # mean starting abundance
n0sd = n0mu/2 # sd starting abundance
Kmu = 1 # mean K
Ksd = 0 # sd K
rmu = 1 # mean r
rsd = 0 # sd r
Amu = -0.8 # mean interaction strength
Asd = 0.2 # sd interaction strength
cmu = 0 # mean dispersal rate
csd = 0.2 # sd dispersal rate (among species)

minval = -999 # log(minval) is lower limit below which abundances are no longer tracked
# note - species only "die" in a cell if their abundance is driven to zero, but the -999 limit lets us avoid
# using lots of computational power to track infintesimally small abundances
M = 5 # number of sites
N0 = 12 # number of species (in global pool)
D = N0 # dimensionality of interaction matrix
steps = 20 # number of simulation steps

# make interaction matrix
Amat = matrix(nrow = N0, ncol = N0)
Amat[row(Amat)!=col(Amat)] = rnorm(N0^2-N0, Amu, Asd)
diag(Amat) = -1

# store prameters in a vector
Pars  <- c(N   = N0,                         # number of species
           M   = M,                          # number of patches
           minval = minval,                  # minimum population size to track
           K   = pmax(0, rnorm(N0, Kmu, Ksd)),   # vector of carrying capacities
           r  = rnorm(N0, rmu, rsd),         # vector of initial growth rates
           A  = c(Amat),                     # interaction coefficients
           d_c = 0.5,                          # proc. noise c
           d_z = 1,                          # proc. noise z
           # >2 means that more abundant species are hit more than less abund.
           d_n = 0,                          # proc. noise nugget
           d_m = 0,                          # proc. noise mean
           d_w = 1,                          # disturbance waiting time
           cv  = pmax(0, rnorm(N0, cmu, csd)),   # dispersal rate
           n0 = abs(rep(n0mu, N0)))          # initial abundance post-colonization

# process noise follows the formula:
#dsd = sqrt(dc*n^dz);
#n_new = n + (dsd+dn)*rnorm(1)+dm;

# make initial abundances
nini  <- abs(rnorm(N0*M, n0mu, n0sd))

# time to record abundances in simulation
times <- seq(0, steps, by = 1)

# log transform initial states - allows for more stable simulations
lnini = log(nini)
lnini[!is.finite(lnini)] = minval

# dummy variable for tracking disturbance times
lnini = c(lnini, 1, runif(1))

# dummy variable for tracking dispersal effects
lnini = c(lnini, rep(1, N0), runif(N0, 0.5, 1))

dyn.load("LVmod_metacom_log.so") # load c code
# run C code
out_C <- ode(y = lnini, 
             times = times, 
             func = "derivs",
             parms = Pars, 
             dllname = "LVmod_metacom_log", 
             initfunc = "initmod",
             events = list(func = "event", root = TRUE), 
             rootfun = "myroot", 
             nout = N0, 
             nroot = N0+1)
dyn.unload("LVmod_metacom_log.so") # unload c code (helps with stablity)

# output: array with 3 dimensions (time, species, sites)
Mout = array(dim =c(nrow(out_C), N0, M))
Mout[] = exp(out_C[,1:(N0*M)+1])


Mout.gamma <- apply(Mout, 1:2, sum)

# -----------------------------------

time1 = 5
time2 = 12

MAT <- Mout.gamma
EXTract(MAT, time1, time2)
PXTract(MAT, time1, time2)


# number of extinctions at gamma scale
P.ext.gamma <- PXTract(Mout.gamma, 
                           time1, 
                           time2)
P.ext.gamma

# mean number of extinctions at alpha scale
P.ext.alpha <- mean(apply(X = Mout,
                          FUN = PXTract, 
                          MARGIN = 3, 
                          time1, 
                          time2))
P.ext.alpha

# plot total biomass of each species at the metacommunity scale
matplot(Mout.gamma, type = "l")
abline(h=0)
abline(v = time1); abline(v=time2)

# plot mean output at gamma scale
par(mfrow=c(2,2), mar=c(4,4,2,2))
matplot(times, apply(Mout, 1:2, mean), lty = 1, type = "l", lwd = 2,
        xlab = "time", ylab = "Species Abundance n, Gamma", col = viridis(N0))
abline(h=0, lty=3)

#  plot alpha-level dynaimcs for each site
par(mfrow=c(2,2), mar=c(4,4,2,2))
for(i in 1:M) {
  matplot(times, Mout[,,i], lty = 1, type = "l", lwd = 2,
          xlab = "time", ylab = "n, alpha", col = viridis(N0))
  abline(h=0, lty=3)
}

# plot fraction of sites occupied over time for each species
Pout = out_C[,ncol(out_C)+1-(N0:1)]
par(mfrow=c(2,2), mar=c(4,4,2,2))
matplot(times, Pout, lty = 1, type = "l", lwd = 2,
        xlab = "time", ylab = "p", col = viridis(N0))
abline(h=0, lty=3)

# show mortality rate vs. abundance
mortality_occurred = (Mout[-1,,] == 0 & Mout[-dim(Mout)[1],,] > 0)
mortality_occurred[Mout[-dim(Mout)[1],,] == 0] = NA # ignore cases where species is already absent

# mortality probability per time-step for each species
pmor = apply(mortality_occurred, c(1,2), function(x) mean(x, na.rm=T))/mean(diff(times))
pmor_N = colMeans(pmor, na.rm=TRUE) # mortality probability for each species
hist(pmor_N)

# mean non-zero abundance per time-step for each species
nmu = apply(Mout, c(1,2), function(x) mean(x[x>0]))
nmu_N = colMeans(nmu, na.rm=T) # mean abundance for each species
hist(nmu_N)

# compare pmor vs. abundance
plot(nmu_N, pmor_N)

# analytical expectation
nsq = seq(min(Mout[Mout>0]), max(Mout), length=1e3)
dsd = sqrt(Pars["d_c"]*nsq^Pars["d_z"])+Pars["d_n"]
pr_mor_analyical = pnorm(0, nsq, dsd)
lines(nsq, pr_mor_analyical, lty = 2)

# try running whole script comparing with d_z < 2 vs. d_z > 2



