# Use all RWC and SWP datsa  fit to standard VG curve
# Hierarchical model to account for 4 populations
# Mean porosity is 0.47

library(tidyverse)
library(udunits2)
library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(broom.mixed)

# Load data
wrc <- read_csv("data_clean/rwc-swp.csv") |> 
  mutate(psi_hm = ud.convert(swp, "MPa", "cmH2O") |> 
           ud.convert("cm", "m"),
         population = factor(population, levels = c("CCR", "JLA", "NRV", "TSZ")))

# Quick plot
ggplot(wrc) +
  geom_point(aes(x = -psi_hm, y = rwc,
                 color = population)) +
  scale_x_log10(name = "Pressure Head (-m)",
                breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                labels = c(0.01, 0.1, 1, 10, 100, 1000))


dat_list <- list(theta = wrc$rwc,
                 Nobs = nrow(wrc),
                 theta.s = 0.47, # mean porosity of 8 samples
                 h.cm = -1*wrc$psi_hm,
                 trt = as.numeric(factor(wrc$population)), # 1 = CCR, 2 = JLA, 3 = NRV, 4 = TSZ
                 Ntrt = length(unique(wrc$population)),
                 Salpha = 0.5,
                 Sn = 0.5,
                 Stheta.r = 0.1) 

# Function for initial values
inits <- function() {
  list(mu.log.alpha = runif(1, -5, 0),
       mu.log.n = runif(1, -1, 2),
       mu.log.theta.r = runif(1, -5, -2),
       tau = runif(1, 0, 1),
       tau.eps.alpha = runif(1, 0, 50),
       tau.eps.n = runif(1, 0, 20),
       tau.eps.theta.r = runif(1, 0, 10000))
}
inits_list <- list(inits(), inits(), inits())

# Or, load previous saved state
load("model/inits/inits.Rdata")


# Intialize model
jm <- jags.model(file = "model/VG_hierarchy.jags",
                 data = dat_list,
                 inits = saved_state[[2]],
                 n.chains = 3)

# Update
# update(jm, 10000)

# Monitor coda samples
jm_coda <- coda.samples(jm, 
                        variable.names = c("deviance", "Dsum", "R2", 
                                           "E.alpha.cm", "E.n", "E.theta.r", # population-level parameters on normal scale
                                           "mu.log.alpha", "mu.log.n", "mu.log.theta.r", # populatio-level parameters on log scale, needed to reinitialize
                                           "sig", "tau", # sample variation, sig needed to reinitialize
                                           "sig.alpha", "sig.n", "sig.theta.r", # sigma on log scale
                                           "tau.eps.alpha", "tau.eps.n", "tau.eps.theta.r", # needed to reinitialize
                                           "alpha.cm", "n", "theta.r"), # treatment-level parameters, normal scale
                        n.iter = 150000, thin = 150)

save(jm_coda, file = "model/coda/jm_coda.Rdata")

# Visualize chains
mcmcplot(jm_coda, parms = c("deviance", "Dsum", "R2",
                            "sig", 
                            "sig.alpha", "sig.n", "sig.theta.r",
                            "E.alpha.cm", "mu.log.alpha",
                            "E.n", "mu.log.n",
                            "E.theta.r", "mu.log.theta.r"))

mcmcplot(jm_coda, parms = c("tau", 
                            "tau.eps.alpha",
                            "tau.eps.n",
                            "tau.eps.theta.r"))
caterplot(jm_coda, parms = c("sig.alpha",
                             "sig.n",
                             "sig.theta.r"), reorder = FALSE)
caterplot(jm_coda, parms = c("mu.log.alpha",
                             "mu.log.n",
                             "mu.log.theta.r"), reorder = FALSE)
caterplot(jm_coda, parms = c("E.alpha.cm",
                             "alpha.cm"), reorder = FALSE)
caterplot(jm_coda, parms = c("E.n",
                             "n"), reorder = FALSE)
caterplot(jm_coda, parms = c("E.theta.r",
                             "theta.r"), reorder = FALSE)
# Save state
# newinits <- initfind(jm_coda, OpenBUGS = FALSE)
# newinits[[1]]
# saved_state <- removevars(initsin = newinits,
#                           variables = c(1:6, 10:14, 19))
# saved_state[[1]]
# save(saved_state, file = "model/inits/inits.Rdata") #for local

# Check convergence
gel <- gelman.diag(jm_coda, multivariate = FALSE)
gel

# Monitor replicated samples
jm_rep <- coda.samples(jm,
                       variable.names = "theta.rep",
                       n.iter = 150000, thin = 150)

save(jm_rep, file = "model/coda/jm_rep.Rdata")
