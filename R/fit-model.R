## fit hierarchical Stan model to coil carbon responses to NutNet experiment ####
# requires installation of Rstan software:
# ** instructions are fussy but need to follow them precisely **
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
# 

# load packages
library(rstan)
library(data.table)
library(shinystan)
library(rstantools)
library(loo)

# stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
corz <- options('mc.cores')$mc.cores

# prepare data
source('R/prepare-data.R')
ls()

## model 5: fully hierarchical parameterization ####
## non-centered Cholesky 
## see Stan manual p 147-148
file.show('model/soilC_responses-cholesky-noncentered.stan')

m5c2 <- stan(file = 'model/soilC_responses-cholesky-noncentered.stan',
             chains = 4,
             cores = corz,
             iter = 2000)

save(m5c2, file = paste0('model/stan-m5c2-', Sys.Date(),'.RDS'))

# pairs plot to look at/for divergent transitions
pairs(m5c2, pars=c('lp__', 'gamma[1,1]', 'gamma[1,6]'))
# summary table by parameter
print(m5c2, pars = 'gamma')
# shinystan for interactive validation & estimation
m5c2shiny <- launch_shinystan(m5c2)

# leave-one-out cross validation for pointwise likelihood
m5_loglik <- extract_log_lik(m5c2)
m5_waic <- waic(m5_loglik)
m5_loo <- loo(m5_loglik)
pareto_k_table(m5_loo)

