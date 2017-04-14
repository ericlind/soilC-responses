## map2stan and stan model exploration
## WARNING: copied from "sandbox" sessions, not guaranteed to work
## 
## For current model see: /model and /R/fit-model.R 
## 
## Model outline:
##  1a. overall model: impact of N,  P,  NP with site-level intercepts (map2stan)
##  1b. intercept predicted by climate & soil & production (map2stan)
##  1c. model with climate var interactions with plot trt (map2stan)
##  2b. model with pH, above- and below-ground biomass at plot scale (stan) 
##  2c. plots level vars are correlated within sites (stan)

## 1a. overall model: impact of N,  P,  NP with site-level intercepts ####
datmod[, site_id := as.integer(factor(site_code))]
datmod[, unique(site_id), site_code][order(V1), site_code]
datmod
plotdat <- unique(datmod[, .(site_code, site_id, block, plot, nit, pho, np, soil_gCm2, live_mass_ln, lnroots, pH_num)])
sitedat <- unique(datmod[, .(map_z, map_var_z, mat_z, temp_var_z, texture_z, site_pH_z, site_above_mass_z, site_below_mass_z, site_id)])
sitedat[, lapply(.SD, as.list), .SDcols=c('site_id')]
as.list(sitedat)

mSoilCX <- map2stan(
  alist(
    soil_gCm2 ~ dlnorm(mu, sigma), 
    mu <- aS[site_id] + bN*nit + bP*pho + bNP*np, 
    aS[site_id] ~ dnorm(a, sigS), 
    c(bN, bP, bNP) ~ dnorm(0, 2), 
    a ~ dnorm(0, 5), 
    sigS ~ dcauchy(0, 2), 
    sigma ~ dcauchy(0, 2)
  ),  data=c(as.list(plotdat), as.list(sitedat)), 
  #chains=1
  cores=4, chains = 4, warmup=1000, iter=2000
)

# check convergence
m1 <- mSoilCX
plot(mSoilCX)
layout(1)
plot(precis(mSoilCX))
plot(precis(mSoilCX), pars=1:3)
plot(precis(mSoilCX, depth=2, pars = c('a', 'aS')))

# posterior predictive check
post1 <- extract.samples(mSoilCX)
names(post1)
lapply(post1, dim) 
#soilC_pred <- data.table(site_id=order(as.character(datmod[, unique(site_id)])))
soilC_pred <- data.table(site_id=1:26)
soilC_pred[, est_mn_soilC:=apply(post1$aS, 2, mean)]
soilC_pred[, est_soilC_lo:=apply(post1$aS, 2, PI)[1, ]]
soilC_pred[, est_soilC_hi:=apply(post1$aS, 2, PI)[2, ]]

p <- ggplot(plotdat, aes(x=factor(site_id), y=log(soil_gCm2)))
p + geom_boxplot() +
  geom_pointrange(data=soilC_pred, aes(y=est_mn_soilC, ymin=est_soilC_lo, ymax=est_soilC_hi), pch=16, col='blue') +
  geom_hline(yintercept = mean(post1$a), lty=1) +
  geom_hline(yintercept = PI(post1$a), lty=3) +
  scale_x_discrete(labels=datmod[, unique(site_code)]) +
  theme(axis.text.x=element_text(angle=90,  hjust=1))

str(data.frame(datmod))
save.image("soil-resp-models.RData")


## 1b. intercept predicted by climate & soil & production ####
str(datmod)
as.list(plotdat)
mSoilCX2 <- map2stan(
  alist(
    soil_gCm2 ~ dlnorm(mu, sigma), 
    mu <- A + bN*nit + bP*pho + bNP*np , 
    A <- aS[site_id] + bmap*map_z + bmapvar*map_var_z + bmat*mat_z + btempvar*temp_var_z + btexture*texture_z + bpH*site_pH_z + bAM*site_above_mass_z + bBM*site_below_mass_z, 
    aS[site_id] ~ dnorm(a, sigS), 
    c(bN, bP, bNP) ~ dnorm(0, 2), 
    c(bmap, bmapvar, bmat, btempvar, btexture, bpH, bAM, bBM) ~ dnorm(0, 2), 
    a ~ dnorm(0, 5), 
    sigS ~ dcauchy(0, 2), 
    sigma ~ dcauchy(0, 2)
  ),  data=c(as.list(datmod)), 
  #chains=1
  cores=4, chains = 4, warmup=1000, iter=2000
)
m2 <- mSoilCX2

# check convergence
m2 <- mSoilCX2
save.image("soil-resp-models.RData")

plot(mSoilCX2)
stancode(mSoilCX2) 
layout(1)
plot(precis(mSoilCX2), pars=1:9)
plot(precis(mSoilCX2, depth=2, pars = c('a', 'aS')))

compare(m2, m1)
str(datmod)


## 1c. model with climate var interactions with plot trt ####
mSoilCX3 <- map2stan(
  alist(
    # likelihood
    soil_gCm2 ~ dlnorm(mu, sigma), 
    
    # linear models
    mu <- A + BN*nit + BP*pho + BNP*np, 
    A <- aS[site_id] + a_map*map_z + a_mapvar*map_var_z + a_mat*mat_z + a_tempvar*temp_var_z + a_texture*texture_z + a_pH*site_pH_z, 
    BN <- bN_S[site_id] + bN_map*map_z + bN_mapvar*map_var_z + bN_mat*mat_z + bN_tempvar*temp_var_z + bN_texture*texture_z + bN_pH*site_pH_z + bN_AM*site_above_mass_z + bN_BM*site_below_mass_z, 
    BP <- bP_S[site_id] + bP_map*map_z + bP_mapvar*map_var_z + bP_mat*mat_z + bP_tempvar*temp_var_z + bP_texture*texture_z + bP_pH*site_pH_z + bP_AM*site_above_mass_z + bP_BM*site_below_mass_z, 
    BNP <- bNP_S[site_id] + bNP_map*map_z + bNP_mapvar*map_var_z + bNP_mat*mat_z + bNP_tempvar*temp_var_z + bNP_texture*texture_z + bNP_pH*site_pH_z + bNP_AM*site_above_mass_z + bNP_BM*site_below_mass_z, 
    
    # mvnorm priors
    c(aS,  bN_S,  bP_S, bNP_S)[site_id] ~ dmvnorm2(c(a, bN, bP, bNP), sigma_site, Rho_site), 
    
    # hyperpriors
    c(bN, bP, bNP) ~ dnorm(0, 2), 
    c(a_map, a_mapvar, a_mat, a_tempvar, a_texture,  a_pH) ~ dnorm(0, 2), 
    c(bN_map, bN_mapvar, bN_mat, bN_tempvar, bN_texture,  bN_pH, bN_AM, bN_BM) ~ dnorm(0, 2), 
    c(bP_map, bP_mapvar, bP_mat, bP_tempvar, bP_texture,  bP_pH, bP_AM, bP_BM) ~ dnorm(0, 2), 
    c(bNP_map, bNP_mapvar, bNP_mat, bNP_tempvar, bNP_texture,  bNP_pH,  bNP_AM, bNP_BM) ~ dnorm(0, 2), 
    a ~ dnorm(0, 5), 
    sigma_site ~ dcauchy(0, 2), 
    Rho_site ~ dlkjcorr(4), 
    sigma ~ dcauchy(0, 2)
  ),  data=data.frame(datmod), 
  #chains=1,
  control = list(adapt_delta = 0.99),
  cores=4, chains = 4, warmup=1000, iter=4000
)

m3 <- mSoilCX3

save.image("soil-resp-models.RData")
sink('Stan-model-code.txt')
stancode(mSoilCX3)
sink()
file.show('Stan-model-code.txt')
list.files("", pattern="RData")
load("soil-resp-models.RData")

# check convergence
plot(mSoilCX3)
pairs(mSoilCX3, pars = c('bN_pH', 'bN_map', 'bN_mapvar', 'bN_mat', 'bN_tempvar', 'bN_AM', 'bN_BM'))
pairs(mSoilCX3, pars = c('bP_pH', 'bP_map', 'bP_mapvar', 'bP_mat', 'bP_tempvar', 'bP_AM', 'bP_BM'))
pairs(mSoilCX3, pars = c('bNP_pH', 'bNP_map', 'bNP_mapvar', 'bNP_mat', 'bNP_tempvar', 'bNP_AM', 'bNP_BM'))
pairs(mSoilCX3, pars = c('bN', 'bP', 'bNP'))

precis(mSoilCX3, depth=2)
layout(1)
plot(precis(mSoilCX3), pars=1:27)
plot(precis(mSoilCX3, depth=2, pars = c('a', 'aS')))
plot(precis(mSoilCX3, depth=2, pars = c('bN', 'bN_S')))
plot(precis(mSoilCX3, depth=2, pars = c('bP', 'bP_S')))
plot(precis(mSoilCX3, depth=2, pars = c('bNP', 'bNP_S')))
plot(precis(mSoilCX3, depth=2, pars = c('a_map', 'a_mapvar', 'a_mat', 'a_tempvar')))
plot(precis(mSoilCX3, depth=2, pars = c('bN_pH', 'bN_map', 'bN_mapvar', 'bN_mat', 'bN_tempvar', 'bNP_pH', 'bNP_AM', 'bNP_BM')))
plot(precis(mSoilCX3, depth=2, pars = c('bP_map', 'bP_mapvar', 'bP_mat', 'bP_tempvar', 'bP_pH', 'bP_AM', 'bP_BM')))
plot(precis(mSoilCX3, depth=2, pars = c('bNP_map', 'bNP_mapvar', 'bNP_mat', 'bNP_tempvar', 'bNP_pH', 'bNP_AM', 'bNP_BM')))
plot(precis(mSoilCX3, depth=2, pars=c('Rho_site')))

compare(m1, m2, m3)

# posterior plots
# posterior predictive check
post3 <- extract.samples(mSoilCX3)
names(post3)
lapply(post3, dim) 

soilC_pred <- unique(datmod[, .(site_id, nit, pho, np, site_code, map_z, MAP, map_var_z, MAP_VAR, mat_z, MAT, temp_var_z, TEMP_VAR, texture_z, texture, site_pH_z, site_pH, site_above_mass_z, site_above_mass, site_below_mass_z, site_below_mass)])
out <- link(mSoilCX3, soilC_pred)  
lapply(out, dim)
soilC_pred[, est_mn_soilC:=apply(out$mu, 2, mean)]
soilC_pred[, est_soilC_lo:=apply(out$mu, 2, PI)[1, ]]
soilC_pred[, est_soilC_hi:=apply(out$mu, 2, PI)[2, ]]
soilC_pred[, A_mn:=apply(out$A, 2, mean)]
soilC_pred[, A_lo:=apply(out$A, 2, PI)[1, ]]
soilC_pred[, A_hi:=apply(out$A, 2, PI)[2, ]]
soilC_pred[, BN_mn:=apply(out$BN, 2, mean)]
soilC_pred[, BN_lo:=apply(out$BN, 2, PI)[1, ]]
soilC_pred[, BN_hi:=apply(out$BN, 2, PI)[2, ]]
soilC_pred[, BP_mn:=apply(out$BP, 2, mean)]
soilC_pred[, BP_lo:=apply(out$BP, 2, PI)[1, ]]
soilC_pred[, BP_hi:=apply(out$BP, 2, PI)[2, ]]
soilC_pred[, BNP_mn:=apply(out$BNP, 2, mean)]
soilC_pred[, BNP_lo:=apply(out$BNP, 2, PI)[1, ]]
soilC_pred[, BNP_hi:=apply(out$BNP, 2, PI)[2, ]]

# overall posterior predictive
soilC_pred[nit==0 & pho==0, trt:='Control']
soilC_pred[nit==1 & pho==0, trt:='N']
soilC_pred[nit==0 & pho==1, trt:='P']
soilC_pred[nit==1 & pho==1, trt:='NP']

p <- ggplot(datmod, aes(x=factor(site_id), y=log(soil_gCm2)))
p + geom_boxplot() +
  geom_pointrange(data=soilC_pred, aes(y=est_mn_soilC, ymin=est_soilC_lo, ymax=est_soilC_hi), pch=16, col='blue') +
  geom_hline(yintercept = mean(out$mu), lty=1) +
  geom_hline(yintercept = PI(out$mu), lty=3) +
  scale_x_discrete(labels=datmod[, unique(site_code)]) +
  theme(axis.text.x=element_text(angle=90,  hjust=1)) + 
  facet_wrap(~trt)


## 1c. plots: Intercepts & Interactions on raw scale ####
p <- ggplot(soilC_pred, aes(x=MAP, y=est_mn_soilC))

p +
  geom_pointrange(aes(y=est_mn_soilC, ymin=est_soilC_lo, ymax=est_soilC_hi, shape=trt)) +
  geom_hline(yintercept = mean(out$mu), lty=1) +
  geom_hline(yintercept = PI(out$mu), lty=3) +
  facet_wrap(~trt)

# MAP 
p <- ggplot(soilC_pred, aes(x=MAP, y=A_mn))
p +
  geom_pointrange(aes(ymin=A_lo, ymax=A_hi)) +
  geom_smooth(method='lm') +
  ylab('site soil C g_m2')

# MAP_VAR 
p <- ggplot(soilC_pred, aes(x=MAP_VAR, y=A_mn))
p +
  geom_pointrange(aes(ymin=A_lo, ymax=A_hi)) +
  geom_smooth(method='lm', lty=3, se=F)

# MAT
p <- ggplot(soilC_pred, aes(x=MAT, y=A_mn))
p +
  geom_pointrange(aes(ymin=A_lo, ymax=A_hi)) +
  geom_smooth(method='lm', lty=3, se=F)

# TEMP_VAR 
p <- ggplot(soilC_pred, aes(x=TEMP_VAR, y=A_mn))
p +
  geom_pointrange(aes(ymin=A_lo, ymax=A_hi)) +
  geom_smooth(method='lm')+
  ylab('site soil C g_m2')

# texture
p <- ggplot(soilC_pred, aes(x=texture, y=A_mn))
p +
  geom_pointrange(aes(ymin=A_lo, ymax=A_hi)) +
  geom_smooth(method='lm')+
  ylab('site soil C g_m2')

# pH
p <- ggplot(soilC_pred, aes(x=site_pH, y=A_mn))
p +
  geom_pointrange(aes(ymin=A_lo, ymax=A_hi)) +
  geom_smooth(method='lm')+
  ylab('site soil C g_m2')

# above ground
p <- ggplot(soilC_pred, aes(x=site_above_mass, y=A_mn))
p +
  geom_pointrange(aes(ymin=A_lo, ymax=A_hi)) +
  geom_smooth(method='lm')+
  ylab('site soil C g_m2')

# below ground
p <- ggplot(soilC_pred, aes(x=site_below_mass, y=A_mn))
p +
  geom_pointrange(aes(ymin=A_lo, ymax=A_hi)) +
  geom_smooth(method='lm')+
  ylab('site soil C g_m2')

# plots: N effect size
# MAP 
p <- ggplot(soilC_pred, aes(x=MAP, y=BN_mn))
p +
  geom_pointrange(aes(ymin=BN_lo, ymax=BN_hi)) +
  #  geom_smooth(method='lm')  +
  ylab('effect of adding N on soil C')

# pH
p <- ggplot(soilC_pred, aes(x=site_pH, y=BN_mn))
p +
  geom_pointrange(aes(ymin=BN_lo, ymax=BN_hi)) +
  #  geom_smooth(method='lm')  +
  ylab('effect of adding N on soil C')


# MAP_VAR 
p <- ggplot(soilC_pred, aes(x=MAP_VAR, y=BN_mn))
p +
  geom_pointrange(aes(ymin=BN_lo, ymax=BN_hi)) +
  geom_smooth(method='lm')

# MAT
p <- ggplot(soilC_pred, aes(x=MAT, y=BN_mn))
p +
  geom_pointrange(aes(ymin=BN_lo, ymax=BN_hi)) +
  geom_smooth(method='lm', lty=3, se=F)

# TEMP_VAR 
p <- ggplot(soilC_pred, aes(x=TEMP_VAR, y=BN_mn))
p +
  geom_pointrange(aes(ymin=BN_lo, ymax=BN_hi)) +
  geom_smooth(method='lm')

# texture
p <- ggplot(soilC_pred, aes(x=texture, y=BN_mn))
p +
  geom_pointrange(aes(ymin=BN_lo, ymax=BN_hi)) +
  geom_smooth(method='lm')

# above ground
p <- ggplot(soilC_pred, aes(x=site_above_mass, y=BN_mn))
p +
  geom_pointrange(aes(ymin=BN_lo, ymax=BN_hi)) +
  geom_smooth(method='lm')

#below ground
p <- ggplot(soilC_pred, aes(x=site_below_mass, y=BN_mn))
p +
  geom_pointrange(aes(ymin=BN_lo, ymax=BN_hi)) +
  geom_smooth(method='lm')

# plots: P effect size
# MAP 
p <- ggplot(soilC_pred, aes(x=MAP, y=BP_mn))
p +
  geom_pointrange(aes(ymin=BP_lo, ymax=BP_hi)) +
  geom_smooth(method='lm')  

# MAP_VAR 
p <- ggplot(soilC_pred, aes(x=MAP_VAR, y=BP_mn))
p +
  geom_pointrange(aes(ymin=BP_lo, ymax=BP_hi)) +
  geom_smooth(method='lm')

# MAT
p <- ggplot(soilC_pred, aes(x=MAT, y=BP_mn))
p +
  geom_pointrange(aes(ymin=BP_lo, ymax=BP_hi)) +
  geom_smooth(method='lm', lty=3, se=F)

# TEMP_VAR 
p <- ggplot(soilC_pred, aes(x=TEMP_VAR, y=BP_mn))
p +
  geom_pointrange(aes(ymin=BP_lo, ymax=BP_hi)) +
  geom_smooth(method='lm')

# plots: NP effect size
# MAP 
p <- ggplot(soilC_pred, aes(x=MAP, y=BNP_mn))
p +
  geom_pointrange(aes(ymin=BNP_lo, ymax=BNP_hi)) +
  geom_smooth(method='lm')  



## 2b. model with pH, above- and below-ground biomass at plot scale ####
# 2b = as covariates, no interactions or pooling
names(datmod)
mSoilCX4 <- map2stan(
  alist(
    soil_gCm2 ~ dlnorm(mu, sigma), 
    mu <- A + bN*nit + bP*pho + bNP*np + bpH*pH_num + bAM*live_mass_ln + bBM*lnroots, 
    A <- aS[site_id] + bmap*map_z + bmapvar*map_var_z + bmat*mat_z + btempvar*temp_var_z + btexture*texture_z, 
    aS[site_id] ~ dnorm(a, sigS), 
    c(bN, bP, bNP, bpH, bAM, bBM) ~ dnorm(0, 2), 
    c(bmap, bmapvar, bmat, btempvar, btexture) ~ dnorm(0, 2), 
    a ~ dnorm(0, 5), 
    sigS ~ dcauchy(0, 2), 
    sigma ~ dcauchy(0, 2)
  ),  data=c(as.list(datmod)), 
  #  chains=1
  cores=4, chains = 4, warmup=1000, iter=2000
)
m4 <- mSoilCX4
save.image("soil-resp-models.RData")

# check convergence
plot(m4)
pairs(m4, pars = c('bN', 'bP', 'bNP', 'bpH', 'bAM', 'bBM', 'lp__'))
layout(1)
plot(precis(m4, pars = c('bN', 'bP', 'bNP', 'bpH', 'bAM', 'bBM')))
plot(precis(m4, pars = c('bmap', 'bmapvar', 'bmat', 'btempvar', 'btexture')))
plot(precis(m4, pars = c('a', 'aS') , depth = 2))
compare(m2, m4)


## 2c. plots level vars are correlated within sites ####
# determined by global parameters
# Stan code from p 144-145 of Stan reference 2.14
file.show('soilC_responses.stan')

## the model matches matrix multiplication approach and notation of original "NutNet generic" model in JAGS
# which in turn was based on Gelman & Hill multilevel model with group-level predictors

# utilize double standardization to allow comparison of factors and continuous

datmod[,  .(nit, nit_z2 = (nit-mean(nit))/(2*sd(nit)))]
datmod[, c('map_z2', 'map_var_z2', 'mat_z2', 'temp_var_z2', 'texture_z2', 'site_pH_z2', 'site_above_mass_z2', 'site_below_mass_z2', 'pH_z2', 'above_mass_z2', 'below_mass_z2', 'nit_z2', 'pho_z2', 'np_z2') := lapply(.SD, function(x) (x-mean(x))/(2*sd(x))), .SDcols=c('MAP', 'MAP_VAR', 'MAT', 'TEMP_VAR', 'texture', 'site_pH', 'site_above_mass', 'site_below_mass', 'pH_num', 'live_mass_ln', 'lnroots', 'nit', 'pho', 'np')]
datmod

# match notation
# parameters & indexes
y <- datmod[, log(soil_gCm2)]
site <- datmod[, site_id]
# N = number observations
N <- datmod[, .N]
# J = number of sites
J <- datmod[, uniqueN(site_id)]
# site-level predictors:
# group parameter matrix including interecept (grand mean)
u <- matrix(c(rep(1, J), unlist(unique(datmod[, .(map_z2, map_var_z2, mat_z2, temp_var_z2, texture_z2)]))), nrow = J)	
unames <- c('Intercept', 'map_z2', 'map_var_z2', 'mat_z2', 'temp_var_z2', 'texture_z2')
# L = number of site-level predictors
L <- ncol(U)
# plot-level predictos:
# plot parameter matrix including intercept (site expectation)
x <- matrix(c(rep(1, N), unlist(datmod[, .(nit_z2, pho_z2, np_z2, above_mass_z2, below_mass_z2, pH_z2)])), nrow = N)
xnames <- c('site_intercept', 'nit_z2', 'pho_z2', 'np_z2', 'above_mass_z2', 'below_mass_z2', 'pH_z2')
# number of plot-level predictors
K <- ncol(x)

# fit!
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options('mc.cores')

m5 <- stan(file = 'soilC_responses.stan',
           chains = 4,
           iter = 5000,
           control = list(adapt_delta = 0.999))
saveRDS(m5, 'stan-m5.RDS')
pairs(m5, pars=c('lp__', 'beta[1,1]', 'gamma[1,1]'))
m5shiny <- launch_shinystan(m5)

# partial pooling N effects
alphas <- paste0('beta[',1:J,',1]')
plot(m5, pars = alphas)
betaNs <- paste0('beta[',1:J,',2]')
plot(m5, pars = betaNs)
betaPs <- paste0('beta[',1:J,',3]')
plot(m5, pars = betaPs)
betaNPs <- paste0('beta[',1:J,',4]')
plot(m5, pars = betaNPs)

# overall N, P, NP
mains <- paste0('gamma[', 2:L, ', 1]')
plot(m5, pars = mains)
muNs <- paste0('gamma[', 1:L, ', 2]')
plot(m5, pars = muNs)
muPs <- paste0('gamma[', 1:L, ', 3]')
plot(m5, pars = muPs)
muNPs <- paste0('gamma[', 1:L, ', 4]')
plot(m5, pars = muNPs)


## non-centered reparameterization with Cholesky decomposition ####
m5c <- stan(file = 'soilC_responses-cholesky.stan',
            chains = 4,
            iter = 10000,
            warmup = 6000,
            control = list(adapt_delta = 0.99))

save(m5c, file = 'stan-m5c.RDS')
pairs(m5c, pars=c('lp__', 'beta[1,1]', 'gamma[1,2]'))
m5cshiny <- launch_shinystan(m5c)
print(m5c, pars = 'gamma')

# partial pooling N effects
alphas <- paste0('beta[',1:J,',1]')
plot(m5c, pars = alphas)
betaNs <- paste0('beta[',1:J,',2]')
plot(m5c, pars = betaNs)
betaPs <- paste0('beta[',1:J,',3]')
plot(m5c, pars = betaPs)
betaNPs <- paste0('beta[',1:J,',4]')
plot(m5c, pars = betaNPs)

# overall N, P, NP
mains <- paste0('gamma[', 2:L, ', 1]')
plot(m5c, pars = mains)
muNs <- paste0('gamma[', 1:L, ', 2]')
plot(m5c, pars = muNs)
muPs <- paste0('gamma[', 1:L, ', 3]')
plot(m5c, pars = muPs)
muNPs <- paste0('gamma[', 1:L, ', 4]')
plot(m5c, pars = muNPs)

## non-centered reparameterization with Cholesky decomposition ####
# back to single standardized scale
file.show('soilC_responses-cholesky.stan')
m5c <- stan(file = 'soilC_responses-cholesky.stan',
            chains = 4,
            iter = 10000,
            warmup = 6000,
            control = list(adapt_delta = 0.99))

save(m5c, file = paste0('stan-m5c-', Sys.Date(),'.RDS'))
pairs(m5c, pars=c('lp__', 'gamma[1,1]', 'gamma[1,6]'))

m5cshiny <- launch_shinystan(m5c)
print(m5c, pars = 'gamma')