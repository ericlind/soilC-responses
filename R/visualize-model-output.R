## process, validate, report, visualize model
library(stan)
library(shinystan)
library(loo)
library(data.table)
library(rethinking)
library(ggplot2)

## load previously fit model (see R/fit-model.R) ####
list.files('model/', pattern = 'RDS')
load('model/stan-m5c2-2017-04-14.RDS')

## call original data and parameter names etc ####
source('R/prepare-data.R')

# interactive validation & estimation
m5c2_shiny <- launch_shinystan(m5c2)

# summary table by parameter
print(m5c2, pars = 'gamma')

## extract MCMC samples ####
post5 <- extract(m5c2)
names(post5)
lapply(post5, dim) # rows = # MCMC samples

## quick-look parameter estimates ####
# partial pooling N effects: coefficients
alphas <- paste0('beta[',1:J,',1]')
plot(m5c2, pars = alphas)
print(m5c2, pars = 'gamma')
betaNs <- paste0('beta[',1:J,',2]')
plot(m5c2, pars = betaNs)
betaPs <- paste0('beta[',1:J,',3]')
plot(m5c2, pars = betaPs)
betaNPs <- paste0('beta[',1:J,',4]')
plot(m5c2, pars = betaNPs)

# overall N, P, NP
mains <- paste0('gamma[', 2:L, ', 1]')
plot(m5c2, pars = mains)
muNs <- paste0('gamma[', 1:L, ', 2]')
plot(m5c2, pars = muNs)
muPs <- paste0('gamma[', 1:L, ', 3]')
plot(m5c2, pars = muPs)
muNPs <- paste0('gamma[', 1:L, ', 4]')
plot(m5c2, pars = muNPs)

## site-level effects ####
betadat <- melt(data.table(as.data.frame(m5c2, pars = c('beta'))))
site_lut <- unique(datmod[, .(site_code, site_id)])
setkey(site_lut, site_id)

betadat[, c("site_id", "var_id") := tstrsplit(variable, ",")]
betadat[, site_id := as.numeric(gsub("beta\\[", "", site_id))]
betadat[, var_id := as.numeric(gsub("\\]", "", var_id))]        
setkey(betadat, site_id)
betadat <- site_lut[betadat]
var_lut <- data.table(var_id = 1:K,
                      predictor = xnames)
setkey(var_lut, var_id)
setkey(betadat, var_id)
betadat <- var_lut[betadat]
siteint <- betadat[predictor == 'site_intercept', 
                   .(soil_gCm2 = mean(value), 
                     ymin = PI(value)[1],
                     ymax = PI(value)[2]), site_id]

# site intercept
pdf('figures/site-intercepts.pdf', width = 10, height = 8)
p <- ggplot(datmod, aes(x=reorder(site_id, log(soil_gCm2)), y=log(soil_gCm2)))
p + geom_boxplot() +
  geom_hline(yintercept = siteint[, mean(soil_gCm2)], lty=1) +
  geom_hline(yintercept = siteint[, PI(soil_gCm2)], lty=3) +
  geom_pointrange(data = siteint, aes(x = reorder(site_id, soil_gCm2), y = soil_gCm2, ymin = ymin, ymax = ymax), col = 'blue') +
  xlab('site')
dev.off()
system('open figures/site-intercepts.pdf')

# effects of treatment
pdf('figures/site-N-effect.pdf', width = 10, height = 6)
p <- ggplot(betadat[predictor == 'nit_z2'], aes(x = reorder(site_code, value), y = value))
p + theme_bw(base_size = 15) + 
  geom_violin(fill = 'lightgray') + 
  geom_hline(yintercept = 0, lty = 2, lwd = 0.5) +
  theme(axis.text.x=element_text(angle=90,  hjust=1)) + 
  ylab('standardized effect size') +
  xlab('site') + 
  annotate(2, 0.8, geom = 'text', label = '+N', size = 8)
dev.off()
system('open figures/site-N-effect.pdf')

pdf('figures/site-P-effect.pdf', width = 10, height = 6)
p <- ggplot(betadat[predictor == 'pho_z2'], aes(x = reorder(site_code, value), y = value))
p + theme_bw(base_size = 15) + 
  geom_violin(fill = 'lightgray') + 
  geom_hline(yintercept = 0, lty = 2, lwd = 0.5) +
  theme(axis.text.x=element_text(angle=90,  hjust=1)) + 
  ylab('standardized effect size') +
  xlab('site') + 
  annotate(2, 0.7, geom = 'text', label = '+P', size = 8)
dev.off()
system('open figures/site-P-effect.pdf')

pdf('figures/site-NP-effect.pdf', width = 10, height = 6)
p <- ggplot(betadat[predictor == 'np_z2'], aes(x = reorder(site_code, value), y = value))
p + theme_bw(base_size = 15) + 
  geom_violin(fill = 'lightgray') + 
  geom_hline(yintercept = 0, lty = 2, lwd = 0.5) +
  theme(axis.text.x=element_text(angle=90,  hjust=1)) + 
  ylab('standardized effect size') +
  xlab('site') + 
  annotate(2, 0.7, geom = 'text', label = '+NP', size = 8)
dev.off()
system('open figures/site-NP-effect.pdf')

## global effects ####
gammadat <- melt(data.table(as.data.frame(m5c2, pars = c('gamma'))))
gammadat
gammadat[, c("u_id", "x_id") := tstrsplit(variable, ",")]
gammadat[, u_id := as.numeric(gsub("gamma\\[", "", u_id))]
gammadat[, x_id := as.numeric(gsub("\\]", "", x_id))]        
u_lut <- data.table(u_id = 1:L,
                    site_predictor = unames)
x_lut <- data.table(x_id = 1:K,
                      plot_predictor = xnames)
setkey(u_lut, u_id)
setkey(gammadat, u_id)
gammadat <- u_lut[gammadat]
setkey(x_lut, x_id)
setkey(gammadat, x_id)
gammadat <- x_lut[gammadat]
gammadat
unames
# global variables plot
# intercept
pdf('figures/grand-mean-soilC-m2.pdf')
p <- ggplot(gammadat[variable == 'gamma[1,1]'], aes(x = value))
p + theme_bw() + 
  geom_density() + 
  xlab(expression(paste('log(soil C ',m^-2,')')))
dev.off()

# overall effects of climate & texture on soil C
pdf('figures/overall-site-predictors.pdf', width = 8, height = 6)
p <- ggplot(gammadat[plot_predictor == 'site_intercept' & site_predictor != 'Intercept'], 
            aes(x = site_predictor,
                          y = value))
p + theme_bw() +
  geom_violin(fill = 'lightgray') + 
  geom_hline(yintercept = 0, lty = 2, lwd = 0.5) 
dev.off()  
system('open figures/overall-site-predictors.pdf')

## "interactions" or influence of global predictors on plot predictors
# 1. plot distribution of N effect per site (y) versus MAP on real scale (x)
siteclim <- unique(datmod[, .(site_id, MAP, map_z2, MAP_VAR, map_var_z2, 
                              MAT, mat_z2, TEMP_VAR, temp_var_z2, 
                              texture, texture_z2)])
setkey(siteclim, site_id)
betaclim <- siteclim[betadat, on = 'site_id']

p <- ggplot(betaclim[predictor == 'nit_z2'], aes(x = MAP, y = value))
p + theme_bw(base_size = 15) + 
  geom_violin(fill = 'lightgray', aes(group = site_code))

#2. summarize (mean, hi, lo) relationship between climate and effect
plot_interact <- function(betavar, gammavar, rep = 1000){
  predseq <- seq(from = min(siteclim[, gammavar, with = F]),
               to = max(siteclim[, gammavar, with = F]),
               length.out = rep)
  fx <- gammadat[plot_predictor == betavar & site_predictor == gammavar, value]
  mn <- sapply(1:rep, function(x) mean(predseq[x] * fx))
  hi <- sapply(1:rep, function(x) quantile(predseq[x] * fx, 0.97))
  lo <- sapply(1:rep, function(x) quantile(predseq[x] * fx, 0.03))
  return(data.table(interaction = paste0(betavar, '_x_', gammavar),
                    predseq = predseq, mn = mn, hi = hi, lo = lo))
}

map_n <- plot_interact(betavar = 'nit_z2', gammavar = 'map_z2')
# add raw scale predictor
map_n[, MAP := (predseq*(2*sd(siteclim$MAP)) + mean(siteclim$MAP))]

pdf('figures/MAP_N-effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'nit_z2'], aes(x = MAP, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = map_n, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = map_n, aes(x = MAP, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('std effect of adding N')
dev.off()
system('open figures/MAP_N-effect.pdf')


map_p <- plot_interact(betavar = 'pho_z2', gammavar = 'map_z2')
# add raw scale predictor
map_p[, MAP := (predseq*(2*sd(siteclim$MAP)) + mean(siteclim$MAP))]

pdf('figures/MAP_P-effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'pho_z2'], aes(x = MAP, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = map_p, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = map_p, aes(x = MAP, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('std effect of adding P')
dev.off()
system('open figures/MAP_P-effect.pdf')


map_np <- plot_interact(betavar = 'np_z2', gammavar = 'map_z2')
# add raw scale predictor
map_np[, MAP := (predseq*(2*sd(siteclim$MAP)) + mean(siteclim$MAP))]

pdf('figures/MAP_NP-effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'np_z2'], aes(x = MAP, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = map_np, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = map_np, aes(x = MAP, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('std effect of adding N & P')
dev.off()
system('open figures/MAP_NP-effect.pdf')

mapvar_n <- plot_interact(betavar = 'nit_z2', gammavar = 'map_var_z2')
# add raw scale predictor
mapvar_n[, MAP_VAR := (predseq*(2*sd(siteclim$MAP_VAR)) + mean(siteclim$MAP_VAR))]

pdf('figures/MAP_VAR_N-effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'nit_z2'], aes(x = MAP_VAR, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = mapvar_n, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = mapvar_n, aes(x = MAP_VAR, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('std effect of adding N') + 
  xlab('precipitation seasonality (MAP variance)')
dev.off()
system('open figures/MAP_VAR_N-effect.pdf')


mapvar_p <- plot_interact(betavar = 'pho_z2', gammavar = 'map_var_z2')
# add raw scale predictor
mapvar_p[, MAP_VAR := (predseq*(2*sd(siteclim$MAP_VAR)) + mean(siteclim$MAP_VAR))]

pdf('figures/MAP_VAR_P-effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'pho_z2'], aes(x = MAP_VAR, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = mapvar_p, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = mapvar_p, aes(x = MAP_VAR, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('std effect of adding P') + 
  xlab('precipitation seasonality (MAP variance)')
dev.off()
system('open figures/MAP_VAR_P-effect.pdf')

mapvar_np <- plot_interact(betavar = 'np_z2', gammavar = 'map_var_z2')
# add raw scale predictor
mapvar_np[, MAP_VAR := (predseq*(2*sd(siteclim$MAP_VAR)) + mean(siteclim$MAP_VAR))]

pdf('figures/MAP_VAR_NP-effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'np_z2'], aes(x = MAP_VAR, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = mapvar_np, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = mapvar_np, aes(x = MAP_VAR, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('std effect of adding N & P') + 
  xlab('precipitation seasonality (MAP variance)')
dev.off()
system('open figures/MAP_VAR_NP-effect.pdf')


# prediction grid with plot-level covar set to site means
soilC_pred <- unique(datmod[, .(site_id, site_intercept = 1, nit_z2, nit, pho_z2, pho, np_z2, np, above_mass_z2 = site_above_mass_z2, below_mass_z2 = site_below_mass_z2, pH_z2 = site_pH_z2, site_code, map_z2, MAP, map_var_z2, MAP_VAR, mat_z2, MAT, temp_var_z2, TEMP_VAR, texture_z2, texture)])
setkey(soilC_pred, site_id, nit, pho, np)
soilC_pred

## prediction on original scale - influence of treatment on soil C ####
# first, create muhat = beta[site] * xbar (x4000)
# then draw from muhat with error sigma
betas <- extract(m5c2, pars = c('beta'))$beta
dim(betas)
head(betas[, 1, ]) # 4000 reps for site 1
xbar <- soilC_pred[site_id == 1, xnames, with = F]
colSums(t(xbar) * betas[1, 1, ])

sigmas <- extract(m5c2, pars = c('sigma'))$sigma

predict_soilC <- function(i = site_id){
  xbar <- soilC_pred[site_id == i, xnames, with = F]
  muhat <- t(sapply(1:nrow(betas[, i, ]), 
                    function(x) colSums(t(xbar) * betas[x, i, ])))
  out <-  t(sapply(1:nrow(muhat), 
                   function(x) rnorm(n = ncol(muhat), 
                                     mean = muhat[x, ], 
                                     sd = sigmas[x])))
  return(out)
}

predlist <- lapply(1:J, predict_soilC)
predlist <- lapply(predlist, as.data.table)
predDT <- rbindlist(predlist, idcol = 'site_id')
setnames(predDT, c('site_id', 'Ctl', 'P', 'N', 'NP'))
predDT
predmelt <- melt(predDT, id.vars = 'site_id',
                variable.name = 'trt', 
                value.name = 'log_soilC_gm2')
predmelt <- site_lut[predmelt]
pdf('figures/trt-posterior-prob-soilC-by-site.pdf', width = 12, 
    height = 8)
p <- ggplot(predmelt, aes(x = factor(trt),
                        y = log_soilC_gm2))
p + geom_violin() + facet_wrap(~ site_code) +
  ylab(expression(paste('posterior prob distribution of soil C', m^-2)))
dev.off()
system('open figures/trt-posterior-prob-soilC-by-site.pdf')

## reformulate as single plot with %C by site, col = trt
pdf('figures/plot-trt-estimates-by-site.pdf', width = 12, height = 9)
p <- ggplot(predmelt, aes(x = reorder(site_code, log_soilC_gm2),
                          y = log_soilC_gm2,
                          fill = factor(trt)))
p + theme_bw(base_size = 12) + 
  geom_violin() + 
  geom_vline(xintercept = c(1:26 + 0.5), col = 'darkgray') + 
  ylab(expression(paste('posterior prob distribution of soil C', m^-2))) + 
  xlab('site_code') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
system('open figures/plot-trt-estimates-by-site.pdf')

## prediction on original scale - influence of climate/soil on soil C ####
# need prediction matrix, set all other values to mean (0 in scaled space)

site_predict <- function(gammavar, rep = 1000){
  predseq <- seq(from = min(siteclim[, gammavar, with = F]),
                 to = max(siteclim[, gammavar, with = F]),
                 length.out = rep)
  fx <- gammadat[plot_predictor == 'site_intercept' & site_predictor == gammavar, value]
  A <- gammadat[variable == 'gamma[1,1]', value]
  mn <- sapply(1:rep, function(x) mean((A + predseq[x] * fx)))
  hi <- sapply(1:rep, function(x) quantile((A + predseq[x] * fx), 0.97))
  lo <- sapply(1:rep, function(x) quantile((A + predseq[x] * fx), 0.03))
  return(data.table(site_predictor = gammavar,
                    predseq = predseq, mn = mn, hi = hi, lo = lo))
}
site_predict(gammavar = 'map_z2')
clims <- unique(gammadat$site_predictor)[2:6]
clim_list <- lapply(clims, function(x) site_predict(gammavar = x))
climfx <- rbindlist(clim_list)
climfx

# add raw scale predictor
map_fx <- climfx[site_predictor == 'map_z2', MAP := (predseq*(2*sd(siteclim$MAP)) + mean(siteclim$MAP))]

pdf('figures/MAP_effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'site_intercept'], aes(x = MAP, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = map_fx, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = map_fx, aes(x = MAP, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('posterior probability of log soil C')
dev.off()
system('open figures/MAP_effect.pdf')

mapvar_fx <- climfx[site_predictor == 'map_var_z2', MAP_VAR := (predseq*(2*sd(siteclim$MAP_VAR)) + mean(siteclim$MAP_VAR))]

pdf('figures/MAP_VAR_effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'site_intercept'], aes(x = MAP_VAR, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = mapvar_fx, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = mapvar_fx, aes(x = MAP_VAR, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('posterior probability of log soil C') + 
  xlab('Seasonality of Precip (MAP_Var)')
dev.off()
system('open figures/MAP_VAR_effect.pdf')

temp_var_fx <- climfx[site_predictor == 'temp_var_z2',] TEMP_VAR := (predseq*(2*sd(siteclim$TEMP_VAR)) + mean(siteclim$TEMP_VAR))

pdf('figures/TEMP_VAR_effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'site_intercept'], aes(x = TEMP_VAR, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = temp_var_fx, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = temp_var_fx, aes(x = TEMP_VAR, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('posterior probability of log soil C') + 
  xlab('SD of Annual Temperature (TEMP_Var)')
dev.off()
system('open figures/TEMP_VAR_effect.pdf')

temp_fx <- climfx[site_predictor == 'mat_z2', MAT := (predseq*(2*sd(siteclim$MAT)) + mean(siteclim$MAT))]

pdf('figures/MAT_effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'site_intercept'], aes(x = MAT, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = temp_fx, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = temp_fx, aes(x = MAT, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('posterior probability of log soil C') + 
  xlab('MAT (C)')
dev.off()
system('open figures/MAT_effect.pdf')

texture_fx <- climfx[site_predictor == 'texture_z2', texture := (predseq*(2*sd(siteclim$texture)) + mean(siteclim$texture))]

pdf('figures/texture_effect.pdf', width = 8, height = 6)
p <- ggplot(betaclim[predictor == 'site_intercept'], aes(x = texture, y = value))
p + theme_bw(base_size = 15) + 
  geom_ribbon(data = texture_fx, aes(y = mn, ymax = hi, ymin = lo), fill = 'lightgray') +
  geom_line(data = texture_fx, aes(x = texture, y = mn)) +
  geom_violin(fill = 'lightgray', aes(group = site_code)) + 
  ylab('posterior probability of log soil C') + 
  xlab('soil texture (% silt + % clay)')
dev.off()
system('open figures/texture_effect.pdf')

