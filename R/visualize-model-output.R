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

## extract MCMC samples
post5 <- extract(m5c2)
names(post5)
lapply(post5, dim) # rows = # MCMC samples

# site-level effects
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

# site intercept plot
pdf('figures/site-intercepts.pdf', width = 10, height = 8)
p <- ggplot(datmod, aes(x=reorder(site_id, log(soil_gCm2)), y=log(soil_gCm2)))
p + geom_boxplot() +
  geom_hline(yintercept = siteint[, mean(soil_gCm2)], lty=1) +
  geom_hline(yintercept = siteint[, PI(soil_gCm2)], lty=3) +
  geom_pointrange(data = siteint, aes(x = reorder(site_id, soil_gCm2), y = soil_gCm2, ymin = ymin, ymax = ymax), col = 'blue') +
  xlab('site')
dev.off()
system('open figures/site-intercepts.pdf')

# vars
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

# global effects
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

# effects of climate & texture on treatment effects
# prediction grid with plot-level covar set to site means
soilC_pred <- unique(datmod[, .(site_id, site_intercept = 1, nit_z2, nit, pho_z2, pho, np_z2, np, above_mass_z2 = site_above_mass_z2, below_mass_z2 = site_below_mass_z2, pH_z2 = site_pH_z2, site_code, map_z2, MAP, map_var_z2, MAP_VAR, mat_z2, MAT, temp_var_z2, TEMP_VAR, texture_z2, texture)])

betas <- data.table(as.data.frame(m5c2, pars = c('beta')))
dim(betas)
sitevarz <- soilC_pred[, xnames, with = F]
predmat <- betas * unlist(sitevarz)
predmat

soilC_pred




# partial pooling N effects
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
