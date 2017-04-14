## Soil C responses (Riggs-Crowther)
library(data.table)

## import and process data ####
dat <- fread('data/data_allplots.csv')
names(dat)

## eliminate saline.us for lack of complete data ####
datmod <- dat[site_code!='saline.us', .(site_code, block, plot, N, P, soil_gCm2, live_mass_ln, lnroots, pH_num, PercentSiltClay, MAP, MAP_VAR, MAT, TEMP_VAR)]
datmod[, unique(PercentSiltClay), site_code]

## add integer treatment effects with usable variable names ####
datmod[, nit:=as.integer(N)]
datmod[, pho:=as.integer(P)]
datmod[, np:=0L]
datmod[nit==1 & pho==1, np:=1L]
datmod[, .(nit, pho, np)]

## mean soil values and production above and below ground ####
datmod[, texture:=mean(PercentSiltClay), site_code]
datmod[, site_pH:=mean(pH_num), site_code]
datmod[, site_above_mass := mean(live_mass_ln), site_code]
datmod[, site_below_mass := mean(lnroots), site_code]

## scale data for model interpretability ####
# center & scale
datmod[, c('map_z', 'map_var_z', 'mat_z', 'temp_var_z', 'texture_z', 'site_pH_z', 'site_above_mass_z', 'site_below_mass_z') := lapply(.SD, scale), .SDcols=c('MAP', 'MAP_VAR', 'MAT', 'TEMP_VAR', 'texture', 'site_pH', 'site_above_mass', 'site_below_mass')]
datmod[, c('map_z', 'map_var_z', 'mat_z', 'temp_var_z', 'texture_z', 'site_pH_z', 'site_above_mass_z', 'site_below_mass_z') := lapply(.SD, as.numeric), .SDcols=c('map_z', 'map_var_z', 'mat_z', 'temp_var_z', 'texture_z', 'site_pH_z', 'site_above_mass_z', 'site_below_mass_z')]

# "Gelman double standardization"
# utilize double standardization to allow comparison of factors and continuous
datmod[, c('map_z2', 'map_var_z2', 'mat_z2', 'temp_var_z2', 'texture_z2', 'site_pH_z2', 'site_above_mass_z2', 'site_below_mass_z2', 'pH_z2', 'above_mass_z2', 'below_mass_z2', 'nit_z2', 'pho_z2', 'np_z2') := lapply(.SD, function(x) (x-mean(x))/(2*sd(x))), .SDcols=c('MAP', 'MAP_VAR', 'MAT', 'TEMP_VAR', 'texture', 'site_pH', 'site_above_mass', 'site_below_mass', 'pH_num', 'live_mass_ln', 'lnroots', 'nit', 'pho', 'np')]

## match hierarchical model notation for Stan input ####
datmod[, site_id := as.integer(factor(site_code))]

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
L <- ncol(u)
# plot-level predictos:
# plot parameter matrix including intercept (site expectation)
x <- matrix(c(rep(1, N), unlist(datmod[, .(nit_z2, pho_z2, np_z2, above_mass_z2, below_mass_z2, pH_z2)])), nrow = N)
xnames <- c('site_intercept', 'nit_z2', 'pho_z2', 'np_z2', 'above_mass_z2', 'below_mass_z2', 'pH_z2')
# number of plot-level predictors
K <- ncol(x)

