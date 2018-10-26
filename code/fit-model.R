# R code to fit bud depths model

# set working directory
setwd("~/Dropbox/research/buds/")

# load packages
library(greta.fda)

# load data@
source("./code/load-data.R")

# prepare data set
data_info$GFORM_COMBINED <- data_info$GFORM
data_info$GFORM_COMBINED <- gsub('tr', 'wd', data_info$GFORM_COMBINED)
data_info$GFORM_COMBINED <- gsub('^sh', 'wd', data_info$GFORM_COMBINED)
data_set <- list(y = data_hist,
                 treat = as.factor(data_info$TREATMENT),
                 gform = data_info$GFORM_COMBINED,
                 maxht = data_info$MAXHTstd,
                 stems = data_info$STEMCATSstd,
                 sla = data_info$SLAstd,
                 year = data_info$YEAR,
                 spp = data_info$SPP)

### look at scatter-matrix of Pr(resprout) and (effect_on_depth)??
#  +,+, +,-, etc. (2 x 2 contingency)

# plots:
#    num buds / mm on y-axis
#    y-axis limits equal across all growth forms

# plot "realisations" (I(p_resprout) * distribution)
# OR hurdle model (need to tweak fda code)
# OR conditional on sprouts

## ADD HURDLE COMPONENT (bernoulli, p ~ traits)
## plot these as a subset/inset of the trait-effects plots??

# plot effects (fitted slopes) and predicted curves | trait values

# fix units on plots

# predicted curves: 4 panels per gform (sequence of observed trait values)
#    ranges don't overlap for any traits, worst for MAXHT
#    - take 4 values in range of observed traits with >= 3 growth forms at each end
#   (could just note gforms on plot)

# put panel within plot(d) of each trait effect, with trees if they go too far beyond
#     range

# define functional response
fda_response <- fda_response(y ~ treat * gform +
                               maxht * treat +
                               stems * treat +
                               sla * treat +
                               (1 | year) + (1 | spp),
                             data = data_set,
                             spline_settings = list(df = 8, degree = 3),
                             priors = list(alpha_sd = 1, beta_sd = 1, sigma_max = 1))

# prepare data for resprout model
resprout$GFORM_COMBINED <- resprout$GFORM
resprout$GFORM_COMBINED <- gsub('tr', 'wd', resprout$GFORM_COMBINED)
resprout$GFORM_COMBINED <- gsub('^sh', 'wd', resprout$GFORM_COMBINED)
resprout$YR_CODE <- as.integer(as.factor(resprout$yr))
resprout$SPP_CODE <- as.integer(as.factor(resprout$spp))
resprout_traits <- resprout[, 13:15]
resprout_traits <- cbind(resprout_traits, model.matrix( ~ resprout$GFORM_COMBINED)[, -1])
resprout_traits <- as.matrix(resprout_traits)

# define prior on resprouting
alpha_resprout <- normal(0, 10)
beta_resprout <- normal(0, 10, dim = ncol(resprout_traits))
sigma_spp <- normal(0, 10, truncation = c(0, Inf))
sigma_yr <- normal(0, 10, truncation = c(0, Inf))
gamma_spp <- normal(0, sigma_spp, dim = max(resprout$SPP_CODE))
gamma_yr <- normal(0, sigma_yr, dim = max(resprout$YR_CODE))
mu_resprout <- alpha_resprout + resprout_traits %*% beta_resprout +
  gamma_spp[resprout$SPP_CODE] + gamma_yr[resprout$YR_CODE]
p_resprout <- ilogit(mu_resprout)

# extract fda_response components
mu <- fda_response$mu
alpha <- fda_response$alpha
beta <- fda_response$beta
gamma <- fda_response$gamma
sigma_main <- fda_response$sigma_main
sigma_bins <- fda_response$sigma_bins
sigma_gamma <- fda_response$sigma_gamma

# set likelihoods
distribution(data_set$y) <- poisson(exp(fda_response$mu))
distribution(resprout$survival) <- binomial(size = 1, p = p_resprout)

# define greta model
mod <- model(mu,
             alpha, beta, gamma,
             sigma_main, sigma_bins, sigma_gamma,
             alpha_resprout, beta_resprout,
             gamma_spp, gamma_yr,
             sigma_spp, sigma_yr)

# sample from model
samples <- mcmc(mod, n_samples = 15000, warmup = 5000)

# summarise fitted
mod_summary <- apply(samples[[1]], 2, mean)
fitted_mean <- exp(matrix(mod_summary[grep('mu', names(mod_summary))], ncol = ncol(data_set$y)))
r2 <- cor(c(fitted_mean), c(data_set$y)) ** 2

alpha_samples <- samples[[1]][, grep('alpha', colnames(samples[[1]]))]
alpha_vals <- alpha_samples %*% as.matrix(fda_response$spline_basis)
beta_samples <- samples[[1]][, grep('beta', colnames(samples[[1]]))]

# summarise observed
y_mean <- matrix(NA, nrow = length(unique(data_set$gform)), ncol = ncol(data_set$y))
rownames(y_mean) <- unique(data_set$gform)
for (i in seq_along(unique(data_set$gform))) {
  y_mean[i, ] <- apply(data_set$y[data_set$gform == unique(data_set$gform)[i], ],
                       2, mean)
}

# extract growth form values to plot
beta_gform_burn <- vector('list', length = 4)
beta_gform_clip <- vector('list', length = 4)
gform_index <- c(NA, 3:5)
gform_index2 <- c(NA, 9, 11, 13)
for (i in seq_along(beta_gform_burn)) {
  clip_tmp <- beta_samples[, grep(paste0('\\[', 1, ','),
                                  colnames(beta_samples))]
  if (i == 1) {
    beta_gform_burn[[i]] <- alpha_vals
    beta_gform_clip[[i]] <- alpha_vals + clip_tmp %*% as.matrix(fda_response$spline_basis)
  } else {
    beta_tmp <- beta_samples[, grep(paste0('\\[', gform_index[i], ','),
                                    colnames(beta_samples))]
    beta_tmp2 <- beta_samples[, grep(paste0('\\[', gform_index2[i], ','),
                                    colnames(beta_samples))]
    beta_gform_burn[[i]] <- alpha_vals + beta_tmp %*% as.matrix(fda_response$spline_basis)
    beta_gform_clip[[i]] <- alpha_vals +
      clip_tmp %*% as.matrix(fda_response$spline_basis) +
      beta_tmp %*% as.matrix(fda_response$spline_basis) +
      beta_tmp2 %*% as.matrix(fda_response$spline_basis)
  }
}

# beta_gform_burn is forb, grass, chenopod, woody (four panels, one plot)
# beta_gform_clip is forb, grass, chenopod, woody (four panels, one plot)

# separate C2 plot for woody plants (1 supp plot)
btmp1 <- beta_samples[, grep(paste0('\\[', 2, ','),
                             colnames(beta_samples))]
btmp2 <- beta_samples[, grep(paste0('\\[', 5, ','),
                             colnames(beta_samples))]
btmp3 <- beta_samples[, grep(paste0('\\[', 14, ','),
                             colnames(beta_samples))]
beta_clip2 <- alpha_vals +
  btmp1 %*% as.matrix(fda_response$spline_basis) +
  btmp2 %*% as.matrix(fda_response$spline_basis) +
  btmp3 %*% as.matrix(fda_response$spline_basis)

# trait depth distributions
trait_seq <- c(-1, -0.5, 0, 0.5, 1, 1.5, 2)
trait_index <- c(6:8)   # MAXHT, STEMS, SLA
trait_index2 <- c(15, 17, 19)
beta_trait_burn <- vector('list', length = 3)
trait_effect_burn <- vector('list', length = 3)
beta_trait_clip <- vector('list', length = 3)
trait_effect_clip <- vector('list', length = 3)
for (i in seq_along(beta_trait_burn)) {
  beta_trait_burn[[i]] <- beta_samples[, grep(paste0('\\[', trait_index[i], ','),
                                              colnames(beta_samples))]
  beta_trait_clip[[i]] <- beta_samples[, grep(paste0('\\[', trait_index2[i], ','),
                                              colnames(beta_samples))]
  trait_effect_burn[[i]] <- vector('list', length = length(trait_seq))
  trait_effect_clip[[i]] <- vector('list', length = length(trait_seq))
  for (j in seq_along(trait_seq)) {
    trait_effect_burn[[i]][[j]] <- alpha_vals +
      trait_seq[j] * (beta_trait_burn[[i]] %*% as.matrix(fda_response$spline_basis))
    trait_effect_clip[[i]][[j]] <- alpha_vals +
      trait_seq[j] * ((beta_trait_burn[[i]] + beta_trait_clip[[i]]) %*% as.matrix(fda_response$spline_basis))
  }
} 
# BURN: add c(3:5) for other growth forms
# CLIP: add c(3:5) + c(9, 11, 13) for other growth forms

pdf(file = './outputs/plots/fitted_burn.pdf', width = 7, height = 7)
plot_name <- c('Forbs', 'Grasses', 'Sub-shrubs', 'Woody plants')
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
for (i in seq_along(beta_gform_burn)) {
  
  # calculate fitted summaries
  gform_mean <- apply(beta_gform_burn[[i]], 2, mean)
  gform_vlow <- apply(beta_gform_burn[[i]], 2, quantile, 0.1)
  gform_low <- apply(beta_gform_burn[[i]], 2, quantile, 0.25)
  gform_high <- apply(beta_gform_burn[[i]], 2, quantile, 0.75)
  gform_vhigh <- apply(beta_gform_burn[[i]], 2, quantile, 0.9)
  
  # plot mean fitted value
  plot(exp(gform_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Count',
       ylim = range(c(exp(gform_mean),
                      exp(gform_vlow), exp(gform_low),
                      exp(gform_vhigh), exp(gform_high))))
  
  # add shaded region to denote sampling limits
  polygon(c(-100, -100, -20, -20),
          c(-20, 20, 20, -20),
          border = NA, col = ggplot2::alpha('gray55', 0.5))
  polygon(c(-20, -20, -5, -5), c(-20, 20, 20, -20),
          border = NA, col = ggplot2::alpha('gray85', 0.5))
  polygon(c(80, 80, 150, 150), c(-20, 20, 20, -20),
          border = NA, col = ggplot2::alpha('gray85', 0.5))
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow), rev(exp(gform_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low), rev(exp(gform_high))),
          border = NA, col = 'gray55')
  lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/fitted_clip.pdf', width = 7, height = 7)
plot_name <- c('Forbs', 'Grasses', 'Sub-shrubs', 'Woody plants')
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
for (i in seq_along(beta_gform_clip)) {
  
  # calculate fitted summaries
  gform_mean <- apply(beta_gform_clip[[i]], 2, mean)
  gform_vlow <- apply(beta_gform_clip[[i]], 2, quantile, 0.1)
  gform_low <- apply(beta_gform_clip[[i]], 2, quantile, 0.25)
  gform_high <- apply(beta_gform_clip[[i]], 2, quantile, 0.75)
  gform_vhigh <- apply(beta_gform_clip[[i]], 2, quantile, 0.9)
  
  # plot mean fitted value
  plot(exp(gform_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Count',
       ylim = range(c(exp(gform_mean),
                      exp(gform_vlow), exp(gform_low),
                      exp(gform_vhigh), exp(gform_high))))
  
  # add shaded region to denote sampling limits
  polygon(c(-100, -100, -20, -20),
          c(-20, 20, 20, -20),
          border = NA, col = ggplot2::alpha('gray55', 0.5))
  polygon(c(-20, -20, -5, -5), c(-20, 20, 20, -20),
          border = NA, col = ggplot2::alpha('gray85', 0.5))
  polygon(c(80, 80, 150, 150), c(-20, 20, 20, -20),
          border = NA, col = ggplot2::alpha('gray85', 0.5))
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow), rev(exp(gform_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low), rev(exp(gform_high))),
          border = NA, col = 'gray55')
  lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/fitted_c2_treatment.pdf', width = 7, height = 7)
par(mfrow = c(1, 1), mar = c(5.1, 5.1, 1.8, 1.1))

# calculate fitted summaries
gform_mean <- apply(beta_clip2, 2, mean)
gform_vlow <- apply(beta_clip2, 2, quantile, 0.1)
gform_low <- apply(beta_clip2, 2, quantile, 0.25)
gform_high <- apply(beta_clip2, 2, quantile, 0.75)
gform_vhigh <- apply(beta_clip2, 2, quantile, 0.9)

# plot mean fitted value
plot(exp(gform_mean) ~ hist_tmp$mids,
     type = 'l', las = 1, bty = 'l',
     lwd = 2, col = 'gray30',
     xlab = 'Depth (mm)', ylab = 'Count',
     ylim = range(c(exp(gform_mean),
                    exp(gform_vlow), exp(gform_low),
                    exp(gform_vhigh), exp(gform_high))))

# add shaded region to denote sampling limits
polygon(c(-100, -100, -20, -20),
        c(-20, 20, 20, -20),
        border = NA, col = ggplot2::alpha('gray55', 0.5))
polygon(c(-20, -20, -5, -5), c(-20, 20, 20, -20),
        border = NA, col = ggplot2::alpha('gray85', 0.5))
polygon(c(80, 80, 150, 150), c(-20, 20, 20, -20),
        border = NA, col = ggplot2::alpha('gray85', 0.5))

# plot shaded regions for credible intervals
polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
        c(exp(gform_vlow), rev(exp(gform_vhigh))),
        border = NA, col = 'gray75')
polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
        c(exp(gform_low), rev(exp(gform_high))),
        border = NA, col = 'gray55')
lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')

dev.off()

pdf(file = './outputs/plots/trait_effect_burn_maxht.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$MAXHTstd, 'scaled:center') + trait_seq * attr(trait_data$MAXHTstd, 'scaled:scale') 
par(mfrow = c(3, 2), mar = c(5.1, 5.1, 1.8, 1.1))
for (i in seq_along(trait_effect_burn[[1]])[-1]) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[1]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.1)
  gform_low <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.25)
  gform_high <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.75)
  gform_vhigh <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.9)
  
  ylim_set <- range(c(exp(gform_mean),
                      exp(gform_vlow), exp(gform_low),
                      exp(gform_vhigh), exp(gform_high)))

  # plot mean fitted value
  plot(exp(gform_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Count',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow), rev(exp(gform_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low), rev(exp(gform_high))),
          border = NA, col = 'gray55')
  lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('Maximum height = ', round(plot_name[i], 1), ' m'), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_burn_stems.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$STEMCATSstd, 'scaled:center') + trait_seq * attr(trait_data$STEMCATSstd, 'scaled:scale') 
par(mfrow = c(3, 2), mar = c(5.1, 5.1, 1.8, 1.1))
for (i in seq_along(trait_effect_burn[[2]])[-length(trait_effect_burn[[2]])]) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[2]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.1)
  gform_low <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.25)
  gform_high <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.75)
  gform_vhigh <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.9)
  
  ylim_set <- range(c(exp(gform_mean),
                      exp(gform_vlow), exp(gform_low),
                      exp(gform_vhigh), exp(gform_high)))
  
  # plot mean fitted value
  plot(exp(gform_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Count',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow), rev(exp(gform_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low), rev(exp(gform_high))),
          border = NA, col = 'gray55')
  lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('Stem counts = ', round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_burn_sla.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$SLAstd, 'scaled:center') + trait_seq * attr(trait_data$SLAstd, 'scaled:scale') 
par(mfrow = c(3, 2), mar = c(5.1, 5.1, 1.8, 1.1))
for (i in seq_along(trait_effect_burn[[3]])[-length(trait_effect_burn[[3]])]) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[3]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.1)
  gform_low <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.25)
  gform_high <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.75)
  gform_vhigh <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.9)
  
  ylim_set <- range(c(exp(gform_mean),
                      exp(gform_vlow), exp(gform_low),
                      exp(gform_vhigh), exp(gform_high)))
  
  # plot mean fitted value
  plot(exp(gform_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Count',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow), rev(exp(gform_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low), rev(exp(gform_high))),
          border = NA, col = 'gray55')
  lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('SLA = ', round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_clip_maxht.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$MAXHTstd, 'scaled:center') + trait_seq * attr(trait_data$MAXHTstd, 'scaled:scale') 
par(mfrow = c(3, 2), mar = c(5.1, 5.1, 1.8, 1.1))
for (i in seq_along(trait_effect_clip[[1]])[-1]) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[1]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.1)
  gform_low <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.25)
  gform_high <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.75)
  gform_vhigh <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.9)
  
  ylim_set <- range(c(exp(gform_mean),
                      exp(gform_vlow), exp(gform_low),
                      exp(gform_vhigh), exp(gform_high)))
  
  # plot mean fitted value
  plot(exp(gform_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Count',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow), rev(exp(gform_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low), rev(exp(gform_high))),
          border = NA, col = 'gray55')
  lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('Maximum height = ', round(plot_name[i], 1), ' m'), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_clip_stems.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$STEMCATSstd, 'scaled:center') + trait_seq * attr(trait_data$STEMCATSstd, 'scaled:scale') 
par(mfrow = c(3, 2), mar = c(5.1, 5.1, 1.8, 1.1))
for (i in seq_along(trait_effect_clip[[2]])[-length(trait_effect_clip[[2]])]) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[2]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.1)
  gform_low <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.25)
  gform_high <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.75)
  gform_vhigh <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.9)
  
  ylim_set <- range(c(exp(gform_mean),
                      exp(gform_vlow), exp(gform_low),
                      exp(gform_vhigh), exp(gform_high)))
  
  # plot mean fitted value
  plot(exp(gform_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Count',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow), rev(exp(gform_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low), rev(exp(gform_high))),
          border = NA, col = 'gray55')
  lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('Stem counts = ', round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_clip_sla.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$SLAstd, 'scaled:center') + trait_seq * attr(trait_data$SLAstd, 'scaled:scale') 
par(mfrow = c(3, 2), mar = c(5.1, 5.1, 1.8, 1.1))
for (i in seq_along(trait_effect_clip[[3]])[-length(trait_effect_clip[[3]])]) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[3]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.1)
  gform_low <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.25)
  gform_high <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.75)
  gform_vhigh <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.9)
  
  ylim_set <- range(c(exp(gform_mean),
                      exp(gform_vlow), exp(gform_low),
                      exp(gform_vhigh), exp(gform_high)))
  
  # plot mean fitted value
  plot(exp(gform_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Count',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow), rev(exp(gform_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low), rev(exp(gform_high))),
          border = NA, col = 'gray55')
  lines(exp(gform_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('SLA = ', round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()
