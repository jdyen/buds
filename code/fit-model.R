# R code to fit bud depths model

# set working directory
setwd("~/Dropbox/research/buds/")

# load packages
library(greta.fda)
library(greta)

# load data@
source("./code/load-data.R")

# prepare data set
data_info$GFORM_COMBINED <- data_info$GFORM
data_info$GFORM_COMBINED <- gsub('tr', 'wd', data_info$GFORM_COMBINED)
data_info$GFORM_COMBINED <- gsub('^sh', 'wd', data_info$GFORM_COMBINED)
data_set <- list(y = data_hist,
                 treat = factor(data_info$TREATMENT, levels = c("c", "b")),
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
resprout_traits <- model.matrix( ~ resprout$treat * resprout$GFORM_COMBINED +
                                   resprout$treat * resprout$MAXHTstd +
                                   resprout$treat * resprout$STEMCATSstd +
                                   resprout$treat * resprout$SLAstd)[, -1]
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
             alpha, beta,
             alpha_resprout, beta_resprout,
             p_resprout)

# sample from model
samples <- mcmc(mod, n_samples = 15000, warmup = 10000, thin = 3)

# summarise fitted
mod_summary <- summary(samples)
samples_averaged <- do.call(abind, list(samples, along = 3))
samples_averaged <- apply(samples_averaged, c(1, 2), mean)
samples_averaged <- samples_averaged[, grep('mu\\[', colnames(samples_averaged), invert = TRUE)]

# summarise fitted
fitted_mean <- exp(matrix(mod_summary$statistics[grep('mu\\[', rownames(mod_summary$statistics)), 'Mean'], ncol = ncol(data_set$y)))
r2 <- cor(c(fitted_mean), c(data_set$y)) ** 2

# pull out means/slopes/etc.
alpha_samples <- samples_averaged[, grep('alpha\\[', colnames(samples_averaged))]
alpha_vals <- alpha_samples %*% as.matrix(fda_response$spline_basis)
beta_samples <- samples_averaged[, grep('beta\\[', colnames(samples_averaged))]
alpha_resprout_samples <- samples_averaged[, grep('alpha_resprout', colnames(samples_averaged))]
beta_resprout_samples <- samples_averaged[, grep('beta_resprout', colnames(samples_averaged))]
p_resprout_samples <- samples_averaged[, grep('p_resprout', colnames(samples_averaged))]

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
gform_index <- c(NA, 2:4)  # grass, subshrub, woody
gform_index2 <- c(NA, 8:10) # burn treatment for grass, subshrub, woody
for (i in seq_along(beta_gform_burn)) {
  burn_tmp <- beta_samples[, grep(paste0('\\[', 1, ','),
                                  colnames(beta_samples))]
  if (i == 1) {
    beta_gform_clip[[i]] <- alpha_vals
    beta_gform_burn[[i]] <- alpha_vals + burn_tmp %*% as.matrix(fda_response$spline_basis)
  } else {
    beta_tmp <- beta_samples[, grep(paste0('\\[', gform_index[i], ','),
                                    colnames(beta_samples))]
    beta_tmp2 <- beta_samples[, grep(paste0('\\[', gform_index2[i], ','),
                                    colnames(beta_samples))]
    beta_gform_clip[[i]] <- alpha_vals + beta_tmp %*% as.matrix(fda_response$spline_basis)
    beta_gform_burn[[i]] <- alpha_vals +
      burn_tmp %*% as.matrix(fda_response$spline_basis) +
      beta_tmp %*% as.matrix(fda_response$spline_basis) +
      beta_tmp2 %*% as.matrix(fda_response$spline_basis)
  }
}
# beta_gform_burn is forb, grass, chenopod, woody (four panels, one plot)
# beta_gform_clip is forb, grass, chenopod, woody (four panels, one plot)

# trait depth distributions
trait_seq <- list(c(-0.5, -0.25, 0, 0.25),    # HT
                  c(-1.0, 0.0, 0.5, 1.0),    # STEMS
                  c(-1.0, -0.5, 0.5, 1.0))    # SLA
trait_index <- c(5:7)   # MAXHT, STEMS, SLA
trait_index2 <- c(11:14)  # burn treatment effect for MAXHT, STEMS, SLA
beta_trait_burn <- vector('list', length = 3)
trait_effect_burn <- vector('list', length = 3)
beta_trait_clip <- vector('list', length = 3)
trait_effect_clip <- vector('list', length = 3)
for (i in seq_along(beta_trait_burn)) {
  beta_trait_clip[[i]] <- beta_samples[, grep(paste0('\\[', trait_index[i], ','),
                                              colnames(beta_samples))]
  beta_trait_burn[[i]] <- beta_samples[, grep(paste0('\\[', trait_index2[i], ','),
                                              colnames(beta_samples))]
  trait_effect_burn[[i]] <- vector('list', length = length(trait_seq[[i]]))
  trait_effect_clip[[i]] <- vector('list', length = length(trait_seq[[i]]))
  for (j in seq_along(trait_seq[[i]])) {
    trait_effect_clip[[i]][[j]] <- alpha_vals +
      trait_seq[[i]][j] * (beta_trait_clip[[i]] %*% as.matrix(fda_response$spline_basis))
    trait_effect_burn[[i]][[j]] <- alpha_vals +
      trait_seq[[i]][j] * ((beta_trait_clip[[i]] + beta_trait_burn[[i]]) %*% as.matrix(fda_response$spline_basis))
  }
} 

# BURN: add c(2:4) for other growth forms
# CLIP: add c(2:4) + c(8:10) for other growth forms

pdf(file = './outputs/plots/fitted_burn.pdf', width = 7, height = 7)
plot_name <- c('Forbs', 'Grasses', 'Sub-shrubs', 'Woody plants')
set.seed(123)
# plot_name <- letters[1:4]
par(mfrow = c(2, 2), mar = c(4.8, 5.1, 2.5, 1.1))
xaxs_lims <- c(-15, 95)
yaxs_lims <- c(0, 8)
n_plot <- 200
for (i in seq_along(beta_gform_burn)) {
  
  # calculate fitted summaries
  gform_mean <- apply(beta_gform_burn[[i]], 2, mean)
  gform_vlow <- apply(beta_gform_burn[[i]], 2, quantile, 0.025)
  gform_low <- apply(beta_gform_burn[[i]], 2, quantile, 0.1)
  gform_high <- apply(beta_gform_burn[[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(beta_gform_burn[[i]], 2, quantile, 0.975)
  
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = 'n', las = 1, bty = 'l',
       xlab = 'Depth (mm)', ylab = 'Buds / mm',
       ylim = yaxs_lims, xlim = xaxs_lims)
  
  sample_plot <- sample(seq_len(nrow(beta_gform_burn[[i]])), size = n_plot, replace = FALSE)
  plot_all <- rep(0, length(gform_mean))
  
  for (j in seq_len(n_plot)) {
    
    # calculate fitted summaries
    plot_tmp <- beta_gform_burn[[i]][sample_plot[j], ]
    preds_tmp <- rep(0, 13)
    if (i > 1)
      preds_tmp[i] <- 1
    p_sprouted <- alpha_resprout_samples[sample_plot[j]] +
      beta_resprout_samples[sample_plot[j], ] %*% preds_tmp
    p_sprouted <- plogis(p_sprouted)
    
    ind_var <- ifelse(runif(1) < p_sprouted, 1, 0)
    plot_tmp <- c(ind_var) * plot_tmp
    # yplot <- rpois(length(plot_tmp), lambda = c(exp(plot_tmp) / 5.5))
    yplot <- c(exp(plot_tmp) / 5.5)
    
    # plot mean fitted value
    lines(yplot ~ hist_tmp$mids,
         lwd = 1, col = ggplot2::alpha('gray50', 0.4))
    
    plot_all <- plot_all + plot_tmp
    
  }
  
  # add shaded region to denote sampling limits
  polygon(c(-20, -20, -5, -5),
          2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          border = NA, col = ggplot2::alpha('gray85', 0.5))
  polygon(c(80, 80, 150, 150),
          2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          border = NA, col = ggplot2::alpha('gray85', 0.5))
  
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'black')
  lines(c(exp(c(plot_all / n_plot)) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'black', lty = 2)
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0.5, adj = -0.02, cex = 1.5)
  
}
dev.off()

pdf(file = './outputs/plots/fitted_clip.pdf', width = 7, height = 7)
plot_name <- c('Forbs', 'Grasses', 'Sub-shrubs', 'Woody plants')
set.seed(123)
# plot_name <- letters[1:4]
par(mfrow = c(2, 2), mar = c(4.8, 5.1, 2.5, 1.1))
xaxs_lims <- c(-15, 95)
yaxs_lims <- c(0, 8)
n_plot <- 200
for (i in seq_along(beta_gform_clip)) {
  
  # calculate fitted summaries
  gform_mean <- apply(beta_gform_clip[[i]], 2, mean)
  gform_vlow <- apply(beta_gform_clip[[i]], 2, quantile, 0.025)
  gform_low <- apply(beta_gform_clip[[i]], 2, quantile, 0.1)
  gform_high <- apply(beta_gform_clip[[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(beta_gform_clip[[i]], 2, quantile, 0.975)
  
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = 'n', las = 1, bty = 'l',
       xlab = 'Depth (mm)', ylab = 'Buds / mm',
       ylim = yaxs_lims, xlim = xaxs_lims)
  
  sample_plot <- sample(seq_len(nrow(beta_gform_clip[[i]])), size = n_plot, replace = FALSE)
  plot_all <- rep(0, length(gform_mean))
  
  for (j in seq_len(n_plot)) {
    
    # calculate fitted summaries
    plot_tmp <- beta_gform_clip[[i]][sample_plot[j], ]
    preds_tmp <- rep(0, 13)
    preds_tmp[1] <- 1
    if (i > 1) {
      preds_tmp[i] <- 1
      preds_tmp[i + 6] <- 1
    }
    
    p_sprouted <- alpha_resprout_samples[sample_plot[j]] +
      beta_resprout_samples[sample_plot[j], ] %*% preds_tmp
    p_sprouted <- plogis(p_sprouted)
    
    ind_var <- ifelse(runif(1) < p_sprouted, 1, 0)
    plot_tmp <- c(ind_var) * plot_tmp
    
    # plot mean fitted value
    lines(c(exp(plot_tmp) / 5.5) ~ hist_tmp$mids,
          lwd = 1, col = ggplot2::alpha('gray50', 0.4))
    
    plot_all <- plot_all + plot_tmp
    
  }
  
  # add shaded region to denote sampling limits
  polygon(c(-20, -20, -5, -5),
          2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          border = NA, col = ggplot2::alpha('gray85', 0.5))
  polygon(c(80, 80, 150, 150),
          2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          border = NA, col = ggplot2::alpha('gray85', 0.5))
  
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'black')
  lines(c(exp(c(plot_all / n_plot)) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'black', lty = 2)
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0.5, adj = -0.02, cex = 1.5)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_burn_maxht.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$MAXHTstd, 'scaled:center') + trait_seq[[1]] * attr(trait_data$MAXHTstd, 'scaled:scale') 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 1.2)
for (i in seq_along(trait_effect_burn[[1]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[1]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Buds / mm',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = 'gray55')
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('Maximum height = ', round(plot_name[i], 0), ' cm'), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_burn_stems.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$STEMCATSstd, 'scaled:center') + trait_seq[[2]] * attr(trait_data$STEMCATSstd, 'scaled:scale') 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 1.2)
for (i in seq_along(trait_effect_burn[[2]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[2]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Buds / mm',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = 'gray55')
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('Stem counts = ', round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_burn_sla.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$SLAstd, 'scaled:center') + trait_seq[[3]] * attr(trait_data$SLAstd, 'scaled:scale') 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 1.2)
for (i in seq_along(trait_effect_burn[[3]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[3]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Buds / mm',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = 'gray55')
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('SLA = ', round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_clip_maxht.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$MAXHTstd, 'scaled:center') + trait_seq[[1]] * attr(trait_data$MAXHTstd, 'scaled:scale') 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 1.2)
for (i in seq_along(trait_effect_clip[[1]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[1]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Buds / mm',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = 'gray55')
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('Maximum height = ', round(plot_name[i], 0), ' cm'), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_clip_stems.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$STEMCATSstd, 'scaled:center') + trait_seq[[2]] * attr(trait_data$STEMCATSstd, 'scaled:scale') 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 1.2)
for (i in seq_along(trait_effect_clip[[2]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[2]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Buds / mm',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = 'gray55')
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('Stem counts = ', round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

pdf(file = './outputs/plots/trait_effect_clip_sla.pdf', width = 7, height = 7)
plot_name <- attr(trait_data$SLAstd, 'scaled:center') + trait_seq[[3]] * attr(trait_data$SLAstd, 'scaled:scale') 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 1.2)
for (i in seq_along(trait_effect_clip[[3]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[3]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Buds / mm',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = 'gray55')
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  
  # add some labels
  mtext(paste0('SLA = ', round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()


# trait effects
pdf(file = './outputs/plots/regression_coefficients.pdf', width = 7, height = 7)
par(mfrow = c(3, 2), mar = c(4.8, 5.1, 2.5, 1.1))
ylim_set <- c(-0.1, 4.5)
# plot_label <- letters[1:6]
plot_label <- rep(c("Max. height", "Stem counts", "SLA"), each = 2)
for (i in seq_len(3)) {
  
  burn_parameters <- ((beta_trait_clip[[i]] + beta_trait_burn[[i]]) %*% as.matrix(fda_response$spline_basis))
  clip_parameters <- (beta_trait_clip[[i]] %*% as.matrix(fda_response$spline_basis))
  
  # fitted params
  plot_mean <- apply(clip_parameters, 2, mean)
  plot_vlow <- apply(clip_parameters, 2, quantile, 0.025)
  plot_low <- apply(clip_parameters, 2, quantile, 0.1)
  plot_high <- apply(clip_parameters, 2, quantile, 0.9)
  plot_vhigh <- apply(clip_parameters, 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(exp(plot_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Effect',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(plot_vlow), rev(exp(plot_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(plot_low), rev(exp(plot_high))),
          border = NA, col = 'gray55')
  
  lines(exp(plot_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  lines(c(min(hist_tmp$mids) - 10, max(hist_tmp$mids) + 10),
        c(1, 1), lty = 2)
  
  #
  mtext(plot_label[(2 * i) - 1], side = 3, line = 0.5, adj = -0.02, cex = 1.5)
  
  # calculate fitted summaries
  plot_mean <- apply(burn_parameters, 2, mean)
  plot_vlow <- apply(burn_parameters, 2, quantile, 0.025)
  plot_low <- apply(burn_parameters, 2, quantile, 0.1)
  plot_high <- apply(burn_parameters, 2, quantile, 0.9)
  plot_vhigh <- apply(burn_parameters, 2, quantile, 0.975)

  # plot mean fitted value
  plot(exp(plot_mean) ~ hist_tmp$mids,
       type = 'l', las = 1, bty = 'l',
       lwd = 2, col = 'gray30',
       xlab = 'Depth (mm)', ylab = 'Effect',
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(plot_vlow), rev(exp(plot_vhigh))),
          border = NA, col = 'gray75')
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(plot_low), rev(exp(plot_high))),
          border = NA, col = 'gray55')
  
  lines(exp(plot_mean) ~ hist_tmp$mids, lwd = 2, col = 'gray30')
  lines(c(min(hist_tmp$mids) - 10, max(hist_tmp$mids) + 10),
        c(1, 1), lty = 2)
  
  #
  mtext(plot_label[(2 * i)], side = 3, line = 0.5, adj = -0.02, cex = 1.5)
  
}
dev.off()

# plot resprouting effects
# combos to plot: clip x gform, burn x gform (= 8 combos inc. forbs)
# slopes to slope: clip x traits, burn x traits (= 6 combos)
colnames(beta_resprout_samples) <- c("Clip", "Grasses", "Sub-shrubs", "Woody plants",
                                     "Max. height", "Stem counts", "SLA",
                                     "Clip:Grasses", "Clip:Sub-shrubs", "Clip:Woody plants",
                                     "Clip:Max. height", "Clip:Stem counts", "Clip:SLA")

pdf(file = './outputs/plots/resprout_trait_effects.pdf', width = 7, height = 6)

par(mfrow = c(1, 2))
burn_int <- cbind(alpha_resprout_samples,
                  sweep(beta_resprout_samples[, 2:4], 1, alpha_resprout_samples, "+"))
clip_int <- cbind(beta_resprout_samples[, 1],
                  sweep(beta_resprout_samples[, 8:10], 1, beta_resprout_samples[, 1], "+"))
int_combos <- plogis(cbind(burn_int + clip_int,
                           burn_int))
colnames(int_combos) <- c("Clip:forbs", "Clip:grasses", "Clip:sub-shrubs", "Clip:woody plants",
                          "Burn:forbs", "Burn:grasses", "Burn:sub-shrubs", "Burn:woody plants")
burn_slope <- beta_resprout_samples[, 5:7]
clip_slope <- burn_slope + beta_resprout_samples[, 11:13]
slope_combos <- cbind(clip_slope, burn_slope)
colnames(slope_combos) <- c("Clip:max. height", "Clip:stem counts", "Clip:SLA",
                          "Burn:max. height", "Burn:stem counts", "Burn:SLA")

int_mean <- apply(int_combos, 2, mean)
int_quant <- apply(int_combos, 2, quantile, p = c(0.025, 0.1, 0.9, 0.975))
par(mar = c(3.1, 8.1, 2.3, 0.2))
plot(seq_len(ncol(int_quant)), type = "n",
     las = 1, bty = "l",
     xlab = "", ylab = "",
     xlim = c(0.0, 1.0),
     xaxt = "n", yaxt = "n")
for (i in seq_along(int_mean)) {
  points(int_mean[i], length(int_mean) + 1 - i, pch = 16, cex = 1)
  lines(c(int_quant[1, i], int_quant[4, i]),
        c(length(int_mean) + 1 - i, length(int_mean) + 1 - i),
        lwd = 1)
  lines(c(int_quant[2, i], int_quant[3, i]),
        c(length(int_mean) + 1 - i, length(int_mean) + 1 - i),
        lwd = 3)
}
lines(c(0.5, 0.5), c(0, length(int_mean) + 1), lty = 2)
axis(1, at = c(0.0, 0.25, 0.5, 0.75, 1.0))
axis(2, at = seq_len(ncol(int_quant)),
     labels = rev(colnames(int_combos)), las = 1)
mtext("a", side = 3, line = 0.5, adj = -0.02, cex = 1.5)

slope_mean <- apply(slope_combos, 2, mean)
slope_quant <- apply(slope_combos, 2, quantile, p = c(0.025, 0.1, 0.9, 0.975))
plot(seq_len(ncol(slope_quant)), type = "n",
     las = 1, bty = "l",
     xlab = "", ylab = "",
     xlim = c(-1.5, 2.0),
     xaxt = "n", yaxt = "n")
for (i in seq_along(slope_mean)) {
  points(slope_mean[i], length(slope_mean) + 1 - i, pch = 16, cex = 1)
  lines(c(slope_quant[1, i], slope_quant[4, i]),
        c(length(slope_mean) + 1 - i, length(slope_mean) + 1 - i),
        lwd = 1)
  lines(c(slope_quant[2, i], slope_quant[3, i]),
        c(length(slope_mean) + 1 - i, length(slope_mean) + 1 - i),
        lwd = 3)
}
lines(c(0, 0), c(0, length(slope_mean) + 1), lty = 2)
axis(1, at = c(-1.5, 0.0, 1.5))
axis(2, at = seq_len(ncol(slope_quant)),
     labels = rev(colnames(slope_combos)), las = 1)
mtext("b", side = 3, line = 0.5, adj = -0.02, cex = 1.5)

# bayesplot::mcmc_intervals(int_combos)
# bayesplot::mcmc_intervals(slope_combos)

dev.off()
