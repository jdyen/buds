# R code to fit sprout depths model

# set working directory
setwd("~/Dropbox/research/sprouts/")

# load packages
library(greta.fda)
library(greta)

# load data (currently removing species with NAs in traits)
source("./code/load-data.R")

# prepare data setg
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

# define functional response
fda_response <- fda_response(y ~ treat * gform +
                               maxht * treat +
                               stems * treat +
                               sla * treat +
                               (1 | year) + (1 | spp),
                             data = data_set,
                             spline_settings = list(df = 8, degree = 3),
                             priors = list(alpha_sd = 1, beta_sd = 1, sigma_sd = 1))

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
mod <- model(mu, alpha, beta,
             alpha_resprout, beta_resprout, p_resprout)

# sample from model
samples <- mcmc(mod, n_samples = 15000, warmup = 10000, thin = 3)

# summarise fitted model
mod_summary <- summary(samples)
samples_averaged <- do.call(abind, list(samples, along = 3))
samples_averaged <- apply(samples_averaged, c(1, 2), mean)
samples_averaged <- samples_averaged[, grep('mu\\[', colnames(samples_averaged), invert = TRUE)]

# extract mean fitted values
fitted_mean <- exp(matrix(mod_summary$statistics[grep('mu\\[', rownames(mod_summary$statistics)), 'Mean'], ncol = ncol(data_set$y)))

# r2 variants
r2 <- cor(c(fitted_mean), c(data_set$y)) ** 2
bayes_r2 <- var(c(fitted_mean)) /
  (var(c(fitted_mean)) + var(c(fitted_mean) - c(data_set$y)))
int_term <- c(data_set$y) / c(fitted_mean)
int_term <- ifelse(int_term == 0, 1, int_term)
int_term2 <- c(data_set$y) / mean(c(data_set$y))
int_term2 <- ifelse(int_term2 == 0, 1, int_term2)
dev_full <- 2 * sum(c(data_set$y) * log(int_term) - (c(data_set$y) - c(fitted_mean)))
dev_null <- 2 * sum(c(data_set$y) * log(int_term2) - (c(data_set$y) - mean(c(data_set$y))))
deviance_r2 <- 1 - dev_full / dev_null

# pull out intercepts and slopes
alpha_samples <- samples_averaged[, grep('alpha\\[', colnames(samples_averaged))]
alpha_vals <- alpha_samples %*% as.matrix(fda_response$spline_basis)
beta_samples <- samples_averaged[, grep('beta\\[', colnames(samples_averaged))]
alpha_resprout_samples <- samples_averaged[, grep('alpha_resprout', colnames(samples_averaged))]
beta_resprout_samples <- samples_averaged[, grep('beta_resprout', colnames(samples_averaged))]
p_resprout_samples <- samples_averaged[, grep('p_resprout', colnames(samples_averaged))]

# calculate resprout mean and r2 for binary part of model
p_resprout_mean <- apply(p_resprout_samples, 2, mean)
prop_resprout <- sum(resprout$survival) / length(resprout$survival)
dev_full_sprout <- - 2 * sum(resprout$survival * log(p_resprout_mean) + (1 - resprout$survival) * log(1 - p_resprout_mean))
dev_null_sprout <- - 2 * sum(resprout$survival * log(prop_resprout) + (1 - resprout$survival) * log(1 - prop_resprout))
deviance_r2_sprout <- 1 - dev_full_sprout / dev_null_sprout

# summarise observed for comparison with fitted
y_mean <- matrix(NA, nrow = length(unique(data_set$gform)), ncol = ncol(data_set$y))
rownames(y_mean) <- unique(data_set$gform)
for (i in seq_along(unique(data_set$gform))) {
  y_mean[i, ] <- apply(data_set$y[data_set$gform == unique(data_set$gform)[i], ],
                       2, mean)
}

# extract growth form effects to plot
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

# compare growth form effects and calculate Pr(x > y)
gform_comparison_burn <- matrix(NA, nrow = 6, ncol = ncol(beta_gform_burn[[1]]))
gform_comparison_clip <- matrix(NA, nrow = 6, ncol = ncol(beta_gform_clip[[1]]))
compare_burn <- vector("list", length = 6)
compare_clip <- vector("list", length = 6)
rownames(gform_comparison_burn) <- rownames(gform_comparison_clip) <- 
  c("forb_gr", "forb_ssh", "forb_wd", "gr_ssh", "gr_wd", "ssh_wd")
names(compare_burn) <- names(compare_clip) <- rownames(gform_comparison_burn)
idx <- 1
for (i in seq_len(length(beta_gform_burn) - 1)) {
  
  for (j in (i + 1):length(beta_gform_burn)) {
    compare_burn[[idx]] <- exp(beta_gform_burn[[i]] / 5.5) - exp(beta_gform_burn[[j]] / 5.5)
    gform_comparison_burn[idx, ] <- apply(compare_burn[[idx]], 2, function(x) sum(x > 0) / length(x))
    compare_clip[[idx]] <- exp(beta_gform_clip[[i]] / 5.5) - exp(beta_gform_clip[[j]] / 5.5)
    gform_comparison_clip[idx, ] <- apply(compare_clip[[idx]], 2, function(x) sum(x > 0) / length(x))
    idx <- idx + 1
  }
  
}

# calculate trait depth distributions
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
