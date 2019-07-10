# create all plots for resprouting manuscript

# work this out:
plot_type <- "jpeg"
file_type <- "jpg"
unit_set <- "in"
width_set <- 7
height_set <- 7
res_set <- 300


get(plot_type)(file = paste0("outputs/plots/fitted_burn.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- c("Forbs", "Grasses", "Sub-shrubs", "Woody plants")
set.seed(123)
# plot_name <- letters[1:4]
par(mfrow = c(2, 2), mar = c(4.8, 5.1, 2.5, 1.1))
xaxs_lims <- c(-15, 95)
yaxs_lims <- c(0, 18)
n_plot <- 200
for (i in seq_along(beta_gform_burn)) {
  
  # calculate fitted summaries
  gform_mean <- apply(beta_gform_burn[[i]], 2, mean)
  gform_vlow <- apply(beta_gform_burn[[i]], 2, quantile, 0.025)
  gform_low <- apply(beta_gform_burn[[i]], 2, quantile, 0.1)
  gform_high <- apply(beta_gform_burn[[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(beta_gform_burn[[i]], 2, quantile, 0.975)
  
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = "n", las = 1, bty = "l",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
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
          lwd = 1, col = ggplot2::alpha("gray50", 0.4))
    
    plot_all <- plot_all + plot_tmp
    
  }
  
  # add shaded region to denote sampling limits
  polygon(c(-20, -20, -5, -5),
          2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  polygon(c(80, 80, 150, 150),
          2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "black")
  lines(c(exp(c(plot_all / n_plot)) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "black", lty = 2)
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0.5, adj = -0.02, cex = 1.5)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/pairwise_comparison_burn.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- c("Forbs, Grasses",
               "Forbs, Sub-shrubs",
               "Forbs, Woody plants",
               "Grasses, Sub-shrubs",
               "Grasses, Woody plants",
               "Sub-shrubs, Woody plants")
par(mfrow = c(3, 2), mar = c(4.8, 5.1, 2.5, 1.1))
xaxs_lims <- c(-15, 95)
yaxs_lims <- c(-2, 2)
for (i in seq_along(compare_burn)) {
  
  # calculate fitted summaries
  gform_mean <- apply(compare_burn[[i]], 2, mean)
  gform_vlow <- apply(compare_burn[[i]], 2, quantile, 0.025)
  gform_low <- apply(compare_burn[[i]], 2, quantile, 0.1)
  gform_high <- apply(compare_burn[[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(compare_burn[[i]], 2, quantile, 0.975)
  
  plot(gform_mean ~ hist_tmp$mids,
       type = "n", las = 1, bty = "l",
       xlab = "Depth (mm)", ylab = "Difference (Sprouts / mm)",
       ylim = yaxs_lims, xlim = xaxs_lims)

  # add shaded region to denote sampling limits
  polygon(c(-20, -20, -5, -5),
          2 * c(yaxs_lims, rev(yaxs_lims)),
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  polygon(c(80, 80, 150, 150),
          2 * c(yaxs_lims, rev(yaxs_lims)),
          border = NA, col = ggplot2::alpha("gray85", 0.5))

  # add shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(gform_vlow, rev(gform_vhigh)),
          border = NA, col = ggplot2::alpha("gray65", 0.5))
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(gform_low, rev(gform_high)),
          border = NA, col = ggplot2::alpha("gray35", 0.75))
  
  lines(gform_mean ~ hist_tmp$mids, lwd = 2, col = "black")
  lines(c(0, 0) ~ c(2 * xaxs_lims), lty = 2)
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0.5, adj = -0.02, cex = 1.1)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/pairwise_comparison_clip.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- c("Forbs, Grasses",
               "Forbs, Sub-shrubs",
               "Forbs, Woody plants",
               "Grasses, Sub-shrubs",
               "Grasses, Woody plants",
               "Sub-shrubs, Woody plants")
par(mfrow = c(3, 2), mar = c(4.8, 5.1, 2.5, 1.1))
xaxs_lims <- c(-15, 95)
yaxs_lims <- c(-2, 2)
for (i in seq_along(compare_clip)) {
  
  # calculate fitted summaries
  gform_mean <- apply(compare_clip[[i]], 2, mean)
  gform_vlow <- apply(compare_clip[[i]], 2, quantile, 0.025)
  gform_low <- apply(compare_clip[[i]], 2, quantile, 0.1)
  gform_high <- apply(compare_clip[[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(compare_clip[[i]], 2, quantile, 0.975)
  
  plot(gform_mean ~ hist_tmp$mids,
       type = "n", las = 1, bty = "l",
       xlab = "Depth (mm)", ylab = "Difference (Sprouts / mm)",
       ylim = yaxs_lims, xlim = xaxs_lims)
  
  # add shaded region to denote sampling limits
  polygon(c(-20, -20, -5, -5),
          2 * c(yaxs_lims, rev(yaxs_lims)),
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  polygon(c(80, 80, 150, 150),
          2 * c(yaxs_lims, rev(yaxs_lims)),
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  
  # add shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(gform_vlow, rev(gform_vhigh)),
          border = NA, col = ggplot2::alpha("gray65", 0.5))
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(gform_low, rev(gform_high)),
          border = NA, col = ggplot2::alpha("gray35", 0.75))
  
  lines(gform_mean ~ hist_tmp$mids, lwd = 2, col = "black")
  lines(c(0, 0) ~ c(2 * xaxs_lims), lty = 2)
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0.5, adj = -0.02, cex = 1.1)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/fitted_intercepts.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- c("Forbs", "Grasses", "Sub-shrubs", "Woody plants")
par(mfrow = c(2, 2), mar = c(4.8, 5.1, 2.5, 1.1))
yaxs_lims <- c(-15, 95)
xaxs_lims <- c(-6.5, 6.5)
eps_set <- 0.2
for (i in seq_along(beta_gform_burn)) {
  
  # calculate fitted summaries
  gform_mean1 <- apply(beta_gform_burn[[i]], 2, mean)
  gform_vlow1 <- apply(beta_gform_burn[[i]], 2, quantile, 0.025)
  gform_low1 <- apply(beta_gform_burn[[i]], 2, quantile, 0.1)
  gform_high1 <- apply(beta_gform_burn[[i]], 2, quantile, 0.9)
  gform_vhigh1 <- apply(beta_gform_burn[[i]], 2, quantile, 0.975)
  gform_mean <- apply(beta_gform_clip[[i]], 2, mean)
  gform_vlow <- apply(beta_gform_clip[[i]], 2, quantile, 0.025)
  gform_low <- apply(beta_gform_clip[[i]], 2, quantile, 0.1)
  gform_high <- apply(beta_gform_clip[[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(beta_gform_clip[[i]], 2, quantile, 0.975)
  
  plot(c(-exp(gform_mean) / 5.5), hist_tmp$mids,
       type = "n", las = 1, bty = "l",
       ylab = "Depth (mm)", xlab = "Sprouts / mm",
       ylim = yaxs_lims, xlim = xaxs_lims,
       xaxt = "n")
  
  axis(1, at = seq(-6, 6, by = 2), labels = abs(seq(-6, 6, by = 2)))
  
  # add shaded region to denote sampling limits
  polygon(2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          c(-20, -20, -5, -5),
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  polygon(2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          c(80, 80, 150, 150),
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  
  # plot shaded regions for credible intervals
  polygon(c(-exp(gform_vlow) / 5.5 - eps_set, rev(-exp(gform_vhigh) / 5.5) - eps_set),
          c(hist_tmp$mids, rev(hist_tmp$mids)),
          border = NA, col = "gray75")
  polygon(c(-exp(gform_low) / 5.5 - eps_set, rev(-exp(gform_high) / 5.5) - eps_set),
          c(hist_tmp$mids, rev(hist_tmp$mids)),
          border = NA, col = "gray55")
  lines(c(-exp(gform_mean) / 5.5) - eps_set, hist_tmp$mids, lwd = 2, col = "black")

  polygon(c(exp(gform_vlow1) / 5.5 + eps_set, rev(exp(gform_vhigh1) / 5.5) + eps_set),
          c(hist_tmp$mids, rev(hist_tmp$mids)),
          border = NA, col = "gray75")
  polygon(c(exp(gform_low1) / 5.5 + eps_set, rev(exp(gform_high1) / 5.5) + eps_set),
          c(hist_tmp$mids, rev(hist_tmp$mids)),
          border = NA, col = "gray55")
  lines(c(exp(gform_mean1) / 5.5) + eps_set, hist_tmp$mids, lwd = 2, col = "black")
  
  lines(c(0, 0),
        c(min(hist_tmp$mids) - 10, max(hist_tmp$mids) + 10),
        lty = 1, lwd = 1, col = "black")
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0.5, adj = -0.02, cex = 1.5)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/burn_effects.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)

plot_name <- c("Forbs", "Grasses", "Sub-shrubs", "Woody plants")
par(mfrow = c(2, 2), mar = c(4.8, 5.1, 2.5, 1.1))
xaxs_lims <- c(-15, 95)
yaxs_lims <- c(-2, 6.5)

for (i in seq_along(beta_gform_burn)) {
  
  # calculate fitted summaries
  burn_int <- (exp(beta_gform_burn[[i]]) - exp(beta_gform_clip[[i]])) / 5.5
  gform_mean <- apply(burn_int, 2, mean)
  gform_vlow <- apply(burn_int, 2, quantile, 0.025)
  gform_low <- apply(burn_int, 2, quantile, 0.1)
  gform_high <- apply(burn_int, 2, quantile, 0.9)
  gform_vhigh <- apply(burn_int, 2, quantile, 0.975)
  
  plot(gform_mean ~ hist_tmp$mids,
       type = "n", las = 1, bty = "l",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
       ylim = yaxs_lims, xlim = xaxs_lims)
  
  # add shaded region to denote sampling limits
  polygon(c(-20, -20, -5, -5),
          2 * c(yaxs_lims, rev(yaxs_lims)),
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  polygon(c(80, 80, 150, 150),
          2 * c(yaxs_lims, rev(yaxs_lims)),
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(gform_vlow, rev(gform_vhigh)),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(gform_low, rev(gform_high)),
          border = NA, col = "gray55")
  
  lines(c(min(hist_tmp$mids) - 10, max(hist_tmp$mids) + 10),
        c(0, 0), lty = 2)
  lines(gform_mean ~ hist_tmp$mids, lwd = 2, col = "black")
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0.5, adj = -0.02, cex = 1.5)
  
}

dev.off()

get(plot_type)(file = paste0("outputs/plots/fitted_clip.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- c("Forbs", "Grasses", "Sub-shrubs", "Woody plants")
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
       type = "n", las = 1, bty = "l",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
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
          lwd = 1, col = ggplot2::alpha("gray50", 0.4))
    
    plot_all <- plot_all + plot_tmp
    
  }
  
  # add shaded region to denote sampling limits
  polygon(c(-20, -20, -5, -5),
          2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  polygon(c(80, 80, 150, 150),
          2 * c(yaxs_lims, rev(yaxs_lims)) - 5,
          border = NA, col = ggplot2::alpha("gray85", 0.5))
  
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "black")
  lines(c(exp(c(plot_all / n_plot)) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "black", lty = 2)
  
  # add some labels
  mtext(plot_name[i], side = 3, line = 0.5, adj = -0.02, cex = 1.5)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/trait_effect_burn_maxht.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- attr(trait_data$MAXHTstd, "scaled:center") + trait_seq[[1]] * attr(trait_data$MAXHTstd, "scaled:scale") 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 3.5)
for (i in seq_along(trait_effect_burn[[1]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[1]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_burn[[1]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = "l", las = 1, bty = "l",
       lwd = 2, col = "gray30",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = "gray55")
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "gray30")
  
  # add some labels
  mtext(paste0("Maximum height = ", round(plot_name[i], 0), " cm"), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/trait_effect_burn_stems.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- attr(trait_data$STEMCATSstd, "scaled:center") + trait_seq[[2]] * attr(trait_data$STEMCATSstd, "scaled:scale") 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 2.5)
for (i in seq_along(trait_effect_burn[[2]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[2]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_burn[[2]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = "l", las = 1, bty = "l",
       lwd = 2, col = "gray30",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = "gray55")
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "gray30")
  
  # add some labels
  mtext(paste0("Stem counts = ", round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/trait_effect_burn_sla.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- attr(trait_data$SLAstd, "scaled:center") + trait_seq[[3]] * attr(trait_data$SLAstd, "scaled:scale") 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 2.5)
for (i in seq_along(trait_effect_burn[[3]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_burn[[3]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_burn[[3]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = "l", las = 1, bty = "l",
       lwd = 2, col = "gray30",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = "gray55")
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "gray30")
  
  # add some labels
  mtext(paste0("SLA = ", round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/trait_effect_clip_maxht.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- attr(trait_data$MAXHTstd, "scaled:center") + trait_seq[[1]] * attr(trait_data$MAXHTstd, "scaled:scale") 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 3.5)
for (i in seq_along(trait_effect_clip[[1]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[1]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_clip[[1]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = "l", las = 1, bty = "l",
       lwd = 2, col = "gray30",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = "gray55")
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "gray30")
  
  # add some labels
  mtext(paste0("Maximum height = ", round(plot_name[i], 0), " cm"), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/trait_effect_clip_stems.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- attr(trait_data$STEMCATSstd, "scaled:center") + trait_seq[[2]] * attr(trait_data$STEMCATSstd, "scaled:scale") 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 2.5)
for (i in seq_along(trait_effect_clip[[2]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[2]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_clip[[2]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = "l", las = 1, bty = "l",
       lwd = 2, col = "gray30",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = "gray55")
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "gray30")
  
  # add some labels
  mtext(paste0("Stem counts = ", round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

get(plot_type)(file = paste0("outputs/plots/trait_effect_clip_sla.", file_type), width = width_set, height = height_set, units = unit_set, res = res_set)
plot_name <- attr(trait_data$SLAstd, "scaled:center") + trait_seq[[3]] * attr(trait_data$SLAstd, "scaled:scale") 
par(mfrow = c(2, 2), mar = c(5.1, 5.1, 1.8, 1.1))
ylim_set <- c(0, 2.5)
for (i in seq_along(trait_effect_clip[[3]])) {
  
  # calculate fitted summaries
  gform_mean <- apply(trait_effect_clip[[3]][[i]], 2, mean)
  gform_vlow <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.025)
  gform_low <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.1)
  gform_high <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.9)
  gform_vhigh <- apply(trait_effect_clip[[3]][[i]], 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids,
       type = "l", las = 1, bty = "l",
       lwd = 2, col = "gray30",
       xlab = "Depth (mm)", ylab = "Sprouts / mm",
       ylim = ylim_set)
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_vlow) / 5.5, rev(exp(gform_vhigh) / 5.5)),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(gform_low) / 5.5, rev(exp(gform_high) / 5.5)),
          border = NA, col = "gray55")
  lines(c(exp(gform_mean) / 5.5) ~ hist_tmp$mids, lwd = 2, col = "gray30")
  
  # add some labels
  mtext(paste0("SLA = ", round(plot_name[i], 1)), side = 3, line = 0, adj = 1, cex = 1.25)
  
}
dev.off()

# trait effects
get(plot_type)(file = paste0("outputs/plots/regression_coefficients.", file_type), width = width_set, height = height_set + 3,
               units = unit_set, res = res_set)
layout(matrix(1:8, ncol = 2, byrow = TRUE), heights = c(0.12, rep(1, 3)))
par(mar = c(0, 0, 0, 0))
plot(1 ~ 1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
text(1, 1, "Clip", cex = 2.8, xpd = TRUE)
plot(1 ~ 1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
text(1, 1, "Burn", cex = 2.8, xpd = TRUE)
par(mar = c(4.8, 5.1, 2.5, 1.1))
ylim_set <- c(0.05, 2.0)
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
       type = "l", las = 1, bty = "l",
       lwd = 2, col = "gray30",
       xlab = "Depth (mm)", ylab = "Effect",
       ylim = ylim_set,
       cex.lab = 1.5,
       log = "y")
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(plot_vlow), rev(exp(plot_vhigh))),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(plot_low), rev(exp(plot_high))),
          border = NA, col = "gray55")
  
  lines(exp(plot_mean) ~ hist_tmp$mids, lwd = 2, col = "gray30")
  lines(c(min(hist_tmp$mids) - 10, max(hist_tmp$mids) + 10),
        c(1, 1), lty = 2)
  
  if (i == 1) {
    data_clip <- data_long$depths_mm[data_long$treat == "c"]
    for (j in seq_along(data_clip)) {
      lines(c(data_clip[j], data_clip[j]), c(0.06, 0.08),
            lty = 1, lwd = 1, col = "gray30")
    }
  }
  
  #
  mtext(plot_label[(2 * i) - 1], side = 3, line = 0.5, adj = -0.02, cex = 1.3)
  
  # calculate fitted summaries
  plot_mean <- apply(burn_parameters, 2, mean)
  plot_vlow <- apply(burn_parameters, 2, quantile, 0.025)
  plot_low <- apply(burn_parameters, 2, quantile, 0.1)
  plot_high <- apply(burn_parameters, 2, quantile, 0.9)
  plot_vhigh <- apply(burn_parameters, 2, quantile, 0.975)
  
  # plot mean fitted value
  plot(exp(plot_mean) ~ hist_tmp$mids,
       type = "l", las = 1, bty = "l",
       lwd = 2, col = "gray30",
       xlab = "Depth (mm)", ylab = "Effect",
       ylim = ylim_set,
       cex.lab = 1.5,
       log = "y")
  
  # plot shaded regions for credible intervals
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(plot_vlow), rev(exp(plot_vhigh))),
          border = NA, col = "gray75")
  polygon(c(hist_tmp$mids, rev(hist_tmp$mids)),
          c(exp(plot_low), rev(exp(plot_high))),
          border = NA, col = "gray55")
  
  lines(exp(plot_mean) ~ hist_tmp$mids, lwd = 2, col = "gray30")
  lines(c(min(hist_tmp$mids) - 10, max(hist_tmp$mids) + 10),
        c(1, 1), lty = 2)
  
  if (i == 1) {
    data_burn <- data_long$depths_mm[data_long$treat == "b"]
    for (j in seq_along(data_burn)) {
      lines(c(data_burn[j], data_burn[j]), c(0.06, 0.08),
            lty = 1, lwd = 1, col = "gray30")
    }
  }

}
dev.off()

# plot resprouting effects
# combos to plot: clip x gform, burn x gform (= 8 combos inc. forbs)
# slopes to slope: clip x traits, burn x traits (= 6 combos)
colnames(beta_resprout_samples) <- c("Clip", "Grasses", "Sub-shrubs", "Woody plants",
                                     "Max. height", "Stem counts", "SLA",
                                     "Clip:Grasses", "Clip:Sub-shrubs", "Clip:Woody plants",
                                     "Clip:Max. height", "Clip:Stem counts", "Clip:SLA")

get(plot_type)(file = paste0("outputs/plots/resprout_trait_effects.", file_type), width = width_set, height = height_set - 1,
               units = unit_set, res = res_set)

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
     xlim = c(-2.0, 2.5),
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
axis(1, at = c(-2.0, -1.0, 0.0, 1.0, 2.0))
axis(2, at = seq_len(ncol(slope_quant)),
     labels = rev(colnames(slope_combos)), las = 1)
mtext("b", side = 3, line = 0.5, adj = -0.02, cex = 1.5)

dev.off()
