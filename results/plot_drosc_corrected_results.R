suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

resource_dir <- "/Users/taehyeon/Downloads/DRoSC-resources"
output_dir <- "/tmp/DRoSC-corrected-figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_levels <- c("Normality", "OBA", "Proposed")
method_colors <- c(Normality = "#F8766D", OBA = "#00BA38", Proposed = "#619CFF")
method_shapes <- c(Normality = 16, OBA = 17, Proposed = 15)

inference_data <- read.csv(
  file.path(resource_dir, "inference", "figure-data", "inference_figure_data.csv"),
  check.names = FALSE
)

safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

collapse_inference_targets <- function(data) {
  data$target_key <- round(data$tau.star, 2)
  group_id <- interaction(
    data$setting, data$config_id, data$phi, data$T0, data$T1,
    data$target_key, data$method,
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(data)), group_id)
  collapsed <- lapply(groups, function(index) {
    group <- data[index, ]
    data.frame(
      setting = group$setting[1],
      config_id = group$config_id[1],
      phi = group$phi[1],
      T0 = group$T0[1],
      T1 = group$T1[1],
      tau.star = group$target_key[1],
      method = as.character(group$method[1]),
      coverage = min(group$coverage, na.rm = TRUE),
      length = max(group$length, na.rm = TRUE),
      union_length = safe_max(group$union_length)
    )
  })
  do.call(rbind, collapsed)
}

inference_data <- collapse_inference_targets(inference_data)
inference_data$method <- factor(inference_data$method, levels = method_levels)

inference_panel <- function(data, outcome, title) {
  plot <- ggplot(data, aes(x = tau.star, y = .data[[outcome]], color = method,
                           shape = method, group = method)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.9) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_shape_manual(values = method_shapes, drop = FALSE) +
    scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5)) +
    labs(x = expression(tau^"*"), y = NULL, title = title, color = NULL, shape = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_line(linewidth = 0.25),
      plot.background = element_rect(fill = "white", color = NA)
    )

  if (outcome == "coverage") {
    plot <- plot + geom_hline(yintercept = 0.95, color = "red", linetype = "dashed",
                              linewidth = 0.55)
  }
  plot
}

save_inference_figure <- function(phi, T0, T1, filename, setting = NULL,
                                  width = 10.7, height = 9.6) {
  data <- inference_data[
    inference_data$phi == phi & inference_data$T0 == T0 & inference_data$T1 == T1,
  ]
  if (!is.null(setting)) data <- data[data$setting == setting, ]

  settings <- if (is.null(setting)) 1:3 else setting
  panels <- list()
  for (current_setting in settings) {
    setting_data <- data[data$setting == current_setting, ]
    panels[[length(panels) + 1L]] <- inference_panel(setting_data, "coverage", "Coverage")
    panels[[length(panels) + 1L]] <- inference_panel(setting_data, "length", "Length")
  }
  figure <- ggarrange(
    plotlist = panels,
    ncol = 2,
    nrow = length(settings),
    common.legend = TRUE,
    legend = "bottom"
  )
  ggsave(file.path(output_dir, filename), figure, width = width, height = height,
         dpi = 200, bg = "white")
}

save_inference_figure(0, 25, 25, "Inference-S2-phi0-pre25-post25.png",
                      setting = 2, height = 3.9)
save_inference_figure(0, 25, 25, "Inference-phi0-pre25-post25.png")
save_inference_figure(0, 25, 50, "Inference-phi0-pre25-post50.png")
save_inference_figure(0, 50, 50, "Inference-phi0-pre50-post50.png")
save_inference_figure(0.5, 25, 25, "Inference-phi0.5-pre25-post25.png")
save_inference_figure(0.5, 25, 50, "Inference-phi0.5-pre25-post50.png")

prop_data <- read.csv(
  file.path(resource_dir, "prop", "figure-data", "prop_figure_data.csv"),
  check.names = FALSE
)

collapse_prop_targets <- function(data) {
  data$target_key <- round(data$tau.star, 2)
  group_id <- interaction(
    data$setting, data$phi, data$T0, data$T1,
    data$target_key, data$threshold,
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(data)), group_id)
  collapsed <- lapply(groups, function(index) {
    group <- data[index, ]
    data.frame(
      setting = group$setting[1],
      phi = group$phi[1],
      T0 = group$T0[1],
      T1 = group$T1[1],
      tau.star = group$target_key[1],
      threshold = as.character(group$threshold[1]),
      coverage = min(group$coverage, na.rm = TRUE),
      length = max(group$length, na.rm = TRUE),
      feasible_proportion = mean(group$feasible_proportion, na.rm = TRUE)
    )
  })
  do.call(rbind, collapsed)
}

prop_data <- collapse_prop_targets(prop_data)
prop_data$threshold <- factor(prop_data$threshold, levels = c("10%", "20%", "30%"))
threshold_colors <- c("10%" = "#F8766D", "20%" = "#00BA38", "30%" = "#619CFF")
threshold_shapes <- c("10%" = 16, "20%" = 17, "30%" = 15)

prop_panel <- function(data, outcome, title) {
  ggplot(data, aes(x = tau.star, y = .data[[outcome]], color = threshold,
                   shape = threshold, group = threshold)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.8) +
    scale_color_manual(values = threshold_colors, drop = FALSE) +
    scale_shape_manual(values = threshold_shapes, drop = FALSE) +
    scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5)) +
    labs(x = expression(tau^"*"), y = NULL, title = title,
         color = "Threshold", shape = "Threshold") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_line(linewidth = 0.25),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

prop_panels <- list()
for (setting in 1:3) {
  setting_data <- prop_data[prop_data$setting == setting, ]
  prop_panels[[length(prop_panels) + 1L]] <- prop_panel(setting_data, "coverage", "Coverage")
  prop_panels[[length(prop_panels) + 1L]] <- prop_panel(setting_data, "length", "Length")
  prop_panels[[length(prop_panels) + 1L]] <- prop_panel(
    setting_data, "feasible_proportion", "Feasible proportion"
  )
}
prop_figure <- ggarrange(
  plotlist = prop_panels,
  ncol = 3,
  nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(file.path(output_dir, "Prop-pre25-post25.png"), prop_figure,
       width = 11, height = 10, dpi = 200, bg = "white")

basque_data <- read.csv(
  file.path(resource_dir, "Basque", "figure-data", "Basque_figure_data.csv"),
  check.names = FALSE
)
basque_environment <- new.env(parent = emptyenv())
load(
  file.path(resource_dir, "Basque", "figure-data", "Basque_figure_resources.RData"),
  envir = basque_environment
)
tau_sc <- unname(basque_environment$tau.SC)
zero_lambda <- min(basque_data$lambda[abs(basque_data$tauHat) < 1e-6])
basque_ci <- basque_data[!is.na(basque_data$CI.tau.1), ]

basque_figure <- ggplot(basque_data, aes(x = lambda)) +
  geom_errorbar(
    data = basque_ci,
    aes(ymin = CI.tau.1, ymax = CI.tau.2),
    color = "orange",
    width = 0.001,
    linewidth = 0.5
  ) +
  geom_line(aes(y = tauHat), color = "black", linewidth = 0.6) +
  geom_hline(yintercept = tau_sc, color = "blue", linetype = "dashed", linewidth = 0.55) +
  geom_vline(xintercept = zero_lambda, color = "red", linetype = "dashed", linewidth = 0.55) +
  scale_x_continuous(breaks = seq(0, 0.06, 0.02)) +
  labs(x = expression(lambda), y = expression(hat(tau))) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_line(linewidth = 0.25)
  )
ggsave(file.path(output_dir, "Basque_lambda_CI5.png"), basque_figure,
       width = 10.2, height = 3.05, dpi = 200, bg = "white")

cat("Corrected figures written to", output_dir, "\n")
