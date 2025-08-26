# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# Function to ensure directories exist or create them
ensure_directory <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Function to calculate sensible y-axis limits
calculate_y_limits <- function(values) {
  range <- max(values, na.rm = TRUE) - min(values, na.rm = TRUE)
  lower <- min(values, na.rm = TRUE) - 0.1 * range
  upper <- max(values, na.rm = TRUE) + 0.1 * range
  c(lower, upper)
}

# Function to generate plots for variables
generate_variable_boxplots <- function(cohort_characteristics, variables, cohort_colors) {
  plots <- list()
  for (var in variables) {
    values <- cohort_characteristics[[var]]
    y_limits <- if (var %in% c("optimal_K", "numberofevents")) calculate_y_limits(values) else c(0, 1)
    plot <- ggplot(cohort_characteristics, aes(x = as.factor(1), y = !!sym(var))) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, fill = "gray", color = "black", lwd = 0.5) +
      geom_jitter(aes(color = cohort, fill = cohort), shape = 21, size = 5, alpha = 0.8, width = 0.2, stroke = 0) +
      scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
      labs(title = paste("Distr. of", var), x = NULL, y = var) +
      scale_fill_manual(values = cohort_colors) +
      scale_color_manual(values = cohort_colors) +
      theme_minimal() +
      theme(
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      coord_cartesian(xlim = c(0, 2))
    plots[[var]] <- plot
  }
  return(plots)
}

# Function to combine and save plots
save_combined_plots <- function(plots, plot_data, filename, script_name) {
  # Combine plots into a grid
  grid_list <- lapply(1:9, function(i) if (i <= length(plots)) plots[[i]] else ggplot() + theme_void())
  combined_plots <- plot_grid(plotlist = grid_list, ncol = 3, align = "hv", labels = "auto")
  
  # Add a header with the script name
  combined_plots_with_header <- ggdraw() +
    draw_label(script_name, size = 16, fontface = 'bold', x = 0.5, y = 1, vjust = 1) +
    draw_plot(combined_plots, 0, -0.05, 1, 1)
  
  # Save plots as PDF and PNG
  ggsave(paste0(filename, ".pdf"), plot = combined_plots_with_header, width = 12, height = 12, units = "in")
  ggsave(paste0(filename, ".png"), plot = combined_plots_with_header, width = 12, height = 12, units = "in")
  
  # Save plot data to CSV
  # Extract the directory path
  directory_path <- dirname(filename)
  # Extract the file name
  file_name <- basename(filename)
  plot_data_file <- file.path(directory_path,paste0("SourceData_", file_name, ".csv"))
  write_csv(plot_data, plot_data_file)
}
