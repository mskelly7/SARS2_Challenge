# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH
# Supplemental Figure S1 – symptom burden and viral load kinetics
# Last updated: 10/07/2026

script_summary <- "This script generates Figure S1 summarizing longitudinal symptom burden and viral load kinetics following SARS-CoV-2 challenge. 
Symptom scores are summarized by low versus high symptom burden, and viral load is summarized separately for nasal (MTS) and throat samples by low 
versus high viral-load groups. Lines show mean values over time, while bars show the percentage of participants with symptoms or detectable viral 
load at each sampling visit."
cat(script_summary)

remove(list=ls())
library(dplyr)
library(ggplot2)
library(cowplot)

project_dir <- "________________________"
data_dir <- file.path(project_dir,"Data")
fig_dir <- file.path(project_dir,"Figures")
suppfig1_dir <- file.path(fig_dir,"FigureS1")
input_file <- "SARS2_challenge_sx_vl.csv"

dir.create(suppfig1_dir, recursive=TRUE, showWarnings=FALSE)

# Plot colors

cols_symptoms_named <- c(
  low = "#C75A5A",
  high = "#5E86C2"
)

cols_vl_named <- c(
  low = "#236343",
  high = "#653875"
)

# Import and validate source data

sx_vl <- read.csv(
  file.path(data_dir,input_file),
  check.names = FALSE,
  na.strings = c("", "NA")
)

required_columns <- c(
  "studyID_random", "timepoint", "sampletype",
  "sx_score", "sx_score_cat", "vl_log10", "vl_cat"
)

missing_columns <- setdiff(required_columns, names(sx_vl))

if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

valid_timepoint <- !is.na(sx_vl$timepoint) &
  grepl("^-?[0-9]+(am|mid|pm)?$", sx_vl$timepoint)

if (any(!valid_timepoint)) {
  stop("Invalid or missing timepoint labels.")
}

# Prepare plotting data

plot_data <- sx_vl %>%
  mutate(
    day = as.numeric(sub("(am|mid|pm)$", "", timepoint)),
    sampling_time = case_when(
      grepl("am$", timepoint) ~ "am",
      grepl("mid$", timepoint) ~ "mid",
      grepl("pm$", timepoint) ~ "pm",
      TRUE ~ "baseline"
    )
  )

sx_data <- plot_data %>%
  filter(
    day >= 0,
    !is.na(studyID_random),
    !is.na(sx_score),
    sx_score_cat %in% c("low", "high")
  ) %>%
  transmute(
    studyID_random,
    day,
    timepoint,
    x = day + case_when(
      sampling_time == "mid" ~ 1 / 3,
      sampling_time == "pm" ~ 2 / 3,
      TRUE ~ 0
    ),
    group = sx_score_cat,
    value = sx_score
  ) %>%
  distinct()

vl_data <- plot_data %>%
  filter(
    day <= 14,
    sampling_time %in% c("am", "pm"),
    !is.na(studyID_random),
    !is.na(vl_log10),
    vl_cat %in% c("low", "high"),
    sampletype %in% c("MTS", "Throat")
  ) %>%
  transmute(
    studyID_random,
    sampletype,
    day,
    timepoint,
    x = day + if_else(sampling_time == "pm", 0.5, 0),
    group = vl_cat,
    value = vl_log10
  ) %>%
  distinct()

# Check that each participant contributes once per visit/sample type.
stopifnot(
  !anyDuplicated(sx_data[c("studyID_random", "timepoint")]),
  !anyDuplicated(
    vl_data[c("studyID_random", "sampletype", "timepoint")]
  )
)
if (
  nrow(sx_data) == 0 ||
  !all(c("MTS", "Throat") %in% vl_data$sampletype)
) {
  stop("Symptom data and both viral-load sample types are required.")
}

# Summarize symptoms and viral load by group and visit

sx_summary <- sx_data %>%
  group_by(group, day, x) %>%
  summarise(
    mean_value = mean(value),
    n_observed = n(),
    n_with_symptoms = sum(value > 0),
    percent = 100 * mean(value > 0),
    .groups = "drop"
  )

vl_summary <- vl_data %>%
  group_by(sampletype, group, day, x) %>%
  summarise(
    mean_value = mean(value),
    n_observed = n(),
    n_positive = sum(value > 0),
    percent = 100 * mean(value > 0),
    .groups = "drop"
  )

vl_y_max <- max(1, ceiling(max(vl_summary$mean_value)))
vl_x_range <- range(vl_summary$x)
vl_day_breaks <- seq(
  min(vl_summary$day),
  max(vl_summary$day),
  by = 1
)

# Plotting theme and panel function

figure_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.title.x = element_text(margin = margin(t = 8, b = 2)),
    axis.text = element_text(size = 11),
    axis.ticks.length = grid::unit(2.5, "mm"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.width = grid::unit(1.1, "cm"),
    legend.margin = margin(t = -6, r = 0, b = 0, l = 0)
  )

make_panel <- function(
    data, title, y_label, percent_label,
    colors, group_labels, visit_spacing,
    y_max = NULL, x_range = NULL, day_breaks = NULL) {
  if (is.null(y_max)) {
    y_max <- max(1, ceiling(max(data$mean_value)))
  }
  if (is.null(x_range)) {
    x_range <- range(data$x)
  }
  if (is.null(day_breaks)) {
    day_breaks <- seq(min(data$day), max(data$day), by = 1)
  }
  # Two bars per visit, centered around the visit position.
  bar_width <- visit_spacing * 0.35
  bar_offset <- visit_spacing * 0.20
  data <- data %>%
    mutate(
      bar_x = x + if_else(group == "low", -bar_offset, bar_offset),
      bar_height = percent / 100 * y_max
    )
  visits <- data %>% distinct(x)
  ggplot(data) +
    # Bars: percentages read against the right axis.
    geom_col(
      aes(x = bar_x, y = bar_height, fill = group),
      width = bar_width,
      alpha = 0.20,
      show.legend = FALSE
    ) +
    # Lines: means read against the left axis.
    geom_line(
      aes(x = x, y = mean_value, color = group, group = group),
      linewidth = 1.2
    ) +
    geom_point(
      aes(x = x, y = mean_value, color = group),
      size = 1.6
    ) +
    # Short unlabeled ticks mark the actual sampling visits.
    geom_rug(
      data = visits,
      aes(x = x),
      inherit.aes = FALSE,
      sides = "b",
      length = grid::unit(0.012, "npc"),
      color = "grey35"
    ) +
    scale_color_manual(
      values = colors,
      limits = c("low", "high"),
      labels = group_labels
    ) +
    scale_fill_manual(
      values = colors,
      limits = c("low", "high"),
      labels = group_labels,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = x_range + c(-visit_spacing / 2, visit_spacing / 2),
      breaks = day_breaks,
      minor_breaks = NULL,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      name = y_label,
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.04)),
      sec.axis = sec_axis(
        ~ . / y_max * 100,
        name = percent_label,
        breaks = seq(0, 100, by = 25)
      )
    ) +
    labs(
      title = title,
      x = "Day relative to inoculation"
    ) +
    figure_theme
}

# Figure S1A – symptom burden

p_sx <- make_panel(
  data = sx_summary,
  title = NULL,
  y_label = "Mean total symptom score",
  percent_label = "Symptoms (%)",
  colors = cols_symptoms_named,
  group_labels = c("Low symptom burden", "High symptom burden"),
  visit_spacing = 1 / 3
) +
  theme(
    plot.margin = margin(t = 10, r = 5.5, b = 5.5, l = 5.5)
  )

# Figure S1B – nasal viral load

p_mts <- make_panel(
  data = filter(vl_summary, sampletype == "MTS"),
  title = "Nasal",
  y_label = expression("Mean"~log[10]~"viral load"),
  percent_label = "Detectable viral load (%)",
  colors = cols_vl_named,
  group_labels = c("Low viral load", "High viral load"),
  visit_spacing = 0.5,
  y_max = vl_y_max,
  x_range = vl_x_range,
  day_breaks = vl_day_breaks
)

# Figure S1C – throat viral load

p_throat <- make_panel(
  data = filter(vl_summary, sampletype == "Throat"),
  title = "Throat",
  y_label = expression("Mean"~log[10]~"viral load"),
  percent_label = "Detectable viral load (%)",
  colors = cols_vl_named,
  group_labels = c("Low viral load", "High viral load"),
  visit_spacing = 0.5,
  y_max = vl_y_max,
  x_range = vl_x_range,
  day_breaks = vl_day_breaks
)

# Assemble Figure S1

bottom_panels <- plot_grid(
  p_mts, p_throat,
  ncol = 2,
  labels = c("B", "C"),
  label_size = 16,
  align = "h",
  axis = "tb"
)

figure_s1 <- plot_grid(
  NULL,
  p_sx,
  bottom_panels,
  ncol = 1,
  labels = c("A", "", ""),
  label_size = 16,
  rel_heights = c(0.05, 1, 1.05)
)

print(figure_s1)

# Save individual panels

ggsave(file.path(suppfig1_dir,"FigureS1a.png"), p_sx, width=15, height=5, dpi=600, bg="white")
ggsave(file.path(suppfig1_dir,"FigureS1b.png"), p_mts, width=7.5, height=5, dpi=600, bg="white")
ggsave(file.path(suppfig1_dir,"FigureS1c.png"), p_throat, width=7.5, height=5, dpi=600, bg="white")

# Save final composite figure

ggsave(
  filename=file.path(fig_dir,"FigureS1.png"),
  plot=figure_s1,
  width=15,
  height=10,
  dpi=600,
  bg="white"
)
