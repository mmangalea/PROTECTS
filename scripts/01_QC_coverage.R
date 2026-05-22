# =============================================================================
# Project:  PROTECTS-Microbiome
# Script:   01_QC_coverage.R
# Analyst:  M. Mangalea
# Date:     2026-05-18
# Purpose:  QC metrics (read counts by subject) and study timeline (Fig 1A).
# =============================================================================

source("scripts/00_packages.R")

options(scipen = 999)
set.seed(42)

# ---- Paths -------------------------------------------------------------------

plots_dir   <- "figures/QC"
results_dir <- "results/QC"

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  read_counts = "tables/read_counts.tsv",
  metadata    = "tables/sample_metadata.tsv",
  timeline    = "data/metadata/sampling_timeline.tsv"
)

# ==============================================================================
# 1. Read Counts QC
# ==============================================================================

read_counts <- read.delim(
  paths$read_counts,
  header = TRUE,
  stringsAsFactors = FALSE
)

metadata <- read.delim(
  paths$metadata,
  header = TRUE,
  stringsAsFactors = FALSE
)

# Extract subject ID from sample name (e.g. "01_A" -> "01")
read_counts <- read_counts %>%
  mutate(subject = str_extract(Samples, "^\\d+"))

readcount_plot <- ggplot(read_counts, aes(x = subject, y = Trimmed_reads, fill = subject)) +
  geom_point(shape = 21, size = 4, alpha = 0.7, position = position_jitter(width = 0.2)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  labs(title = NULL,
       x = NULL,
       y = "Total Reads (in millions)",
       fill = "Subject") +
  coord_cartesian() +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  theme_classic() +
  theme(
    axis.line.x = element_line(size = 1, colour = "black", linetype = 1),
    axis.line.y = element_line(size = 1, colour = "black", linetype = 1),
    axis.ticks = element_line(size = 1, color = "black"),
    axis.ticks.length = unit(.25, "cm"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

print(readcount_plot)

ggsave(
  filename = file.path(plots_dir, "SUPP_FIG_1A_readcounts.svg"),
  plot = readcount_plot,
  scale = 1.5, width = 4, height = 4, units = "in", dpi = 300,
  device = grDevices::svg
)

# ==============================================================================
# 2. Study Timeline — Figure 1A
# ==============================================================================

timeline_data <- read.delim(paths$timeline, stringsAsFactors = FALSE)

# Format participant IDs as zero-padded two-digit strings
timeline_data$Participant_ID <- sprintf("%02d", timeline_data$Participant_ID)

# Convert to long format (values are days relative to enrollment)
timeline_data_long <- timeline_data %>%
  pivot_longer(cols = -Participant_ID, names_to = "Sampling_Event", values_to = "Day") %>%
  filter(!is.na(Day) & Day != "") %>%
  mutate(Day = as.numeric(Day)) %>%
  mutate(Event_Type = case_when(
    Sampling_Event == "Enrollment" ~ "Enrollment",
    Sampling_Event %in% c("Serology_1", "Serology_2", "Serology_3") ~ "Serology",
    Sampling_Event %in% c("Vaccine_1", "Vaccine_2") ~ "Vaccination",
    TRUE ~ "Microbiome"
  ))

timeline_data_long$Participant_ID <- factor(
  timeline_data_long$Participant_ID,
  levels = rev(sort(unique(timeline_data_long$Participant_ID)))
)

##### Plot for Publication - Figure 1A ######

timeline_plot <- ggplot(timeline_data_long, aes(x = Day, y = Participant_ID, group = Participant_ID)) +
  geom_line(color = "black", size = 1) +
  geom_point(aes(fill = Sampling_Event, shape = Event_Type), color = "black", size = 5, stroke = 0.5) +
  labs(x = "Days Relative to Enrollment", y = "Participant", title = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.spacing.y = unit(0.5, "cm"),
    legend.margin = margin(10, 10, 10, 10),
    legend.key.width = unit(1.5, "lines"),
    legend.key.height = unit(1.5, "lines"),
    legend.justification = c(1, 0.5),
    legend.box.just = "right",
    legend.box.background = element_rect(fill = "White", colour = NA)
  ) +
  scale_fill_manual(name = "Sample Type",
                    values = c("Vaccine_1" = "purple1",
                               "Vaccine_2" = "purple4",
                               "Serology_1" = "red1",
                               "Serology_2" = "red3",
                               "Serology_3" = "red4",
                               "Stool_Sample_1" = "skyblue",
                               "Stool_Sample_2" = "dodgerblue",
                               "Stool_Sample_3" = "dodgerblue4",
                               "Enrollment" = "forestgreen")) +
  scale_shape_manual(name = "Event Type",
                     values = c("Microbiome" = 21,
                                "Enrollment" = 22,
                                "Serology" = 25,
                                "Vaccination" = 24)) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21), title.position = "top", label.hjust = .5),
    shape = guide_legend(title.position = "top", label.hjust = .5)
  )

print(timeline_plot)

ggsave(
  filename = file.path(plots_dir, "FIG_1A_timeline.svg"),
  plot = timeline_plot,
  width = 5, height = 4, scale = 1.5,
  device = grDevices::svg
)

# ------------------------------------------------------------------------------
# 3. Nonpareil coverage curves
# ------------------------------------------------------------------------------

# library(Nonpareil)

options(scipen = 0)  # Use scientific notation (e.g. 1e+07) for axis labels

npo_dir <- "data/processed/nonpareil"
samples <- read.delim(file.path(npo_dir, "samples.tsv"), stringsAsFactors = FALSE)

# Build full paths to .npo files
npo_files <- file.path(npo_dir, samples$File)
stopifnot(all(file.exists(npo_files)))

# Generate Nonpareil curves
svg(file.path(plots_dir, "SUPP_FIG_1B_nonpareil_curves.svg"), width = 12, height = 9)
nps <- Nonpareil.set(
  npo_files,
  col    = samples$Col,
  labels = samples$Name,
  plot.opts = list(plot.observed = FALSE, model.lwd = 2,
                   legend.opts = list(cex = 0.75))
)
dev.off()

# Extract summary metrics
nps_summary <- summary(nps)

coverage <- data.frame(
  sample   = rownames(nps_summary),
  coverage = nps_summary[, "C"] * 100
)

diversity <- data.frame(
  sample    = rownames(nps_summary),
  diversity = nps_summary[, "diversity"]
)

seq_effort <- data.frame(
  sample     = rownames(nps_summary),
  seq_effort = nps_summary[, "LRstar"] / 1e9
)

# Write results
qc_dir <- "results/QC"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(coverage,   file.path(qc_dir, "nonpareil_coverage.csv"),   row.names = FALSE)
write.csv(diversity,  file.path(qc_dir, "nonpareil_diversity.csv"),  row.names = FALSE)
write.csv(seq_effort, file.path(qc_dir, "nonpareil_seq_effort.csv"), row.names = FALSE)

message("01_QC_coverage.R complete.")
