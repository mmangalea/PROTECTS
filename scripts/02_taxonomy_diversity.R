# =============================================================================
# Project:  PROTECTS-Microbiome
# Script:   02_taxonomy_diversity.R
# Analyst:  M. Mangalea
# Date:     2026-02-25
# Purpose:  Run taxonomy and diversity analyses for PROTECTS-Microbiome, including
#           plotting and outputting manuscript figures and updated summary statistics.
# =============================================================================

source("scripts/00_packages.R")

options(scipen = 999)
set.seed(42)

#---- Paths --------------------------------------------------------------------

input_dir  <- "data/processed"
metadata_dir <- "data/metadata"
plots_dir  <- "figures/taxonomy_diversity"
results_dir <- "results"

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
    metadata   = file.path(metadata_dir, "protects_metadata_heatmap.txt"),
    phyla      = file.path(input_dir, "merged_phyla.tsv"),
    shannon    = file.path(input_dir, "diversity_shannon.tsv"),
    richness   = file.path(input_dir, "richness_observed.tsv"),
    gini       = file.path(input_dir, "dominance_gini.tsv"),
    bray       = file.path(input_dir, "bray_curtis.tsv"),
    species    = file.path(input_dir, "merged_species.tsv")
)

# Read metadata once because the diversity sections all depend on it.
sample_metadata <- read.delim(
  paths$metadata,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  dplyr::rename(Sample = ID) |>
  dplyr::mutate(site = as.factor(site))

# ---- Input validation --------------------------------------------------------
invisible(lapply(paths, function(f) stopifnot(file.exists(f))))

# ==============================================================================
# Figure 1B: Relative Abundance by Individual
# ==============================================================================

merged_phyla <- read.delim(
  paths$phyla,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

col_sums <- colSums(merged_phyla, na.rm = TRUE)

shortfall <- 100 - col_sums

bacteria_unclassified_values <- merged_phyla["Bacteria_unclassified", , drop = FALSE]
unclassified_combined <- bacteria_unclassified_values + shortfall
rownames(unclassified_combined) <- "Unclassified"

merged_phyla <- merged_phyla[!grepl("^Bacteria_unclassified", rownames(merged_phyla)), ]

archaea_rows <- grep("^(Candidatus_|Euryarchaeota)", rownames(merged_phyla))
archaea_combined <- merged_phyla[archaea_rows, ] %>%
  summarise(across(everything(), sum, na.rm = TRUE))
rownames(archaea_combined) <- "Archaea"
merged_phyla <- merged_phyla[-archaea_rows, ]

merged_phyla_combined <- rbind(merged_phyla, archaea_combined, unclassified_combined)

ordered_taxa_names <- rownames(merged_phyla_combined)
merged_phyla_combined_long <- merged_phyla_combined %>%
  rownames_to_column(var = "Taxonomy") %>%
  pivot_longer(-Taxonomy, names_to = "Sample", values_to = "Value") %>%
  mutate(Taxonomy=factor(Taxonomy, levels=ordered_taxa_names))

merged_phyla_combined_long <- merged_phyla_combined_long %>%
  separate(Sample, into = c("Individual", "SampleTime"), sep = "_", remove = FALSE) %>%
  mutate(SampleTime = recode(SampleTime, "A" = "S1", "B" = "S2", "C" = "S3"))

# making custom color palette for the all taxa #
library(Polychrome)

num_taxa <- unique(merged_phyla_combined_long$Taxonomy)
num_colors <- length(num_taxa)

colors <- green.armytage.colors(num_colors)

taxa_palette <- setNames(colors, num_taxa)

stacked_phyla_plot <- ggplot(merged_phyla_combined_long, aes(x = SampleTime, y = Value, fill = Taxonomy)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = taxa_palette) +
  labs(x = "Sample", y = "MetaPhlAn4 Relative Abundance") +
  labs(fill = "Taxa") +
  theme_classic() +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), labels = paste0(seq(0, 100, by = 10), "%"))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "grey", colour = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle=90, size=14, colour="black", vjust = 1, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
  facet_wrap(~Individual, scales="free_x", ncol=10)

stacked_phyla_plot

ggsave(filename = file.path(plots_dir, "FIG_1B.svg"),
       plot = stacked_phyla_plot,
       width = 5,
       height = 3,
       units = 'in',
       dpi = 300,
       scale = 2,
       device = grDevices::svg)

write.csv(
  merged_phyla_combined_long,
  file.path(results_dir, "stacked_phyla_long_2026.csv"),
  row.names = FALSE
)


#---- Diversity: Shannon -------------------------------------------------------

shannon <- read.delim(
  paths$shannon,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  col.names = c("Sample", "Value"),
  stringsAsFactors = FALSE
)

shannon <- shannon %>%
  mutate(Sample = as.character(Sample))

merged_meta <- merge(
  shannon,
  sample_metadata,
  by = "Sample",
  all = TRUE
) %>%
  separate(
    Sample,
    into = c("Individual", "SampleTime"),
    sep = "_",
    remove = FALSE
  )

shannon_summary <- summary(merged_meta$Value)
kable(data.frame(
  Statistic = names(shannon_summary),
  Value = as.numeric(shannon_summary)
), caption = "Shannon Diversity Summary Statistics", digits = 3)

# library(ggsignif)

# Create Shannon diversity plot comparing spike inhibition groups
shannon_plot <- ggplot(merged_meta, 
                       aes(x = spike_inhibition, y = Value, 
                           color = spike_inhibition, fill = spike_inhibition)) +
  stat_boxplot(geom='errorbar', color="black", width=0.25, lwd=0.25) + 
  geom_signif(comparisons = list(c("no","yes")), 
              map_signif_level = FALSE, textsize=4, 
              test = "wilcox.test",
              y_position = max(merged_meta$Value, na.rm = TRUE) * 1)+
  geom_boxplot(width=0.5, color="black", outlier.shape = NA , fatten=1, lwd=0.25, alpha=0.25) + 
  geom_point(aes(color = spike_inhibition),size=5, shape=21, alpha=0.5)+
  geom_line(aes(group=Individual), lwd = 0.25)+
  scale_color_manual(name=NULL,
                     values=c("black", "black"),
                     breaks=c("no","yes"),
                     labels=c("Low", "High"))+
  scale_shape_manual(name=NULL,
                     values=c(21, 22, 23),
                     breaks=c("no","yes"),
                     labels=c("Low", "High"))+
  scale_fill_manual(name=NULL,
                    values=c("orange" ,"purple4"),
                    breaks=c("no","yes"),
                    labels=c("No", "Yes"))+
  scale_y_continuous(limits = c(3,4.5))+
  scale_x_discrete(limits=c("no","yes"),
                   labels=c("Low", "High")) +
  labs(title="Shannon",
       X=NULL,
       y="Shannon diversity")+
  stat_summary(fun=mean,
               geom="point")+
  coord_cartesian()+
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
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  )

print(shannon_plot)

ggsave(
  filename = file.path(plots_dir, "FIG_1_shannon.svg"),
  plot = shannon_plot,
  scale = 1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300,
  device = grDevices::svg
)

#---- Diversity: Richness ------------------------------------------------------

richness <- read.delim(
  paths$richness,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  col.names = c("Sample", "Value"),
  stringsAsFactors = FALSE
)

richness <- richness %>%
  mutate(Sample = as.character(Sample))

merged_rich <- merge(
  richness,
  sample_metadata,
  by = "Sample",
  all = TRUE
) %>%
  separate(
    Sample,
    into = c("Individual", "SampleTime"),
    sep = "_",
    remove = FALSE
)

richness_plot <- ggplot(merged_rich, 
                        aes(x = spike_inhibition, y = Value, 
                            color = spike_inhibition, fill = spike_inhibition)) +
  stat_boxplot(geom='errorbar', color="black", width=0.25, lwd=0.25) + 
  geom_signif(comparisons = list(c("no","yes")), 
              map_signif_level = FALSE, textsize=4, 
              test = "wilcox.test",
              y_position = max(merged_rich$Value, na.rm = TRUE) * 1)+
  geom_boxplot(width=0.5, color="black", outlier.shape = NA , fatten=1, lwd=0.25, alpha=0.25) + 
  geom_point(aes(color = spike_inhibition),size=6, shape=21, alpha=0.5)+
  geom_line(aes(group=Individual), lwd = 0.25)+
  scale_color_manual(name=NULL,
                     values=c("black", "black"),
                     breaks=c("no","yes"),
                     labels=c("< 80%", "> 80%"))+
  scale_fill_manual(name=NULL,
                    values=c("orange" ,"purple4"),
                    breaks=c("no","yes"),
                    labels=c("No", "Yes"))+
  scale_x_discrete(limits=c("no","yes"),
                   labels=c("< 80%", "> 80%")) +
  labs(title = NULL,
       x = "Spike Inhibition",
       y = "Species Observed") +
  stat_summary(fun=mean,
               geom="point")+
  coord_cartesian()+
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
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  )

print(richness_plot)

ggsave(
  filename = file.path(plots_dir, "FIG_1_richness.svg"),
  plot = richness_plot,
  scale = 1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300,
  device = grDevices::svg
)

#----Diversity: Gini -----------------------------------------------------------

gini <- read.delim(
  paths$gini,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  col.names = c("Sample", "Value"),
  stringsAsFactors = FALSE
)

gini <- gini %>%
  mutate(Sample = as.character(Sample))

merged_gini <- merge(
  gini,
  sample_metadata,
  by = "Sample",
  all = TRUE
) %>%
  separate(
    Sample,
    into = c("Individual", "SampleTime"),
    sep = "_",
    remove = FALSE
  )

ginidominance_plot <- ggplot(merged_gini, 
                             aes(x = spike_inhibition, y = Value, 
                                 color = spike_inhibition, fill = spike_inhibition)) +
  stat_boxplot(geom='errorbar', color="black", width=0.25, lwd=0.25) + 
  geom_signif(comparisons = list(c("no","yes")), 
              map_signif_level = FALSE, textsize=4, 
              test = "wilcox.test",
              y_position = max(merged_gini$Value, na.rm = TRUE) * 1)+
  geom_boxplot(width=0.5, color="black", outlier.shape = NA , fatten=1, lwd=0.25, alpha=0.25) + 
  geom_point(aes(color = spike_inhibition),size=6, shape=21, alpha=0.5)+
  geom_line(aes(group=Individual), lwd = 0.25)+
  scale_color_manual(name=NULL,
                     values=c("black", "black"),
                     breaks=c("no","yes"),
                     labels=c("< 80%", "> 80%"))+
  scale_fill_manual(name=NULL,
                    values=c("orange" ,"purple4"),
                    breaks=c("no","yes"),
                    labels=c("No", "Yes"))+
  scale_x_discrete(limits=c("no","yes"),
                   labels=c("< 80%", "> 80%")) +
  labs(title = NULL,
       x = "Spike Inhibition",
       y = "Gini Dominance") +
  stat_summary(fun=mean,
               geom="point")+
  coord_cartesian()+
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
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  )

print(ginidominance_plot)

ggsave(
  filename = file.path(plots_dir, "FIG_1_ginidominance.svg"),
  plot = ginidominance_plot,
  scale = 1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300,
  device = grDevices::svg
)

# ==============================================================================
# Figure 1C: Combined Alpha Diversity plots
# ==============================================================================

# Create combined dataset for all three diversity metrics
# Prepare Shannon data
shannon_combined <- merged_meta %>%
  select(Individual, SampleTime, spike_inhibition, Value) %>%
  mutate(Analysis = "Shannon Diversity")

# Prepare richness data
richness_combined <- merged_rich %>%
  select(Individual, SampleTime, spike_inhibition, Value) %>%
  mutate(Analysis = "Species Richness")

# Prepare Gini dominance data
gini_combined <- merged_gini %>%
  select(Individual, SampleTime, spike_inhibition, Value) %>%
  mutate(Analysis = "Gini Dominance")

# Combine all datasets
combined_diversity <- bind_rows(shannon_combined, richness_combined, gini_combined)

combined_diversity$Analysis <- factor(combined_diversity$Analysis, 
                                      levels = c("Shannon Diversity", "Species Richness", "Gini Dominance"))

# Create combined faceted plot
combined_plot <- ggplot(combined_diversity, 
                        aes(x = spike_inhibition, y = Value, 
                            color = spike_inhibition, fill = spike_inhibition)) +
  stat_boxplot(geom = 'errorbar', color = "black", 
               width = 0.25, lwd = 0.25) + 
  geom_signif(comparisons = list(c("no", "yes")), 
              map_signif_level = FALSE, textsize = 3, 
              test = "wilcox.test") +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, 
               fatten = 1, lwd = 0.25, alpha = 0.25) + 
  geom_point(aes(color = spike_inhibition), size = 3, 
             shape = 21, alpha = 0.6) +
  geom_line(aes(group = Individual), lwd = 0.25, alpha = 0.7) +
  scale_color_manual(name = NULL,
                     values = c("black", "black"),
                     breaks = c("no", "yes"),
                     labels = c("Low", "High")) +
  scale_fill_manual(name = NULL,
                    values = c("orange", "purple4"),
                    breaks = c("no", "yes"),
                    labels = c("No", "Yes")) +
  scale_x_discrete(limits = c("no", "yes"),
                   labels = c("Low", "High")) +
  labs(title = NULL,
       x = "Spike Inhibition",
       y = "Metric Value") +
  stat_summary(fun = mean, geom = "point", size = 8, 
               color = "black", shape = 21, alpha = 0.5) +
  facet_wrap(~ Analysis, scales = "free_y", ncol = 3,
             strip.position = "top",
             labeller = as_labeller(c(
               "Shannon Diversity" = "Alpha Diversity\n(Shannon)",
               "Species Richness" = "Richness\n(Species Observed)",
               "Gini Dominance" = "Species Dominance\n(Gini)"
             ))) +
  theme_classic() +
  coord_cartesian() +
  theme(
    axis.line.x = element_line(size = 0.5, colour = "black", linetype = 1),
    axis.line.y = element_line(size = 0.5, colour = "black", linetype = 1),
    axis.ticks = element_line(size = 0.5, color = "black"),
    axis.ticks.length = unit(.25, "cm"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black", size = 0.5),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

print(combined_plot)

ggsave(
  filename = file.path(plots_dir, "FIG_1C.svg"),
  plot = combined_plot,
  scale = 1,
  width = 7,
  height = 4,
  units = "in",
  dpi = 300,
  device = grDevices::svg
)

# ==============================================================================
# Beta Diversity: Bray-Curtis NMDS + PERMANOVA — Figure 1D
# ==============================================================================

bray <- read.delim(paths$bray, header = TRUE, row.names = 1, check.names = FALSE)

nmds_bray <- metaMDS(bray)

nmds_scores <- as.data.frame(nmds_bray$points)
nmds_scores$Sample <- row.names(nmds_scores)

merged_nmds <- merge(nmds_scores, sample_metadata, by = "Sample", all = TRUE) %>%
  separate(Sample, into = c("Individual", "SampleTime"), sep = "_", remove = FALSE) %>%
  mutate(site = as.character(site))
merged_nmds$spike_inhibition <- as.factor(merged_nmds$spike_inhibition)
merged_nmds$Individual <- as.factor(merged_nmds$Individual)
merged_nmds$site <- as.factor(merged_nmds$site)

# PERMANOVA — single variable
permanova_result <- adonis2(bray ~ spike_inhibition, data = merged_nmds)
print(permanova_result)

# PERMANOVA — individual
permanova_individual <- adonis2(bray ~ Individual, data = merged_nmds)
print(permanova_individual)

# Marginal PERMANOVA — all variables
permanova_marginal <- adonis2(
  bray ~ antibiotics + status + site + spike_inhibition,
  data = merged_nmds,
  by = "margin"
)
print(permanova_marginal)

# Centroids
centroids <- merged_nmds %>%
  group_by(spike_inhibition) %>%
  summarise(MDS1 = mean(MDS1, na.rm = TRUE),
            MDS2 = mean(MDS2, na.rm = TRUE))

# Extract R² and P for Figure 1D annotation
r2_spike       <- round(permanova_result$R2[1], 4)
p_spike        <- round(permanova_result$`Pr(>F)`[1], 3)
r2_spike_marg  <- round(permanova_marginal["spike_inhibition", "R2"], 4)
p_spike_marg   <- round(permanova_marginal["spike_inhibition", "Pr(>F)"], 3)

##### Plot for Publication - Figure 1D ######

nmds_plot <- ggplot(merged_nmds, aes(x = MDS1, y = MDS2, color = spike_inhibition)) +
  geom_point(size = 3, shape = 16, alpha = 0.95) +
  geom_point(data = centroids, aes(x = MDS1, y = MDS2, color = spike_inhibition),
             size = 8, shape = 16, alpha = 0.55) +
  scale_color_manual(name = NULL,
                     values = c("orange", "purple4"),
                     breaks = c("no", "yes"),
                     labels = c("Low (<80%)", "High (\u226580%)")) +
  labs(x = "NMDS1", y = "NMDS2", title = NULL) +
  scale_y_continuous(limits = c(-0.8, 0.8)) +
  scale_x_continuous(limits = c(-0.8, 0.8)) +
  stat_ellipse(type = "t", level = 0.95) +
  annotate("text", x = -0.8, y = 0.75,
           label = paste0("PERMANOVA\nR\u00B2 = ", r2_spike, "\nP = ", p_spike),
           size = 3, hjust = 0) +
  annotate("text", x = 0.8, y = 0.75,
           label = paste0("Marginal PERMANOVA\nR\u00B2 = ", r2_spike_marg, "\nP = ", p_spike_marg),
           size = 3, hjust = 1) +
  theme_classic() +
  coord_cartesian() +
  theme(
    axis.line.x = element_line(size = 0.5, colour = "black", linetype = 1),
    axis.line.y = element_line(size = 0.5, colour = "black", linetype = 1),
    axis.ticks = element_line(size = 0.5, color = "black"),
    axis.ticks.length = unit(.25, "cm"),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

print(nmds_plot)

ggsave(
  filename = file.path(plots_dir, "FIG_1D_nmds.svg"),
  plot = nmds_plot,
  scale = 1, width = 5, height = 4, units = "in", dpi = 300,
  device = grDevices::svg
)

# ==============================================================================
# Supplemental Figure 3A: NMDS by Antibiotics, Status, Site
# ==============================================================================

# library(patchwork)

# Shared theme for supplemental panels
nmds_theme <- theme_classic() +
  theme(
    axis.line.x = element_line(size = 0.5, colour = "black", linetype = 1),
    axis.line.y = element_line(size = 0.5, colour = "black", linetype = 1),
    axis.ticks = element_line(size = 0.5, color = "black"),
    axis.ticks.length = unit(.25, "cm"),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Extract marginal R² and P values
r2_abx_multi   <- round(permanova_marginal["antibiotics", "R2"], 4)
p_abx_multi    <- round(permanova_marginal["antibiotics", "Pr(>F)"], 3)
r2_inf_multi   <- round(permanova_marginal["status", "R2"], 4)
p_inf_multi    <- round(permanova_marginal["status", "Pr(>F)"], 3)
r2_site_multi  <- round(permanova_marginal["site", "R2"], 4)
p_site_multi   <- round(permanova_marginal["site", "Pr(>F)"], 3)

# --- Panel: Antibiotics ---
centroids_abx <- merged_nmds %>%
  group_by(antibiotics) %>%
  summarise(MDS1 = mean(MDS1, na.rm = TRUE),
            MDS2 = mean(MDS2, na.rm = TRUE))

supp_abx <- ggplot(merged_nmds, aes(x = MDS1, y = MDS2, color = antibiotics)) +
  geom_point(size = 3, shape = 16, alpha = 0.95) +
  geom_point(data = centroids_abx, aes(x = MDS1, y = MDS2, color = antibiotics),
             size = 8, shape = 16, alpha = 0.55) +
  scale_color_manual(name = "Recent Antibiotics",
                     values = c("darkkhaki", "brown"),
                     breaks = c("no", "yes"),
                     labels = c("No", "Yes")) +
  labs(x = "NMDS1", y = "NMDS2", title = "Recent Antibiotic Use") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  stat_ellipse(type = "t", level = 0.95) +
  annotate("text", x = -1, y = 0.9,
           label = paste0("Marginal PERMANOVA\nR\u00B2 = ", r2_abx_multi, "\nP = ", p_abx_multi),
           size = 3, hjust = 0) +
  coord_cartesian() +
  nmds_theme +
  theme(legend.position = "bottom")

# --- Panel: Infection Status ---
centroids_infection <- merged_nmds %>%
  group_by(status) %>%
  summarise(MDS1 = mean(MDS1, na.rm = TRUE),
            MDS2 = mean(MDS2, na.rm = TRUE))

supp_infection <- ggplot(merged_nmds, aes(x = MDS1, y = MDS2, color = status)) +
  geom_point(size = 3, shape = 16, alpha = 0.95) +
  geom_point(data = centroids_infection, aes(x = MDS1, y = MDS2, color = status),
             size = 8, shape = 16, alpha = 0.55) +
  scale_color_manual(name = "Prior Infection",
                     values = c("forestgreen", "orange"),
                     breaks = c("uninfected", "infected"),
                     labels = c("No", "Yes")) +
  labs(x = "NMDS1", y = "NMDS2", title = "Prior SARS-CoV-2 Infection") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  stat_ellipse(type = "t", level = 0.95) +
  annotate("text", x = -1, y = 0.9,
           label = paste0("Marginal PERMANOVA\nR\u00B2 = ", r2_inf_multi, "\nP = ", p_inf_multi),
           size = 3, hjust = 0) +
  coord_cartesian() +
  nmds_theme +
  theme(legend.position = "bottom")

# --- Panel: Site ---
centroids_site <- merged_nmds %>%
  mutate(site = factor(site)) %>%
  group_by(site) %>%
  summarise(MDS1 = mean(MDS1, na.rm = TRUE),
            MDS2 = mean(MDS2, na.rm = TRUE))

supp_site <- ggplot(merged_nmds, aes(x = MDS1, y = MDS2, color = factor(site))) +
  geom_point(size = 3, shape = 16, alpha = 0.95) +
  geom_point(data = centroids_site, aes(x = MDS1, y = MDS2, color = factor(site)),
             size = 8, shape = 16, alpha = 0.55) +
  scale_color_manual(name = "NH Location",
                     values = c("navy", "skyblue"),
                     breaks = c("1", "2"),
                     labels = c("NH1", "NH2")) +
  labs(x = "NMDS1", y = "NMDS2", title = "Nursing Home Location") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  stat_ellipse(type = "t", level = 0.95) +
  annotate("text", x = -1, y = 0.9,
           label = paste0("Marginal PERMANOVA\nR\u00B2 = ", r2_site_multi, "\nP = ", p_site_multi),
           size = 3, hjust = 0) +
  coord_cartesian() +
  nmds_theme +
  theme(legend.position = "bottom")

# --- Combine Supplemental Figure 3A ---
supp_fig3a <- supp_abx + supp_infection + supp_site +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

print(supp_fig3a)

ggsave(
  filename = file.path(plots_dir, "SUPP_FIG_3A_nmds_categorical.svg"),
  plot = supp_fig3a,
  scale = 1.5, width = 12, height = 4, units = "in", dpi = 300,
  device = grDevices::svg
)

# ==============================================================================
# Supplemental Figure 2A: Species Heatmap
# ==============================================================================

merged_species <- read.delim(
  paths$species,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

sample_metadata_hm <- read.delim(
  paths$metadata,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

sample_metadata_hm <- sample_metadata_hm[colnames(merged_species), , drop = FALSE]
if (any(is.na(rownames(sample_metadata_hm)))) {
  stop("Metadata alignment produced NA rownames; species column names and metadata rownames still don't match.")
}

species_mat <- as.matrix(merged_species)
species_mat[species_mat > 0] <- log10(species_mat[species_mat > 0])

species_mat_scaled <- t(scale(t(species_mat)))

species_mat_scaled_t <- t(species_mat_scaled)

sample_labels <- rownames(species_mat_scaled_t)
sample_labels <- sub("_A$", "_S1", sample_labels)
sample_labels <- sub("_B$", "_S2", sample_labels)
sample_labels <- sub("_C$", "_S3", sample_labels)
rownames(species_mat_scaled_t) <- sample_labels

response_labels <- c(
  yes = "High Responder (>=80% Spike Inhibition)",
  no  = "Low Responder (<80% Spike Inhibition)"
)

sample_metadata_hm$response <- factor(
  response_labels[as.character(sample_metadata_hm$spike_inhibition)],
  levels = response_labels
)

ha_right <- HeatmapAnnotation(
  Response = sample_metadata_hm$response,
  col = list(Response = c(
    "High Responder (>=80% Spike Inhibition)" = "purple4",
    "Low Responder (<80% Spike Inhibition)"  = "orange"
  )),
  which = "row",
  show_annotation_name = FALSE,
  annotation_legend_param = list(Response = list(title = ""))
)

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}

svglite::svglite(file.path(plots_dir, "SUPP_FIG_2A.svg"), 
                 width = 8, height = 5)

ht <- Heatmap(
  species_mat_scaled_t,
  name                        = "Scaled Abundance\n(Log-Transformed)",
  # clustering
  cluster_rows                = TRUE,
  cluster_columns             = TRUE,
  clustering_distance_rows    = "spearman",
  clustering_method_rows      = "complete",
  # splits (adjust or remove as desired)
  row_split                   = 10,
  # labels & orientation
  show_row_names              = TRUE,   # de-identified sample IDs are the rownames now
  row_names_side              = "right",
  row_names_gp                = gpar(fontsize = 12),
  row_dend_side               = "left",
  show_column_names           = FALSE,  # species names hidden (often too dense)
  # titles
  row_title = "Samples",
  column_title = "Species",
  # appearance
  border                      = TRUE,
  # annotation (right side, next to sample labels)
  right_annotation            = ha_right,
  heatmap_legend_param        = list(direction = "horizontal")
)

draw(
  ht,
  heatmap_legend_side    = "bottom",
  annotation_legend_side = "bottom",
  merge_legend           = TRUE
)

dev.off()

# ==============================================================================
# Supplemental Figure 2B: NMDS by Individual
# ==============================================================================

r2_ind <- round(permanova_individual$R2[1], 4)
p_ind  <- round(permanova_individual$`Pr(>F)`[1], 3)

num_ID <- sort(unique(merged_nmds$Individual))
individual_palette <- setNames(dark.colors(length(num_ID)), num_ID)

centroids_individual <- merged_nmds %>%
  group_by(Individual) %>%
  summarise(MDS1 = mean(MDS1, na.rm = TRUE),
            MDS2 = mean(MDS2, na.rm = TRUE))

supp_individual <- ggplot(merged_nmds, aes(x = MDS1, y = MDS2, color = Individual)) +
  geom_point(size = 3, shape = 16, alpha = 0.95) +
  geom_point(data = centroids_individual, aes(x = MDS1, y = MDS2, color = Individual),
             size = 8, shape = 16, alpha = 0.55) +
  scale_color_manual(name = "NH Resident", values = individual_palette) +
  labs(x = "NMDS1", y = "NMDS2", title = NULL) +
  scale_y_continuous(limits = c(-0.8, 0.8)) +
  scale_x_continuous(limits = c(-0.8, 0.8)) +
  annotate("text", x = -0.8, y = 0.6,
           label = paste0("PERMANOVA\nR\u00B2 = ", r2_ind, "\nP = ", p_ind),
           size = 3, hjust = 0) +
  coord_cartesian() +
  nmds_theme +
  theme(legend.position = "bottom")

print(supp_individual)

ggsave(
  filename = file.path(plots_dir, "SUPP_FIG_2B.svg"),
  plot = supp_individual,
  scale = 1, width = 4, height = 4, units = "in", dpi = 300,
  device = grDevices::svg
)

message("02_taxonomy_divercompsity.R complete.")
