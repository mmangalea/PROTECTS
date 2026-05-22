# =============================================================================
# Project:  PROTECTS-Microbiome
# Script:   03_MaAsLin2_metabolites.R
# Analyst:  Mihnea (Mike) Mangalea
# Date:     2026-05-18
# Purpose:  MaAsLin2 multivariable association analyses (metabolic pathways,
#           GO terms, species) and metabolite visualizations (heatmap, pathway
#           boxplots, volcano/effect-size plots).
# =============================================================================

source("scripts/00_packages.R")

options(scipen = 999)
set.seed(42)

# ---- Paths -------------------------------------------------------------------

input_dir    <- "data/processed"
metadata_dir <- "data/metadata"
plots_dir    <- "figures/metabolites"
results_dir  <- "results/MaAsLin2"

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  metadata       = file.path(metadata_dir, "protects_metadata_heatmap.txt"),
  pathabund_cpm  = file.path(input_dir, "pathabundance_cpm.tsv"),
  genefamilies   = file.path(input_dir, "genefamilies_GO_rename_cpm.tsv"),
  species        = file.path(input_dir, "merged_species.tsv")
)

# ---- Helper: un-stratify nested pathway tables -------------------------------

unstrat_pathways <- function(dat_path) {
  dat_path[!grepl("\\|", rownames(dat_path)), ]
}

# ---- Helper: plot pathways and test by spike inhibition ----------------------


make_pathway_plot <- function(
    data,
    pathway_ids,
    facet_var = c("feature_description", "feature_name"),
    ncol = 3,
    wrap_width = 24,
    y_multiplier = 0.9
) {
  facet_var <- match.arg(facet_var)
  
  plot_data <- data %>%
    filter(
      feature_name %in% pathway_ids,
      !is.na(spike_inhibition),
      spike_inhibition %in% c("no", "yes")
    ) %>%
    mutate(
      spike_inhibition = factor(spike_inhibition, levels = c("no", "yes"))
    ) %>%
    droplevels()
  
  ggplot(
    plot_data,
    aes(x = spike_inhibition, y = cpm, fill = spike_inhibition)
  ) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(
      shape = 21,
      size = 4,
      alpha = 0.7,
      position = position_jitter(width = 0.2)
    ) +
    geom_signif(
      comparisons = list(c("no", "yes")),
      map_signif_level = FALSE,
      textsize = 4,
      test = "wilcox.test",
      y_position = max(plot_data$cpm, na.rm = TRUE) * y_multiplier
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("no" = "orange", "yes" = "purple4"),
      breaks = c("no", "yes"),
      labels = c("Low Responder", "High Responder")
    ) +
    theme_minimal() +
    labs(title = NULL, x = NULL, y = "Copies per Million") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(fill = "lightblue", colour = "black", linewidth = 1),
      strip.text = element_text(
        size = 9,
        face = "bold",
        lineheight = 0.9,
        margin = margin(t = 4, r = 4, b = 4, l = 4)
      ),
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14)
    ) +
    facet_wrap(
      as.formula(paste("~", facet_var)),
      ncol = ncol,
      labeller = labeller(.default = label_wrap_gen(width = wrap_width))
    )
}

# ==============================================================================
# 1. MaAsLin2: Metabolic Pathway Abundances (CPM)
# ==============================================================================

df_input_protectspathcpm <- read.delim(
  paths$pathabund_cpm,
  row.names = 1,
  stringsAsFactors = FALSE
)

df_input_protectsmetadata <- read.delim(
  paths$metadata,
  row.names = 1,
  stringsAsFactors = FALSE
)

df_input_protectspathcpm <- unstrat_pathways(df_input_protectspathcpm)

write.table(
  df_input_protectspathcpm,
  file = file.path(results_dir, "protects_pathabund_cpm_unstratified.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

fit_pathways <- Maaslin2(
  input_data      = df_input_protectspathcpm,
  input_metadata  = df_input_protectsmetadata,
  output          = file.path(results_dir, "pathways_cpm"),
  min_prevalence  = 0,
  normalization   = "NONE",
  analysis_method = "CPLM",
  fixed_effects   = c("antibiotics", "spike_inhibition", "age", "status", "stay_years"),
  random_effects  = c("subject", "site"),
  reference       = c("spike_inhibition,no")
)

# ==============================================================================
# 2. MaAsLin2: GO Terms (gene families regrouped)
# ==============================================================================

df_input_genefamiliescpm <- read.delim(
  paths$genefamilies,
  row.names = 1,
  stringsAsFactors = FALSE
)

df_input_protectsmetadata <- read.delim(
  paths$metadata,
  row.names = 1,
  stringsAsFactors = FALSE
)

df_input_genefamiliescpm <- unstrat_pathways(df_input_genefamiliescpm)

fit_go <- Maaslin2(
  input_data      = df_input_genefamiliescpm,
  input_metadata  = df_input_protectsmetadata,
  output          = file.path(results_dir, "goterms_unstrat"),
  min_prevalence  = 0.99,
  normalization   = "NONE",
  analysis_method = "CPLM",
  fixed_effects   = c("antibiotics", "spike_inhibition", "age", "status", "stay_years"),
  random_effects  = c("subject", "site"),
  reference       = c("spike_inhibition,no")
)

# ==============================================================================
# 3. MaAsLin2: Species
# ==============================================================================

df_input_species <- read.delim(
  paths$species,
  row.names = 1,
  stringsAsFactors = FALSE
)

df_input_protectsmetadata <- read.delim(
  paths$metadata,
  row.names = 1,
  stringsAsFactors = FALSE
)

fit_species <- Maaslin2(
  input_data      = df_input_species,
  input_metadata  = df_input_protectsmetadata,
  min_prevalence  = 0,
  normalization   = "TSS",
  analysis_method = "LM",
  output          = file.path(results_dir, "species"),
  fixed_effects   = c("antibiotics", "spike_inhibition", "age", "status", "stay_years"),
  random_effects  = c("subject", "site"),
  reference       = c("spike_inhibition,no")
)

# ==============================================================================
# 4. Pathway Effect Size Plot (significant features)
# ==============================================================================

# Read MaAsLin2 all_results from pathway analysis
all_features <- read.delim(
  file.path(results_dir, "pathways_cpm", "all_results.tsv"),
  stringsAsFactors = FALSE
)

metadata_colors <- c(
  ">80% Spike Inhibition" = "purple4",
  "Recent Antibiotics"    = "brown",
  "Prior Infection"       = "forestgreen",
  "Age"                   = "blue",
  "Stay Years"            = "orange"
)

sig_features <- all_features %>%
  filter(qval < 0.05) %>%
  mutate(metadata = str_replace_all(metadata,
                                    c("status" = "Prior Infection",
                                      "antibiotics" = "Recent Antibiotics",
                                      "spike_inhibition" = ">80% Spike Inhibition",
                                      "stay_years" = "Stay Years",
                                      "age" = "Age"))) %>%
  arrange(coef) %>%
  mutate(neg_log_qval = -log10(qval))

effectsize_plot <- ggplot(sig_features, aes(x = coef, y = feature, fill = metadata)) +
  geom_point(size = 5, shape = 21, color = "black", stroke = 0.5, alpha = 0.6) +
  geom_errorbar(aes(xmin = coef - stderr, xmax = coef + stderr), width = 0.2) +
  geom_vline(xintercept = 0, colour = "red", linetype = 2) +
  scale_fill_manual(values = metadata_colors) +
  theme_classic() +
  labs(title = NULL, x = "Coefficient", y = NULL, fill = "Metadata") +
  theme(
    axis.line.x = element_line(size = 0.75, colour = "black", linetype = 1),
    axis.line.y = element_line(size = 0.75, colour = "black", linetype = 1),
    axis.ticks = element_line(size = 0.5, color = "black"),
    axis.ticks.length = unit(.25, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )

print(effectsize_plot)

ggsave(
  file.path(plots_dir, "FIG_2C.svg"),
  plot = effectsize_plot,
  scale = 1.7, width = 8, height = 4, units = "in", dpi = 300,
  device = grDevices::svg
)

# ==============================================================================
# 5. Volcano Plot (all pathway features)
# ==============================================================================

all_features_labeled <- all_features %>%
  mutate(metadata = str_replace_all(metadata,
                                    c("status" = "Prior Infection",
                                      "antibiotics" = "Recent Antibiotics",
                                      "spike_inhibition" = ">80% Spike Inhibition",
                                      "stay_years" = "Stay Years",
                                      "age" = "Age")))

keyvals <- rep("black", nrow(all_features_labeled))
names(keyvals) <- rep("Unchanged", nrow(all_features_labeled))

keyvals[all_features_labeled$metadata == ">80% Spike Inhibition" & all_features_labeled$qval < 0.05] <- "purple4"
names(keyvals)[all_features_labeled$metadata == ">80% Spike Inhibition" & all_features_labeled$qval < 0.05] <- ">80% Spike Inhibition"

keyvals[all_features_labeled$metadata == "Recent Antibiotics" & all_features_labeled$qval < 0.05] <- "brown"
names(keyvals)[all_features_labeled$metadata == "Recent Antibiotics" & all_features_labeled$qval < 0.05] <- "Recent Antibiotics"

keyvals[all_features_labeled$metadata == "Prior Infection" & all_features_labeled$qval < 0.05] <- "forestgreen"
names(keyvals)[all_features_labeled$metadata == "Prior Infection" & all_features_labeled$qval < 0.05] <- "Prior Infection"

keyvals[all_features_labeled$metadata == "Age" & all_features_labeled$qval < 0.05] <- "blue"
names(keyvals)[all_features_labeled$metadata == "Age" & all_features_labeled$qval < 0.05] <- "Age"

keyvals[all_features_labeled$metadata == "Stay Years" & all_features_labeled$qval < 0.05] <- "orange"
names(keyvals)[all_features_labeled$metadata == "Stay Years" & all_features_labeled$qval < 0.05] <- "Stay Years"

svg(file.path(plots_dir, "FIG_2B.svg"), width = 8, height = 8)

EnhancedVolcano(all_features_labeled,
                lab = NA,
                x = "coef",
                y = "qval",
                ylim = c(0, 3),
                xlim = c(-40, 40),
                pCutoff = 0.05,
                FCcutoff = 0,
                pointSize = 8.0,
                labSize = 3.0,
                colCustom = keyvals,
                xlab = "Coefficient",
                ylab = "-log10(q-value)",
                legendPosition = "top",
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                lengthConnectors = unit(2.0, "mm"),
                colAlpha = 0.6,
                border = "full")

dev.off()

# ==============================================================================
# 6. Metabolite Heatmap (all unstratified pathways, ComplexHeatmap)
# ==============================================================================

sample_metadata_hm <- read.delim(
  paths$metadata,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

collapsed_pathways <- unstrat_pathways(
  read.delim(
    paths$pathabund_cpm,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
)

sample_metadata_hm <- sample_metadata_hm[colnames(collapsed_pathways), ]

metabolites_log <- collapsed_pathways
metabolites_log[metabolites_log > 0] <- log(metabolites_log[metabolites_log > 0], base = 10)
metabolites_log_scaled <- t(scale(t(metabolites_log)))

metabolites_log_scaled <- t(metabolites_log_scaled)

sample_labels <- rownames(metabolites_log_scaled)
sample_labels <- sub("_A$", "_S1", sample_labels)
sample_labels <- sub("_B$", "_S2", sample_labels)
sample_labels <- sub("_C$", "_S3", sample_labels)
rownames(metabolites_log_scaled) <- sample_labels

annotation_colors <- list(
  spike_inhibition = c(
    "Low Responder (<80% Spike Inhibition)" = "orange",
    "High Responder (≥80% Spike Inhibition)" = "purple4"
  )
)

sample_metadata_hm$spike_inhibition <- ifelse(
  sample_metadata_hm$spike_inhibition == "yes",
  "High Responder (≥80% Spike Inhibition)",
  "Low Responder (<80% Spike Inhibition)"
)

sample_metadata_hm$spike_inhibition <- factor(
  sample_metadata_hm$spike_inhibition,
  levels = c(
    "Low Responder (<80% Spike Inhibition)",
    "High Responder (≥80% Spike Inhibition)"
  )
)

row_annotation <- rowAnnotation(
  spike_inhibition = sample_metadata_hm$spike_inhibition,
  col = annotation_colors,
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    spike_inhibition = list(title = NULL)
  )
)

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}

svglite::svglite(file.path(plots_dir, "FIG_2A.svg"), 
                 width = 11, height = 8)

ht <- Heatmap(
  metabolites_log_scaled,
  name = "Expression",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  row_dend_side = "left",
  row_names_side = "right",
  row_split = 4,
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Scaled Log-Transformed Values",
    title_position = "topleft",
    direction = "horizontal"
  ),
  row_title = "Samples",
  column_title = "Metabolites",
  border = TRUE,
  right_annotation = row_annotation,
  heatmap_width = unit(7, "in"),
  heatmap_height = unit(4, "in")
)

draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)

dev.off()
# ==============================================================================
# 7. Key Pathway Boxplots (SCFA, TCA, Biosynthesis, and more)
# ==============================================================================

pathabund_cpm <- unstrat_pathways(df_input_protectspathcpm) %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols = -feature, names_to = "sample", values_to = "cpm") %>%
  mutate(
    # Remove artificial leading X added to numeric-starting sample IDs
    sample = sub("^X(?=\\d)", "", sample, perl = TRUE)
  ) %>%
  separate(
    feature,
    into = c("feature_name", "feature_description"),
    sep = ":",
    extra = "merge",
    fill = "right"
  )

metadata_long <- df_input_protectsmetadata %>%
  rownames_to_column(var = "sample") %>%
  mutate(
    # Apply the same cleanup to metadata sample IDs
    sample = sub("^X(?=\\d)", "", sample, perl = TRUE),
    spike_inhibition = tolower(trimws(spike_inhibition))
  )

merged_data <- pathabund_cpm %>%
  left_join(metadata_long, by = "sample") %>%
  mutate(
    feature_name = gsub("-", "_", feature_name),
    spike_inhibition = factor(spike_inhibition, levels = c("no", "yes"))
  )

# ------------------------------------------------------------------------------
# A. SCFA PATHWAYS (Figure 2D)
# ------------------------------------------------------------------------------

SCFA_main <- c(
  "P161_PWY",         # acetylene degradation (anaerobic) — acetate end-product
  "P124_PWY",         # Bifidobacterium shunt
  "P461_PWY",         # hexitol fermentation to lactate, formate, ethanol and acetate
  "PWY_5022",         # 4-aminobutanoate degradation V
  "PWY_5676",         # acetyl-CoA fermentation to butanoate II
  "CENTFERM_PWY",     # pyruvate fermentation to butanoate
  "P122_PWY",         # heterolactic fermentation
  "ANAEROFRUCAT_PWY", # homolactic fermentation
  "PWY_7013"          # (S)-propane-1,2-diol degradation
)

scfamain_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = SCFA_main,
  facet_var = "feature_description",
  ncol = 9,
  wrap_width = 12,
  y_multiplier = 0.9
)

print(scfamain_plot)

ggsave(
  file.path(plots_dir, "FIG_2D.svg"),
  plot = scfamain_plot,
  scale = 1.5, width = 8, height = 4, units = "in", dpi = 300,
  device = grDevices::svg
)

# ------------------------------------------------------------------------------
# B. BIOSYNTHESIS PATHWAYS
# ------------------------------------------------------------------------------

biosynthesis_pathways <- c(
  "P108_PWY",    # pyruvate fermentation to propanoate I
  "PWY_5723",    # Rubisco shunt
  "PRPP_PWY",    # superpathway of histidine, purine, and pyrimidine biosynthesis
  "PWY_6285",    # superpathway of fatty acids biosynthesis (E. coli)
  "PWY0_1337",   # oleate beta-oxidation
  "PWY0_301"     # L-ascorbate degradation I (bacterial, anaerobic)
)

biosynthesis_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = biosynthesis_pathways,
  facet_var = "feature_description",
  ncol = 6,
  wrap_width = 18,
  y_multiplier = 0.9
)

print(biosynthesis_plot)

ggsave(
  file.path(plots_dir, "SUPP_FIG_4A.svg"),
  plot = biosynthesis_plot,
  scale = 1.5, width = 6, height = 4, units = "in", dpi = 300,
  device = grDevices::svg
)

# ------------------------------------------------------------------------------
# C. SCFA PATHWAYS (extended)
# ------------------------------------------------------------------------------

SCFA_extended <- c(
  "FERMENTATION_PWY", # mixed acid fermentation
  "P41_PWY",          # pyruvate fermentation to acetate and (S)-lactate I
  "PWY_5100",         # pyruvate fermentation to acetate and lactate II
  "PWY_5677",         # succinate fermentation to butanoate
  "PWY_5494",         # pyruvate fermentation to propanoate II (acrylate pathway)
  "P163_PWY",         # L-lysine fermentation to acetate and butanoate
  "PWY_6590",         # superpathway of Clostridium acetobutylicum acidogenic fermentation
  "PWY_7383"          # anaerobic energy metabolism
)

scfaext_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = SCFA_extended,
  facet_var = "feature_description",
  ncol = 8,
  wrap_width = 18,
  y_multiplier = 0.9
)

print(scfaext_plot)

# ------------------------------------------------------------------------------
# D. TCA PATHWAYS
# ------------------------------------------------------------------------------

TCA_pathways <- c(
  "P23_PWY",           # reductive TCA cycle I
  "TCA",               # TCA cycle I (prokaryotic)
  "REDCITCYC",         # TCA cycle VI (Helicobacter)
  "P105_PWY",          # TCA cycle IV (2-oxoglutarate decarboxylase)
  "PWY_5690",          # TCA cycle II (plants and fungi)
  "PWY_6969",          # TCA cycle V (2-oxoglutarate synthase)
  "PWY_7254",          # TCA cycle VII (acetate-producers)
  "P42_PWY",           # incomplete reductive TCA cycle
  "PWY_5392",          # reductive TCA cycle II
  "PWY_5913",          # partial TCA cycle
  "GLYOXYLATE_BYPASS", # glyoxylate cycle
  "TCA_GLYOX_BYPASS"   # superpathway of glyoxylate bypass and TCA
)

tca_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = TCA_pathways,
  facet_var = "feature_description",
  ncol = 6,
  wrap_width = 18,
  y_multiplier = 0.9
)

print(tca_plot)

ggsave(
  file.path(plots_dir, "SUPP_FIG_4B.svg"),
  plot = tca_plot,
  scale = 1.5, width = 6, height = 4, units = "in", dpi = 300,
  device = grDevices::svg
)

# ------------------------------------------------------------------------------
# E. PYRUVATE PATHWAYS
# ------------------------------------------------------------------------------

pyruvate <- c(
  "P41_PWY",                     # pyruvate fermentation to acetate and (S)-lactate I
  "PWY_5100",                    # pyruvate fermentation to acetate and lactate II
  "PWY_6588",                    # pyruvate fermentation to acetone
  "FERMENTATION_PWY",            # mixed acid fermentation
  "PWY_5464",                    # superpathway of cytosolic glycolysis, PDH and TCA
  "GLYCOLYSIS_TCA_GLYOX_BYPASS"  # superpathway of glycolysis, PDH, TCA and glyoxylate bypass
)

pyruvate_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = pyruvate,
  facet_var = "feature_description",
  ncol = 6,
  wrap_width = 20,
  y_multiplier = 0.8
)

print(pyruvate_plot)

# ------------------------------------------------------------------------------
# F. GLYCOLYSIS & GLUCONEOGENESIS PATHWAYS
# ------------------------------------------------------------------------------

glycolysis_gluconeogenesis <- c(
  "GLYCOLYSIS",        # glycolysis I
  "PWY_1042",          # glycolysis IV
  "ANAGLYCOLYSIS_PWY", # glycolysis III
  "PWY_5484",          # glycolysis II
  "GLYCOLYSIS_E_D",    # superpathway of glycolysis and Entner-Doudoroff pathway
  "PWY_8004",          # Entner-Doudoroff pathway I
  "GLUCONEO_PWY",      # gluconeogenesis I
  "PWY66_399"          # gluconeogenesis III
)

glygluc_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = glycolysis_gluconeogenesis,
  facet_var = "feature_description",
  ncol = 4,
  wrap_width = 18,
  y_multiplier = 0.9
)

print(glygluc_plot)

# ------------------------------------------------------------------------------
# G. AMINO ACID BIOSYNTHESIS PATHWAYS
# ------------------------------------------------------------------------------

aa_biosynthesis <- c(
  "ARGSYN_PWY",       # L-arginine biosynthesis I
  "HISTSYN_PWY",      # L-histidine biosynthesis
  "VALSYN_PWY",       # L-valine biosynthesis
  "ILEUSYN_PWY",      # L-isoleucine biosynthesis I
  "TRPSYN_PWY",       # L-tryptophan biosynthesis
  "DAPLYSINESYN_PWY", # L-lysine biosynthesis I
  "THRESYN_PWY",      # superpathway of L-threonine biosynthesis
  "GLUTORN_PWY"       # L-ornithine biosynthesis I
)

aabio_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = aa_biosynthesis,
  facet_var = "feature_description",
  ncol = 4,
  wrap_width = 18,
  y_multiplier = 0.90
)

print(aabio_plot)

# ------------------------------------------------------------------------------
# H. METHANOGENESIS PATHWAYS
# ------------------------------------------------------------------------------

methanogenesis <- c(
  "METH_ACETATE_PWY",  # methanogenesis from acetate
  "METHANOGENESIS_PWY" # methanogenesis from H2 and CO2
)

methano_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = methanogenesis,
  facet_var = "feature_description",
  ncol = 2,
  wrap_width = 22,
  y_multiplier = 0.9
)

print(methano_plot)

# ------------------------------------------------------------------------------
# I. FERMENTATION, BROADER CONTEXT
# ------------------------------------------------------------------------------

fermentation_broad <- c(
  "PWY_6590",         # Clostridium acetobutylicum acidogenic fermentation
  "PWY_6588",         # pyruvate fermentation to acetone
  "P125_PWY",         # superpathway of (R,R)-butanediol biosynthesis
  "PWY_6396",         # superpathway of 2,3-butanediol biosynthesis
  "PWY4LZ_257",       # superpathway of fermentation
  "FERMENTATION_PWY"  # mixed acid fermentation
)

fermbroad_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = fermentation_broad,
  facet_var = "feature_description",
  ncol = 6,
  wrap_width = 18,
  y_multiplier = 0.9
)

print(fermbroad_plot)

# ------------------------------------------------------------------------------
# J. LIPID METABOLISM
# ------------------------------------------------------------------------------

lipid_biosynthesis <- c(
  "FASYN_ELONG_PWY",   # fatty acid elongation, saturated
  "FASYN_INITIAL_PWY", # fatty acid biosynthesis initiation
  "PWY_5971",          # palmitate biosynthesis
  "PWY_5989",          # stearate biosynthesis II
  "PWY_6282",          # palmitoleate biosynthesis I
  "PWY_6284",          # unsaturated fatty acids biosynthesis
  "PWY_6285",          # fatty acids biosynthesis
  "PWY_7663",          # gondoate biosynthesis
  "PWY_7664"           # oleate biosynthesis IV
)

lipid_oxidation <- c(
  "FAO_PWY",   # fatty acid beta-oxidation I
  "PWY_5136",  # fatty acid beta-oxidation II
  "PWY_5138",  # fatty acid beta-oxidation IV
  "PWY0_1337"  # oleate beta-oxidation
)

lipid_metabolism <- c(lipid_biosynthesis, lipid_oxidation)

lipid_plot <- make_pathway_plot(
  data = merged_data,
  pathway_ids = lipid_metabolism,
  facet_var = "feature_description",
  ncol = 7,
  wrap_width = 10,
  y_multiplier = 0.8
)

print(lipid_plot)

message("03_MaAsLin2_metabolites.R complete.")
