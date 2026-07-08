############################################################
## Descriptive supernatant metabolomics analysis
## BMC + adipocyte vs BMC + fibroblast co-culture experiment
##
## Primary biological contrast:
##   (BMC + Adipocytes - Adipocytes) vs
##   (BMC + 3T3-L1 - 3T3-L1)
##
## This is a descriptive analysis for n=3/group.
## P-values are intentionally not used as the main readout.
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(pheatmap)
  library(writexl)
})


## -----------------------------
## User settings
## -----------------------------

lipid_file <- "Input_files/lipid_supernatants.xlsx"
water_file <- "Input_files/watersoluble_supernatants.xlsx"
outdir <- "supernatant_metabolomics_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

quality_cutoff <- 3
pseudocount <- 1e-9

condition_order <- c(
  "DMEM",
  "3T3-L1",
  "Adipocytes",
  "BMC + 3T3-L1",
  "BMC + Adipocytes"
)

focus_groups <- c(
  "Acylcarnitine",
  "Lysophosphatidylcholine",
  "Lysophospholipid / PAF",
  "Sphingomyelin",
  "Plasmalogen",
  "Phosphatidylethanolamine",
  "Diacylglycerol",
  "Fatty acid / Carboxylic acid"
)

## -----------------------------
## Metabolite/lipid grouping
## -----------------------------

get_group_new <- function(x) {
  x_upper <- toupper(x)

  if (grepl("^AC-\\(", x_upper) || x_upper == "CARNITINE") {
    return("Acylcarnitine")
  } else if (grepl("^SM-\\(", x_upper)) {
    return("Sphingomyelin")
  } else if (grepl("^PC-\\(", x_upper)) {
    return("Phosphatidylcholine")
  } else if (grepl("^PEA-\\(", x_upper) || grepl("^PE-\\(", x_upper)) {
    return("Phosphatidylethanolamine")
  } else if (grepl("^PI-\\(", x_upper)) {
    return("Phosphatidylinositol")
  } else if (grepl("^PS-\\(", x_upper)) {
    return("Phosphatidylserine")
  } else if (grepl("^PA-\\(", x_upper)) {
    return("Phosphatidic acid")
  } else if (grepl("^LPA-\\(", x_upper)) {
    return("Lysophosphatidic acid")
  } else if (grepl("^LPEA-\\(", x_upper) || grepl("^LPE-\\(", x_upper)) {
    return("Lysophosphatidylethanolamine")
  } else if (grepl("^LPS-\\(", x_upper)) {
    return("Lysophosphatidylserine")
  } else if (grepl("^LPC-\\(", x_upper)) {
    return("Lysophosphatidylcholine")
  } else if (grepl("^LPI-\\(", x_upper)) {
    return("Lysophosphatidylinositol")
  } else if (grepl("^LYSOPAF-", x_upper)) {
    return("Lysophospholipid / PAF")
  } else if (grepl("^PLASC-\\(", x_upper) || grepl("^PLASEA-\\(", x_upper)) {
    return("Plasmalogen")
  } else if (grepl("^TAG-", x_upper)) {
    return("Triacylglycerol")
  } else if (grepl("^DAG-\\(", x_upper)) {
    return("Diacylglycerol")
  } else if (grepl("^MAG-\\(", x_upper)) {
    return("Monoacylglycerol")
  } else if (grepl("^CA-\\(", x_upper) || grepl("^FA-\\(", x_upper)) {
    return("Fatty acid / Carboxylic acid")
  } else if (x_upper %in% c("CREATINE", "CREATININE", "CYSTINE", "CYSTATHIONINE",
                            "AMINOADIPIC ACID", "GLUARG", "GLYCINE", "ACETYLGLYCINE",
                            "DIMETHYLARGININE", "HYPOTAURINE", "N-AC-MET")) {
    return("Amino acid & derivatives")
  } else if (x_upper %in% c("A-KETOGLUTARATE", "ACONITATE", "PYRUVATE")) {
    return("TCA cycle")
  } else if (x_upper == "ACETOACETATE") {
    return("Ketone body")
  } else if (x_upper %in% c("3-HYDROXYMETHYLGLUTARATE")) {
    return("Organic acid")
  } else if (x_upper %in% c("DEHYDROASCORBATE", "PYRIDOXAL", "MNAM")) {
    return("Vitamin derivative")
  } else if (x_upper %in% c("ALDOHEXOSE", "GLUCOSAMINE", "FUCOSE", "DRIBOSE-P")) {
    return("Carbohydrate / Sugar")
  } else if (x_upper %in% c("CDP-ETHANOLAMINE", "PHOSPHOETHANOLAMINE")) {
    return("Phospho intermediate")
  } else {
    return("Other/Unclassified")
  }
}

classify_chain_length <- function(name) {
  ## Extract the first carbon number from strings such as AC-(16:0)
  carbon <- str_match(name, "\\((\\d+):")[, 2] |> as.numeric()
  case_when(
    is.na(carbon) ~ NA_character_,
    carbon <= 5 ~ "Short-chain",
    carbon <= 12 ~ "Medium-chain",
    carbon <= 20 ~ "Long-chain",
    carbon > 20 ~ "Very-long-chain"
  )
}

## -----------------------------
## Import function
## -----------------------------

read_supernatant_file <- function(file, assay_type) {
  raw <- readxl::read_excel(file, col_names = FALSE)

  h1 <- raw[1, ] |> unlist(use.names = FALSE) |> as.character()
  h2 <- raw[2, ] |> unlist(use.names = FALSE) |> as.character()

  ## Fill merged condition headers to the right.
  h1_fill <- h1
  current <- NA_character_
  for (i in seq_along(h1_fill)) {
    if (!is.na(h1_fill[i]) && h1_fill[i] != "") current <- h1_fill[i]
    h1_fill[i] <- current
  }

  sample_cols <- which(str_detect(h2, "^Sample_"))
  sample_map <- tibble(
    sample = h2[sample_cols],
    condition = h1_fill[sample_cols]
  ) |>
    mutate(
      condition = factor(condition, levels = condition_order),
      replicate = row_number()
    )

  idx_class <- which(h2 == "Class")[1]
  idx_name  <- which(h2 == "Name")[1]
  idx_qual  <- which(h2 == "Qual.")[1]
  idx_method <- which(h2 == "Method")[1]

  dat <- raw[-c(1, 2), ]

  features <- tibble(
    feature_id = seq_len(nrow(dat)),
    assay_type = assay_type,
    original_class = as.character(dat[[idx_class]]),
    name = as.character(dat[[idx_name]]),
    quality = suppressWarnings(as.numeric(dat[[idx_qual]])),
    method = if (!is.na(idx_method)) as.character(dat[[idx_method]]) else NA_character_,
    group = map_chr(name, get_group_new),
    chain_length = classify_chain_length(name)
  )

  values <- dat[, sample_cols]
  names(values) <- h2[sample_cols]

  long <- bind_cols(features, values) |>
    pivot_longer(
      cols = all_of(sample_map$sample),
      names_to = "sample",
      values_to = "value"
    ) |>
    left_join(sample_map, by = "sample") |>
    mutate(
      value = suppressWarnings(as.numeric(value)),
      condition = factor(condition, levels = condition_order)
    )

  long
}

lipid_long <- read_supernatant_file(lipid_file, "Lipid")
water_long <- read_supernatant_file(water_file, "Water-soluble")

## Main analysis objects.
lipid_qc <- lipid_long |>
  filter(!is.na(quality), quality <= quality_cutoff)

water_qc <- water_long |>
  filter(!is.na(quality), quality <= quality_cutoff, original_class != "Stds.")

all_qc <- bind_rows(lipid_qc, water_qc)

## -----------------------------
## Summaries and contrasts
## -----------------------------

make_condition_summary <- function(dat) {
  dat |>
    group_by(assay_type, feature_id, original_class, name, quality, method, group, chain_length, condition) |>
    summarise(
      n = sum(!is.na(value)),
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      mean_value = if_else(is.nan(mean_value), NA_real_, mean_value),
      median_value = if_else(is.nan(median_value), NA_real_, median_value)
    )
}

make_contrast_table <- function(dat) {
  summary_wide <- make_condition_summary(dat) |>
    select(assay_type, feature_id, original_class, name, quality, method, group, chain_length,
           condition, mean_value, median_value) |>
    pivot_wider(
      names_from = condition,
      values_from = c(mean_value, median_value),
      names_glue = "{.value}__{condition}"
    )

  summary_wide |>
    mutate(
      ## Simple contrasts. Useful but not the main biological contrast.
      log2FC_BMC_Adipo_vs_BMC_3T3 = log2((`median_value__BMC + Adipocytes` + pseudocount) /
                                           (`median_value__BMC + 3T3-L1` + pseudocount)),
      log2FC_Adipo_vs_3T3 = log2((median_value__Adipocytes + pseudocount) /
                                  (`median_value__3T3-L1` + pseudocount)),
      log2FC_BMC_Adipo_vs_Adipo = log2((`median_value__BMC + Adipocytes` + pseudocount) /
                                        (median_value__Adipocytes + pseudocount)),
      log2FC_BMC_3T3_vs_3T3 = log2((`median_value__BMC + 3T3-L1` + pseudocount) /
                                    (`median_value__3T3-L1` + pseudocount)),

      ## Primary feeder-background-adjusted BMC effect.
      bmc_effect_adipocyte_abs = `median_value__BMC + Adipocytes` - median_value__Adipocytes,
      bmc_effect_3T3_abs = `median_value__BMC + 3T3-L1` - `median_value__3T3-L1`,
      delta_bmc_effect_abs = bmc_effect_adipocyte_abs - bmc_effect_3T3_abs,

      ## Log-ratio version of the same concept.
      bmc_effect_adipocyte_log2 = log2((`median_value__BMC + Adipocytes` + pseudocount) /
                                        (median_value__Adipocytes + pseudocount)),
      bmc_effect_3T3_log2 = log2((`median_value__BMC + 3T3-L1` + pseudocount) /
                                  (`median_value__3T3-L1` + pseudocount)),
      delta_bmc_effect_log2 = bmc_effect_adipocyte_log2 - bmc_effect_3T3_log2,

      biological_pattern = case_when(
        log2FC_Adipo_vs_3T3 > 0.5 & bmc_effect_adipocyte_log2 < -0.5 ~
          "Adipocyte-enriched; reduced after BMC co-culture",
        log2FC_Adipo_vs_3T3 > 0.5 & bmc_effect_adipocyte_log2 > 0.5 ~
          "Adipocyte-enriched; further increased after BMC co-culture",
        delta_bmc_effect_log2 > 0.5 ~
          "Stronger BMC-associated accumulation in adipocyte context",
        delta_bmc_effect_log2 < -0.5 ~
          "Stronger BMC-associated depletion in adipocyte context",
        TRUE ~ "No strong descriptive pattern"
      )
    ) |>
    arrange(desc(abs(delta_bmc_effect_log2)))
}

lipid_contrasts <- make_contrast_table(lipid_qc)
water_contrasts <- make_contrast_table(water_qc)

focus_lipid_contrasts <- lipid_contrasts |>
  filter(group %in% focus_groups)

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95", color = NA),
      legend.position = "right"
    )
}

## -----------------------------
## Figure 1: class-level primary BMC effect (Fig. S13)
## -----------------------------

p_class_primary <- focus_lipid_contrasts |>
  mutate(group = fct_reorder(group, delta_bmc_effect_log2, .fun = median, na.rm = TRUE)) |>
  ggplot(aes(x = group, y = delta_bmc_effect_log2)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.6) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 1.7) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Delta BMC effect, log2 scale",
    title = "Feeder-background-adjusted BMC effect across lipid classes",
    subtitle = "Primary contrast: log2[(BMC + adipocytes)/adipocytes] - log2[(BMC + 3T3-L1)/3T3-L1]"
  ) +
  theme_pub()

p_class_primary
#save_plot(p_class_primary, "figure_1_class_level_primary_BMC_effect.pdf", width = 8.5, height = 5.5)

label_df <- focus_lipid_contrasts |>
  group_by(group) |>
  slice_max(delta_bmc_effect_log2, n = 3, with_ties = FALSE) |>
  bind_rows(
    focus_lipid_contrasts |>
      group_by(group) |>
      slice_min(delta_bmc_effect_log2, n = 3, with_ties = FALSE)
  ) |>
  ungroup() |>
  distinct(name, group, .keep_all = TRUE)

p_class_primary <- focus_lipid_contrasts |>
  mutate(
    group = fct_reorder(group, delta_bmc_effect_log2, .fun = median, na.rm = TRUE)
  ) |>
  ggplot(aes(x = group, y = delta_bmc_effect_log2)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.6) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 1.7) +
  geom_text_repel(
    data = label_df |>
      mutate(group = fct_reorder(group, delta_bmc_effect_log2, .fun = median, na.rm = TRUE)),
    aes(label = name),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Delta BMC effect, log2 scale",
    title = "Feeder-background-adjusted BMC effect across lipid classes",
    subtitle = "Primary contrast: log2[(BMC + adipocytes)/adipocytes] - log2[(BMC + 3T3-L1)/3T3-L1]"
  ) +
  theme_pub()

p_class_primary


plot_df <- focus_lipid_contrasts |>
  mutate(
    group = fct_reorder(group, delta_bmc_effect_log2, .fun = median, na.rm = TRUE),
    label_name = gsub(" \\(\\+\\)", "", name)
  )

p_class_primary_labeled_all <- ggplot(
  plot_df,
  aes(x = delta_bmc_effect_log2, y = group)
) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_boxplot(
    aes(group = group),
    orientation = "y",
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.6
  ) +
  geom_jitter(
    height = 0.15,
    width = 0,
    alpha = 0.7,
    size = 1.7
  ) +
  geom_text_repel(
    aes(label = label_name),
    size = 2.4,
    max.overlaps = Inf,
    max.time = 20,
    max.iter = 200000,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    force = 3,
    force_pull = 0.05,
    seed = 123
  ) +
  labs(
    x = "Delta BMC effect, log2 scale",
    y = NULL,
    title = "Feeder-background-adjusted BMC effect across lipid classes",
    subtitle = "Primary contrast: log2[(BMC + adipocytes)/adipocytes] - log2[(BMC + 3T3-L1)/3T3-L1]"
  ) +
  theme_pub() +
  theme(
    plot.margin = ggplot2::margin(10, 120, 10, 10)
  )

p_class_primary_labeled_all

plot_df <- focus_lipid_contrasts |>
  mutate(
    group = fct_reorder(group, delta_bmc_effect_log2, .fun = median, na.rm = TRUE),
    label_name = gsub(" \\(\\+\\)", "", name),
    
    # extract first chain-length number before ":"
    chain_length = stringr::str_match(name, "\\((\\d+):")[,2],
    chain_length = as.numeric(chain_length),
    
    chain_group = case_when(
      is.na(chain_length) ~ "Unknown",
      chain_length <= 18 ~ "≤18",
      chain_length <= 22 ~ "19–22",
      chain_length >= 23 ~ "≥23"
    )
  )

label_df <- plot_df |>
  group_by(group) |>
  slice_max(delta_bmc_effect_log2, n = 2, with_ties = FALSE) |>
  bind_rows(
    plot_df |>
      group_by(group) |>
      slice_min(delta_bmc_effect_log2, n = 2, with_ties = FALSE)
  ) |>
  ungroup() |>
  distinct(name, group, .keep_all = TRUE)

p_class_primary_clean <- ggplot(
  plot_df,
  aes(x = delta_bmc_effect_log2, y = group)
) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_boxplot(
    aes(group = group),
    orientation = "y",
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.45
  ) +
  geom_jitter(
    height = 0.15,
    width = 0,
    alpha = 0.55,
    size = 1.9,
    color = "black"
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = label_name),
    size = 2.7,
    max.overlaps = Inf,
    max.time = 20,
    max.iter = 200000,
    box.padding = 0.7,
    point.padding = 0.6,
    force = 8,
    force_pull = 0.01,
    seed = 123
  ) +
  labs(
    x = "Delta BMC effect, log2 scale",
    y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  theme_pub() +
  theme(
    plot.margin = ggplot2::margin(10, 80, 10, 10)
  )

p_class_primary_clean
ggsave("Figures_Manuscript/p_class_primary_clean.svg", plot = p_class_primary_clean, width = 10, height = 6)

