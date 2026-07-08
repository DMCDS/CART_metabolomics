## Reanalysis of invitro experiments

############################################################
# In vitro CAR-T / metabolite co-culture analysis
# Author: David / ChatGPT
# Input: invitro_combined.xlsx
############################################################

# =========================
# 0. Packages
# =========================

packages <- c(
  "readxl", "dplyr", "tidyr", "stringr", "purrr",
  "ggplot2", "emmeans", "broom", "forcats", "readr"
)

installed <- rownames(installed.packages())
for (p in packages) {
  if (!p %in% installed) install.packages(p)
}

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(emmeans)
library(broom)
library(forcats)
library(readr)

# =========================
# 1. User settings
# =========================

input_file <- "input_files/invitro_combined.xlsx"
sheet_name <- "combined_long_for_R_updated_v3"

outdir <- "invitro_analysis_results"
dir.create(outdir, showWarnings = FALSE)

# Conditions to exclude globally
exclude_analytes <- c("Lyso16:1", "Palmitate")

# AC-focused analytes.
# I treat ALCAR and Combi as AC-related conditions.
# Remove "Combi" here if you want only individual AC species.
ac_analytes <- c("AC10:0", "AC14:0", "AC18:1", "ALCAR")

control_analytes <- c("HPLM", "DMSO")

cart_groups <- c("CAR 74", "CAR 77")

# =========================
# 2. Load and clean data
# =========================

df_raw <- read_excel(input_file, sheet = sheet_name)
df <- df_raw %>%
  mutate(
    screenshot = as.character(screenshot),
    concentration = as.character(concentration),
    analyte = as.character(analyte),
    analyte = recode(analyte, "ALCAR" = "AC2:0"),
    group = as.character(group),
    sample = as.character(sample),
    value = as.numeric(value)
  ) %>%
  filter(!analyte %in% exclude_analytes) %>%
  mutate(
    concentration = factor(concentration, levels = c("phys", "supra")),
    group = factor(group),
    analyte = factor(analyte),
    treatment_class = case_when(
      analyte %in% ac_analytes ~ "AC",
      analyte %in% control_analytes ~ as.character(analyte),
      TRUE ~ "Other"
    ),
    treatment_class = factor(
      treatment_class,
      levels = c("HPLM", "DMSO", "AC", "Other")
    ),
    is_cart = group %in% cart_groups
  )

# Keep only ACs and controls for the main analyses
df_main <- df %>%
  filter(analyte %in% c(ac_analytes, control_analytes))

# Basic QC table
qc_counts <- df_main %>%
  count(screenshot, concentration, analyte, group, name = "n")

#write_csv(qc_counts, file.path(outdir, "qc_counts_by_endpoint_condition_group.csv"))

# =========================
# 3. Split endpoints
# =========================

viab <- df_main %>%
  filter(screenshot == "viability") %>%
  filter(group %in% c(cart_groups, "Target Only", "UT TCells")) %>%
  mutate(
    residual_viability = value,
    killing_score = 100 - value
  )

markers <- df_main %>%
  filter(screenshot != "viability") %>%
  filter(group %in% cart_groups) %>%
  mutate(
    log_value = log2(value + 1)
  )

# Marker annotation
markers <- markers %>%
  mutate(
    cell_type = case_when(
      str_starts(screenshot, "CD4_") ~ "CD4",
      str_starts(screenshot, "CD8_") ~ "CD8",
      TRUE ~ NA_character_
    ),
    marker = str_remove(screenshot, "^CD4_"),
    marker = str_remove(marker, "^CD8_"),
    marker_category = case_when(
      marker %in% c("CD25", "CD69", "PD1") ~ "activation",
      marker %in% c("GSH", "cROS", "mSOX", "MitoRed", "MitoGreen") ~ "metabolic",
      TRUE ~ "other"
    )
  )

#write_csv(viab, file.path(outdir, "viability_cleaned_for_analysis.csv"))
#write_csv(markers, file.path(outdir, "markers_cleaned_for_analysis.csv"))

# =========================
# 4. Summary statistics
# =========================

viability_summary <- viab %>%
  filter(group %in% cart_groups) %>%
  group_by(concentration, analyte, treatment_class, group) %>%
  summarise(
    n = n(),
    mean_viability = mean(residual_viability, na.rm = TRUE),
    sd_viability = sd(residual_viability, na.rm = TRUE),
    median_viability = median(residual_viability, na.rm = TRUE),
    iqr_viability = IQR(residual_viability, na.rm = TRUE),
    mean_killing_score = mean(killing_score, na.rm = TRUE),
    sd_killing_score = sd(killing_score, na.rm = TRUE),
    .groups = "drop"
  )

marker_summary <- markers %>%
  group_by(screenshot, cell_type, marker, marker_category, concentration, analyte, treatment_class, group) %>%
  summarise(
    n = n(),
    mean_mfi = mean(value, na.rm = TRUE),
    sd_mfi = sd(value, na.rm = TRUE),
    median_mfi = median(value, na.rm = TRUE),
    iqr_mfi = IQR(value, na.rm = TRUE),
    mean_log2_mfi = mean(log_value, na.rm = TRUE),
    sd_log2_mfi = sd(log_value, na.rm = TRUE),
    .groups = "drop"
  )

#write_csv(viability_summary, file.path(outdir, "viability_summary.csv"))
#write_csv(marker_summary, file.path(outdir, "marker_summary.csv"))

# =========================
# 5. Viability / killing analysis
# =========================
# Main question:
# Do ACs improve killing compared with HPLM or DMSO?
#
# Since lower residual viability = better killing,
# negative estimates for AC - control indicate improved killing.

viab_cart <- viab %>%
  filter(group %in% cart_groups) %>%
  filter(treatment_class %in% c("HPLM", "DMSO", "AC")) %>%
  droplevels()

# 5A. Overall model across CAR constructs
# Adjusts for construct and concentration.
fit_viab_overall <- lm(
  residual_viability ~ treatment_class + group + concentration,
  data = viab_cart
)

emm_viab_overall <- emmeans(fit_viab_overall, ~ treatment_class)

contrast_viab_overall <- contrast(
  emm_viab_overall,
  method = list(
    "AC vs HPLM" = c(-1, 0, 1),
    "AC vs DMSO" = c(0, -1, 1)
  ),
  adjust = "BH"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(
    endpoint = "viability",
    model = "overall_adjusted_for_construct_and_concentration",
    interpretation = ifelse(estimate < 0, "lower viability / better killing with AC", "higher viability / less killing with AC")
  )

#write_csv(
#  contrast_viab_overall,
#  file.path(outdir, "viability_AC_vs_controls_overall.csv")
#)

# 5B. Interaction model: does AC effect differ by construct?
fit_viab_interaction <- lm(
  residual_viability ~ treatment_class * group + concentration,
  data = viab_cart
)

anova_viab_interaction <- anova(fit_viab_interaction) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term")

#write_csv(
#  anova_viab_interaction,
#  file.path(outdir, "viability_interaction_model_anova.csv")
#)

# 5C. AC vs controls within each CAR-T construct
emm_viab_by_construct <- emmeans(
  fit_viab_interaction,
  ~ treatment_class | group
)

contrast_viab_by_construct <- contrast(
  emm_viab_by_construct,
  method = list(
    "AC vs HPLM" = c(-1, 0, 1),
    "AC vs DMSO" = c(0, -1, 1)
  ),
  adjust = "BH"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(
    endpoint = "viability",
    model = "stratified_by_construct",
    interpretation = ifelse(estimate < 0, "lower viability / better killing with AC", "higher viability / less killing with AC")
  )

#write_csv(
#  contrast_viab_by_construct,
#  file.path(outdir, "viability_AC_vs_controls_by_construct.csv")
#)

# 5D. Individual AC species vs HPLM/DMSO
viab_species <- viab %>%
  filter(group %in% c("CAR 74", "CAR 77", "UT TCells")) %>%
  filter(analyte %in% c(ac_analytes, control_analytes)) %>%
  mutate(
    analyte = factor(analyte, levels = c("HPLM", "DMSO", ac_analytes))
  ) %>%
  droplevels()

fit_viab_species <- lm(
  residual_viability ~ analyte * group + concentration,
  data = viab_species
)

fit_viab_species2 <- lm(
  residual_viability ~ analyte * group,
  data = viab_species |> filter(concentration == "supra")
)
summary(fit_viab_species2)

emm_viab_species_by_construct <- emmeans(
  fit_viab_species,
  ~ analyte | group
)

#### Plotting viability

# Define conditions to keep
plot_analytes <- c("Untransduced", "DMSO", "AC10:0", "AC14:0", "AC18:1", "ALCAR")

# Prepare pooled plotting data
plot_viab_pooled <- viab %>%
  filter(
    group %in% c("UT TCells", "CAR 74", "CAR 77"),
    analyte %in% c("HPLM", "DMSO", "AC10:0", "AC14:0", "AC18:1", "ALCAR")
  ) %>%
  mutate(
    analyte = if_else(group == "UT TCells", "Untransduced", as.character(analyte)),
    analyte = factor(analyte, levels = plot_analytes),
    concentration = factor(concentration, levels = c("phys", "supra")),
    group = factor(group, levels = c("UT TCells", "CAR 74", "CAR 77"))
  )

plot_phys <- plot_viab_pooled |>
  filter(
    analyte == "HPLM" |
      analyte == "DMSO" |
      (analyte %in% c("AC10:0", "AC14:0", "AC18:1", "ALCAR") &
         concentration == "phys")
  ) |>
  mutate(
    plot_concentration = "phys",
    analyte = factor(
      analyte,
      levels = c("HPLM", "DMSO", "AC10:0", "AC14:0", "AC18:1", "ALCAR")
    )
  )

p_phys <- ggplot(plot_phys, aes(x = analyte, y = residual_viability)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6
  ) +
  geom_jitter(
    aes(shape = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(ref.group = 'DMSO', method = 't.test', label = 'p.format')+
  facet_wrap(~group)+
  theme_bw(base_size = 13) +
  labs(
    title = "Residual viability under physiologic AC concentrations",
    subtitle = "CAR 74 and CAR 77 pooled; individual points show independent replicates",
    x = NULL,
    y = "Residual viability",
    shape = "CAR-T construct"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_phys

## supra
plot_supra <- plot_viab_pooled %>%
  filter(concentration == "supra")

p_supra <- ggplot(plot_supra, aes(x = analyte, y = residual_viability)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6
  ) +
  geom_jitter(
    aes(shape = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(ref.group = 'DMSO', method = "t.test", label = 'p.format')+
  facet_wrap(~group)+
  theme_bw(base_size = 13) +
  labs(
    title = "Residual viability under supra-physiologic AC concentrations",
    subtitle = "CAR 74 and CAR 77 pooled; individual points show independent replicates",
    x = NULL,
    y = "Residual viability",
    shape = "CAR-T construct"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_supra

wilcox.test(residual_viability ~ analyte, plot_supra |> filter(analyte %in% c('AC14:0', 'DMSO'), group == 'CAR 74'))

#### Comparing constructs
plot_phys_construct <- plot_viab_pooled %>%
  filter(
    analyte %in% c("Untransduced", "DMSO") |
      (analyte %in% c("AC10:0", "AC14:0", "AC18:1", "ALCAR") &
         concentration == "phys")
  ) %>%
  mutate(
    plot_concentration = "phys",
    group = case_when(
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("41BB costim", "CD28z costim"))
  )

p_phys_construct <- ggplot(
  plot_phys_construct,
  aes(x = analyte, y = residual_viability, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(shape = group, color = group),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    aes(group = group),
    method = "wilcoxon",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    tip.length = 0.01,
    hide.ns = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "41BB costim" = "grey60",
      "CD28z costim" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.20))
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Residual viability",
    fill = "CAR-T construct",
    shape = "CAR-T construct",
    color = "CAR-T construct"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

p_phys_construct


plot_supra_construct <- plot_viab_pooled %>%
  filter(
    analyte == "DMSO" |
      (analyte %in% c("AC10:0", "AC14:0", "AC18:1", "ALCAR") &
         concentration == "supra")
  ) %>%
  mutate(
    plot_concentration = "supra",
    analyte = factor(
      analyte,
      levels = c("DMSO", "AC10:0", "AC14:0", "AC18:1", "ALCAR")
    ),
    group = case_when(
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("41BB costim", "CD28z costim"))
  )

p_supra_construct <- ggplot(
  plot_supra_construct,
  aes(x = analyte, y = residual_viability, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(shape = group, color = group),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    aes(group = group),
    method = "wilcoxon",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    tip.length = 0.01,
    hide.ns = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "41BB costim" = "grey60",
      "CD28z costim" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.20))
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Residual viability",
    fill = "CAR-T construct",
    shape = "CAR-T construct",
    color = "CAR-T construct"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

p_supra_construct



# Conditions to show on x-axis
plot_analytes <- c(
  "Untransduced",
  "DMSO",
  "AC10:0",
  "AC14:0",
  "AC18:1",
  "ALCAR"
)

ac_plot_analytes <- c("AC10:0", "AC14:0", "AC18:1", "ALCAR")
cart_groups <- c("CAR 74", "CAR 77")

# Make sure group names do not contain accidental spaces
viab_plot_base <- viab %>%
  mutate(
    group = trimws(as.character(group)),
    analyte = as.character(analyte),
    concentration = as.character(concentration)
  ) %>%
  filter(
    group %in% c("UT TCells", "CAR 74", "CAR 77")
  )

# ----------------------------------------------------------
# Helper function:
# Build plotting dataset for either phys or supra
# ----------------------------------------------------------

make_viability_plot_data <- function(conc_use) {
  
  # Generic negative control:
  # UT TCells are only available in HPLM, but we include them
  # in both phys and supra plots as "Untransduced".
  ut_control <- viab_plot_base %>%
    filter(
      group == "UT TCells",
      analyte == "HPLM"
    ) %>%
    mutate(
      analyte = "Untransduced",
      plot_concentration = conc_use
    )
  
  # DMSO control:
  # DMSO may be saved as supra, but should be shown in both plots.
  dmso_control <- viab_plot_base %>%
    filter(
      group %in% cart_groups,
      analyte == "DMSO"
    ) %>%
    mutate(
      plot_concentration = conc_use
    )
  
  # AC species:
  # These are filtered by actual phys or supra concentration.
  ac_data <- viab_plot_base %>%
    filter(
      group %in% cart_groups,
      analyte %in% ac_plot_analytes,
      concentration == conc_use
    ) %>%
    mutate(
      plot_concentration = conc_use
    )
  
  bind_rows(
    ut_control,
    dmso_control,
    ac_data
  ) %>%
    mutate(
      analyte = factor(analyte, levels = plot_analytes),
      group = factor(group, levels = c("UT TCells", "CAR 74", "CAR 77")),
      plot_concentration = factor(plot_concentration, levels = c("phys", "supra"))
    )
}

plot_phys <- make_viability_plot_data("phys")
plot_supra <- make_viability_plot_data("supra")

# Conditions to show on x-axis
plot_analytes <- c(
  "Untransduced",
  "DMSO",
  "AC10:0",
  "AC14:0",
  "AC18:1",
  "ALCAR"
)

ac_plot_analytes <- c("AC10:0", "AC14:0", "AC18:1", "ALCAR")
cart_groups <- c("CAR 74", "CAR 77")

# Make sure group names do not contain accidental spaces
viab_plot_base <- viab %>%
  mutate(
    group = trimws(as.character(group)),
    analyte = as.character(analyte),
    concentration = as.character(concentration)
  ) %>%
  filter(
    group %in% c("UT TCells", "CAR 74", "CAR 77")
  )

# ----------------------------------------------------------
# Helper function:
# Build plotting dataset for either phys or supra
# ----------------------------------------------------------

make_viability_plot_data <- function(conc_use) {
  
  # Generic negative control:
  # UT TCells are only available in HPLM, but we include them
  # in both phys and supra plots as "Untransduced".
  ut_control <- viab_plot_base %>%
    filter(
      group == "UT TCells",
      analyte == "HPLM"
    ) %>%
    mutate(
      analyte = "Untransduced",
      plot_concentration = conc_use
    )
  
  # DMSO control:
  # DMSO may be saved as supra, but should be shown in both plots.
  dmso_control <- viab_plot_base %>%
    filter(
      group %in% cart_groups,
      analyte == "DMSO"
    ) %>%
    mutate(
      plot_concentration = conc_use
    )
  
  # AC species:
  # These are filtered by actual phys or supra concentration.
  ac_data <- viab_plot_base %>%
    filter(
      group %in% cart_groups,
      analyte %in% ac_plot_analytes,
      concentration == conc_use
    ) %>%
    mutate(
      plot_concentration = conc_use
    )
  
  bind_rows(
    ut_control,
    dmso_control,
    ac_data
  ) %>%
    mutate(
      analyte = factor(analyte, levels = plot_analytes),
      group = factor(group, levels = c("UT TCells", "CAR 74", "CAR 77")),
      plot_concentration = factor(plot_concentration, levels = c("phys", "supra"))
    )
}

plot_phys <- make_viability_plot_data("phys")
plot_supra <- make_viability_plot_data("supra")


p_phys <- ggplot(
  plot_phys,
  aes(x = analyte, y = residual_viability)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6
  ) +
  geom_jitter(
    aes(shape = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_phys %>%
      filter(
        group %in% cart_groups,
        analyte != "Untransduced"
      ),
    aes(group = analyte),
    ref.group = "DMSO",
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  facet_wrap(~ group) +
  theme_bw(base_size = 13) +
  labs(
    title = "Residual viability under physiologic AC concentrations",
    subtitle = "Untransduced T cells are shown as a generic negative control; AC species are compared with DMSO",
    x = NULL,
    y = "Residual viability",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_phys


### Supplementary Figure: comparison killing phys
p_phys_supp <- ggplot(
  plot_phys,
  aes(x = analyte, y = residual_viability)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6
  ) +
  geom_jitter(
    aes(shape = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_phys %>%
      filter(
        group %in% cart_groups,
        analyte != "Untransduced"
      ),
    aes(group = analyte),
    ref.group = "DMSO",
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Residual viability under physiologic AC concentrations",
    subtitle = "CAR 74 and CAR 77 pooled; untransduced T cells shown as a negative control",
    x = NULL,
    y = "Residual viability",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_phys_supp

plot_phys_supp_killing <- plot_phys %>%
  mutate(
    group = case_when(
      analyte == "Untransduced" ~ "Untransduced",
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(
      group,
      levels = c("Untransduced", "41BB costim", "CD28z costim")
    ),
    analyte = factor(
      analyte,
      levels = c("Untransduced", "DMSO", "AC10:0", "AC14:0", "AC18:1", "ALCAR")
    )
  ) %>%
  filter(!is.na(group))

p_phys_supp_killing <- ggplot(
  plot_phys_supp_killing,
  aes(x = analyte, y = killing_score)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6,
    color = "black",
    fill = "grey75"
  ) +
  geom_jitter(
    aes(shape = group, color = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(ref.group = 'DMSO')+
  scale_color_manual(
    values = c(
      "Untransduced" = "grey60",
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Untransduced" = 16,
      "41BB costim" = 17,
      "CD28z costim" = 15
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    shape = "Cell product",
    color = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_phys_supp_killing

### 
p_supra <- ggplot(
  plot_supra,
  aes(x = analyte, y = residual_viability)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6
  ) +
  geom_jitter(
    aes(shape = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_supra %>%
      filter(
        group %in% cart_groups,
        analyte != "Untransduced"
      ),
    aes(group = analyte),
    ref.group = "DMSO",
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  facet_wrap(~ group) +
  theme_bw(base_size = 13) +
  labs(
    title = "Residual viability under supra-physiologic AC concentrations",
    subtitle = "Untransduced T cells are shown as a generic negative control; AC species are compared with DMSO",
    x = NULL,
    y = "Residual viability",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_supra


### Supra supplementary figure
p_supra_supp <- ggplot(
  plot_supra,
  aes(x = analyte, y = residual_viability)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6
  ) +
  geom_jitter(
    aes(shape = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_supra %>%
      filter(
        group %in% cart_groups,
        analyte != "Untransduced"
      ),
    aes(group = analyte),
    ref.group = "DMSO",
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Residual viability under supra-physiologic AC concentrations",
    subtitle = "CAR 74 and CAR 77 pooled; untransduced T cells shown as a negative control",
    x = NULL,
    y = "Residual viability",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_supra_supp

plot_supra_supp_killing <- plot_supra %>%
  mutate(
    group = case_when(
      analyte == "Untransduced" ~ "Untransduced",
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(
      group,
      levels = c("Untransduced", "41BB costim", "CD28z costim")
    ),
    analyte = factor(
      analyte,
      levels = c("Untransduced", "DMSO", "AC10:0", "AC14:0", "AC18:1", "ALCAR")
    )
  ) %>%
  filter(!is.na(group))


p_supra_supp_killing <- ggplot(
  plot_supra_supp_killing,
  aes(x = analyte, y = killing_score)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6,
    color = "black",
    fill = "grey75"
  ) +
  geom_jitter(
    aes(shape = group, color = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(ref.group = 'DMSO')+
  scale_color_manual(
    values = c(
      "Untransduced" = "grey60",
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Untransduced" = 16,
      "41BB costim" = 17,
      "CD28z costim" = 15
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    shape = "Cell product",
    color = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_supra_supp_killing


#### Comparison of concentrations
ac_concentration_data <- viab_plot_base %>%
  filter(
    group %in% cart_groups,
    analyte %in% c("AC10:0", "AC14:0", "AC18:1", "ALCAR"),
    concentration %in% c("phys", "supra")
  ) %>%
  mutate(
    analyte = factor(analyte, levels = c("AC10:0", "AC14:0", "AC18:1", "ALCAR")),
    concentration = factor(concentration, levels = c("phys", "supra")),
    group = factor(group, levels = c("CAR 74", "CAR 77"))
  )


p_ac_phys_vs_supra_killing <- ggplot(
  ac_concentration_data,
  aes(x = analyte, y = killing_score, fill = concentration)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(,
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    aes(group = concentration),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    fill = "AC concentration",
    shape = "CAR-T construct"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_ac_phys_vs_supra_killing

## concentration by product
ac_concentration_data_plot <- ac_concentration_data %>%
  mutate(
    group = case_when(
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("41BB costim", "CD28z costim"))
  )

p_ac_phys_vs_supra_by_construct <- ggplot(
  ac_concentration_data_plot,
  aes(x = analyte, y = killing_score, fill = concentration)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(color = concentration),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    aes(group = concentration),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  facet_wrap(~ group) +
  scale_fill_manual(
    values = c(
      "phys" = "darkblue",
      "supra" = "darkred"
    )
  ) +
  scale_color_manual(
    values = c(
      "phys" = "darkblue",
      "supra" = "darkred"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    fill = "AC concentration",
    color = "AC concentration"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_ac_phys_vs_supra_by_construct

## Construct comparison
p_phys_construct <- ggplot(
  plot_phys %>%
    filter(
      group %in% cart_groups,
      analyte != "Untransduced"
    ),
  aes(x = analyte, y = killing_score, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(shape = group),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_phys %>%
      filter(
        group %in% cart_groups,
        analyte != "Untransduced"
      ),
    aes(group = group),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  theme_bw(base_size = 13) +
  labs(x = NULL,
    y = "Killing score",
    fill = "Cell product",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_phys_construct


p_supra_construct <- ggplot(
  plot_supra %>%
    filter(
      group %in% cart_groups,
      analyte != "Untransduced"
    ),
  aes(x = analyte, y = killing_score, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(shape = group),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_supra %>%
      filter(
        group %in% cart_groups,
        analyte != "Untransduced"
      ),
    aes(group = group),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  theme_bw(base_size = 13) +
  labs(x = NULL,
    y = "Killing score",
    fill = "Cell product",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_supra_construct





# Pairwise comparisons of each AC species against HPLM and DMSO
contrast_viab_species_by_construct <- contrast(
  emm_viab_species_by_construct,
  method = "trt.vs.ctrl",
  ref = which(levels(viab_species$analyte) == "HPLM"),
  adjust = "BH"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(control = "HPLM")

contrast_viab_species_vs_dmso <- contrast(
  emm_viab_species_by_construct,
  method = "trt.vs.ctrl",
  ref = which(levels(viab_species$analyte) == "DMSO"),
  adjust = "BH"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(control = "DMSO")

contrast_viab_species_all <- bind_rows(
  contrast_viab_species_by_construct,
  contrast_viab_species_vs_dmso
) %>%
  filter(str_detect(contrast, "AC|ALCAR|Combi")) %>%
  mutate(
    endpoint = "viability",
    interpretation = ifelse(estimate < 0, "lower viability / better killing with metabolite", "higher viability / less killing with metabolite")
  )

#write_csv(
#  contrast_viab_species_all,
#  file.path(outdir, "viability_individual_AC_species_vs_controls_by_construct.csv")
#)

# =========================
# 6. Activation/metabolic marker analysis
# =========================
# Main question:
# Do added metabolites increase marker MFI vs HPLM/DMSO?
#
# Uses log2(MFI + 1) because MFI values are usually right-skewed.
# Positive estimates for AC - control indicate higher marker expression.

# 6A. Pooled AC vs HPLM/DMSO per marker, adjusted for construct and concentration

run_marker_overall_model <- function(dat) {
  dat <- dat %>%
    filter(treatment_class %in% c("HPLM", "DMSO", "AC")) %>%
    droplevels()
  
  if (n_distinct(dat$treatment_class) < 3 || nrow(dat) < 6) return(NULL)
  
  fit <- lm(log_value ~ treatment_class + group + concentration, data = dat)
  emm <- emmeans(fit, ~ treatment_class)
  
  contrast(
    emm,
    method = list(
      "AC vs HPLM" = c(-1, 0, 1),
      "AC vs DMSO" = c(0, -1, 1)
    ),
    adjust = "BH"
  ) %>%
    summary(infer = TRUE) %>%
    as.data.frame()
}

marker_contrasts_overall <- markers %>%
  group_by(screenshot, cell_type, marker, marker_category) %>%
  group_modify(~ {
    res <- run_marker_overall_model(.x)
    if (is.null(res)) tibble()
    else as_tibble(res)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    interpretation = ifelse(estimate > 0, "higher MFI with AC", "lower MFI with AC")
  )

write_csv(
  marker_contrasts_overall,
  file.path(outdir, "marker_AC_vs_controls_overall.csv")
)

# 6B. Pooled AC vs controls per marker and per CAR-T construct

run_marker_by_construct_model <- function(dat) {
  dat <- dat %>%
    filter(treatment_class %in% c("HPLM", "DMSO", "AC")) %>%
    droplevels()
  
  if (n_distinct(dat$treatment_class) < 3 || n_distinct(dat$group) < 2 || nrow(dat) < 12) return(NULL)
  
  fit <- lm(log_value ~ treatment_class * group + concentration, data = dat)
  emm <- emmeans(fit, ~ treatment_class | group)
  
  contrast(
    emm,
    method = list(
      "AC vs HPLM" = c(-1, 0, 1),
      "AC vs DMSO" = c(0, -1, 1)
    ),
    adjust = "BH"
  ) %>%
    summary(infer = TRUE) %>%
    as.data.frame()
}

marker_contrasts_by_construct <- markers %>%
  group_by(screenshot, cell_type, marker, marker_category) %>%
  group_modify(~ {
    res <- run_marker_by_construct_model(.x)
    if (is.null(res)) tibble()
    else as_tibble(res)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    interpretation = ifelse(estimate > 0, "higher MFI with AC", "lower MFI with AC")
  )

write_csv(
  marker_contrasts_by_construct,
  file.path(outdir, "marker_AC_vs_controls_by_construct.csv")
)

# 6C. Individual AC species vs HPLM/DMSO per marker and construct

run_marker_species_model <- function(dat) {
  dat <- dat %>%
    filter(analyte %in% c(control_analytes, ac_analytes)) %>%
    mutate(analyte = factor(analyte, levels = c("HPLM", "DMSO", ac_analytes))) %>%
    droplevels()
  
  if (n_distinct(dat$analyte) < 4 || nrow(dat) < 12) return(NULL)
  
  fit <- lm(log_value ~ analyte * group + concentration, data = dat)
  emm <- emmeans(fit, ~ analyte | group)
  
  vs_hplm <- contrast(
    emm,
    method = "trt.vs.ctrl",
    ref = which(levels(dat$analyte) == "HPLM"),
    adjust = "BH"
  ) %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    mutate(control = "HPLM")
  
  vs_dmso <- contrast(
    emm,
    method = "trt.vs.ctrl",
    ref = which(levels(dat$analyte) == "DMSO"),
    adjust = "BH"
  ) %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    mutate(control = "DMSO")
  
  bind_rows(vs_hplm, vs_dmso) %>%
    filter(str_detect(contrast, "AC|ALCAR|Combi"))
}

marker_species_contrasts <- markers %>%
  group_by(screenshot, cell_type, marker, marker_category) %>%
  group_modify(~ {
    res <- run_marker_species_model(.x)
    if (is.null(res)) tibble()
    else as_tibble(res)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    interpretation = ifelse(estimate > 0, "higher MFI with metabolite", "lower MFI with metabolite")
  )

write_csv(
  marker_species_contrasts,
  file.path(outdir, "marker_individual_AC_species_vs_controls_by_construct.csv")
)

# =========================
# 7. Construct comparison
# =========================
# Compare CAR 74 vs CAR 77 under the same analyte/concentration.

# 7A. Viability: CAR 74 vs CAR 77 within each analyte/concentration
fit_construct_viab <- lm(
  residual_viability ~ group * analyte * concentration,
  data = viab_species
)

emm_construct_viab <- emmeans(
  fit_construct_viab,
  ~ group | analyte + concentration
)

contrast_construct_viab <- pairs(
  emm_construct_viab,
  adjust = "BH"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(
    endpoint = "viability",
    interpretation = ifelse(estimate < 0, "CAR 74 lower viability than CAR 77", "CAR 74 higher viability than CAR 77")
  )

write_csv(
  contrast_construct_viab,
  file.path(outdir, "viability_CAR74_vs_CAR77_within_condition.csv")
)

# 7B. Markers: CAR 74 vs CAR 77 within marker/analyte/concentration

run_construct_marker_model <- function(dat) {
  dat <- dat %>%
    filter(group %in% cart_groups) %>%
    filter(analyte %in% c(control_analytes, ac_analytes)) %>%
    droplevels()
  
  if (n_distinct(dat$group) < 2 || nrow(dat) < 8) return(NULL)
  
  fit <- lm(log_value ~ group * analyte * concentration, data = dat)
  emm <- emmeans(fit, ~ group | analyte + concentration)
  
  pairs(emm, adjust = "BH") %>%
    summary(infer = TRUE) %>%
    as.data.frame()
}

marker_construct_contrasts <- markers %>%
  group_by(screenshot, cell_type, marker, marker_category) %>%
  group_modify(~ {
    res <- run_construct_marker_model(.x)
    if (is.null(res)) tibble()
    else as_tibble(res)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    interpretation = ifelse(estimate > 0, "higher log2 MFI in first construct listed", "lower log2 MFI in first construct listed")
  )

write_csv(
  marker_construct_contrasts,
  file.path(outdir, "marker_CAR74_vs_CAR77_within_condition.csv")
)

# =========================
# 8. Plots
# =========================

# 8A. Viability plot: residual viability
p_viab <- viab_species %>%
  ggplot(aes(x = analyte, y = residual_viability, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65, position = position_dodge(width = 0.8)) +
  geom_point(
    aes(color = group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.75,
    size = 1.8
  ) +
  facet_wrap(~ concentration) +
  theme_bw(base_size = 12) +
  labs(
    title = "Residual target viability by metabolite condition and CAR-T construct",
    x = "Condition",
    y = "Residual viability",
    fill = "Group",
    color = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_viab

#ggsave(
#  filename = file.path(outdir, "plot_viability_by_condition_construct.pdf"),
#  plot = p_viab,
#  width = 10,
#  height = 5
#)

# 8B. Killing score plot
p_kill <- viab_species %>%
  ggplot(aes(x = analyte, y = killing_score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65, position = position_dodge(width = 0.8)) +
  geom_point(
    aes(color = group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.75,
    size = 1.8
  ) +
  facet_wrap(~ concentration) +
  theme_bw(base_size = 12) +
  labs(
    title = "Killing score by metabolite condition and CAR-T construct",
    subtitle = "Killing score defined as 100 - residual viability",
    x = "Condition",
    y = "Killing score",
    fill = "Group",
    color = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_kill

#ggsave(
#  filename = file.path(outdir, "plot_killing_score_by_condition_construct.pdf"),
#  plot = p_kill,
#  width = 10,
#  height = 5
#)

# 8C. Marker overview heatmap: AC vs HPLM overall
marker_heatmap_data <- marker_contrasts_overall %>%
  filter(contrast == "AC vs DMSO") %>%
  mutate(
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

p_marker_heatmap <- marker_heatmap_data %>%
  ggplot(aes(x = marker, y = cell_type, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 5) +
  facet_wrap(~ marker_category, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "Effect of ACs vs HPLM on marker MFI",
    subtitle = "Estimate is difference in log2(MFI + 1); positive = higher marker expression with ACs",
    x = "Marker",
    y = "Cell type",
    fill = "log2 effect"
  )

p_marker_heatmap

#ggsave(
#  filename = file.path(outdir, "plot_marker_heatmap_AC_vs_HPLM.pdf"),
#  plot = p_marker_heatmap,
#  width = 9,
#  height = 4
#)

# 8D. Individual marker plots
# One PDF per marker endpoint
marker_plot_dir <- file.path(outdir, "marker_plots")
dir.create(marker_plot_dir, showWarnings = FALSE)

unique_markers <- unique(markers$screenshot)

for (m in unique_markers) {
  p <- markers %>%
    filter(screenshot == m) %>%
    mutate(analyte = factor(analyte, levels = c("HPLM", "DMSO", ac_analytes))) %>%
    ggplot(aes(x = analyte, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, position = position_dodge(width = 0.8)) +
    geom_point(
      aes(color = group),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      alpha = 0.75,
      size = 1.5
    ) +
    facet_wrap(~ concentration) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(m, " MFI by condition and CAR-T construct"),
      x = "Condition",
      y = "MFI",
      fill = "Group",
      color = "Group"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    filename = file.path(marker_plot_dir, paste0("plot_", m, ".pdf")),
    plot = p,
    width = 10,
    height = 5
  )
}

# =========================
# 9. Optional: normalized viability/killing against Target Only
# =========================
# This only makes sense if Target Only exists for the same analyte/concentration.
# It estimates killing relative to target-only viability/confluency.

target_ref <- viab %>%
  filter(group == "Target Only") %>%
  group_by(concentration, analyte) %>%
  summarise(
    target_only_mean = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

viab_normalized <- viab %>%
  filter(group %in% cart_groups) %>%
  left_join(target_ref, by = c("concentration", "analyte")) %>%
  mutate(
    viability_fraction_of_target = value / target_only_mean,
    percent_killing_vs_target = 100 * (1 - viability_fraction_of_target)
  )

write_csv(
  viab_normalized,
  file.path(outdir, "viability_normalized_to_target_only.csv")
)

# =========================
# 10. Session info
# =========================

sink(file.path(outdir, "session_info.txt"))
sessionInfo()
sink()

message("Analysis complete. Results written to: ", outdir)



# -------------------------------------------------------------------------
# Required input columns assumed:
# marker_long:
#   sample_id / donor_id optional
#   construct          e.g. "CD28z", "41BB"
#   condition          e.g. "DMSO", "AC", "HPLM", etc.
#   cell_type          e.g. "CD4", "CD8"
#   marker
#   marker_category
#   MFI
# -------------------------------------------------------------------------
marker_long <- markers %>%
  filter(is_cart) %>%
  mutate(
    condition = factor(treatment_class),
    condition = relevel(condition, ref = "DMSO"),
    construct = factor(group),
    cell_type = factor(cell_type),
    marker = factor(marker),
    marker_category = factor(marker_category),
    MFI = value,
    log2_mfi = log_value
  )

levels(marker_long$condition)
levels(marker_long$construct)
unique(marker_long$cell_type)
unique(marker_long$marker)

marker_long <- marker_long %>%
  mutate(
    construct_type = case_when(
      group %in% c("CAR 74", "CAR 77") ~ "CD28z",
      group %in% c("CAR 78", "CAR 80") ~ "41BB",
      TRUE ~ NA_character_
    ),
    construct_type = factor(construct_type)
  )

baseline_check <- marker_long %>%
  filter(condition == "DMSO", !is.na(construct_type)) %>%
  group_by(marker_category, marker, cell_type) %>%
  summarise(
    n_constructs = n_distinct(construct_type),
    constructs_present = paste(sort(unique(as.character(construct_type))), collapse = ", "),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(n_constructs, marker_category, marker, cell_type)

baseline_check

baseline_group_celltype <- marker_long %>%
  filter(condition == "DMSO") %>%
  group_by(marker_category, marker, cell_type) %>%
  filter(n_distinct(construct) >= 2) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ construct, data = .x)),
    tidy = map(model, broom::tidy)
  ) %>%
  unnest(tidy) %>%
  filter(str_detect(term, "^construct")) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

p_baseline_group_celltype <- baseline_group_celltype %>%
  ggplot(aes(x = marker, y = cell_type, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 5) +
  facet_wrap(~ marker_category, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "Baseline marker expression differences between CAR groups",
    subtitle = "DMSO only; estimate is CAR 77 vs CAR 74 if CAR 74 is the reference",
    x = "Marker",
    y = "Cell type",
    fill = "log2 effect"
  )

p_baseline_group_celltype


marker_contrasts_by_group_overall <- marker_long %>%
  filter(condition %in% c("DMSO", "AC")) %>%
  droplevels() %>%
  group_by(marker_category, marker, construct) %>%
  filter(
    n_distinct(condition) >= 2,
    n_distinct(cell_type) >= 2
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition + cell_type, data = .x)),
    emm = map(model, ~ emmeans::emmeans(.x, ~ condition)),
    contrasts = map(
      emm,
      ~ contrast(.x, method = "trt.vs.ctrl", ref = "DMSO") %>%
        as.data.frame()
    )
  ) %>%
  select(marker_category, marker, construct, contrasts) %>%
  unnest(contrasts) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

marker_ac_vs_dmso_by_group_overall <- marker_contrasts_by_group_overall %>%
  filter(contrast == "AC - DMSO")

p_ac_vs_dmso_group_overall <- marker_ac_vs_dmso_by_group_overall %>%
  ggplot(aes(x = marker, y = construct, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 5) +
  facet_wrap(~ marker_category, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "Effect of ACs vs DMSO on marker MFI by CAR group",
    subtitle = "Adjusted for CD4/CD8 cell type; positive = higher marker expression with ACs",
    x = "Marker",
    y = "CAR group",
    fill = "log2 effect"
  )

p_ac_vs_dmso_group_overall


marker_contrasts_by_group_celltype <- marker_long %>%
  filter(condition %in% c("DMSO", "AC")) %>%
  droplevels() %>%
  group_by(marker_category, marker, construct, cell_type) %>%
  filter(n_distinct(condition) >= 2) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition, data = .x)),
    emm = map(model, ~ emmeans::emmeans(.x, ~ condition)),
    contrasts = map(
      emm,
      ~ contrast(.x, method = "trt.vs.ctrl", ref = "DMSO") %>%
        as.data.frame()
    )
  ) %>%
  select(marker_category, marker, construct, cell_type, contrasts) %>%
  unnest(contrasts) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

marker_ac_vs_dmso_by_group_celltype <- marker_contrasts_by_group_celltype %>%
  filter(contrast == "AC - DMSO")


p_ac_vs_dmso_group_celltype <- marker_ac_vs_dmso_by_group_celltype %>%
  ggplot(aes(x = marker, y = cell_type, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 5) +
  facet_grid(construct ~ marker_category, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "Effect of ACs vs DMSO on marker MFI by CAR group and T-cell subset",
    subtitle = "Positive estimate = higher marker expression with ACs",
    x = "Marker",
    y = "Cell type",
    fill = "log2 effect"
  )

p_ac_vs_dmso_group_celltype


marker_group_interaction_overall <- marker_long %>%
  filter(condition %in% c("DMSO", "AC")) %>%
  droplevels() %>%
  group_by(marker_category, marker) %>%
  filter(
    n_distinct(condition) >= 2,
    n_distinct(construct) >= 2,
    n_distinct(cell_type) >= 2
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition * construct + cell_type, data = .x)),
    tidy = map(model, broom::tidy)
  ) %>%
  unnest(tidy) %>%
  filter(str_detect(term, "conditionAC:construct|construct.*:conditionAC")) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

p_group_interaction_overall <- marker_group_interaction_overall %>%
  mutate(y_label = "AC effect difference") %>%
  ggplot(aes(x = marker, y = y_label, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 5) +
  facet_wrap(~ marker_category, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "CAR group-specific sensitivity to ACs",
    subtitle = "Interaction term: does AC-vs-DMSO effect differ between CAR 74 and CAR 77?",
    x = "Marker",
    y = NULL,
    fill = "interaction effect"
  )

p_group_interaction_overall


marker_group_interaction_celltype <- marker_long %>%
  filter(condition %in% c("DMSO", "AC")) %>%
  droplevels() %>%
  group_by(marker_category, marker, cell_type) %>%
  filter(
    n_distinct(condition) >= 2,
    n_distinct(construct) >= 2
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition * construct, data = .x)),
    tidy = map(model, broom::tidy)
  ) %>%
  unnest(tidy) %>%
  filter(str_detect(term, "conditionAC:construct|construct.*:conditionAC")) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

p_group_interaction_celltype <- marker_group_interaction_celltype %>%
  ggplot(aes(x = marker, y = cell_type, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 5) +
  facet_wrap(~ marker_category, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "CAR group-specific AC sensitivity by CD4/CD8 subset",
    subtitle = "Interaction term: does AC-vs-DMSO effect differ between CAR 74 and CAR 77?",
    x = "Marker",
    y = "Cell type",
    fill = "interaction effect"
  )

p_group_interaction_celltype


#### Baseline marker check

# library(dplyr)
# library(tidyr)
# library(purrr)
# library(stringr)
# library(ggplot2)
# library(broom)
 library(rstatix)
# library(ggpubr)

# Combine CD4 and CD8 cell types for DMSO baseline comparison
baseline_cd4_cd8 <- marker_long %>%
  filter(condition == "DMSO") %>%
  filter(str_detect(cell_type, regex("CD4|CD8", ignore_case = TRUE))) %>%
  mutate(
    cell_type_combined = "CD4/CD8 T cells",
    construct = case_when(
      construct == "CAR 74" ~ "41BB",
      construct == "CAR 77" ~ "CD28z",
      TRUE ~ as.character(construct)
    ),
    construct = factor(construct, levels = c("41BB", "CD28z"))
  )

# Optional: define black/grey colors dynamically based on construct levels
construct_cols <- setNames(
  c("grey70", "black")[seq_along(levels(baseline_cd4_cd8$construct))],
  levels(baseline_cd4_cd8$construct)
)

# Wilcoxon test for plotting p-values on boxplots
stat_cd4_cd8 <- baseline_cd4_cd8 %>%
  group_by(marker_category, marker, cell_type_combined) %>%
  filter(n_distinct(construct) >= 2) %>%
  wilcox_test(log2_mfi ~ construct) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p, method = "BH"),
    p_adj_label = paste0("FDR=", formatC(p_adj_global, digits = 3, format = "f"))
  ) %>%
  add_xy_position(
    x = "construct",
    fun = "max",
    data = baseline_cd4_cd8
  )

# Lower brackets slightly for each marker/facet
y_ranges <- baseline_cd4_cd8 %>%
  group_by(marker_category, marker) %>%
  summarise(
    y_range = max(log2_mfi, na.rm = TRUE) - min(log2_mfi, na.rm = TRUE),
    .groups = "drop"
  )

stat_cd4_cd8 <- stat_cd4_cd8 %>%
  left_join(y_ranges, by = c("marker_category", "marker")) %>%
  mutate(
    y.position = y.position - 0.10 * y_range
  )

# Boxplot visualization
p_baseline_cd4_cd8_box <- baseline_cd4_cd8 %>%
  ggplot(aes(x = construct, y = log2_mfi, fill = construct)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black"
  ) +
  geom_jitter(
    aes(color = construct),
    width = 0.12,
    size = 2,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_pwc()+
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.20))
  ) +
  facet_wrap(
    ~ marker,
    scales = "free_y",
    ncol = 4
  ) +
  scale_fill_manual(values = construct_cols) +
  scale_color_manual(values = construct_cols) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "none",
    plot.title = element_blank()
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "log2 MFI"
  )

p_baseline_cd4_cd8_box

p_baseline_cd4_cd8_box_adj <- ggadjust_pvalue(
  p_baseline_cd4_cd8_box,
  p.adjust.method = "BH",
  label = "p.adj.format"
)

p_baseline_cd4_cd8_box_adj



# library(dplyr)
# library(tidyr)
# library(purrr)
# library(ggplot2)
# library(stringr)
# library(emmeans)

str(marker_long)

marker_long |> filter(condition == "DMSO")

marker_contrasts_by_group_overall <- marker_long %>%
  filter(concentration == 'phys')|>
  filter(condition %in% c("DMSO", "AC")) %>%
  droplevels() %>%
  mutate(
    construct = case_when(
      construct == "CAR 74" ~ "41BB",
      construct == "CAR 77" ~ "CD28z",
      TRUE ~ as.character(construct)
    ),
    construct = factor(construct, levels = c("41BB", "CD28z"))
  ) %>%
  group_by(marker_category, marker, construct) %>%
  filter(
    n_distinct(condition) >= 2,
    n_distinct(cell_type) >= 2
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition + cell_type, data = .x)),
    emm = map(model, ~ emmeans::emmeans(.x, ~ condition)),
    contrasts = map(
      emm,
      ~ contrast(.x, method = "trt.vs.ctrl", ref = "DMSO") %>%
        as.data.frame()
    )
  ) %>%
  select(marker_category, marker, construct, contrasts) %>%
  unnest(contrasts) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

marker_ac_vs_dmso_by_group_overall <- marker_contrasts_by_group_overall %>%
  filter(contrast == "AC - DMSO")

max_abs <- max(abs(marker_ac_vs_dmso_by_group_overall$estimate), na.rm = TRUE)

p_ac_vs_dmso_group_overall <- marker_ac_vs_dmso_by_group_overall %>%
  ggplot(aes(x = marker, y = construct, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 5) +
  facet_grid(
    ~ marker_category,
    scales = "free_x",
    space = "free_x",
    labeller = as_labeller(c(
      activation = "Activation",
      metabolic = "Metabolic"
    ))
  )+
  scale_fill_gradient2(
    low = "steelblue3",
    mid = "lightgrey",
    high = "darkred",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "log2 effect"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank()
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "",
    y = ""
  )

p_ac_vs_dmso_group_overall


marker_contrasts_by_group_overall_phys <- marker_long %>%
  filter(
    condition == "DMSO" |
      (condition == "AC" & concentration == "phys")
  ) %>%
  droplevels() %>%
  mutate(
    construct = case_when(
      as.character(construct) == "CAR 74" ~ "41BB",
      as.character(construct) == "CAR 77" ~ "CD28z",
      TRUE ~ as.character(construct)
    ),
    construct = factor(construct, levels = c("41BB", "CD28z")),
    condition = factor(condition, levels = c("DMSO", "AC"))
  ) %>%
  group_by(marker_category, marker, construct) %>%
  filter(
    n_distinct(condition) >= 2,
    n_distinct(cell_type) >= 2
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition + cell_type, data = .x)),
    emm = map(model, ~ emmeans::emmeans(.x, ~ condition)),
    contrasts = map(
      emm,
      ~ contrast(.x, method = "trt.vs.ctrl", ref = "DMSO") %>%
        as.data.frame()
    )
  ) %>%
  select(marker_category, marker, construct, contrasts) %>%
  unnest(contrasts) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

marker_ac_vs_dmso_by_group_overall_phys <- marker_contrasts_by_group_overall_phys %>%
  filter(contrast == "AC - DMSO")


max_abs <- max(abs(marker_ac_vs_dmso_by_group_overall_phys$estimate), na.rm = TRUE)

p_ac_vs_dmso_group_overall_phys <- marker_ac_vs_dmso_by_group_overall_phys %>%
  ggplot(aes(x = marker, y = construct, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 5) +
  facet_grid(
    ~ marker_category,
    scales = "free_x",
    space = "free_x",
    labeller = as_labeller(c(
      activation = "Activation markers",
      metabolic = "Mitochondria & Redox markers"
    ))
  ) +
  scale_fill_gradient2(
    low = "steelblue3",
    mid = "lightgrey",
    high = "darkred",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "log2 effect"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "")

p_ac_vs_dmso_group_overall_phys



marker_contrasts_by_group_overall_supra <- marker_long %>%
  filter(
    condition == "DMSO" |
      (condition == "AC" & concentration == "supra")
  ) %>%
  droplevels() %>%
  mutate(
    construct = case_when(
      as.character(construct) == "CAR 74" ~ "41BB",
      as.character(construct) == "CAR 77" ~ "CD28z",
      TRUE ~ as.character(construct)
    ),
    construct = factor(construct, levels = c("41BB", "CD28z")),
    condition = factor(condition, levels = c("DMSO", "AC"))
  ) %>%
  group_by(marker_category, marker, construct) %>%
  filter(
    n_distinct(condition) >= 2,
    n_distinct(cell_type) >= 2
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition + cell_type, data = .x)),
    emm = map(model, ~ emmeans::emmeans(.x, ~ condition)),
    contrasts = map(
      emm,
      ~ contrast(.x, method = "trt.vs.ctrl", ref = "DMSO") %>%
        as.data.frame()
    )
  ) %>%
  select(marker_category, marker, construct, contrasts) %>%
  unnest(contrasts) %>%
  ungroup() %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

marker_ac_vs_dmso_by_group_overall_supra <- marker_contrasts_by_group_overall_supra %>%
  filter(contrast == "AC - DMSO")


max_abs <- max(abs(marker_ac_vs_dmso_by_group_overall_supra$estimate), na.rm = TRUE)

p_ac_vs_dmso_group_overall_supra <- marker_ac_vs_dmso_by_group_overall_supra %>%
  ggplot(aes(x = marker, y = construct, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 5) +
  facet_grid(
    ~ marker_category,
    scales = "free_x",
    space = "free_x",
    labeller = as_labeller(c(
      activation = "Activation markers",
      metabolic = "Mitochondria & Redox markers"
    ))
  ) +
  scale_fill_gradient2(
    low = "steelblue3",
    mid = "lightgrey",
    high = "darkred",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "log2 effect"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "")

p_ac_vs_dmso_group_overall_supra



## Marker specific analysis for supplement

# Physiologic
ac_levels_phys <- marker_long %>%
  filter(condition == "AC", concentration == "phys") %>%
  pull(analyte) %>%
  unique() %>%
  as.character()

dmso_phys_expanded <- marker_long %>%
  filter(condition == "DMSO") %>%
  tidyr::crossing(analyte_compare = ac_levels_phys)

ac_phys <- marker_long %>%
  filter(condition == "AC", concentration == "phys") %>%
  mutate(analyte_compare = as.character(analyte))

marker_specific_phys <- bind_rows(
  ac_phys,
  dmso_phys_expanded
) %>%
  mutate(
    analyte_compare = factor(analyte_compare, levels = ac_levels_phys),
    construct = case_when(
      as.character(construct) == "CAR 74" ~ "41BB",
      as.character(construct) == "CAR 77" ~ "CD28z",
      TRUE ~ as.character(construct)
    ),
    construct = factor(construct, levels = c("41BB", "CD28z")),
    condition = factor(condition, levels = c("DMSO", "AC"))
  ) %>%
  droplevels()

marker_contrasts_by_ac_phys <- marker_specific_phys %>%
  group_by(marker_category, marker, construct, analyte_compare) %>%
  filter(
    n_distinct(condition) >= 2,
    n_distinct(cell_type) >= 2
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition + cell_type, data = .x)),
    emm = map(model, ~ emmeans::emmeans(.x, ~ condition)),
    contrasts = map(
      emm,
      ~ contrast(.x, method = "trt.vs.ctrl", ref = "DMSO") %>%
        as.data.frame()
    )
  ) %>%
  select(marker_category, marker, construct, analyte_compare, contrasts) %>%
  unnest(contrasts) %>%
  ungroup() %>%
  filter(contrast == "AC - DMSO") %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

marker_contrasts_by_ac_phys %>%
  filter(marker == "GSH") %>%
  arrange(p.value)

max_abs_ac <- max(abs(marker_contrasts_by_ac_phys$estimate), na.rm = TRUE)

marker_contrasts_by_ac_phys_plot <- marker_contrasts_by_ac_phys %>%
  mutate(
    construct = recode(as.character(construct),
                       "41BB" = "41BB",
                       "CD28z" = "CD28"
    ),
    construct = factor(construct, levels = c("CD28", "41BB")),
    analyte_compare = recode(as.character(analyte_compare),
                             "ALCAR" = "AC2:0"
    ),
    analyte_compare = factor(
      analyte_compare,
      levels = c("AC2:0", "AC10:0", "AC14:0", "AC18:1")
    )
  )


max_abs_ac_phys <- max(abs(marker_contrasts_by_ac_phys_plot$estimate), na.rm = TRUE)

p_ac_specific_phys <- marker_contrasts_by_ac_phys %>%
  mutate(
    construct = recode(as.character(construct),
                       "41BB" = "41BB",
                       "CD28z" = "CD28ζ"
    ),
    construct = factor(construct, levels = c("CD28ζ", "41BB")),
    analyte_compare = recode(as.character(analyte_compare),
                             "ALCAR" = "AC2:0"
    ),
    analyte_compare = factor(
      analyte_compare,
      levels = c("AC2:0", "AC10:0", "AC14:0", "AC18:1")
    )
  ) %>%
  ggplot(aes(x = analyte_compare, y = marker, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 4) +
  facet_grid(
    marker_category ~ construct,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      marker_category = as_labeller(c(
        activation = "Activation markers",
        metabolic = "Mitochondria & Redox markers"
      )),
      construct = label_value
    )
  ) +
  scale_fill_gradient2(
    low = "steelblue3",
    mid = "lightgrey",
    high = "darkred",
    midpoint = 0,
    limits = c(-max_abs_ac_phys, max_abs_ac_phys),
    name = "log2 effect"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "")

p_ac_specific_phys



# Supraphysiologic
ac_levels_supra <- marker_long %>%
  filter(condition == "AC", concentration == "supra") %>%
  pull(analyte) %>%
  unique() %>%
  as.character()

dmso_supra_expanded <- marker_long %>%
  filter(condition == "DMSO") %>%
  tidyr::crossing(analyte_compare = ac_levels_supra)

ac_supra <- marker_long %>%
  filter(condition == "AC", concentration == "supra") %>%
  mutate(analyte_compare = as.character(analyte))

marker_specific_supra <- bind_rows(
  ac_supra,
  dmso_supra_expanded
) %>%
  mutate(
    analyte_compare = factor(analyte_compare, levels = ac_levels_supra),
    construct = case_when(
      as.character(construct) == "CAR 74" ~ "41BB",
      as.character(construct) == "CAR 77" ~ "CD28z",
      TRUE ~ as.character(construct)
    ),
    construct = factor(construct, levels = c("41BB", "CD28z")),
    condition = factor(condition, levels = c("DMSO", "AC"))
  ) %>%
  droplevels()

marker_contrasts_by_ac_supra <- marker_specific_supra %>%
  group_by(marker_category, marker, construct, analyte_compare) %>%
  filter(
    n_distinct(condition) >= 2,
    n_distinct(cell_type) >= 2
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2_mfi ~ condition + cell_type, data = .x)),
    emm = map(model, ~ emmeans::emmeans(.x, ~ condition)),
    contrasts = map(
      emm,
      ~ contrast(.x, method = "trt.vs.ctrl", ref = "DMSO") %>%
        as.data.frame()
    )
  ) %>%
  select(marker_category, marker, construct, analyte_compare, contrasts) %>%
  unnest(contrasts) %>%
  ungroup() %>%
  filter(contrast == "AC - DMSO") %>%
  mutate(
    p_adj_global = p.adjust(p.value, method = "BH"),
    significant = p_adj_global < 0.05,
    label = ifelse(significant, "*", "")
  )

marker_contrasts_by_ac_supra %>%
  filter(marker == "GSH") %>%
  arrange(p.value)

max_abs_ac <- max(abs(marker_contrasts_by_ac_supra$estimate), na.rm = TRUE)

marker_contrasts_by_ac_supra_plot <- marker_contrasts_by_ac_supra %>%
  mutate(
    construct = recode(as.character(construct),
                       "41BB" = "41BB",
                       "CD28z" = "CD28"
    ),
    construct = factor(construct, levels = c("CD28", "41BB")),
    analyte_compare = recode(as.character(analyte_compare),
                             "ALCAR" = "AC2:0"
    ),
    analyte_compare = factor(
      analyte_compare,
      levels = c("AC2:0", "AC10:0", "AC14:0", "AC18:1")
    )
  )


max_abs_ac_supra <- max(abs(marker_contrasts_by_ac_supra_plot$estimate), na.rm = TRUE)

p_ac_specific_supra <- marker_contrasts_by_ac_supra_plot %>%
  ggplot(aes(x = analyte_compare, y = marker, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 4) +
  facet_grid(
    marker_category ~ construct,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      marker_category = as_labeller(c(
        activation = "Activation markers",
        metabolic = "Mitochondria & Redox markers"
      )),
      construct = label_value
    )
  ) +
  scale_fill_gradient2(
    low = "steelblue3",
    mid = "lightgrey",
    high = "darkred",
    midpoint = 0,
    limits = c(-max_abs_ac_supra, max_abs_ac_supra),
    name = "log2 effect"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "")

p_ac_specific_supra


### Different marker status under DMSO control conditions

plot_phys_construct <- plot_phys %>%
  filter(
    group %in% cart_groups,
    analyte != "Untransduced"
  ) %>%
  mutate(
    group = case_when(
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("41BB costim", "CD28z costim"))
  )

p_phys_construct <- ggplot(
  plot_phys_construct,
  aes(x = analyte, y = killing_score, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(shape = group, color = group),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_phys_construct,
    aes(group = group),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = c(
      "41BB costim" = "grey60",
      "CD28z costim" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Killing score",
    fill = "Cell product",
    color = "Cell product",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

p_phys_construct


plot_supra_construct <- plot_supra %>%
  filter(
    group %in% cart_groups,
    analyte != "Untransduced"
  ) %>%
  mutate(
    group = case_when(
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("41BB costim", "CD28z costim"))
  )

p_supra_construct <- ggplot(
  plot_supra_construct,
  aes(x = analyte, y = killing_score, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(shape = group, color = group),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_supra_construct,
    aes(group = group),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = c(
      "41BB costim" = "grey60",
      "CD28z costim" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Killing score",
    fill = "Cell product",
    color = "Cell product",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

p_supra_construct



#####
##### ANALYSIS OF SUPERNATANTS
#####


# Water-soluble metabolites: 519 features total; 258 usable with quality score ≤3.
# Lipids: 690 features total; 241 usable with quality score ≤3.
# Both files contain 5 conditions with 3 replicates each:
#   DMEM
# 3T3-L1
# Adipocytes
# BMC + 3T3-L1
# BMC + Adipocytes


#### New Klilings

# ============================================================
# Patch: relabel ALCAR as AC2:0 and rebuild only final figures
# ============================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

cart_groups <- c("CAR 74", "CAR 77")

ac_plot_analytes <- c("AC2:0", "AC10:0", "AC14:0", "AC18:1")

plot_analytes <- c(
  "Untransduced",
  "DMSO",
  "AC2:0",
  "AC10:0",
  "AC14:0",
  "AC18:1"
)

# Start from viab, not df_main, so we do not accidentally drop ALCAR rows
viab_plot_base <- viab %>%
  mutate(
    group = str_trim(as.character(group)),
    analyte = str_trim(as.character(analyte)),
    concentration = str_trim(as.character(concentration)),
    analyte = case_when(
      analyte == "ALCAR" ~ "AC2:0",
      TRUE ~ analyte
    )
  ) %>%
  filter(
    group %in% c("UT TCells", "CAR 74", "CAR 77")
  )


make_viability_plot_data <- function(conc_use) {
  
  ut_control <- viab_plot_base %>%
    filter(
      group == "UT TCells",
      analyte == "HPLM"
    ) %>%
    mutate(
      analyte = "Untransduced",
      plot_concentration = conc_use
    )
  
  dmso_control <- viab_plot_base %>%
    filter(
      group %in% cart_groups,
      analyte == "DMSO"
    ) %>%
    mutate(
      plot_concentration = conc_use
    )
  
  ac_data <- viab_plot_base %>%
    filter(
      group %in% cart_groups,
      analyte %in% ac_plot_analytes,
      concentration == conc_use
    ) %>%
    mutate(
      plot_concentration = conc_use
    )
  
  bind_rows(ut_control, dmso_control, ac_data) %>%
    mutate(
      analyte = factor(analyte, levels = plot_analytes),
      group = factor(group, levels = c("UT TCells", "CAR 74", "CAR 77")),
      plot_concentration = factor(plot_concentration, levels = c("phys", "supra"))
    )
}

plot_phys <- make_viability_plot_data("phys")
plot_supra <- make_viability_plot_data("supra")


# ============================================================
# 1. Physiologic construct comparison
# ============================================================

plot_phys_construct <- plot_phys %>%
  filter(
    group %in% cart_groups,
    analyte != "Untransduced"
  ) %>%
  mutate(
    group = case_when(
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("41BB costim", "CD28z costim")),
    analyte = factor(
      analyte,
      levels = c("DMSO", "AC2:0", "AC10:0", "AC14:0", "AC18:1")
    )
  )

p_phys_construct <- ggplot(
  plot_phys_construct,
  aes(x = analyte, y = killing_score, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(shape = group, color = group),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_phys_construct,
    aes(group = group),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = c(
      "41BB costim" = "grey60",
      "CD28z costim" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    fill = "Cell product",
    color = "Cell product",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

p_phys_construct
ggsave("Figures_Manuscript/p_phys_construct.svg", plot = p_phys_construct, width = 5, height = 3.5)


# ============================================================
# 2. Supra construct comparison
# ============================================================

plot_supra_construct <- plot_supra %>%
  filter(
    group %in% cart_groups,
    analyte != "Untransduced"
  ) %>%
  mutate(
    group = case_when(
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("41BB costim", "CD28z costim")),
    analyte = factor(
      analyte,
      levels = c("DMSO", "AC2:0", "AC10:0", "AC14:0", "AC18:1")
    )
  )

p_supra_construct <- ggplot(
  plot_supra_construct,
  aes(x = analyte, y = killing_score, fill = group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(shape = group, color = group),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    data = plot_supra_construct,
    aes(group = group),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = c(
      "41BB costim" = "grey60",
      "CD28z costim" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    fill = "Cell product",
    color = "Cell product",
    shape = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

p_supra_construct
ggsave("Figures_Manuscript/p_supra_construct.svg", plot = p_supra_construct, width = 5, height = 3.5)



# ============================================================
# 3. Physiologic supplementary killing plot
# ============================================================

plot_phys_supp_killing <- plot_phys %>%
  mutate(
    group = case_when(
      analyte == "Untransduced" ~ "Untransduced",
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(
      group,
      levels = c("Untransduced", "41BB costim", "CD28z costim")
    ),
    analyte = factor(
      analyte,
      levels = c("Untransduced", "DMSO", "AC2:0", "AC10:0", "AC14:0", "AC18:1")
    )
  ) %>%
  filter(!is.na(group))

p_phys_supp_killing <- ggplot(
  plot_phys_supp_killing,
  aes(x = analyte, y = killing_score)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6,
    color = "black",
    fill = "grey75"
  ) +
  geom_jitter(
    aes(shape = group, color = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(ref.group = 'DMSO')+
  scale_color_manual(
    values = c(
      "Untransduced" = "grey60",
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Untransduced" = 16,
      "41BB costim" = 17,
      "CD28z costim" = 15
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    shape = "Cell product",
    color = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_phys_supp_killing
ggsave("Figures_Manuscript/p_phys_supp_killing.svg", plot = p_phys_supp_killing, width = 5, height = 4)



# ============================================================
# 4. Supra supplementary killing plot
# ============================================================

plot_supra_supp_killing <- plot_supra %>%
  mutate(
    group = case_when(
      analyte == "Untransduced" ~ "Untransduced",
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(
      group,
      levels = c("Untransduced", "41BB costim", "CD28z costim")
    ),
    analyte = factor(
      analyte,
      levels = c("Untransduced", "DMSO", "AC2:0", "AC10:0", "AC14:0", "AC18:1")
    )
  ) %>%
  filter(!is.na(group))

p_supra_supp_killing <- ggplot(
  plot_supra_supp_killing,
  aes(x = analyte, y = killing_score)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.6,
    color = "black",
    fill = "grey75"
  ) +
  geom_jitter(
    aes(shape = group, color = group),
    width = 0.12,
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(ref.group = 'DMSO')+
  scale_color_manual(
    values = c(
      "Untransduced" = "grey60",
      "41BB costim" = "grey35",
      "CD28z costim" = "black"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Untransduced" = 16,
      "41BB costim" = 17,
      "CD28z costim" = 15
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    shape = "Cell product",
    color = "Cell product"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_supra_supp_killing
ggsave("Figures_Manuscript/p_supra_supp_killing.svg", plot = p_supra_supp_killing, width = 5, height = 4)


# ============================================================
# 5. Phys vs supra by construct
# ============================================================

ac_concentration_data <- viab_plot_base %>%
  filter(
    group %in% cart_groups,
    analyte %in% ac_plot_analytes,
    concentration %in% c("phys", "supra")
  ) %>%
  mutate(
    analyte = factor(
      analyte,
      levels = c("AC2:0", "AC10:0", "AC14:0", "AC18:1")
    ),
    concentration = factor(concentration, levels = c("phys", "supra")),
    group = factor(group, levels = c("CAR 74", "CAR 77"))
  )

ac_concentration_data_plot <- ac_concentration_data %>%
  mutate(
    group = case_when(
      group == "CAR 74" ~ "41BB costim",
      group == "CAR 77" ~ "CD28z costim",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("41BB costim", "CD28z costim"))
  )

p_ac_phys_vs_supra_by_construct <- ggplot(
  ac_concentration_data_plot,
  aes(x = analyte, y = killing_score, fill = concentration)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(color = concentration),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  geom_pwc(
    aes(group = concentration),
    method = "wilcox_test",
    label = "p.format",
    p.adjust.method = "BH",
    group.by = "x.var",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  facet_wrap(~ group) +
  scale_fill_manual(
    values = c(
      "phys" = "darkblue",
      "supra" = "darkred"
    )
  ) +
  scale_color_manual(
    values = c(
      "phys" = "darkblue",
      "supra" = "darkred"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Killing score",
    fill = "AC concentration",
    color = "AC concentration"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_ac_phys_vs_supra_by_construct

ggsave("Figures_Manuscript/p_ac_phys_vs_supra_by_construct.svg", plot = p_ac_phys_vs_supra_by_construct, width = 8, height = 3.5)

