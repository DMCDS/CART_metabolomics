# ==================================================================================================
# 5_invitro.R
# Purpose: Revision-added in vitro CAR-T lipid supplementation analyses.
# Main input:
#   - input_files/invitro_combined.xlsx, sheet "combined_long_for_R_updated_v3"
# Main outputs:
#   - Statistical result tables in invitro_analysis_results/
#   - Manuscript/revision figures in Figures_Manuscript/
# Notes for reviewers/readers:
#   - The script analyzes BCMA CAR-T constructs CAR 74 and CAR 77, untransduced controls, and NCI-H929 targets.
#   - Main readouts include killing/viability and flow-cytometry markers related to activation, redox state,
#     mitochondrial biology, and construct-specific metabolic remodeling.
#   - This copy keeps the original all-in-one structure but adds comments to clarify the analysis flow.
# ==================================================================================================

###
### Data loading and variable definition ----
###

input_file <- "input_files/invitro_combined.xlsx"
sheet_name <- "combined_long_for_R_updated_v3"

outdir <- "invitro_analysis_results"
dir.create(outdir, showWarnings = FALSE)

# Conditions to exclude globally
exclude_analytes <- c("Lyso16:1", "Palmitate")
ac_analytes <- c("AC2:0", "AC10:0", "AC14:0", "AC18:1")
control_analytes <- c("HPLM", "DMSO")

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

###
### Clean data ----
### 

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

###
#### Expression of markers at baseline (Fig. 6C) ----
###

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

# define black/grey colors dynamically based on construct levels
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


###
### Marker analysis phys vs supra (Fig. 6D) ----
###

str(marker_long)

marker_long |> filter(condition == "DMSO")

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


###
### Marker specific analysis for supplement (Fig. S16) ----
###

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


###
### Killing analysis (Fig. 6B, S14) -----
###

plot_analytes <- c(
  "Untransduced",
  "DMSO",
  "AC2:0",
  "AC10:0",
  "AC14:0",
  "AC18:1"
)

# 
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



# 1. Physiologic construct comparison (Fig. 6B)

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



# 2. Supra construct comparison (Fig. 6B)


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




# 3. Physiologic supplementary killing plot (Fig. S14A)


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




# 4. Supra supplementary killing plot (Fig. S14A)


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



# 5. Phys vs supra by construct (Fig. S14B)


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

