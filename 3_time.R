###
### Time course analysis of interesting metabolites across survival/CRS
###

## Creating a data frame with metabolites from all time points - not normalized
cohort_1_all_timepoints <- read_excel("Input_files/1_CART_Metabolomics_all_samples_includingNA.xlsx")
cohort_1_all_timepoints$Sample

cohort_1_all_timepoints <- cohort_1_all_timepoints %>%
  mutate(Sample.ID.Number = str_extract(Sample, "\\d+$"))

cohort_1_all_timepoints$Sample.ID.Number

bl_all_master <- all_master |>
  slice(1:(n() - 20))

bl_cohort_1_all_timepoints <- left_join(cohort_1_all_timepoints, bl_all_master, by = "Sample.ID.Number")

bl_cohort_1_all_timepoints_long <- bl_cohort_1_all_timepoints |>
  pivot_longer(cols = c(Alanine:`Cer-(24:01)`), names_to = "metabolite", values_to = "level")

bl_cohort_1_all_timepoints_long$level <- as.numeric(bl_cohort_1_all_timepoints_long$level)
bl_cohort_1_all_timepoints_long$Timepoint <- as.character(bl_cohort_1_all_timepoints_long$Timepoint)


## Unifying interesting vectors for both outcomes
t2_all_crs_corr_filter
bl_all_surv_corr_filter

time_all_metabolites <- sort(unique(c(
  t2_all_crs_corr_filter,
  bl_all_surv_corr_filter)))

time_all_metabolites

## Unifying interesting vectors for both outcomes
t2_all_crs_corr_filter
bl_all_surv_corr_filter

time_all_metabolites <- sort(unique(c(
  t2_all_crs_corr_filter,
  bl_all_surv_corr_filter)))

time_all_metabolites

## Looking at time course of ACs
# Extracting all ACs in one vector
moi_all_ac <- get_ACs(time_all_metabolites)

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_ac,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level,
             color = TAT_survival,
             fill  = TAT_survival)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.3) +
  stat_summary(fun = median, geom = "line", aes(group = TAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_color_manual(values = c("darkblue", "orange"))+
  scale_fill_manual(values = c("darkblue", "orange"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_ac,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format",
           bracket.nudge.y =-0.2)+
  stat_summary(aes(group = 1),
               fun = median, geom = "line",
               size = 0.8)+
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Day of/after CAR-T infusion",
       y     = "Serum levels")

p_time_ac_vat_high <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_ac,
    !is.na(level),
    VAT_survival == "high"
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0"))+
  scale_fill_manual(values = c("#6e51a0"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_ac_vat_high
ggsave("Figures_Manuscript/time/time_ac_vat_high.svg", plot = p_time_ac_vat_high, height = 8, width = 8)


p_time_ac_vat_low <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_ac,
    !is.na(level),
    VAT_survival == "low"
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#c63735"))+
  scale_fill_manual(values = c("#c63735"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_ac_vat_low
ggsave("Figures_Manuscript/time/time_ac_vat_low.svg", plot = p_time_ac_vat_low, height = 8, width = 8)


# Subset of ACs for the figure
moi_all_ac_subset <- c("AC-(02:0)", "AC-(06:0)", "AC-(10:0)","AC-(22:1)")

p_time_ac_vat_subset <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_ac_subset,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0", "#c63735"))+
  scale_fill_manual(values = c("#6e51a0", "#c63735"), name = c("VAT"), labels = c("High", "Low"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", nrow = 1) +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_ac_vat_subset
ggsave("Figures_Manuscript/time/time_ac_vat_subset.svg", plot = p_time_ac_vat_subset, height = 2.5, width = 7)


## Looking at time course of LPCs
# Extracting all LPCs in one vector
moi_all_lpc <- get_LPCs(time_all_metabolites)

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lpc,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level,
             color = TAT_survival,
             fill  = TAT_survival)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.3) +
  stat_summary(fun = median, geom = "line", aes(group = TAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_color_manual(values = c("darkblue", "orange"))+
  scale_fill_manual(values = c("darkblue", "orange"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lpc,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format",
           bracket.nudge.y =-0.2)+
  stat_summary(aes(group = 1),
               fun = median, geom = "line",
               size = 0.8)+
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Day of/after CAR-T infusion",
       y     = "Serum levels")


p_time_lpc_vat_high <-bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lpc,
    !is.na(level),
    VAT_survival == "high"
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0"))+
  scale_fill_manual(values = c("#6e51a0"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_lpc_vat_high
ggsave("Figures_Manuscript/time/time_lpc_vat_high.svg", plot = p_time_lpc_vat_high, height = 6, width = 8)


p_time_lpc_vat_low <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lpc,
    !is.na(level),
    VAT_survival == "low"
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#c63735"))+
  scale_fill_manual(values = c("#c63735"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_lpc_vat_low
ggsave("Figures_Manuscript/time/time_lpc_vat_low.svg", plot = p_time_lpc_vat_low, height = 6, width = 8)

# LPC subset
moi_all_lpc_subset <- c("LPC-(16:01)", "LPC-(18:01)", "LPC-(20:00)", "LPC-(24:01)")

p_time_lpc_vat_subset <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lpc_subset,
    !is.na(level),
    !is.na(VAT_survival),
    level < 300
  ) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0", "#c63735"))+
  scale_fill_manual(values = c("#6e51a0", "#c63735"), name = c("VAT"), labels = c("High", "Low"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", nrow = 1) +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_lpc_vat_subset
ggsave("Figures_Manuscript/time/time_lpc_vat_subset.svg", plot = p_time_lpc_vat_subset, height = 2.5, width = 7)

## Looking at time course of PEAs
# Extracting all PEAs in one vector
moi_all_pea <- get_PEAs(time_all_metabolites)

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_pea,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level,
             color = VAT_survival,
             fill  = VAT_survival)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.3) +
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_color_manual(values = c("darkblue", "orange"))+
  scale_fill_manual(values = c("darkblue", "orange"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_pea,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format",
           bracket.nudge.y =-0.2)+
  stat_summary(aes(group = 1),
               fun = median, geom = "line",
               size = 0.8)+
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Day of/after CAR-T infusion",
       y     = "Serum levels")


p_time_pea_vat_high <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_pea,
    !is.na(level),
    VAT_survival == "high"
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0"))+
  scale_fill_manual(values = c("#6e51a0"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_pea_vat_high
ggsave("Figures_Manuscript/time/time_pea_vat_high.svg", plot = p_time_pea_vat_high, height = 8, width = 8)

p_time_pea_vat_low <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_pea,
    !is.na(level),
    VAT_survival == "low"
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#c63735"))+
  scale_fill_manual(values = c("#c63735"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_pea_vat_low
ggsave("Figures_Manuscript/time/time_pea_vat_low.svg", plot = p_time_pea_vat_low, height = 8, width = 8)

# PEA subset
moi_all_pea_subset <- c("PEA-(28:00)", "PEA-(32:02)", "PEA-(38:06)", "PEA-(40:06)")

p_time_pea_vat_subset <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_pea_subset,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0", "#c63735"))+
  scale_fill_manual(values = c("#6e51a0", "#c63735"), name = c("VAT"), labels = c("High", "Low"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", nrow = 1) +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_pea_vat_subset
ggsave("Figures_Manuscript/time/time_pea_vat_subset.svg", plot = p_time_pea_vat_subset, height = 2.5, width = 7)


bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_pea_subset,
    !is.na(level),
    !is.na(VAT_survival)
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level,
             color = SAT_survival,
             fill  = SAT_survival)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.3) +
  stat_summary(fun = median, geom = "line", aes(group = SAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_color_manual(values = c("darkblue", "orange"))+
  scale_fill_manual(values = c("darkblue", "orange"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")




## Looking at time course of SMs
# Extracting all SMs in one vector
moi_all_SM <- get_SMs(time_all_metabolites)

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_SM,
    !is.na(level),
    !is.na(SAT_survival),
    level < 300,
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level,
             color = SAT_survival,
             fill  = SAT_survival)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.3) +
  stat_summary(fun = median, geom = "line", aes(group = SAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_color_manual(values = c("darkblue", "orange"))+
  scale_fill_manual(values = c("darkblue", "orange"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")


bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_SM,
    !is.na(level),
    !is.na(VAT_survival),
    level < 300,
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format",
           bracket.nudge.y =-0.2)+
  stat_summary(aes(group = 1),
               fun = median, geom = "line",
               size = 0.8)+
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Day of/after CAR-T infusion",
       y     = "Serum levels")

p_time_sm_vat_high <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_SM,
    !is.na(level),
    VAT_survival == "high",
    level < 300,
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0"))+
  scale_fill_manual(values = c("#6e51a0"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 3) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_sm_vat_high
ggsave("Figures_Manuscript/time/time_sm_vat_high.svg", plot = p_time_sm_vat_high, height = 4, width = 6)

p_time_sm_vat_low <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_SM,
    !is.na(level),
    VAT_survival == "low",
    level < 300,
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#c63735"))+
  scale_fill_manual(values = c("#c63735"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 3) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_sm_vat_low
ggsave("Figures_Manuscript/time/time_sm_vat_low.svg", plot = p_time_sm_vat_low, height = 4, width = 6)

# SM subset 
moi_all_SM_subset <- c("SM-(d18:1/16:01)")

p_time_sm_vat_subset <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_SM_subset,
    !is.na(level),
    !is.na(SAT_survival),
    level < 300,
  ) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0", "#c63735"))+
  scale_fill_manual(values = c("#6e51a0", "#c63735"), name = c("VAT"), labels = c("High", "Low"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", nrow = 1) +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_sm_vat_subset
ggsave("Figures_Manuscript/time/time_sm_vat_subset.svg", plot = p_time_sm_vat_subset, height = 2.5, width = 2.6)


## Looking at time course of Plasmalogene
# Extracting all Plasmalogens in one vector
moi_all_Plas <- get_Plass(time_all_metabolites)

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_Plas,
    !is.na(level),
    !is.na(SAT_survival),
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level,
             color = SAT_survival,
             fill  = SAT_survival)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.3) +
  stat_summary(fun = median, geom = "line", aes(group = SAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_color_manual(values = c("darkblue", "orange"))+
  scale_fill_manual(values = c("darkblue", "orange"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")


bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_Plas,
    !is.na(level),
    !is.na(VAT_survival),
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format",
           bracket.nudge.y =-0.2)+
  stat_summary(aes(group = 1),
               fun = median, geom = "line",
               size = 0.8)+
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Day of/after CAR-T infusion",
       y     = "Serum levels")

p_time_plas_vat_high <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_Plas,
    !is.na(level),
    VAT_survival == "high",
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0"))+
  scale_fill_manual(values = c("#6e51a0"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 3) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_plas_vat_high
ggsave("Figures_Manuscript/time/time_plas_vat_high.svg", plot = p_time_plas_vat_high, height = 4, width = 6)

p_time_plas_vat_low <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_Plas,
    !is.na(level),
    VAT_survival == "low",
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#c63735"))+
  scale_fill_manual(values = c("#c63735"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 3) +
  theme_classic() +
  labs(x  = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_plas_vat_low
ggsave("Figures_Manuscript/time/time_plas_vat_low.svg", plot = p_time_plas_vat_low, height = 4, width = 6)

# SM subset 
moi_all_Plas_subset<- c("PlasC-(40:05)", "PlasEA-(38:03)")

p_time_plas_vat_subset <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_Plas_subset,
    !is.na(level),
    !is.na(VAT_survival),
  ) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0", "#c63735"))+
  scale_fill_manual(values = c("#6e51a0", "#c63735"), name = c("VAT"), labels = c("High", "Low"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", nrow = 1) +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_plas_vat_subset
ggsave("Figures_Manuscript/time/time_plas_vat_subset.svg", plot = p_time_plas_vat_subset, height = 2.5, width = 4.1)



###
###

## Looking at time course of Plasmalogene
# Extracting all Plasmalogens in one vector
moi_all_lysoPAF <- get_lysoPAF(time_all_metabolites)

bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lysoPAF,
    !is.na(level),
    !is.na(VAT_survival),
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level,
             color = SAT_survival,
             fill  = SAT_survival)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.3) +
  stat_summary(fun = median, geom = "line", aes(group = SAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_color_manual(values = c("darkblue", "orange"))+
  scale_fill_manual(values = c("darkblue", "orange"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")


bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lysoPAF,
    !is.na(level),
    !is.na(VAT_survival),
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint),
             y     = level)) +
  geom_boxplot(width    = 0.6, alpha    = 0.4) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format",
           bracket.nudge.y =-0.2)+
  stat_summary(aes(group = 1),
               fun = median, geom = "line",
               size = 0.8)+
  # stat_summary(fun = median, geom = "point", aes(group = VAT_survival), position= position_dodge(0.8),
  #              size = 2.5) +
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x     = "Day of/after CAR-T infusion",
       y     = "Serum levels")

p_time_lysopaf_vat_high <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lysoPAF,
    !is.na(level),
    VAT_survival == "high",
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0"))+
  scale_fill_manual(values = c("#6e51a0"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 2) +
  theme_classic() +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_lysopaf_vat_high
ggsave("Figures_Manuscript/time/time_lysopaf_vat_high.svg", plot = p_time_lysopaf_vat_high, height = 4, width = 4)

p_time_lysopaf_vat_low <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lysoPAF,
    !is.na(level),
    VAT_survival == "low",
  ) %>%
  ggplot(aes(x     = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  geom_pwc(method = "wilcox_test", ref.group = 1, p.adjust.method = "fdr", label = "p.adj.format", label.size = 3, bracket.nudge.y = -0.15, step.increase = 0.08)+
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#c63735"))+
  scale_fill_manual(values = c("#c63735"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", ncol = 2) +
  theme_classic() +
  labs(x  = "Timepoint",y = "Serum level")+
  guides(color = "none", fill = "none")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_lysopaf_vat_low
ggsave("Figures_Manuscript/time/time_lysopaf_vat_low.svg", plot = p_time_lysopaf_vat_low, height = 4, width = 4)

#  subset 
moi_all_lysopaf_subset<- c("lysoPAF-(16:0)", "lysoPAF-(18:0)")

p_time_lysopaf_vat_subset <- bl_cohort_1_all_timepoints_long %>%
  filter(
    metabolite %in% moi_all_lysopaf_subset,
    !is.na(level),
    !is.na(SAT_survival),
  ) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = VAT_survival, fill  = VAT_survival)) +
  geom_boxplot(width = 0.6, alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width  = 0.65), 
              size = 2, alpha = 0.2) +
  stat_summary(fun = median, geom = "line", aes(group = VAT_survival), position= position_dodge(0.8),
               size = 0.8) +
  scale_color_manual(values = c("#6e51a0", "#c63735"))+
  scale_fill_manual(values = c("#6e51a0", "#c63735"), name = c("VAT"), labels = c("High", "Low"))+
  scale_x_discrete(labels = c("0", "3-5", "14"))+
  facet_wrap(~ metabolite, scales = "free_y", nrow = 1) +
  labs(x     = "Timepoint",y = "Serum level")+
  guides(color = "none")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9))

p_time_lysopaf_vat_subset
ggsave("Figures_Manuscript/time/time_lysopaf_vat_subset.svg", plot = p_time_lysopaf_vat_subset, height = 2.5, width = 4.1)

