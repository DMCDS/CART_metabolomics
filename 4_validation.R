# ==================================================================================================
# 4_validation.R
# Purpose: Validation cohort analyses for metabolite signatures, adiposity groups, survival, CRS,
#          and selected correlation analyses.
# Main inputs:
#   - Input_files/1_validation_metabolomics_master.xlsx
#   - Input_files/all_master.xlsx
# Main outputs:
#   - Validation cohort figures and intermediate files in Figures_Manuscript/ and Input_files/
# Dependencies:
#   - Run source("0_packages.R") first.
#   - Some plots/statistics reuse discovery-derived metabolite lists/cutoffs from earlier scripts.
# Notes for reviewers/readers:
#   - This script applies the discovery-cohort logic to an independent validation dataset.
#   - The script also generates comparison plots between training and validation cohorts.
# ==================================================================================================

###
### Data wrangling and preparation of data frame ----
###

validation_met_master <- read_excel("Input_files/1_validation_metabolomics_master.xlsx")
str(validation_met_master)
all_master_new <- read.xlsx("Input_files/all_master.xlsx")

all_master_new$Responder <- as.character(all_master_new$Responder)
all_master_new$CRS_high <- as.character(all_master_new$CRS_high)
all_master_new$Maximaler.CRS.Grad <- as.character(all_master_new$Maximaler.CRS.Grad)
all_master_new$ICANS_high <- as.character(all_master_new$ICANS_high)
all_master_new <- all_master_new |>
  mutate(VAT_survival = ifelse(VAT > median(VAT, na.rm =T), "high", "low"),
         VAT_CRS = ifelse(VAT > 161.8, "high", "low"),
         SAT_survival = ifelse(SAT > median(SAT, na.rm =T), "high", "low"),
         SAT_CRS = ifelse(SAT > 209, "high", "low"),
         TAT_survival = ifelse(TAT > median(TAT, na.rm =T), "high", "low"),
         TAT_CRS = ifelse(TAT > 310, "high", "low"))

#Removing controls 
validation_met_master <- validation_met_master |> filter(Sample != "Control")

validation_met_master <- validation_met_master %>%
  rename(!!!rename_map)

# Combining metabolites with clinical data
str(all_master)
str(validation_met_master)
validation_master <- left_join(validation_met_master, all_master_new, by=join_by("Sample_ID" == "Sample"))
str(validation_master)

## Validation of association with body composition and metabolite levels

val_t1_filter <- validation_master |> filter(Time == "Day 0", !is.na(VAT)) |> select(Sample_ID, any_of(qc_mb_cohort_combined))
val_t1_sample_vector <- val_t1_filter |> select(Sample_ID) |> unlist() |> as.vector()

val_t1_VAT_vector <- all_master_new |> filter(Sample %in% c(val_t1_sample_vector)) |> select(VAT_CRS) 

val_t1_filter$VAT_survival <- val_t1_VAT_vector
val_t1_filter <- val_t1_filter|>
  rename("Sample" = "Sample_ID")|>
  select(Sample, VAT_survival, everything())|>
  mutate(VAT_survival = ifelse(VAT_survival == "high", 1,0))|>
  filter(!is.na(Alanine))
write.csv(val_t1_filter, "Input_files/val_t1_filter.csv", row.names = F)


###
### Normalization of validation cohort metabolites at day 0 for survival analyses ----
###

## 1.1 Loading of data for normalization and analyses
mSet_valt1norm<-InitDataObjects("pktable", "stat", FALSE)
mSet_valt1norm<-Read.TextData(mSet_valt1norm, "Input_files/val_t1_filter.csv", "rowu", "disc");
mSet_valt1norm<-SanityCheckData(mSet_valt1norm)
mSet_valt1norm<-ReplaceMin(mSet_valt1norm);
mSet_valt1norm<-SanityCheckData(mSet_valt1norm)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_valt1norm<-FilterVariable(mSet_valt1norm, "median", 0, "F")
mSet_valt1norm<-PreparePrenormData(mSet_valt1norm)

## Normalization by sum and data scaling based on auto-scaling
mSet_valt1norm<-Normalization(mSet_valt1norm, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

# Extraction and saving of the normalized data into a new tibble
val_t1_norm <- as.data.frame(mSet_valt1norm[["dataSet"]][["norm"]])

# Load original table and cbind sample name and time point label
val_t1_original <- read.csv("Input_files/val_t1_filter.csv", na = "NA")

## Left join original data and normalized to link vat groups
val_t1_norm$Sample <- row.names(val_t1_norm)
val_t1_norm <- val_t1_norm %>% select(Sample, everything()) %>%
  arrange(Sample)

val_t1_sample_id <- val_t1_original %>% select(Sample, VAT_survival)

val_t1_norm <- left_join(val_t1_norm, val_t1_sample_id, by = "Sample")
val_t1_norm <- val_t1_norm %>%
  select(Sample, VAT_survival, everything())

val_t1_norm$VAT_survival <- as.character(val_t1_norm$VAT_survival)


val_t1_norm_bc <- left_join(val_t1_norm, all_master, by = "Sample")
val_t1_norm_bc <- val_t1_norm_bc |> filter(!is.na(Entity))


validation_master_t1 <- validation_master |> filter(Time == "Day 0")

###
### Normalization of validation cohort metabolites at day 2-5 for CRS analysis  ---- 
###

## Validation of association with body composition and metabolite levels

val_t2_filter <- validation_master |> filter(Time == "Day 2-5", !is.na(VAT)) |> select(Sample_ID, any_of(qc_mb_cohort_combined))
val_t2_sample_vector <- val_t2_filter |> select(Sample_ID) |> unlist() |> as.vector()

val_t2_VAT_vector <- all_master_new |> filter(Sample %in% c(val_t2_sample_vector)) |> select(VAT_CRS) 

val_t2_filter$VAT_CRS <- val_t2_VAT_vector
val_t2_filter <- val_t2_filter|>
  rename("Sample" = "Sample_ID")|>
  select(Sample, VAT_CRS, everything())|>
  mutate(VAT_CRS = ifelse(VAT_CRS == "high", 1,0))|>
  filter(!is.na(Alanine))
write.csv(val_t2_filter, "Input_files/val_t2_filter.csv", row.names = F)


## 1.1 Loading of data for normalization and analyses
mSet_valt2norm<-InitDataObjects("pktable", "stat", FALSE)
mSet_valt2norm<-Read.TextData(mSet_valt2norm, "Input_files/val_t2_filter.csv", "rowu", "disc");
mSet_valt2norm<-SanityCheckData(mSet_valt2norm)
mSet_valt2norm<-ReplaceMin(mSet_valt2norm);
mSet_valt2norm<-SanityCheckData(mSet_valt2norm)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_valt2norm<-FilterVariable(mSet_valt2norm, "median", 0, "F")
mSet_valt2norm<-PreparePrenormData(mSet_valt2norm)

## Normalization by sum and data scaling based on auto-scaling
mSet_valt2norm<-Normalization(mSet_valt2norm, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

# Extraction and saving of the normalized data into a new tibble
val_t2_norm <- as.data.frame(mSet_valt2norm[["dataSet"]][["norm"]])

# Load original table and cbind sample name and time point label
val_t2_original <- read.csv("Input_files/val_t2_filter.csv", na = "NA")

## Left join original data and normalized to link vat groups
val_t2_norm$Sample <- row.names(val_t2_norm)
val_t2_norm <- val_t2_norm %>% select(Sample, everything()) %>%
  arrange(Sample)

val_t2_sample_id <- val_t2_original %>% select(Sample, VAT_CRS)

val_t2_norm <- left_join(val_t2_norm, val_t2_sample_id, by = "Sample")
val_t2_norm <- val_t2_norm %>%
  select(Sample, VAT_CRS, everything())

val_t2_norm$VAT_CRS <- as.character(val_t2_norm$VAT_CRS)


val_t2_norm_bc <- left_join(val_t2_norm, all_master, by = "Sample")


validation_master_t2 <- validation_master |> filter(Time == "Day 2-5")


###
### Comparison of training and validation cohort plus basic CRS analyses (Fig. 4B-E) ----
###

p_vat_training_vs_validation <- all_master |>
  ggplot(aes(x=cohort, y=VAT))+
  geom_jitter(aes(fill = factor(cohort), color = factor(cohort)), width = 0.2, size = 3, alpha = 0.9) +
  geom_boxplot(aes(fill =  factor(cohort)), outlier.shape = NA, alpha = 0.7) +
  geom_pwc() +
  scale_x_discrete(labels = c("Training", "Validation"))+
  labs(x = "", y = "VAT [cm²]") +
  guides(color = "none", fill = "none")+
  scale_fill_manual(values = c("grey", "darkblue")) +
  scale_color_manual(values = c("grey", "darkblue")) +
  ylim(0,500)+
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )

p_vat_training_vs_validation
ggsave("Figures_Manuscript/p_vat_training_vs_validation.svg", plot = p_vat_training_vs_validation, width = 2, height = 3)


p_sat_training_vs_validation <- all_master |>
  ggplot(aes(x=cohort, y=SAT))+
  geom_jitter(aes(fill = factor(cohort), color = factor(cohort)), width = 0.2, size = 3, alpha = 0.9) +
  geom_boxplot(aes(fill =  factor(cohort)), outlier.shape = NA, alpha = 0.7) +
  geom_pwc() +
  scale_x_discrete(labels = c("Training", "Validation"))+
  labs(x = "", y = "SAT [cm²]") +
  guides(color = "none", fill = "none")+
  scale_fill_manual(values = c("grey", "lightblue")) +
  scale_color_manual(values = c("grey", "lightblue")) +
  ylim(0,650)+
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )

p_sat_training_vs_validation
ggsave("Figures_Manuscript/p_sat_training_vs_validation.svg", plot = p_sat_training_vs_validation, width = 2, height = 3)

p_crs_validation_pie <- all_master |> 
  filter(cohort == "validation") |> 
  ggplot(aes(x = "", fill = factor(Maximaler.CRS.Grad, levels = c(0, 1, 2, 3)))) +
  geom_bar(width = 1) +
  coord_polar(theta = "y", start = 0, direction = -1) +  # start at 12 o'clock and go clockwise
  scale_fill_manual(values = c("#FFF5E1", "#FFDAB9", "#FFC04C", "#FF8C00"), name = "CRS grade") +
  theme_void()  # clean theme for pie chart
#theme(legend.title = element_blank())

ggsave("Figures_Manuscript/p_crs_validation_pie.svg", plot = p_crs_validation_pie, width = 2, height = 3)


all_master |> 
  filter(cohort == "validation") |> 
  group_by(Maximaler.CRS.Grad) |>
  count() |>
  ungroup()|>
  mutate(perc = n/sum(n))


p_crs_validation_vat <- all_master |> 
  filter(cohort == "validation") |> 
  ggplot(aes(x = factor(CRS_high), y = VAT)) +
  geom_jitter(aes(fill = factor(CRS_high), color = factor(CRS_high)), width = 0.2, size = 3, alpha = 0.9) +
  geom_boxplot(aes(fill =  factor(CRS_high)), outlier.shape = NA, alpha = 0.7) +
  geom_pwc(ref.group = "0", method = "wilcox.test", label.size = 4
           #tip.length = 0.01,
  ) +
  scale_fill_manual(values = c("#FFDAB9", "#FF8C00")) +
  scale_color_manual(values = c("#FFDAB9", "#FF8C00")) +
  scale_x_discrete(labels = c("0-1", "≥ 2"))+
  labs(x = "CRS grade", y = "VAT [cm²]") +
  guides(color = "none", fill = "none")+
  ylim(0,500)+
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )
p_crs_validation_vat
ggsave("Figures_Manuscript/p_crs_validation_vat.svg", plot = p_crs_validation_vat, width = 2, height = 3)

p_crs_validation_sat <- all_master |> 
  filter(cohort == "validation") |> 
  ggplot(aes(x = factor(CRS_high), y = SAT)) +
  geom_jitter(aes(fill = factor(CRS_high), color = factor(CRS_high)), width = 0.2, size = 3, alpha = 0.9) +
  geom_boxplot(aes(fill =  factor(CRS_high)), outlier.shape = NA, alpha = 0.7) +
  geom_pwc(ref.group = "0", method = "wilcox.test", label.size = 4
           #tip.length = 0.01,
  ) +
  scale_fill_manual(values = c("#FFDAB9", "#FF8C00")) +
  scale_color_manual(values = c("#FFDAB9", "#FF8C00")) +
  scale_x_discrete(labels = c("0-1", "≥ 2"))+
  labs(x = "CRS grade", y = "SAT [cm²]") +
  guides(color = "none", fill = "none")+
  ylim(0,550)+
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )
p_crs_validation_sat
ggsave("Figures_Manuscript/p_crs_validation_sat.svg", plot = p_crs_validation_sat, width = 2, height = 3)


crs_ci_val_vat <- survfit2(Surv(CRS_onset_high, as.numeric(`CRS_high`)) ~ VAT_CRS, data = all_master |> filter(cohort == "validation"))
p_crs_ci_val_vat <- ggsurvplot(crs_ci_val_vat,
                           fun = "event",
                           ylim= c(0,1), xlim=c(0, 19), break.x.by = 3, ylab="CRS grade ≥ 2", xlab="Days after CAR-T infusion",
                           pval= TRUE, pval.coord = c(1, 0.92), pval.size = 4,
                           size = 1.15,
                           axes.offset = FALSE,
                           risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                           tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                           conf.int = FALSE,
                           ggtheme = theme_classic2(10),
                           font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                           fontsize=3,
                           legend.labs = c("High", "Low"),
                           legend.title = c("VAT"),
                           palette = c("darkblue", "darkgrey")
)
p_crs_ci_val_vat
ggsave(filename = "Figures_Manuscript/p_crs_ci_val_vat.svg", plot = p_crs_ci_val_vat$plot,
       width = 3, height = 2.5)

crs_validation_moi <- c(crs_lysopaf_moi, crs_ac_moi, crs_pea_moi, crs_lpc_moi)

val_t2_norm_bc <- val_t2_norm_bc |> filter(!is.na(Entity))

###
### Validation of metabolites of interest (moi) in the validation cohort (Fig. 4F) ----
###

t2_val_glm_adjusted <- data.frame(marker = character(),
                                  coefficient = numeric(),
                                  std_error = numeric(),
                                  p_value = numeric(),
                                  lower95 = numeric(),
                                  upper95 = numeric(),
                                  stringsAsFactors = FALSE)

for(i in crs_validation_moi) {
  # Fit logistic regression model
  formula_str <- paste0("as.factor(CRS_high) ~ `",i,"` + Geschlecht + Costim")
  
  #print(formula_str)
  model <- glm(formula_str, data = val_t2_norm_bc, family = binomial)
  
  # Extract coefficients, standard errors, and p-values, and confidence intervals
  coef_summary <- summary(model)$coefficients[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
  coinf <- exp(confint(model))
  
  # Create a data frame with results for the current marker
  marker_results <- data.frame(
    marker = i,
    coefficient = coef_summary["Estimate"],
    std_error = coef_summary["Std. Error"],
    p_value = coef_summary["Pr(>|z|)"],
    lower95 = coinf[2,1],
    upper95 = coinf[2,2])
  
  # Append results to the main data frame
  t2_val_glm_adjusted <- rbind(t2_val_glm_adjusted, marker_results)
}


t2_val_glm_adjusted <- t2_val_glm_adjusted|>
  mutate(OR = exp(coefficient), FDR = p.adjust(p_value, method = "BH"))

t2_val_glm_adjusted$group <- sapply(as.character(t2_val_glm_adjusted$marker), get_group_new)

t2_val_glm_adjusted |>
  #filter(FDR <= 0.1)|>
  filter(!is.na(lower95))|>
  ggplot()+
  geom_point(aes(x = marker, y = OR, color = group),
             size = 3, shape = 19, alpha = 0.6)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = upper95))+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #
  ylim(-1, 10) +
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  labs(color = "Group")+
  #scale_color_manual(values = c("darkgrey", "darkblue", "orange", "pink"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )


# Now create the columns needed for meta-analysis
meta_val_crs <- t2_val_glm_adjusted %>%
  mutate(
    logOR    = log(OR),
    logLower = log(lower95),
    logUpper = log(upper95),
    # approximate standard error for log(HR)
    SE_logOR = (logUpper - logLower) / (2 * 1.96)
  )


# Double-check that meta_df:
meta_val_crs <- meta_val_crs |> filter(!(marker %in% c("PI-(38:07)", "PI-(40:03)","PI-(40:08)","PI-(40:09)")))

meta_results_val_crs <- meta_val_crs %>%
  filter(!is.na(lower95), upper95 != "Inf")|>
  group_by(group) %>%
  nest() %>%
  # n_mets is the number of metabolites in each group
  mutate(
    n_mets = map_int(data, nrow),
    # Fit a random-effects meta-analysis for each group
    fit = map(data, ~ rma.uni(
      yi  = .x$logOR,
      sei = .x$SE_logOR,
      method = "REML"
    )),
    combined_logOR = map_dbl(fit, ~ as.numeric(.x$b)),
    ci.lb          = map_dbl(fit, ~ .x$ci.lb),
    ci.ub          = map_dbl(fit, ~ .x$ci.ub),
    pval           = map_dbl(fit, ~ .x$pval)
  ) %>%
  # Exponentiate to get back to HR scale
  mutate(
    combined_OR = exp(combined_logOR),
    lower95     = exp(ci.lb),
    upper95     = exp(ci.ub)
  ) %>%
  select(
    group,
    n_mets,
    combined_logOR,
    ci.lb,
    ci.ub,
    pval,
    combined_OR,
    lower95,
    upper95
  )

meta_crs_val_significant <- c("Acylcarnitine", "Lysophosphatidylcholine")

p_meta_crs_val <- meta_results_val_crs |>
  mutate(col_group = if_else(group %in% meta_crs_val_significant, "highlight", "other")) %>% 
  filter(group != "Other/Unclassified")|>
  ggplot()+
  geom_linerange(aes(x = reorder(group, combined_OR), ymin = lower95, 
                     ymax = upper95), size=0.8, alpha = 0.8)+
  geom_point(aes(x = reorder(group, combined_OR), y = combined_OR, size = n_mets, color=col_group),
             shape = 19, alpha = 0.8)+
  coord_flip()+ 
  labs(y = "Combined OR (95%CI)", x = "")+
  scale_color_manual(values = c("highlight" = "#CC0000", "other" = "black"), guide = "none") +
  scale_size_continuous(breaks = c(3,11), range = c(2,7))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  guides(color = "none")+
  labs(size="Metabolites")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

print(p_meta_crs_val)
ggsave("Figures_Manuscript/p_meta_crs_val.svg", plot = p_meta_crs_val, width = 5.5, height = 3)


###
### Generating an AC signature to test association with CRS and survival (Fig. 4G-H) ----
###

## Creating a acetylcarnitine signature and testing for CRS
validation_master_t2 <- validation_master |> filter(Time == "Day 2-5")

# 1. calculate medians for each acetylcarnitine
val_ac_med_vec <- validation_master_t2 %>% 
  summarise(across(all_of(crs_ac_moi), ~ median(.x, na.rm = TRUE))) %>% 
  unlist()

# 2. createvalidation_master_t2# 2. create the AC_signature: count of ACs above median
validation_master_t2 <- validation_master_t2 %>%
  rowwise() %>%
  mutate(
    AC_signature = sum(c_across(all_of(crs_ac_moi)) > val_ac_med_vec)
  ) %>%
  ungroup()

# 3. high vs low
validation_master_t2 <- validation_master_t2 %>%
  mutate(
    AC_sig_median = if_else(
      AC_signature > median(AC_signature), 
      "High", 
      "Low"
    )
  )

validation_master_t2 <- validation_master_t2 %>%
  mutate(
    AC_sig_group = case_when(
      AC_signature %in% c(0) ~ "0",
      AC_signature %in% c(1,2) ~ "1-2",
      AC_signature %in% c(3,4) ~ "3-4"
    )
  )

# 4. Cox models for PFS and OS
crs_ci_val_ac_sig <- survfit(Surv(CRS_onset_high, as.numeric(CRS_high)) ~ AC_sig_group, 
                             data = validation_master_t2)
p_crs_ci_val_ac_sig <- ggsurvplot(crs_ci_val_ac_sig,
                               fun = "event",
                               ylim= c(0,1), xlim=c(0, 19), break.x.by = 3, ylab="CRS grade ≥ 2", xlab="Days after CAR-T infusion",
                               pval= TRUE, pval.coord = c(1, 0.92), pval.size = 4,
                               size = 1.15,
                               axes.offset = T,
                               risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                               tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                               conf.int = FALSE,
                               ggtheme = theme_classic2(10),
                               font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                               fontsize=3,
                               legend.labs = c("0", "1-2", "3-4"),
                               legend.title = c("AC signature"),
                               palette = c("#88DC88", "#22A822", "#006600")
)

ggsave(filename = "Figures_Manuscript/p_crs_ci_val_ac_sig.svg", plot = p_crs_ci_val_ac_sig$plot,
       width = 3, height = 2.5)

p_crs_validation_ac_sig <- validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  ggplot(aes(x = AC_sig_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  # scale_x_discrete(limits = c("LPC_low", "LPC_high"))+
  scale_fill_manual(values=c("#FFF5E1", "#FFDAB9", "#FFC04C", "#FF8C00"))+
  labs(x="AC signature", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

p_crs_validation_ac_sig
ggsave("Figures_Manuscript/p_crs_validation_ac_sig.svg", plot = p_crs_validation_ac_sig, width = 4, height = 3)

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  with(table(AC_sig_group, Maximaler.CRS.Grad)) |>
  chisq.test(.)


## Testing AC signature for survival in the validation cohort 
cox_pfs <- coxph(
  Surv(PFS_days, PFS_event) ~ AC_sig_group, 
  data = validation_master_t2
)
cox_os  <- coxph(
  Surv(OS_days, OS_event) ~ AC_sig_group, 
  data = validation_master_t2
)

pfs_ac_signature <- survfit(Surv(PFS_days/30.44, PFS_event) ~ AC_sig_group, 
                            data = validation_master_t2)
ggsurvplot(pfs_ac_signature,
           pval = T)
p_km_pfs_ac_signature <- ggsurvplot(pfs_ac_signature,
                              ylab = "Estimated PFS",
                              xlab = "Months after CAR-T infusion",
                              break.time.by = 3,
                              xlim = c(0,18),
                              censor.size = 5,
                              pval = TRUE,
                              pval.coord = c(0.3, 0.1),
                              pval.size = 4,
                              size = 1.5,
                              axes.offset = F,
                              risk.table = F,
                              risk.table.title = "No. at risk",
                              risk.table.heigbcma.ht = .2,
                              survival.plot.heigbcma.ht = 0.9,
                              tables.y.text = FALSE,
                              tables.theme = theme_cleantable(base_size = 2),
                              conf.int = F,
                              ggtheme = theme_classic2(10),
                              font.title = c(9, "bold"),
                              font.tickslab = c(10),
                              font.legend.labs = c(10),
                              font.x = c(10, "bold"),
                              font.y = c(10, "bold"),
                              fontsize = 3,
                              legend.title = "AC signature",
                              legend.labs= c("0", "1-2", "3-4"),
                              palette = c("#88DC88", "#22A822", "#006600")
)
p_km_pfs_ac_signature
ggsave(filename = "Figures_Manuscript/km/p_km_pfs_ac_signature.svg", plot = p_km_pfs_ac_signature$plot,
       width = 3, height = 2.5)

os_ac_signature <- survfit(Surv(OS_days/30.44, OS_event) ~ AC_sig_group, 
                           data = validation_master_t2)
ggsurvplot(os_ac_signature,
           pval = T)
p_km_os_ac_signature <- ggsurvplot(os_ac_signature,
                                    ylab = "Estimated OS",
                                    xlab = "Months after CAR-T infusion",
                                    break.time.by = 3,
                                    xlim = c(0,18),
                                    censor.size = 5,
                                    pval = TRUE,
                                    pval.coord = c(0.3, 0.1),
                                    pval.size = 4,
                                    size = 1.5,
                                    axes.offset = F,
                                    risk.table = F,
                                    risk.table.title = "No. at risk",
                                    risk.table.heigbcma.ht = .2,
                                    survival.plot.heigbcma.ht = 0.9,
                                    tables.y.text = FALSE,
                                    tables.theme = theme_cleantable(base_size = 2),
                                    conf.int = F,
                                    ggtheme = theme_classic2(10),
                                    font.title = c(9, "bold"),
                                    font.tickslab = c(10),
                                    font.legend.labs = c(10),
                                    font.x = c(10, "bold"),
                                    font.y = c(10, "bold"),
                                    fontsize = 3,
                                    legend.title = "AC signature",
                                    legend.labs= c("0", "1-2", "3-4"),
                                    palette = c("#88DC88", "#22A822", "#006600")
)
p_km_os_ac_signature
ggsave(filename = "Figures_Manuscript/km/p_km_os_ac_signature.svg", plot = p_km_os_ac_signature$plot,
       width = 3, height = 2.5)

