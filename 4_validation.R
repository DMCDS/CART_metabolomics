###
### Validation analysis
###

validation_met_master <- read_excel("Input_files/1_validation_metabolomics_master.xlsx")
str(validation_met_master)
all_master_new <- read.xlsx("Input_files/all_master.xlsx")


#

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

## Combining all metabolites of interest

moi_all_val <- sort(unique(c(moi_all_ac, moi_all_lpc, moi_all_lysoPAF, moi_all_pea)))
surv_validation_moi <- c(surv_plas_moi, surv_ac_moi, surv_pea_moi, surv_sm_moi)

bl_val_all_cox_adjusted <- data.frame(marker = character(),
                                  HR = numeric(),
                                  lower95 = numeric(),
                                  higher95 = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (i in surv_validation_moi) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(PFS_days, PFS_event)  ~ `",i,"` + Geschlecht + Costim")
  model_cox <- coxph(as.formula(formula_str), data = val_t1_norm_bc)
  
  # Extract coefficients and their standard errors
  summary_model_cox <- summary(model_cox)
  
  # Create a data frame with results for the current marker
  marker_results_cox <- data.frame(
    marker = i,
    HR = summary_model_cox$conf.int[1,"exp(coef)"],
    lower95 = summary_model_cox$conf.int[1,"lower .95"],
    higher95 = summary_model_cox$conf.int[1,"upper .95"],
    p_value = summary_model_cox$coefficients[1,"Pr(>|z|)"]
  )
  
  # Append results to the main data frame
  bl_val_all_cox_adjusted  <- rbind(bl_val_all_cox_adjusted, marker_results_cox)
}


bl_val_all_cox_adjusted <- bl_val_all_cox_adjusted %>%
  mutate(group = sapply(marker, get_group_new)) |>
  filter(higher95 != "Inf")# your custom get_group() function

# Now create the columns needed for meta-analysis
meta_val_survival <- bl_val_all_cox_adjusted %>%
  mutate(
    logHR    = log(HR),
    logLower = log(lower95),
    logUpper = log(higher95),
    # approximate standard error for log(HR)
    SE_logHR = (logUpper - logLower) / (2 * 1.96)
  )


# Double-check that meta_df has the columns you expect:
meta_val_survival <- meta_val_survival |> filter(!(marker %in% c("PI-(38:07)", "PI-(40:03)","PI-(40:08)","PI-(40:09)")))

meta_results_val_survival <- meta_val_survival %>%
  group_by(group) %>%
  nest() %>%
  # n_mets is the number of metabolites in each group
  mutate(
    n_mets = map_int(data, nrow),
    # Fit a random-effects meta-analysis for each group
    fit = map(data, ~ rma.uni(
      yi  = .x$logHR,
      sei = .x$SE_logHR,
      method = "REML"
    )),
    combined_logHR = map_dbl(fit, ~ as.numeric(.x$b)),
    ci.lb          = map_dbl(fit, ~ .x$ci.lb),
    ci.ub          = map_dbl(fit, ~ .x$ci.ub),
    pval           = map_dbl(fit, ~ .x$pval)
  ) %>%
  # Exponentiate to get back to HR scale
  mutate(
    combined_HR = exp(combined_logHR),
    lower95     = exp(ci.lb),
    upper95     = exp(ci.ub)
  ) %>%
  # Select and arrange columns as desired
  select(
    group,
    n_mets,
    combined_logHR,
    ci.lb,
    ci.ub,
    pval,
    combined_HR,
    lower95,
    upper95
  )

p_surv_val_group <- meta_results_val_survival |>
  #filter(FDR <= 0.1)|>
  ggplot()+
  geom_point(aes(x = reorder(group, combined_HR), y = combined_HR, size = n_mets, color=group),
             shape = 19, alpha = 0.6)+
  geom_linerange(aes(x = reorder(group, combined_HR), ymin = lower95, 
                     ymax = upper95), size=0.8)+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #ylim = c(-1, 4)
  labs(y = "Combined Hazard Ratio (95%CI)", x = "")+
  #scale_color_manual(values = c("darkgrey", "darkblue", "orange", "pink"))+
  scale_size_continuous(breaks = c(1,5,10), range = c(3,8))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  guides(color = "none")+
  labs(size="Metabolites")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

print(p_surv_val_group)

ggsave("Figures_Manuscript/meta_results_val_survival.svg", plot = p_surv_val_group, width = 7, height = 3)

pfs_val_PEA3403 <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`PEA-(34:03)` > median(`PEA-(34:03)`), "PEA_high", "PEA_low"), data=val_t1_norm_bc)
summary(pfs_val_PEA3403)
summary(pfs_val_PEA3403, times = 365)
p_km_pfsval_PEA3403 <- ggsurvplot(pfs_val_PEA3403,
           ylab = "Estimated PFS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
           censor.size = 5,
           pval = TRUE,
           pval.coord = c(0.3, 0.05),
           pval.size = 3,
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
           legend.title = "PEA-(34:03)",
           legend.labs= c("High", "Low"),
           palette = c("black", "darkgrey")
)

os_val_PEA3403 <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`AC-(13:0)` > median(`AC-(13:0)`), "PEA_high", "PEA_low"), data=val_t1_norm_bc)
summary(os_val_PEA3403)
summary(os_val_PEA3403, times = 365)
p_km_osval_PEA3403 <- ggsurvplot(os_val_PEA3403,
           ylab = "Estimated OS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
           censor.size = 5,
           pval = TRUE,
           pval.coord = c(0.3, 0.05),
           pval.size = 3,
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
           legend.title = "PEA-(34:03)",
           legend.labs= c("High", "Low"),
           palette = c("black", "darkgrey")
)


pfs_val_AC221 <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`AC-(22:1)` > median(`AC-(22:1)`), "AC_high", "AC_low"), data=val_t1_norm_bc)
summary(pfs_val_AC221)
summary(pfs_val_AC221, times = 365)
p_km_pfsval_AC221 <- ggsurvplot(pfs_val_AC221,
           ylab = "Estimated PFS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
           censor.size = 5,
           pval = TRUE,
           pval.coord = c(0.3, 0.05),
           pval.size = 3,
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
           legend.title = "AC-(22:1)",
           legend.labs= c("High", "Low"),
           palette = c("black", "darkgrey")
)

os_val_AC221 <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`AC-(22:1)` > median(`AC-(22:1)`), "AC_high", "AC_low"), data=val_t1_norm_bc)
summary(os_val_AC221)
summary(os_val_AC221, times = 365)
p_km_osval_AC221 <- ggsurvplot(os_val_AC221,
           ylab = "Estimated OS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
           censor.size = 5,
           pval = TRUE,
           pval.coord = c(0.3, 0.05),
           pval.size = 3,
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
           legend.title = "AC-(22:1)",
           legend.labs= c("High", "Low"),
           palette = c("black", "darkgrey")
)

pfs_val_SMI181 <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`SM-(d18:1/18:02)` > median(`SM-(d18:1/18:02)`), "AC_high", "AC_low"), data=validation_master |> filter(Time == "Day 0"))
summary(pfs_val_SMI181)
summary(pfs_val_SMI181, times = 365)
p_km_pfsval_SMI181 <- ggsurvplot(pfs_val_SMI181,
                                ylab = "Estimated PFS",
                                xlab = "Months after CAR-T infusion",
                                break.time.by = 3,
                                xlim = c(0,26),
                                censor.size = 5,
                                pval = TRUE,
                                pval.coord = c(0.3, 0.05),
                                pval.size = 3,
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
                                legend.title = "SM-(d18:1/18:02)",
                                legend.labs= c("High", "Low"),
                                palette = c("black", "darkgrey")
)

os_val_SMI181 <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`SM-(d18:1/18:02)` > median(`SM-(d18:1/18:02)`), "AC_high", "AC_low"), data=val_t1_norm_bc)
summary(os_val_SMI181)
summary(os_val_SMI181, times = 365)
p_km_osval_SMI181 <- ggsurvplot(os_val_SMI181,
                               ylab = "Estimated OS",
                               xlab = "Months after CAR-T infusion",
                               break.time.by = 3,
                               xlim = c(0,26),
                               censor.size = 5,
                               pval = TRUE,
                               pval.coord = c(0.3, 0.05),
                               pval.size = 3,
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
                               legend.title = "SM-(d18:1/18:02)",
                               legend.labs= c("High", "Low"),
                               palette = c("black", "darkgrey")
)

pfs_val_Plas <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`PlasEA-(38:06)` > median(`PlasEA-(38:06)`), "AC_high", "AC_low"), data=validation_master |> filter(Time == "Day 0"))
summary(pfs_val_Plas)
summary(pfs_val_Plas, times = 365)
p_km_pfsval_Plas <- ggsurvplot(pfs_val_Plas,
                                ylab = "Estimated PFS",
                                xlab = "Months after CAR-T infusion",
                                break.time.by = 3,
                                xlim = c(0,26),
                                censor.size = 5,
                                pval = TRUE,
                                pval.coord = c(0.3, 0.05),
                                pval.size = 3,
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
                                legend.title = "PlasEA-(38:06)",
                                legend.labs= c("High", "Low"),
                                palette = c("black", "darkgrey")
)

os_val_Plas <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`PlasEA-(38:06)` > median(`PlasEA-(38:06)`), "AC_high", "AC_low"), data=val_t1_norm_bc)
summary(os_val_Plas)
summary(os_val_Plas, times = 365)
p_km_osval_Plas <- ggsurvplot(os_val_Plas,
                               ylab = "Estimated OS",
                               xlab = "Months after CAR-T infusion",
                               break.time.by = 3,
                               xlim = c(0,26),
                               censor.size = 5,
                               pval = TRUE,
                               pval.coord = c(0.3, 0.05),
                               pval.size = 3,
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
                               legend.title = "PlasEA-(38:06)",
                               legend.labs= c("High", "Low"),
                               palette = c("black", "darkgrey")
)

### MAYBE no validation because of differente time point

### Testing if CRS validation works

validation_master$CRS_high
validation_master %>%
  filter(Time == "Day 2-5", !is.na(Maximaler.CRS.Grad), Entity != "MM")|>
   mutate(lysoPAF_group = ifelse(`lysoPAF-(18:1)` >= median(`lysoPAF-(18:1)`), "lysoPAF_high", "lysoPAF_low")) %>%
  ggplot(aes(x = lysoPAF_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  scale_x_discrete(limits = c("lysoPAF_low", "lysoPAF_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

validation_master %>%
  filter(Time == "Day 2-5", !is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  mutate(lysoPAF_group = ifelse(`lysoPAF-(18:1)` >= median(`lysoPAF-(18:1)`), "lysoPAF_high", "lysoPAF_low")) %>%
  with(table(lysoPAF_group, Maximaler.CRS.Grad)) |>
  fisher.test(.)

ggplot(validation_master, aes(x=`lysoPAF-(18:1)`, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

cor.test(validation_master$`lysoPAF-(18:1)`, validation_master$STLV, method = "spearman")


### Testing if CRS validation works

validation_master$CRS_high
validation_master %>%
  filter(Time == "Day 6-8", !is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  mutate(LPC_group = ifelse(`LPC-(20:01)` >= median(`LPC-(20:01)`), "LPC_high", "LPC_low")) %>%
  ggplot(aes(x = LPC_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  scale_x_discrete(limits = c("LPC_low", "LPC_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

validation_master %>%
  filter(Time == "Day 6-8", !is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  mutate(lysoPAF_group = ifelse(`lysoPAF-(18:1)` >= median(`lysoPAF-(18:1)`), "lysoPAF_high", "lysoPAF_low")) %>%
  with(table(lysoPAF_group, Maximaler.CRS.Grad)) |>
  fisher.test(.)

ggplot(validation_master |> filter(Time == "Day 2-5"), aes(x=`lysoPAF-(18:1)`, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

cor.test(validation_master$`lysoPAF-(18:1)`, validation_master$STLV, method = "spearman")


validation_master <- validation_master %>%
  mutate(
    CRS_group_new = case_when(
      Maximaler.CRS.Grad <= 1 ~ "0-1",
      Maximaler.CRS.Grad == 2 ~ "2",
      Maximaler.CRS.Grad == 3 ~ "3")
  )

validation_master$CRS_high
validation_master %>%
  filter(Time == "Day 6-8", !is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  mutate(AC_group = ifelse(`AC-(10:0)` >= median(`AC-(10:0)`), "AC_high", "AC_low")) %>%
  ggplot(aes(x = AC_group, fill = factor(CRS_group_new))) +
  geom_bar(position = "fill") +
  scale_x_discrete(limits = c("AC_low", "AC_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()


validation_master$CRS_group_new

validation_master %>%
  filter(Time == "Day 2-5", !is.na(CRS_group_new), Entity != "MM")|>
  mutate(AC_group = ifelse(`AC-(10:0)` >= median(`AC-(10:0)`), "AC_high", "AC_low")) %>%
  with(table(AC_group, CRS_group_new)) |>
  chisq.test(.)

ggplot(validation_master |> filter(Time == "Day 2-5"), aes(x=`AC-(10:0)`, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

cor.test(validation_master$`AC-(10:0)`, validation_master$VAT, method = "spearman")


#AC
t2_vat_crs_norm_bc %>%
  mutate(AC_group = ifelse(`AC-(10:0)` >= median(`AC-(10:0)`), "AC_high", "AC_low")) %>%
  ggplot(aes(x = AC_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  scale_x_discrete(limits = c("AC_low", "AC_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

t2_vat_crs_norm_bc %>%
  mutate(AC_group = ifelse(`AC-(10:0)` >= median(`AC-(10:0)`), "AC_high", "AC_low")) %>%
  with(table(AC_group, Maximaler.CRS.Grad)) |>
  fisher.test(.)

ggplot(t2_vat_crs_norm_bc, aes(x=`AC-(10:0)`, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

ggplot(t2_vat_crs_norm_bc, aes(x=`AC-(10:0)`, y=VAT))+
  geom_point()+
  geom_smooth(method = "lm")



###
###
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


t2_val_glm_adjusted <- data.frame(marker = character(),
                                  coefficient = numeric(),
                                  std_error = numeric(),
                                  p_value = numeric(),
                                  lower95 = numeric(),
                                  upper95 = numeric(),
                                  stringsAsFactors = FALSE)


for(i in t2_all_crs_corr_filter) {
  # Fit logistic regression model
  formula_str <- paste0("as.factor(CRS_high) ~ `",i,"` + Costim + STLV + Geschlecht")
  
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

str(t2_val_glm_adjusted)
t2_val_glm_adjusted |>
  #filter(FDR <= 0.1)|>
  ggplot()+
  geom_point(aes(x = marker, y = OR, color = group),
             size = 3, shape = 19, alpha = 0.6)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = upper95))+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #ylim = c(-1, 4)
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

## Meta-analysis of group effects

# Now create the columns needed for meta-analysis
meta_val_crs <- t2_val_glm_adjusted %>%
  mutate(
    logOR    = log(OR),
    logLower = log(lower95),
    logUpper = log(upper95),
    # approximate standard error for log(HR)
    SE_logOR = (logUpper - logLower) / (2 * 1.96)
  )

## !! Multiple values are the same for different metabolites - need to be removed
str(meta_val_crs)

# Double-check that meta_df has the columns you expect:
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
  # Select and arrange columns as desired
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

meta_results_val_crs |>
  filter(!group %in% c("Lysophosphatidylinositol", "Carbohydrate / Sugar"))|>
  #filter(FDR <= 0.1)|>
  ggplot()+
  geom_point(aes(x = reorder(group, combined_OR), y = combined_OR, size = n_mets, color=group),
             shape = 19, alpha = 0.6)+
  geom_linerange(aes(x = reorder(group, combined_OR), ymin = lower95, 
                     ymax = upper95), size=0.8)+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #ylim = c(-1, 4)
  labs(y = "Combined Odds Ratio (95%CI)", x = "")+
  #scale_color_manual(values = c("darkgrey", "darkblue", "orange", "pink"))+
  scale_size_continuous(breaks = c(1,5,10), range = c(3,8))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  guides(color = "none")+
  labs(size="Metabolites")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )




## Neuer Validierungsversuch
## Generierung von Metabolitenscores basierend auf den CRS Metaboliten
crs_ac_moi <- get_ACs(t2_all_crs_corr_filter)
surv_ac_moi <- get_ACs(bl_all_surv_corr_filter)
surv_sm_moi <- get_SMs(bl_all_surv_corr_filter)
surv_plas_moi <- get_Plass(bl_all_surv_corr_filter)
surv_pea_moi <- get_PEAs(bl_all_surv_corr_filter)
str(validation_master_t2)



## Generierung von Metabolitenscores basierend auf den CRS Metaboliten
crs_pea_moi <- get_PEAs(t2_all_crs_corr_filter)
surv_ac_moi <- get_ACs(bl_all_surv_corr_filter)
str(validation_master_t2)

# 1. calculate medians for each acetylcarnitine
validation_master_t2$`AC-(12:1)`
med_vec_pea <- validation_master_t2 %>% 
  summarise(across(all_of(crs_pea_moi), ~ median(.x, na.rm = TRUE))) %>% 
  unlist()

# med_vec_ac <- validation_master_t1 |>
#   summarise(across(all_of(crs_ac_moi), ~ median(.x, na.rm = TRUE))) %>% 
#   unlist()

# 2. createvalidation_master_t2# 2. create the AC_signature: count of ACs above median
validation_master_t2 <- validation_master_t2 %>%
  rowwise() %>%
  mutate(
    PEA_signature = sum(c_across(all_of(crs_pea_moi)) > med_vec_pea)
  ) %>%
  ungroup()

# validation_master_t1 <- validation_master_t1 %>%
#   rowwise() %>%
#   mutate(
#     AC_signature = sum(c_across(all_of(surv_ac_moi)) > med_vec_ac)
#   ) %>%
#   ungroup()

# 3. optionally dichotomize into high vs low
validation_master_t2 <- validation_master_t2 %>%
  mutate(
    PEA_sig_median = if_else(
      PEA_signature > median(PEA_signature), 
      "High", 
      "Low"
    )
  )

validation_master_t2 <- validation_master_t2 %>%
  mutate(
    PEA_sig_group = case_when(
      AC_signature %in% c(0,1) ~ "0-1",
      AC_signature %in% c(2,3,4) ~ "2-4"
    )
  )

# validation_master_t1 <- validation_master_t1 %>%
#   mutate(
#     AC_sig_median = if_else(
#       AC_signature > median(AC_signature), 
#       "High", 
#       "Low"
#     )
#   )
# 
# validation_master_t1 <- validation_master_t1 %>%
#   mutate(
#     AC_sig_group = case_when(
#       AC_signature %in% c(0,1) ~ "0-1",
#       AC_signature %in% c(2,3,4) ~ "2-4"
#     )
#   )


# 4. Cox models for PFS and OS
#    – swap in your actual column names for time/event below
coxph(
  Surv(PFS_days, PFS_event) ~ PEA_sig_median, 
  data = validation_master_t2
)
coxph(
  Surv(OS_days, OS_event) ~ PEA_sig_median, 
  data = validation_master_t2
)

pfs_pea_signature <- survfit(Surv(PFS_days, PFS_event) ~ PEA_sig_median, 
                            data = validation_master_t2)
ggsurvplot(pfs_pea_signature,
           pval = T)

os_pea_signature <- survfit(Surv(OS_days, OS_event) ~ PEA_sig_median, 
                           data = validation_master_t2)
ggsurvplot(os_pea_signature,
           pval = T)

validation_master_t2$CRS_onset_high
crs_pea_signature <- survfit(Surv(CRS_onset_high, as.numeric(CRS_high)) ~ PEA_sig_median, 
                            data = validation_master_t2)
ggsurvplot(crs_pea_signature,
           xlim = c(0,19),
           fun = "event",
           pval = T)


validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  ggplot(aes(x = PEA_signature, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  # scale_x_discrete(limits = c("LPC_low", "LPC_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  with(table(PEA_sig_median, Maximaler.CRS.Grad)) |>
  chisq.test(.)

ggplot(validation_master_t2 |> filter(Time == "Day 2-5"), aes(x=PEA_sig_median, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

cor.test(validation_master_t2$PEA_signature, validation_master_t2$STLV, method = "spearman")




crs_ac_moi <- get_ACs(t2_all_crs_corr_filter)
crs_pea_moi <- get_PEAs(t2_all_crs_corr_filter)
crs_lysopaf_moi <- get_lysoPAF(t2_all_crs_corr_filter)
crs_lpc_moi <- get_LPCs(t2_all_crs_corr_filter)

# 1) define your metabolite sets
pos_mets <- c(crs_ac_moi, crs_pea_moi)            # ACs + PEAs ↑ score
neg_mets <- c(crs_lysopaf_moi, crs_lpc_moi)        # LPCs + lysoPAFs ↓ score


pos_mets <- c("AC-(10:0)", "AC-(13:0)", "PEA-(36:02)", "PEA-(38:06)")            # ACs + PEAs ↑ score
neg_mets <- c("LPC-(20:03)", "LPC-(22:00)", "lysoPAF-(18:0)", "lysoPAF-(16:0)")  


all_mets <- c(pos_mets, neg_mets)

# 2) compute medians for each
medians <- validation_master_t2 %>% 
  summarise(across(all_of(all_mets), ~ median(.x, na.rm = TRUE))) %>% 
  unlist()

# 3) rowwise compute the signature
validation_master_t2 <- validation_master_t2 %>%
  rowwise() %>%
  mutate(
    # count how many pos_mets exceed their medians
    pos_count = sum(c_across(all_of(pos_mets)) > medians[pos_mets], na.rm = TRUE),
    # count how many neg_mets exceed their medians
    neg_count = sum(c_across(all_of(neg_mets)) > medians[neg_mets], na.rm = TRUE),
    # directional signature: + for pos, – for neg
    CRS_signature = pos_count - neg_count
  ) %>%
  ungroup()

# 4) peek at your new signature
validation_master_t2 %>% 
  select(Sample_ID, pos_count, neg_count, CRS_signature, Maximaler.CRS.Grad, Entity)

validation_master_t2 %>%
 # filter(!is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  ggplot(aes(x = CRS_signature, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  # scale_x_discrete(limits = c("LPC_low", "LPC_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()


###

validation_master_t2_long <- validation_master_t2 |>
  pivot_longer(cols = all_of(all_mets), names_to = "metabolite", values_to = "level")

validation_master_t2_long|>
  ggplot(aes(x=factor(Maximaler.CRS.Grad), y=level)) +
  geom_point()+
  geom_boxplot()+
  facet_wrap(vars(metabolite), scales = "free_y")+
  theme_classic()

validation_master_t2_long|>
  ggplot(aes(x=VAT, y=level)) +
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(vars(metabolite), scales = "free_y")+
  theme_classic()

validation_master_t2_long|>
  ggplot(aes(x=factor(CRS_high), y=level)) +
  geom_point()+
  geom_boxplot()+
  facet_wrap(vars(metabolite), scales = "free_y")

validation_master_t2|>
  ggplot(aes(x=factor(Maximaler.CRS.Grad), y=`PEA-(40:07)`))+
  geom_point() +
  geom_boxplot()

### time course in validation cohort

validation_master |>
  pivot_longer(cols = all_of(moi_all_ac), names_to = "acs", values_to = "level")|>
  ggplot(aes(x     = as.factor(Time),
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
  scale_x_discrete(labels = c("-6", "2-5", "6-8"))+
  facet_wrap(~ acs, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")

validation_master |>
  pivot_longer(cols = all_of(moi_all_ac_subset), names_to = "acs", values_to = "level")|>
  ggplot(aes(x     = as.factor(Time),
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
  scale_x_discrete(labels = c("-6", "2-5", "6-8"))+
  facet_wrap(~ acs, scales = "free_y") +
  theme_classic() +
  labs(x     = "Time point [day]",
       y     = "Serum level",
       color = "VAT survival",
       fill  = "VAT survival")



####

crs_ac_moi <- get_ACs(t2_all_crs_corr_filter)
crs_pea_moi <- get_PEAs(t2_all_crs_corr_filter)
crs_lysopaf_moi <- get_lysoPAF(t2_all_crs_corr_filter)
crs_lpc_moi <- get_LPCs(t2_all_crs_corr_filter)

validation_master_t2 <- validation_master_t2 %>%
  rowwise() %>%
  mutate(
    AC_sum      = sum(c_across(all_of(crs_ac_moi)),      na.rm = TRUE),
    PEA_sum     = sum(c_across(all_of(crs_pea_moi)),      na.rm = TRUE),
    LPC_sum     = sum(c_across(all_of(crs_lpc_moi)),      na.rm = TRUE),
    lysoPAF_sum = sum(c_across(all_of(crs_lysopaf_moi)),  na.rm = TRUE)
  ) %>%
  ungroup()

# 2. Compute medians of those sums
med_sums <- validation_master_t2 %>%
  summarise(
    AC_sum_med      = median(AC_sum,      na.rm = TRUE),
    PEA_sum_med     = median(PEA_sum,     na.rm = TRUE),
    LPC_sum_med     = median(LPC_sum,     na.rm = TRUE),
    lysoPAF_sum_med = median(lysoPAF_sum, na.rm = TRUE)
  ) %>%
  unlist()

# 3. Assign points based on median cut‑points
#    +1 for AC and PEA above median; –1 for LPC and lysoPAF above median
validation_master_t2 <- validation_master_t2 %>%
  mutate(
    AC_point      = case_when(AC_sum      > med_sums["AC_sum_med"]      ~  1,
                              TRUE                                       ~  0),
    PEA_point     = case_when(PEA_sum     > med_sums["PEA_sum_med"]     ~  1,
                              TRUE                                       ~  0),
    LPC_point     = case_when(LPC_sum     > med_sums["LPC_sum_med"]     ~ -1,
                              TRUE                                       ~  0),
    lysoPAF_point = case_when(lysoPAF_sum > med_sums["lysoPAF_sum_med"] ~ -1,
                              TRUE                                       ~  0),
    summed_signature = AC_point + PEA_point + LPC_point + lysoPAF_point
  )

validation_master_t2$summed_signature

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  mutate(summed_signature_group = ifelse(summed_signature <0, "low", "high"))|>
  ggplot(aes(x = summed_signature_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  # scale_x_discrete(limits = c("LPC_low", "LPC_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  mutate(summed_signature_group = ifelse(summed_signature <0, "low", "high"))|>
  with(table(summed_signature_group, Maximaler.CRS.Grad)) |>
  chisq.test(.)

ggplot(validation_master_t2 |> filter(Time == "Day 2-5"), aes(x=PEA_sig_median, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

cor.test(validation_master_t2$PEA_signature, validation_master_t2$STLV, method = "spearman")


#####
##### Just looking at correlations between MOIs and body composition in the validation cohort
#####

moi_all_ac

val_t2_norm_bc
moi_all_ac_subset
crs_ac_moi
surv_ac_moi

crs_pea_moi

val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_ac_moi),
               names_to = "metabolite",
               values_to = "level")|>
  ggplot(aes(x = SAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

cor_results <- val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_ac_moi),
               names_to = "metabolite",
               values_to = "level") %>%
  group_by(metabolite) %>%
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "pearson"),
    p_value = cor.test(VAT, level, method = "pearson")$p.value,
    n = sum(!is.na(level) & !is.na(VAT)),
    .groups = "drop"
  )

cor_results


val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_pea_moi),
               names_to = "metabolite",
               values_to = "level")|>
  ggplot(aes(x = SAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

cor_results <- val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_ac_moi),
               names_to = "metabolite",
               values_to = "level") %>%
  group_by(metabolite) %>%
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "pearson"),
    p_value = cor.test(VAT, level, method = "pearson")$p.value,
    n = sum(!is.na(level) & !is.na(VAT)),
    .groups = "drop"
  )

cor_results

val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_lpc_moi),
               names_to = "metabolite",
               values_to = "level")|>
  ggplot(aes(x = VAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

t2_vat_crs_norm_bc %>%
  pivot_longer(cols = all_of(crs_lpc_moi),
               names_to = "metabolite",
               values_to = "level")|>
  #filter(Timepoint == "2", metabolite %in% crs_ac_moi) |>
  ggplot(aes(x = SAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "SAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

cor_results <- t2_vat_crs_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_lpc_moi),
               names_to = "metabolite",
               values_to = "level") %>%
  group_by(metabolite) %>%
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "pearson"),
    p_value = cor.test(SAT, level, method = "pearson")$p.value,
    n = sum(!is.na(level) & !is.na(VAT)),
    .groups = "drop"
  )

cor_results


cor_results <- val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_lpc_moi),
               names_to = "metabolite",
               values_to = "level") %>%
  group_by(metabolite) %>%
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "pearson"),
    p_value = cor.test(SAT, level, method = "pearson")$p.value,
    n = sum(!is.na(level) & !is.na(VAT)),
    .groups = "drop"
  )

cor_results

val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_lysopaf_moi),
               names_to = "metabolite",
               values_to = "level")|>
  ggplot(aes(x = SAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

bl_cohort_1_all_timepoints_long %>%
  filter(Timepoint == "2", metabolite %in% crs_lysopaf_moi) |>
  ggplot(aes(x = SAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "SAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))


cor_results <- val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_lpc_moi),
               names_to = "metabolite",
               values_to = "level") %>%
  group_by(metabolite) %>%
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "pearson"),
    p_value = cor.test(SAT, level, method = "pearson")$p.value,
    n = sum(!is.na(level) & !is.na(VAT)),
    .groups = "drop"
  )

cor_results


val_t2_norm_bc_dlbcl <- val_t2_norm_bc |> filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))

cor.test(val_t2_norm_bc_dlbcl$`AC-(02:0)`, val_t2_norm_bc_dlbcl$VAT, method = "pearson")

val_t2_norm_bc %>%
  pivot_longer(cols = all_of(crs_pea_moi),
               names_to = "metabolite",
               values_to = "level")|>
  ggplot(aes(x = TAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

val_t1_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(crs_lpc_moi),
               names_to = "metabolite",
               values_to = "level")|>
  ggplot(aes(x = VAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

cor_results <- val_t1_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL"))|>
  pivot_longer(cols = all_of(moi_all_lpc),
               names_to = "metabolite",
               values_to = "level") %>%
  group_by(metabolite) %>%
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

# View the result# View the rcrs_lpc_moiesult
cor_results

cor.test(val_t2_norm_bc$`LPC-(24:01)`, val_t2_norm_bc$VAT, method = "spearman")

val_t2_norm_bc %>%
  #filter(!Sample %in% c("A71", "A76"))|>
  filter(!is.na(Maximaler.CRS.Grad), Entity == "DLBCL")|>
  pivot_longer(cols = all_of(crs_lysopaf_moi),
               names_to = "metabolite",
               values_to = "level")|>
  ggplot(aes(x = VAT, y = level, color = factor(metabolite))) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

cor_results <- val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity == "DLBCL")|>
  pivot_longer(cols = all_of(crs_lysopaf_moi),
               names_to = "metabolite",
               values_to = "level") %>%
  group_by(metabolite) %>%
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "pearson")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

# View the result
cor_results

val_t2_norm_bc |> select(`lysoPAF-(18:1)`, VAT, Sample, Entity)

val_t2_norm_bc %>%
  ggplot(aes(x = VAT, y = SAT)) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))

bl_all_master %>%
  ggplot(aes(x = VAT, y = SAT)) +
  geom_jitter(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~ metabolite, scales = "free_y") +
  theme_classic() +
  labs(x = "VAT", y = "Metabolite level", color = "VAT group or value") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9))


validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  mutate(summed_signature_group = ifelse(summed_signature <0, "low", "high"))|>
  ggplot(aes(x = summed_signature_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  # scale_x_discrete(limits = c("LPC_low", "LPC_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  ggplot(aes(x = CRS_high, y = `AC-(08:1)`)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc()

val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  ggplot(aes(x = CRS_high, y = `AC-(02:0)`)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc()


val_t2_norm_bc <- val_t2_norm_bc %>%
  mutate(
    CRS_group_new = case_when(
      Maximaler.CRS.Grad <= 1 ~ "0-1",
      Maximaler.CRS.Grad == 2 ~ "2",
      Maximaler.CRS.Grad == 3 ~ "3")
  )



val_t2_norm_bc %>%
  filter(!is.na(Maximaler.CRS.Grad), Entity != "MM")|>
  ggplot(aes(x = CRS_high, y = `PEA-(40:07)`)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc()

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  ggplot(aes(x = CRS_high, y = `LPC-(20:02)`)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc()

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  ggplot(aes(x = CRS_high, y = `AC-(08:1)`)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc()

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  pivot_longer(cols = crs_ac_moi, names_to = "metabolite", values_to = "level")|>
  ggplot(aes(x = CRS_high, y = `level`)) +
  geom_boxplot(aes(group = metabolite), position = position_dodge()) +
  geom_jitter(aes(color = metabolite)) +
  geom_pwc()


validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  pivot_longer(cols = crs_lysopaf_moi, names_to = "metabolite", values_to = "level")|>
  ggplot(aes(x = CRS_high, y = `level`)) +
  geom_boxplot() +
  geom_jitter(aes(color = metabolite)) +
  geom_pwc()

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  ggplot(aes(x = CRS_high, y = `lysoPAF-(18:0)`)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc()


validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  ggplot(aes(x = CRS_high, y = VAT)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc()

validation_master_t2 %>%
  filter(!is.na(Maximaler.CRS.Grad))|>
  ggplot(aes(x = CRS_high, y = SAT)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc()



validation_master_t2$`AC-(08:1)`
  
  geom_bar(position = "fill") +
  # scale_x_discrete(limits = c("LPC_low", "LPC_high"))+
  scale_fill_manual(values=c("lightgrey","#FFCCCC", "#FF5151" ,"#990000"))+
  labs(x="", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()
  
)

  
####
#### Looking at correlations with AT and metabolites from baseline cohort 
  
## --- 1.  Correlation tables -------------------------------------------
t2_lpc_cor_vat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


t2_ac_cor_vat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
    summarise(
      cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
      p_value = cor.test(VAT, level, method = "spearman")$p.value,
      n = sum(!is.na(level) & !is.na(TAT)),
      .groups = "drop"
    )

t2_pea_cor_vat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


t2_lysopaf_cor_vat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

## --- 2.  Long data with correlation value attached --------------------
t2_lpc_corr_vat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_lpc_cor_vat_tbl, by = "metabolite")

t2_ac_corr_vat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_ac_cor_vat_tbl, by = "metabolite") 

t2_pea_corr_vat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_pea_cor_vat_tbl, by = "metabolite")     

t2_lysopaf_corr_vat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_lysopaf_cor_vat_tbl, by = "metabolite")


## --- 3.  Scatter plot coloured by correlation coefficient -------------
p_t2_lpc_corr_vat <- ggplot(t2_lpc_corr_vat_plot_data,
       aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "LPC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_t2_ac_corr_vat <- ggplot(t2_ac_corr_vat_plot_data,
       aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "AC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_t2_pea_corr_vat <- ggplot(t2_pea_corr_vat_plot_data,
       aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "PEA\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))


p_t2_lysopaf_corr_vat <- ggplot(t2_lysopaf_corr_vat_plot_data,
       aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "LysoPAF\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))


## --- 1.  Correlation tables -------------------------------------------
t2_lpc_cor_sat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

t2_ac_cor_sat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


t2_pea_cor_sat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


t2_lysopaf_cor_sat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

## --- 2.  Long data with correlation value attached --------------------
t2_lpc_corr_sat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_lpc_cor_sat_tbl, by = "metabolite")

t2_ac_corr_sat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_ac_cor_sat_tbl, by = "metabolite")    

t2_pea_corr_sat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_pea_cor_sat_tbl, by = "metabolite")     

t2_lysopaf_corr_sat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_lysopaf_cor_sat_tbl, by = "metabolite")


## --- 3.  Scatter plot coloured by correlation coefficient -------------
p_t2_lpc_corr_sat <- ggplot(t2_lpc_corr_sat_plot_data,
       aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "LPC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_t2_ac_corr_sat <- ggplot(t2_ac_corr_sat_plot_data,
       aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "AC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_t2_pea_corr_sat <-ggplot(t2_pea_corr_sat_plot_data,
       aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "PEA\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))


p_t2_lysopaf_corr_sat <- ggplot(t2_lysopaf_corr_sat_plot_data,
       aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "LysoPAF\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

###
### Baseline

## --- 1.  Correlation tables -------------------------------------------
bl_ac_cor_vat_tbl <- bl_vat_surv_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(surv_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

bl_pea_cor_vat_tbl <- bl_vat_surv_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(surv_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

bl_sm_cor_vat_tbl <- bl_vat_surv_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(surv_sm_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

bl_plas_cor_vat_tbl <- bl_vat_surv_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(surv_plas_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

## --- 2.  Long data with correlation value attached --------------------
bl_ac_corr_vat_plot_data <- bl_vat_surv_norm_bc %>% 
  pivot_longer(cols  = all_of(surv_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(bl_ac_cor_vat_tbl, by = "metabolite") 

bl_pea_corr_vat_plot_data <- bl_vat_surv_norm_bc %>% 
  pivot_longer(cols  = all_of(surv_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(bl_pea_cor_vat_tbl, by = "metabolite")     

bl_sm_corr_vat_plot_data <- bl_vat_surv_norm_bc %>% 
  pivot_longer(cols  = all_of(surv_sm_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(bl_sm_cor_vat_tbl, by = "metabolite")

bl_plas_corr_vat_plot_data <- bl_vat_surv_norm_bc %>% 
  pivot_longer(cols  = all_of(surv_plas_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(bl_plas_cor_vat_tbl, by = "metabolite")


## --- 3.  Scatter plot coloured by correlation coefficient -------------
p_bl_ac_corr_vat <- ggplot(bl_ac_corr_vat_plot_data,
       aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "AC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_bl_pea_corr_vat <- ggplot(bl_pea_corr_vat_plot_data,
       aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "PEA\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_bl_sm_corr_vat <- ggplot(bl_sm_corr_vat_plot_data,
       aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "SM\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_bl_plas_corr_vat <- ggplot(bl_plas_corr_vat_plot_data,
       aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "Plasmalogen\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

## --- 1.  Correlation tables -------------------------------------------
bl_ac_cor_sat_tbl <- bl_vat_surv_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(surv_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

bl_pea_cor_sat_tbl <- bl_vat_surv_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(surv_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

bl_sm_cor_sat_tbl <- bl_vat_surv_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(surv_sm_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

bl_plas_cor_sat_tbl <- bl_vat_surv_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(surv_plas_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

## --- 2.  Long data with correlation value attached --------------------
bl_ac_corr_sat_plot_data <- bl_vat_surv_norm_bc %>% 
  pivot_longer(cols  = all_of(surv_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(bl_ac_cor_sat_tbl, by = "metabolite") 

bl_pea_corr_sat_plot_data <- bl_vat_surv_norm_bc %>% 
  pivot_longer(cols  = all_of(surv_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(bl_pea_cor_sat_tbl, by = "metabolite")     

bl_sm_corr_sat_plot_data <- bl_vat_surv_norm_bc %>% 
  pivot_longer(cols  = all_of(surv_sm_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(bl_sm_cor_sat_tbl, by = "metabolite")

bl_plas_corr_sat_plot_data <- bl_vat_surv_norm_bc %>% 
  pivot_longer(cols  = all_of(surv_plas_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(bl_plas_cor_sat_tbl, by = "metabolite")


## --- 3.  Scatter plot coloured by correlation coefficient -------------
p_bl_ac_corr_sat <- ggplot(bl_ac_corr_sat_plot_data,
       aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "AC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_bl_pea_corr_sat <- ggplot(bl_pea_corr_sat_plot_data,
       aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "PEA\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_bl_sm_corr_sat <- ggplot(bl_sm_corr_sat_plot_data,
       aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "SM\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_bl_plas_corr_sat <- ggplot(bl_plas_corr_sat_plot_data,
       aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "Plasmalogen\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))


####
#### Looking at correlations with AT and metabolites in the validation cohort

## --- 1.  Correlation tables -------------------------------------------
val_t2_lpc_cor_vat_tbl <- val_t2_norm_bc %>% 
 # filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

val_t2_ac_cor_vat_tbl <- val_t2_norm_bc %>% 
  #filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

val_t2_pea_cor_vat_tbl <- val_t2_norm_bc %>% 
 # filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


val_t2_lysopaf_cor_vat_tbl <- val_t2_norm_bc %>% 
  #filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

## --- 2.  Long data with correlation value attached --------------------
t2_lpc_corr_vat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_lpc_cor_vat_tbl, by = "metabolite")

t2_ac_corr_vat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_ac_cor_vat_tbl, by = "metabolite") 

t2_pea_corr_vat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_pea_cor_vat_tbl, by = "metabolite")     

t2_lysopaf_corr_vat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_lysopaf_cor_vat_tbl, by = "metabolite")


## --- 3.  Scatter plot coloured by correlation coefficient -------------
p_t2_lpc_corr_vat <- ggplot(t2_lpc_corr_vat_plot_data,
                            aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "LPC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_t2_ac_corr_vat <- ggplot(t2_ac_corr_vat_plot_data,
                           aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "AC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_t2_pea_corr_vat <- ggplot(t2_pea_corr_vat_plot_data,
                            aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "PEA\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))


p_t2_lysopaf_corr_vat <- ggplot(t2_lysopaf_corr_vat_plot_data,
                                aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "LysoPAF\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))


## --- 1.  Correlation tables -------------------------------------------
t2_lpc_cor_sat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

t2_ac_cor_sat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


t2_pea_cor_sat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


t2_lysopaf_cor_sat_tbl <- t2_vat_crs_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad), Entity %in% c("DLBCL", "MCL")) %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(SAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(SAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

## --- 2.  Long data with correlation value attached --------------------
t2_lpc_corr_sat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_lpc_cor_sat_tbl, by = "metabolite")

t2_ac_corr_sat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_ac_cor_sat_tbl, by = "metabolite")    

t2_pea_corr_sat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_pea_cor_sat_tbl, by = "metabolite")     

t2_lysopaf_corr_sat_plot_data <- t2_vat_crs_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(t2_lysopaf_cor_sat_tbl, by = "metabolite")


## --- 3.  Scatter plot coloured by correlation coefficient -------------
p_t2_lpc_corr_sat <- ggplot(t2_lpc_corr_sat_plot_data,
                            aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "LPC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_t2_ac_corr_sat <- ggplot(t2_ac_corr_sat_plot_data,
                           aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "AC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_t2_pea_corr_sat <-ggplot(t2_pea_corr_sat_plot_data,
                           aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "PEA\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))


p_t2_lysopaf_corr_sat <- ggplot(t2_lysopaf_corr_sat_plot_data,
                                aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "LysoPAF\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9),
          axis.title = element_text(size = 9)
        )

p_corr_legend <- ggplot(t2_lysopaf_corr_sat_plot_data,
                                aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Spearman r") +
  #guides(color="")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "LysoPAF\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9)
  )
ggsave("Figures_Manuscript/correlation/corr_legend.svg", plot = p_corr_legend, width = 3, height = 3)

### Function to save plots ----
# Save plots function with output directory
save_plot <- function(plot_list, output_dir, width = 5, height = 6) {
  for (name in names(plot_list)) {
    plot <- plot_list[[name]]

    # Check if the object is a ggplot
    if (inherits(plot, "ggplot")) {
      ggsave(filename = file.path(output_dir, paste0(name, ".svg")),
             plot = plot, width = width, height = height)
    } else {
      warning(paste("Object", name, "is not a ggplot. Skipping."))
    }
  }
}

# Retrieve the plots using the pattern
pattern_corr <- "p_[a-zA-Z0-9]+_[a-zA-Z0-9]+_corr_[a-zA-Z0-9]"
plots_corr <- mget(ls(pattern = pattern_corr))

# Define the directory to save the plots
output_dir <- "Figures_Manuscript/correlation"
dir.create(output_dir, showWarnings = FALSE)

# Save the plots to the specified directory
save_plot(plots_corr, output_dir, width = 3, height = 3)



###
### Comparison of training and validation cohort


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

## !! Multiple values are the same for different metabolites - need to be removed
str(meta_val_crs)

# Double-check that meta_df has the columns you expect:
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
  # Select and arrange columns as desired
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

## Validation of correlations
val_t2_lpc_cor_vat_tbl <- val_t2_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad)) %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

val_t2_ac_cor_vat_tbl <- val_t2_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad)) %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


val_t2_pea_cor_vat_tbl <- val_t2_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad)) %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )


val_t2_lysopaf_cor_vat_tbl <- val_t2_norm_bc %>% 
  filter(!is.na(Maximaler.CRS.Grad)) %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  group_by(metabolite) %>% 
  summarise(
    cor = cor(VAT, level, use = "complete.obs", method = "spearman"),
    p_value = cor.test(VAT, level, method = "spearman")$p.value,
    n = sum(!is.na(level) & !is.na(TAT)),
    .groups = "drop"
  )

## --- 2.  Long data with correlation value attached --------------------
val_t2_lpc_corr_vat_plot_data <- val_t2_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lpc_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(val_t2_lpc_cor_vat_tbl, by = "metabolite")

val_t2_ac_corr_vat_plot_data <- val_t2_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_ac_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(val_t2_ac_cor_vat_tbl, by = "metabolite")    

val_t2_pea_corr_vat_plot_data <- val_t2_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_pea_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(val_t2_pea_cor_vat_tbl, by = "metabolite")     

val_t2_lysopaf_corr_vat_plot_data <- val_t2_norm_bc %>% 
  pivot_longer(cols  = all_of(crs_lysopaf_moi),
               names_to  = "metabolite",
               values_to = "level") %>% 
  left_join(val_t2_lysopaf_cor_vat_tbl, by = "metabolite")


## --- 3.  Scatter plot coloured by correlation coefficient -------------
p_val_t2_lpc_corr_vat <- ggplot(val_t2_lpc_corr_vat_plot_data,
                            aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "LPC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_val_t2_ac_corr_vat <- ggplot(val_t2_ac_corr_vat_plot_data,
                           aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "AC\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))

p_val_t2_pea_corr_vat <-ggplot(val_t2_pea_corr_vat_plot_data,
                           aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "PEA\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))


p_val_t2_lysopaf_corr_vat <- ggplot(val_t2_lysopaf_corr_vat_plot_data,
                                aes(x = VAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Pearson r") +
  guides(color="none")+
  theme_classic() +
  labs(x = expression(VAT~"["*cm^2*"]"), y = "LysoPAF\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9)
  )

p_corr_legend <- ggplot(t2_lysopaf_corr_sat_plot_data,
                        aes(x = SAT, y = level, group = factor(metabolite), color = cor)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.4, 0.4),
                         name = "Spearman r") +
  #guides(color="")+
  theme_classic() +
  labs(x = expression(SAT~"["*cm^2*"]"), y = "LysoPAF\nserum level [norm.]") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9)
  )



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

# 3. optionally dichotomize into high vs low
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
#    – swap in your actual column names for time/event below
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

## Testing AC survival in validation cohort 
##

cox_pfs <- coxph(
  Surv(PFS_days, PFS_event) ~ AC_sig_group, 
  data = validation_master_t2
)
cox_os  <- coxph(
  Surv(OS_days, OS_event) ~ AC_sig_group, 
  data = validation_master_t2
)

pfs_ac_signature <- survfit(Surv(PFS_days/30.44, PFS_event) ~ AC_sig_group, 
                            data = validation_master_t2 |> filter(Entity != "MM"))
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

validation_master_t2$`M3.Response.(Best.ORR)`
validation_master_t2 %>%
  ggplot(aes(x = AC_sig_group, fill = factor(`M3.Response.(Best.ORR)`))) +
  geom_bar(position = "fill") +
  # scale_x_discrete(limits = c("LPC_low", "LPC_high"))+
  scale_fill_manual(values=c("#FFF5E1", "#FFDAB9", "#FFC04C", "#FF8C00"))+
  labs(x="AC signature", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

# ggplot(validation_master_t2 |> filter(Time == "Day 2-5"), aes(x=AC_signature, y=VAT))+
#   geom_point()+
#   geom_smooth(method = "lm")
# 
# cor.test(validation_master_t2$AC_signature, validation_master_t2$VAT, method = "spearman")
# 
# 
