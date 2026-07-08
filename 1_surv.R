# ==================================================================================================
# 1_surv.R
# Purpose: Discovery cohort baseline analyses focused on patient/body-composition characteristics,
#          baseline metabolomics, and survival-associated lipid signatures.
# Main inputs:
#   - Input_files/all_master.xlsx
#   - Input_files/0_CART_Metabolomics_Data_MASTER.xlsx
#   - Input_files/1_validation_metabolomics_master.xlsx
# Main outputs:
#   - Intermediate normalized/filtered baseline files written to Input_files/
#   - Manuscript figures written to Figures_Manuscript/
# Dependencies:
#   - Run source("0_packages.R") first.
#   - Downstream scripts reuse objects created here, especially all_master and baseline metabolite vectors.
# Notes for reviewers/readers:
#   - "BC" denotes body-composition variables, including BMI, waist, WtHR, VAT, SAT, and TAT.
#   - Metabolomics normalization is performed with MetaboAnalystR.
#   - Survival models use Cox proportional hazards regression; several later plots summarize metabolite classes.
# ==================================================================================================

### Data wrangling and patient characteristics ----
all_master <- read.xlsx("Input_files/all_master.xlsx")
all_master <- all_master |> filter(Entity != "ALL")

## Analyzing numbers
str(all_master)

all_master$Responder <- as.character(all_master$Responder)
all_master$CRS_high <- as.character(all_master$CRS_high)
all_master$Maximaler.CRS.Grad <- as.character(all_master$Maximaler.CRS.Grad)
all_master$ICANS_high <- as.character(all_master$ICANS_high)

all_master |> count(Entity)

all_master |> count(Responder)

## PFS and OS of entire cohort (Fig. S1) as well as with previously published VAT cut-off (Rejeski et al)
km_pfs_cohort <- survfit(Surv(PFS_days/30.44, PFS_event) ~ 1, data=all_master |> filter(cohort == 'training'))
p_kml_pfs_cohort <- ggsurvplot(km_pfs_cohort,
                           ylab = "Estimated PFS",
                           xlab = "Months after CAR-T infusion",
                           break.time.by = 3,
                           xlim = c(0,26),
                           censor.size = 5,
                           pval = F,
                           size = 1.5,
                           axes.offset = F,
                           risk.table = TRUE,
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
                           legend.title = "PFS",
                           palette = c("Black")
)
p_kml_pfs_cohort

km_os_cohort <- survfit(Surv(OS_days/30.44, OS_event) ~ 1, data=all_master|> filter(cohort == 'training'))
p_kml_os_cohort <- ggsurvplot(km_os_cohort,
                              ylab = "Estimated OS",
                              xlab = "Months after CAR-T infusion",
                              break.time.by = 3,
                              xlim = c(0,26),
                              censor.size = 5,
                              pval = F,
                              size = 1.5,
                              axes.offset = F,
                              risk.table = TRUE,
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
                              legend.title = "PFS",
                              palette = c("Black")
)

km_pfs_vat <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(VAT >= 61.04, "high", "low"), data=all_master |> filter(cohort == 'training'))
p_kml_pfs_vat <- ggsurvplot(km_pfs_vat,
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
                            risk.table = TRUE,
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
                            legend.title = "VAT",
                            legend.labs= c("High", "Low"),
                            palette = c("#003366", "#CC0000" )
)
p_kml_pfs_vat

km_os_vat <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(VAT >= 61.04, "high", "low"), data=all_master |> filter(cohort == 'training'))
p_kml_os_vat <- ggsurvplot(km_os_vat,
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
                           risk.table = TRUE,
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
                           legend.title = "VAT",
                           legend.labs= c("High", "Low"),
                           palette = c("#003366", "#CC0000" )
)
p_kml_os_vat

### Defining cut-offs for adipose tissue compartments based on CRS development of grade 2 higher ----
# Prepare CRS dataset
crs_data <- all_master %>%
  filter(Entity != "MM") %>%
  mutate(
    CRS_any = ifelse(Maximaler.CRS.Grad >= 1, 1, 0),
    CRS_grade2p = ifelse(Maximaler.CRS.Grad >= 2, 1, 0),
    CRS_grade3p = ifelse(Maximaler.CRS.Grad >= 3, 1, 0)
  )

# Function to calculate ROC and Youden cutoff
get_youden_cutoff <- function(data, marker, outcome) {
  
  df <- data %>%
    filter(
      !is.na(.data[[marker]]),
      !is.na(.data[[outcome]])
    )
  
  # Skip if endpoint has only one class
  if (length(unique(df[[outcome]])) < 2) {
    return(
      tibble(
        marker = marker,
        outcome = outcome,
        n = nrow(df),
        events = sum(df[[outcome]] == 1),
        auc = NA_real_,
        threshold = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        youden = NA_real_
      )
    )
  }
  
  roc_obj <- pROC::roc(
    response = df[[outcome]],
    predictor = df[[marker]],
    direction = "<",
    quiet = TRUE
  )
  
  cutoff <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  tibble(
    marker = marker,
    outcome = outcome,
    n = nrow(df),
    events = sum(df[[outcome]] == 1),
    auc = as.numeric(pROC::auc(roc_obj)),
    threshold = as.numeric(cutoff["threshold"]),
    sensitivity = as.numeric(cutoff["sensitivity"]),
    specificity = as.numeric(cutoff["specificity"]),
    youden = sensitivity + specificity - 1
  )
}

# Run for all markers and CRS definitions
youden_results <- expand.grid(
  marker = c("VAT", "SAT", "TAT"),
  outcome = c("CRS_any", "CRS_grade2p", "CRS_grade3p"),
  stringsAsFactors = FALSE
) %>%
  purrr::pmap_dfr(
    function(marker, outcome) {
      get_youden_cutoff(crs_data, marker = marker, outcome = outcome)
    }
  )

youden_results

# Adding cut-offs to all_master table
all_master <- all_master |>
  mutate(VAT_survival = ifelse(VAT > median(VAT, na.rm =T), "high", "low"),
         VAT_CRS = ifelse(VAT > 161.8, "high", "low"),
         SAT_survival = ifelse(SAT > median(SAT, na.rm =T), "high", "low"),
         SAT_CRS = ifelse(SAT > 209, "high", "low"),
         TAT_survival = ifelse(TAT > median(TAT, na.rm =T), "high", "low"),
         TAT_CRS = ifelse(TAT > 310, "high", "low"))

### Loading and cleaning of metabolite data tables ----
cohort_1 <- read_excel("Input_files/0_CART_Metabolomics_Data_MASTER.xlsx") # Discovery cohort
cohort_2 <- read_excel("Input_files/1_validation_metabolomics_master.xlsx") # Validation cohort

cohort_1 <- cohort_1%>%
  mutate(
    across(`Alanine`:`Cer-(24:01)`, as.double)
  )

cohort_2 <- cohort_2%>%
  mutate(
    across(`Alanine`:`Cer-(24:01)`, as.double)
  )


## Renaming and unifying metabolite names 
rename_map <- c(
  # One-to-one renames
  "Choline"                = "Choline (+)",
  "Phosphocholine"         = "Phosphocholine (+)",
  "Acetylcholine"          = "Acetylcholine (+)",
  "CDP-Choline"            = "CDP-Choline (+)",
  "Betaine"                = "Betaine (+)",
  "SAM"                    = "SAM (+)",
  "Decarboxy-SAM"          = "Decarboxy-SAM (+)",
  "Carnitine"             = "Carnitine (+)",
  "Methylguanosine"        = "Methylguanosine (+)",
  "NAD"                    = "NAD (+)",
  "NADP"                   = "NADP (+)",
  "Thiamine"               = "Thiamine (+)",
  "Thiamine phosphate"     = "Thiamine phosphate (+)",
  "Ergosterol"             = "Ergosterol (-H2O)",
  "Cholesterol"            = "Cholesterol (-H2O)",
  
  # Acylcarnitines: drop " (+)" from the name
  "AC-(10:0)"              = "AC-(10:0) (+)",
  "AC-(11:0)"              = "AC-(11:0) (+)",
  "AC-(12:0)"              = "AC-(12:0) (+)",
  "AC-(12:1)"              = "AC-(12:1) (+)",
  "AC-(13:0)"              = "AC-(13:0) (+)",
  "AC-(14:0)"              = "AC-(14:0) (+)",
  "AC-(14:1)"              = "AC-(14:1) (+)",
  "AC-(15:0)"              = "AC-(15:0) (+)",
  "AC-(16:0)"              = "AC-(16:0) (+)",
  "AC-(16:1)"              = "AC-(16:1) (+)",
  "AC-(17:0)"              = "AC-(17:0) (+)",
  "AC-(18:0)"              = "AC-(18:0) (+)",
  "AC-(18:1)"              = "AC-(18:1) (+)",
  "AC-(18:2)"              = "AC-(18:2) (+)",
  "AC-(18:3)"              = "AC-(18:3) (+)",
  "AC-(19:0)"              = "AC-(19:0) (+)",
  "AC-(20:0)"              = "AC-(20:0) (+)",
  "AC-(20:1)"              = "AC-(20:1) (+)",
  "AC-(20:2)"              = "AC-(20:2) (+)",
  "AC-(20:3)"              = "AC-(20:3) (+)",
  "AC-(20:4)"              = "AC-(20:4) (+)",
  "AC-(20:5)"              = "AC-(20:5) (+)",
  "AC-(21:0)"              = "AC-(21:0) (+)",
  "AC-(22:0)"              = "AC-(22:0) (+)",
  "AC-(22:1)"              = "AC-(22:1) (+)",
  "AC-(22:2)"              = "AC-(22:2) (+)",
  "AC-(22:3)"              = "AC-(22:3) (+)",
  "AC-(22:4)"              = "AC-(22:4) (+)",
  "AC-(22:5)"              = "AC-(22:5) (+)",
  "AC-(22:6)"              = "AC-(22:6) (+)",
  "AC-(23:0)"              = "AC-(23:0) (+)",
  "AC-(24:0)"              = "AC-(24:0) (+)",
  "AC-(24:1)"              = "AC-(24:1) (+)",
  "AC-(25:0)"              = "AC-(25:0) (+)",
  "AC-(26:0)"              = "AC-(26:0) (+)",
  "AC-(26:1)"              = "AC-(26:1) (+)",
  
  # MAG: remove " -H2O"
  "MAG-(14:00)"            = "MAG-(14:00) -H2O",
  "MAG-(16:00)"            = "MAG-(16:00) -H2O",
  "MAG-(16:01)"            = "MAG-(16:01) -H2O",
  "MAG-(18:00)"            = "MAG-(18:00) -H2O",
  "MAG-(18:01)"            = "MAG-(18:01) -H2O",
  "MAG-(18:02)"            = "MAG-(18:02) -H2O",
  "MAG-(18:03)"            = "MAG-(18:03) -H2O",
  "MAG-(20:00)"            = "MAG-(20:00) -H2O",
  "MAG-(20:03)"            = "MAG-(20:03) -H2O",
  "MAG-(20:04)"            = "MAG-(20:04) -H2O",
  "MAG-(20:05)"            = "MAG-(20:05) -H2O",
  "MAG-(22:00)"            = "MAG-(22:00) -H2O",
  "MAG-(22:05)"            = "MAG-(22:05) -H2O",
  "MAG-(22:06)"            = "MAG-(22:06) -H2O",
  "MAG-(24:00)"            = "MAG-(24:00) -H2O",
  "MAG-(24:01)"            = "MAG-(24:01) -H2O",
  
  # DAG: remove " -H2O"
  "DAG-(32:00)"            = "DAG-(32:00) -H2O",
  "DAG-(34:00)"            = "DAG-(34:00) -H2O",
  "DAG-(34:01)"            = "DAG-(34:01) -H2O",
  "DAG-(36:00)"            = "DAG-(36:00) -H2O",
  "DAG-(36:01)"            = "DAG-(36:01) -H2O",
  "DAG-(36:02)"            = "DAG-(36:02) -H2O",
  "DAG-(36:04)"            = "DAG-(36:04) -H2O",
  "DAG-(38:04)"            = "DAG-(38:04) -H2O",
  "DAG-(38:05)"            = "DAG-(38:05) -H2O",
  "DAG-(38:06)"            = "DAG-(38:06) -H2O",
  "DAG-(40:06)"            = "DAG-(40:06) -H2O",
  "DAG-(40:07)"            = "DAG-(40:07) -H2O",
  "DAG-(40:08)"            = "DAG-(40:08) -H2O"
)

rename_map_2 <- c(
  # One-to-one renames
  "Choline"                = "Choline.(+)",
  "Phosphocholine"         = "Phosphocholine.(+)",
  "Acetylcholine"          = "Acetylcholine.(+)",
  "CDP-Choline"            = "CDP-Choline.(+)",
  "Betaine"                = "Betaine.(+)",
  "SAM"                    = "SAM.(+)",
  "Decarboxy-SAM"          = "Decarboxy-SAM.(+)",
  "Carnitine"             = "Carnitine.(+)",
  "Methylguanosine"        = "Methylguanosine.(+)",
  "NAD"                    = "NAD.(+)",
  "NADP"                   = "NADP.(+)",
  "Thiamine"               = "Thiamine.(+)",
  "Thiamine phosphate"     = "Thiamine phosphate.(+)",
  "Ergosterol"             = "Ergosterol.(-H2O)",
  "Cholesterol"            = "Cholesterol.(-H2O)",
  
  # Acylcarnitines: drop ".(+)" from the name
  "AC-(10:0)"              = "AC-(10:0).(+)",
  "AC-(11:0)"              = "AC-(11:0).(+)",
  "AC-(12:0)"              = "AC-(12:0).(+)",
  "AC-(12:1)"              = "AC-(12:1).(+)",
  "AC-(13:0)"              = "AC-(13:0).(+)",
  "AC-(14:0)"              = "AC-(14:0).(+)",
  "AC-(14:1)"              = "AC-(14:1).(+)",
  "AC-(15:0)"              = "AC-(15:0).(+)",
  "AC-(16:0)"              = "AC-(16:0).(+)",
  "AC-(16:1)"              = "AC-(16:1).(+)",
  "AC-(17:0)"              = "AC-(17:0).(+)",
  "AC-(18:0)"              = "AC-(18:0).(+)",
  "AC-(18:1)"              = "AC-(18:1).(+)",
  "AC-(18:2)"              = "AC-(18:2).(+)",
  "AC-(18:3)"              = "AC-(18:3).(+)",
  "AC-(19:0)"              = "AC-(19:0).(+)",
  "AC-(20:0)"              = "AC-(20:0).(+)",
  "AC-(20:1)"              = "AC-(20:1).(+)",
  "AC-(20:2)"              = "AC-(20:2).(+)",
  "AC-(20:3)"              = "AC-(20:3).(+)",
  "AC-(20:4)"              = "AC-(20:4).(+)",
  "AC-(20:5)"              = "AC-(20:5).(+)",
  "AC-(21:0)"              = "AC-(21:0).(+)",
  "AC-(22:0)"              = "AC-(22:0).(+)",
  "AC-(22:1)"              = "AC-(22:1).(+)",
  "AC-(22:2)"              = "AC-(22:2).(+)",
  "AC-(22:3)"              = "AC-(22:3).(+)",
  "AC-(22:4)"              = "AC-(22:4).(+)",
  "AC-(22:5)"              = "AC-(22:5).(+)",
  "AC-(22:6)"              = "AC-(22:6).(+)",
  "AC-(23:0)"              = "AC-(23:0).(+)",
  "AC-(24:0)"              = "AC-(24:0).(+)",
  "AC-(24:1)"              = "AC-(24:1).(+)",
  "AC-(25:0)"              = "AC-(25:0).(+)",
  "AC-(26:0)"              = "AC-(26:0).(+)",
  "AC-(26:1)"              = "AC-(26:1).(+)",
  
  # MAG: remove ".-H2O"
  "MAG-(14:00)"            = "MAG-(14:00).-H2O",
  "MAG-(16:00)"            = "MAG-(16:00).-H2O",
  "MAG-(16:01)"            = "MAG-(16:01).-H2O",
  "MAG-(18:00)"            = "MAG-(18:00).-H2O",
  "MAG-(18:01)"            = "MAG-(18:01).-H2O",
  "MAG-(18:02)"            = "MAG-(18:02).-H2O",
  "MAG-(18:03)"            = "MAG-(18:03).-H2O",
  "MAG-(20:00)"            = "MAG-(20:00).-H2O",
  "MAG-(20:03)"            = "MAG-(20:03).-H2O",
  "MAG-(20:04)"            = "MAG-(20:04).-H2O",
  "MAG-(20:05)"            = "MAG-(20:05).-H2O",
  "MAG-(22:00)"            = "MAG-(22:00).-H2O",
  "MAG-(22:05)"            = "MAG-(22:05).-H2O",
  "MAG-(22:06)"            = "MAG-(22:06).-H2O",
  "MAG-(24:00)"            = "MAG-(24:00).-H2O",
  "MAG-(24:01)"            = "MAG-(24:01).-H2O",
  
  # DAG: remove ".-H2O"
  "DAG-(32:00)"            = "DAG-(32:00).-H2O",
  "DAG-(34:00)"            = "DAG-(34:00).-H2O",
  "DAG-(34:01)"            = "DAG-(34:01).-H2O",
  "DAG-(36:00)"            = "DAG-(36:00).-H2O",
  "DAG-(36:01)"            = "DAG-(36:01).-H2O",
  "DAG-(36:02)"            = "DAG-(36:02).-H2O",
  "DAG-(36:04)"            = "DAG-(36:04).-H2O",
  "DAG-(38:04)"            = "DAG-(38:04).-H2O",
  "DAG-(38:05)"            = "DAG-(38:05).-H2O",
  "DAG-(38:06)"            = "DAG-(38:06).-H2O",
  "DAG-(40:06)"            = "DAG-(40:06).-H2O",
  "DAG-(40:07)"            = "DAG-(40:07).-H2O",
  "DAG-(40:08)"            = "DAG-(40:08).-H2O"
)

cohort_2 <- cohort_2 %>%
  rename(!!!rename_map)

qc_cohort_1 <- cohort_1 |> 
  filter(is.na(Sample)) |>
  select(-Sample, -Timepoint)|>
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "quality")

qc_cohort_2 <- cohort_2 |> 
  filter(Sample == "NA") |>
  select(-Sample, -Time, -Kombi)|>
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "quality")

## Removing metabolites with quality score of 4 or higher
qc_mb_cohort_1 <- qc_cohort_1 |> filter(quality <= 3) |> select(metabolite) |> unlist() |> as.vector()
qc_mb_cohort_2 <- qc_cohort_2 |> filter(quality <= 3) |> select(metabolite) |> unlist() |> as.vector()

qc_mb_cohort_combined <- union(qc_mb_cohort_1, qc_mb_cohort_2) ## These are the columns to keep in a unified data frame

###
### Creating a unified data frame containing the data from day 0 only ----
### 
bl_cohort_1 <- cohort_1 |>
  filter(!is.na(Sample))|>
  filter(Timepoint == 1)|>
  select(-Timepoint)

bl_cohort_2 <- cohort_2 |>
  filter(Sample != "NA") |>
  filter(Time == "Day 0") |>
  select(-Sample, -Time) |>
  rename("Sample" = "Kombi")

common_cols <- intersect(colnames(bl_cohort_1), colnames(bl_cohort_2))

bl_cohort_1 <- bl_cohort_1[ , common_cols]
bl_cohort_2 <- bl_cohort_2[ , common_cols]

# missing_in_2 <- setdiff(names(bl_cohort_1), names(bl_cohort_2))
# missing_in_1 <- setdiff(names(bl_cohort_2), names(bl_cohort_1))
# bl_cohort_2[missing_in_2] <- NA
# bl_cohort_1[missing_in_1] <- NA

bl_cohort_merged <- bind_rows(bl_cohort_1, bl_cohort_2)

bl_cohort_filter <- bl_cohort_merged |> select(Sample, any_of(qc_mb_cohort_combined))

## Adding Sample Name from Metadata table
# Removing two sample ids without day 0 measurements
bl_sample_id <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(Sample) |> unlist() |> as.vector() 

bl_cohort_filter$Sample_new <- bl_sample_id
bl_cohort_filter <- bl_cohort_filter |> select(Sample, Sample_new, everything())

VAT_survival_vector <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(VAT_survival) |> unlist() |> as.vector() 
VAT_crs_vector <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(VAT_CRS) |> unlist() |> as.vector() 
SAT_survival_vector <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(SAT_survival) |> unlist() |> as.vector() 
SAT_crs_vector <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(SAT_CRS) |> unlist() |> as.vector() 
TAT_survival_vector <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(TAT_survival) |> unlist() |> as.vector() 
TAT_crs_vector <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(TAT_CRS) |> unlist() |> as.vector() 

###
### Creating a combined data frame for day 3-5 (2nd time point) for CRS analysis ----
###

## Creating a unified data frame containing the data from day 0 only
t2_cohort_1 <- cohort_1 |>
  filter(!is.na(Sample))|>
  filter(Timepoint == 2)|>
  select(-Timepoint)

t2_cohort_2 <- cohort_2 |>
  filter(Sample != "NA") |>
  filter(Time == "Day 2-5") |>
  select(-Sample, -Time) |>
  rename("Sample" = "Kombi")

## Keeping only shared metabolites
t2_cohort_1 <- t2_cohort_1[ , common_cols]
t2_cohort_2 <- t2_cohort_2[ , common_cols]

## Merging data sets
t2_cohort_merged <- bind_rows(t2_cohort_1, t2_cohort_2)

##Excluding ALL sample
t2_cohort_merged <- t2_cohort_merged |> filter(Sample != "VSample_B_9")

## Selecting only metabolites passing Lamivudin QC
t2_cohort_filter <- t2_cohort_merged |> select(Sample, any_of(qc_mb_cohort_combined))

## Adding Sample Name from Metadata table
t2_sample_id <- all_master |> filter(!(Sample %in% c("A60", "A66", "A72", "A75", "A78"))) |> select(Sample) |> unlist() |> as.vector() #Removing two sample ids without time point 2 

t2_cohort_filter$Sample_new <- t2_sample_id
t2_cohort_filter <- t2_cohort_filter |> select(Sample, Sample_new, everything())


## Creating new columns with cutoffs for VAT_survival and VAT_CRS
all_master <- all_master |>
  mutate(VAT_survival = ifelse(VAT > 61.04, "high", "low"),
         VAT_CRS = ifelse(VAT > 161.8, "high", "low"))

VAT_survival_vector <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(VAT_survival) |> unlist() |> as.vector() 
VAT_crs_vector <- all_master |> filter(Sample != "A60") |> filter(Sample != "A73") |> select(VAT_CRS) |> unlist() |> as.vector() 

# Removing patients with missing metabolite measurements
t2_VAT_crs_vector <- all_master |> filter(!(Sample %in% c("A60", "A66", "A72", "A75", "A78"))) |> select(VAT_CRS) |> unlist() |> as.vector() 
t2_SAT_crs_vector <- all_master |> filter(!(Sample %in% c("A60", "A66", "A72", "A75", "A78"))) |> select(SAT_CRS) |> unlist() |> as.vector() 
t2_TAT_crs_vector <- all_master |> filter(!(Sample %in% c("A60", "A66", "A72", "A75", "A78"))) |> select(TAT_CRS) |> unlist() |> as.vector() 

## Creating data frames including different VAT, SAT and TAT groups for metabolomic analyses at baseline and saving as cvs
bl_cohort_filter$VAT_survival <- VAT_survival_vector
bl_cohort_filter$VAT_CRS <- VAT_crs_vector
bl_cohort_filter$SAT_survival <- SAT_survival_vector
bl_cohort_filter$SAT_CRS <- SAT_crs_vector
bl_cohort_filter$TAT_survival <- TAT_survival_vector
bl_cohort_filter$TAT_CRS <- TAT_crs_vector

bl_cohort_VAT_survival <- bl_cohort_filter |> 
  dplyr::select(Sample_new, VAT_survival, everything(), -Sample, -VAT_CRS, -SAT_survival, -SAT_CRS, -TAT_survival, -TAT_CRS) |>
  rename("Sample" = "Sample_new")|>
  mutate(VAT_survival = ifelse(VAT_survival == "high", 1, 0))|>
  filter(!is.na(Alanine))|>
  filter(!is.na(VAT_survival))|>
  slice(1:(n() - 18))
write.csv(bl_cohort_VAT_survival, "Input_files/bl_cohort_VAT_survival.csv", row.names = F)

bl_cohort_VAT_crs <- bl_cohort_filter |> 
  dplyr::select(Sample_new, VAT_CRS, everything(), -Sample, -VAT_survival, -SAT_survival, -SAT_CRS, -TAT_survival, -TAT_CRS) |>
  rename("Sample" = "Sample_new")|>
  mutate(VAT_CRS = ifelse(VAT_CRS == "high", 1, 0))|>
  filter(!is.na(Alanine))|>
  filter(!is.na(VAT_CRS))|>
  slice(1:(n() - 18))
write.csv(bl_cohort_VAT_crs, "Input_files/bl_cohort_VAT_crs.csv", row.names = F)

bl_cohort_SAT_survival <- bl_cohort_filter |> 
  dplyr::select(Sample_new, SAT_survival, everything(), -Sample, -VAT_survival, -VAT_CRS, -SAT_CRS, -TAT_survival, -TAT_CRS) |>
  rename("Sample" = "Sample_new")|>
  mutate(SAT_survival = ifelse(SAT_survival == "high", 1, 0))|>
  filter(!is.na(Alanine))|>
  filter(!is.na(SAT_survival))|>
  slice(1:(n() - 18))
write.csv(bl_cohort_SAT_survival, "Input_files/bl_cohort_SAT_survival.csv", row.names = F)

bl_cohort_SAT_crs <- bl_cohort_filter |> 
  dplyr::select(Sample_new, SAT_CRS, everything(), -Sample, -VAT_survival, -VAT_CRS, -SAT_survival, -TAT_survival, -TAT_CRS) |>
  rename("Sample" = "Sample_new")|>
  mutate(SAT_CRS = ifelse(SAT_CRS == "high", 1, 0))|>
  filter(!is.na(Alanine))|>
  filter(!is.na(SAT_CRS))|>
  slice(1:(n() - 18))
write.csv(bl_cohort_SAT_crs, "Input_files/bl_cohort_SAT_crs.csv", row.names = F)

bl_cohort_TAT_survival <- bl_cohort_filter |> 
  dplyr::select(Sample_new, TAT_survival, everything(), -Sample, -VAT_survival, -VAT_CRS,-SAT_survival, -SAT_CRS, -TAT_CRS) |>
  rename("Sample" = "Sample_new")|>
  mutate(TAT_survival = ifelse(TAT_survival == "high", 1, 0))|>
  filter(!is.na(Alanine))|>
  filter(!is.na(TAT_survival))|>
  slice(1:(n() - 18))
write.csv(bl_cohort_TAT_survival, "Input_files/bl_cohort_TAT_survival.csv", row.names = F)

bl_cohort_TAT_crs <- bl_cohort_filter |> 
  dplyr::select(Sample_new, VAT_CRS, everything(), -Sample, -VAT_survival, -VAT_CRS,-SAT_survival, -SAT_CRS, -TAT_survival) |>
  rename("Sample" = "Sample_new")|>
  mutate(TAT_CRS = ifelse(TAT_CRS == "high", 1, 0))|>
  filter(!is.na(Alanine))|>
  filter(!is.na(TAT_CRS))|>
  slice(1:(n() - 18))
write.csv(bl_cohort_TAT_crs, "Input_files/bl_cohort_TAT_crs.csv", row.names = F)


## Creating data frames including different VAT, SAT and TAT groups for metabolomic analyses at timepoint 2 and saving as cvs
t2_cohort_filter$VAT_CRS <- t2_VAT_crs_vector
t2_cohort_filter$SAT_CRS <- t2_SAT_crs_vector
t2_cohort_filter$TAT_CRS <- t2_TAT_crs_vector

t2_cohort_VAT_crs <- t2_cohort_filter |>
  dplyr::select(Sample_new, VAT_CRS, everything(), -Sample, -SAT_CRS, -TAT_CRS) |>
    rename("Sample" = "Sample_new")|>
    mutate(VAT_CRS = ifelse(VAT_CRS == "high", 1, 0))|>
    filter(!is.na(Alanine))|>
    filter(!is.na(VAT_CRS))|>
    slice(1:(n() - 15))
write.csv(t2_cohort_VAT_crs, "Input_files/t2_cohort_VAT_crs.csv", row.names = F)

t2_cohort_SAT_crs <- t2_cohort_filter |>
  dplyr::select(Sample_new, SAT_CRS, everything(), -Sample, -VAT_CRS, -TAT_CRS) |>
  rename("Sample" = "Sample_new")|>
  mutate(SAT_CRS = ifelse(SAT_CRS == "high", 1, 0))|>
  filter(!is.na(Alanine))|>
  filter(!is.na(SAT_CRS))|>
  slice(1:(n() - 15))
write.csv(t2_cohort_SAT_crs, "Input_files/t2_cohort_SAT_crs.csv", row.names = F)


t2_cohort_TAT_crs <- t2_cohort_filter |>
  dplyr::select(Sample_new, TAT_CRS, everything(), -Sample, -VAT_CRS, -SAT_CRS) |>
  rename("Sample" = "Sample_new")|>
  mutate(TAT_CRS = ifelse(TAT_CRS == "high", 1, 0))|>
  filter(!is.na(Alanine))|>
  filter(!is.na(TAT_CRS))|>
  slice(1:(n() - 15))
write.csv(t2_cohort_TAT_crs, "Input_files/t2_cohort_TAT_crs.csv", row.names = F)


###
### Identification of metabolites from VAT differences for survival and feature selection ----
### 

## 1. Analysis of day 0 differences to extract metabolites

## 1.1 Loading of data for normalization and analyses
mSet_blvatsurv<-InitDataObjects("pktable", "stat", FALSE)
mSet_blvatsurv<-Read.TextData(mSet_blvatsurv, "Input_files/bl_cohort_VAT_survival.csv", "rowu", "disc");
mSet_blvatsurv<-SanityCheckData(mSet_blvatsurv)
mSet_blvatsurv<-ReplaceMin(mSet_blvatsurv);
mSet_blvatsurv<-SanityCheckData(mSet_blvatsurv)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_blvatsurv<-FilterVariable(mSet_blvatsurv, "median", 0, "F")
mSet_blvatsurv<-PreparePrenormData(mSet_blvatsurv)

## Normalization by sum and data scaling based on auto-scaling
mSet_blvatsurv<-Normalization(mSet_blvatsurv, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

# Extraction and saving of the normalized data into a new tibble
bl_vat_surv <- as.data.frame(mSet_blvatsurv[["dataSet"]][["norm"]])

# Load original table and cbind sample name and time point label
bl_vat_surv_original <- read.csv("Input_files/bl_cohort_VAT_survival.csv", na = "NA")

## Left join original data and normalized to link vat groups
bl_vat_surv$Sample <- row.names(bl_vat_surv)
bl_vat_surv <- bl_vat_surv %>% select(Sample, everything()) %>%
  arrange(Sample)

bl_vat_surv_sample_id <- bl_vat_surv_original %>% select(Sample, VAT_survival)

bl_vat_surv_norm <- left_join(bl_vat_surv, bl_vat_surv_sample_id, by = "Sample")
bl_vat_surv_norm <- bl_vat_surv_norm %>%
  select(Sample, VAT_survival, everything())

bl_vat_surv_norm$VAT_survival <- as.character(bl_vat_surv_norm$VAT_survival)

## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 0.05, and FC > 1.5
mSet_blvatsurv<-Volcano.Anal(mSet_blvatsurv, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

bl_vat_surv_volcano <- as.data.frame(mSet_blvatsurv[["analSet"]][["volcano"]][["fc.log"]])
bl_vat_surv_volcano$metabolite <- rownames(bl_vat_surv_volcano)
bl_vat_surv_volcano$log_p <- mSet_blvatsurv[["analSet"]][["volcano"]][["p.log"]]
bl_vat_surv_volcano$log_fc <- mSet_blvatsurv[["analSet"]][["volcano"]][["fc.log"]]
bl_vat_surv_volcano$inx.up <- mSet_blvatsurv[["analSet"]][["volcano"]][["inx.up"]]
bl_vat_surv_volcano$inx.down <- mSet_blvatsurv[["analSet"]][["volcano"]][["inx.down"]]
bl_vat_surv_volcano$inx.p <- mSet_blvatsurv[["analSet"]][["volcano"]][["inx.p"]]

bl_vat_surv_volcano <- bl_vat_surv_volcano %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

p_volcano_bl_vat_surv <- ggplot(bl_vat_surv_volcano, aes(x = log_fc, y = log_p))+
  geom_point(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 30)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p-value"))+
  guides(color = "none") +
  # xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

p_volcano_bl_vat_surv
ggsave("Figures_Manuscript/volcano_example.svg", plot = p_volcano_bl_vat_surv, width = 3, height =3)

bl_vat_surv_volcano_metabolites <- bl_vat_surv_volcano %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()

### Adding PCA analysis
mSet_blvatsurv<-PCA.Anal(mSet_blvatsurv)

## Extraction of PCA component values
bl_vat_surv_PCA <- as.data.frame(mSet_blvatsurv[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_blvatsurv order for samples (used for PLSDA)
bl_vat_surv_PCA_sample_order <- bl_vat_surv_PCA %>%
  mutate(Sample = rownames(bl_vat_surv_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
bl_vat_surv_PCA$Sample <- rownames(bl_vat_surv_PCA)
bl_vat_surv_PCA <- bl_vat_surv_PCA %>%
  select(Sample, everything())

bl_vat_surv_PCA <- left_join(bl_vat_surv_PCA, bl_vat_surv_sample_id, by = "Sample")

bl_vat_surv_PCA <- bl_vat_surv_PCA %>%
  select(Sample, VAT_survival, everything())

bl_vat_surv_PCA$VAT_survival <- as.character(bl_vat_surv_PCA$VAT_survival)

## Visualization of PCA comp1 vs comp2
p_pca_bl_vat_surv <- bl_vat_surv_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = VAT_survival)) +
  geom_jitter() +  
  geom_text_repel(aes(label = Sample), max.overlaps = 30)+
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  #  coord_fixed(ratio = 1)+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

p_pca_bl_vat_surv


## Perform PLSDA
mSet_blvatsurv<-PLSR.Anal(mSet_blvatsurv, reg=TRUE)

## Extraction of PLSDA component values
bl_vat_surv_PLSDA <- as.matrix.data.frame(bl_vat_surv_PLSDA <- mSet_blvatsurv[["analSet"]][["plsr"]][["scores"]])
bl_vat_surv_PLSDA <- as.data.frame(bl_vat_surv_PLSDA)
# Sample information lost in mSet_blvatsurv upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
bl_vat_surv_PLSDA <- bl_vat_surv_PLSDA %>%
  mutate(Sample = bl_vat_surv_PCA_sample_order)

bl_vat_surv_PLSDA <- left_join(bl_vat_surv_PLSDA, bl_vat_surv_sample_id, by = "Sample")

bl_vat_surv_PLSDA$VAT_survival <- as.character(bl_vat_surv_PLSDA$VAT_survival)

## Visualization of PLSDA comp1 vs comp2
p_plsda_bl_vat_surv <- bl_vat_surv_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = VAT_survival)) +
  geom_jitter() +
  geom_text_repel(aes(label = Sample), max.overlaps = 30)+
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  #scale_color_manual(values = c("black", "darkmagenta"))+
  scale_color_manual(
    values = c("black", "darkmagenta"),
    labels = c("Low", "High"),
    name = "VAT Survival"
  ) +
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

p_plsda_bl_vat_surv
ggsave("Figures_Manuscript/plsda_example.svg", plot = p_plsda_bl_vat_surv, width = 5, height =5)


bl_vat_surv_PLSDA_VIP <- as.data.frame(mSet_blvatsurv[["analSet"]][["plsr"]][["vip.mat"]])
bl_vat_surv_PLSDA_VIP <- bl_vat_surv_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

p_plasda_vip_bl_vat_surv <- bl_vat_surv_PLSDA_VIP %>%
  arrange(desc(`Comp. 1`)) %>%
  filter(`Comp. 1` > 2) %>%
  ggplot(aes(x = reorder(metabolite,`Comp. 1`), y = `Comp. 1`, fill = `Comp. 1`)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "PLS-DA Comp.1 VIP Score")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill = "none")+
  coord_flip() +
  scale_fill_gradient(low = "lightgrey", high = "darkmagenta")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

p_plasda_vip_bl_vat_surv
ggsave("Figures_Manuscript/plsda_vip_example.svg", plot = p_plasda_vip_bl_vat_surv, width = 3, height =3)


bl_vat_surv_plsda_metabolites <- bl_vat_surv_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) %>%
  select(metabolite)|>
  unlist() %>%
  as.vector()

## Combining volcano and plsda metabolites

bl_vat_surv_metabolites <- union(bl_vat_surv_plsda_metabolites, bl_vat_surv_volcano_metabolites)

###

###
### Identification of metabolites from SAT differences for survival and feature selection ----
### 

## 1. Analysis of day 0 differences to extract metabolites

## 1.1 Loading of data for normalization and analyses
mSet_blsatsurv<-InitDataObjects("pktable", "stat", FALSE)
mSet_blsatsurv<-Read.TextData(mSet_blsatsurv, "Input_files/bl_cohort_SAT_survival.csv", "rowu", "disc");
mSet_blsatsurv<-SanityCheckData(mSet_blsatsurv)
mSet_blsatsurv<-ReplaceMin(mSet_blsatsurv);
mSet_blsatsurv<-SanityCheckData(mSet_blsatsurv)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_blsatsurv<-FilterVariable(mSet_blsatsurv, "median", 0, "F")
mSet_blsatsurv<-PreparePrenormData(mSet_blsatsurv)

## Normalization by sum and data scaling based on auto-scaling
mSet_blsatsurv<-Normalization(mSet_blsatsurv, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

# Extraction and saving of the normalized data into a new tibble
bl_sat_surv <- as.data.frame(mSet_blsatsurv[["dataSet"]][["norm"]])

# Load original table and cbind sample name and time point label
bl_sat_surv_original <- read.csv("Input_files/bl_cohort_SAT_survival.csv", na = "NA")

## Left join original data and normalized to link vat groups
bl_sat_surv$Sample <- row.names(bl_sat_surv)
bl_sat_surv <- bl_sat_surv %>% select(Sample, everything()) %>%
  arrange(Sample)

bl_sat_surv_sample_id <- bl_sat_surv_original %>% select(Sample, SAT_survival)

bl_sat_surv_norm <- left_join(bl_sat_surv, bl_sat_surv_sample_id, by = "Sample")
bl_sat_surv_norm <- bl_sat_surv_norm %>%
  select(Sample, SAT_survival, everything())

bl_sat_surv_norm$SAT_survival <- as.character(bl_sat_surv_norm$SAT_survival)

## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 0.05, and FC > 1.5
mSet_blsatsurv<-Volcano.Anal(mSet_blsatsurv, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

bl_sat_surv_volcano <- as.data.frame(mSet_blsatsurv[["analSet"]][["volcano"]][["fc.log"]])
bl_sat_surv_volcano$metabolite <- rownames(bl_sat_surv_volcano)
bl_sat_surv_volcano$log_p <- mSet_blsatsurv[["analSet"]][["volcano"]][["p.log"]]
bl_sat_surv_volcano$log_fc <- mSet_blsatsurv[["analSet"]][["volcano"]][["fc.log"]]
bl_sat_surv_volcano$inx.up <- mSet_blsatsurv[["analSet"]][["volcano"]][["inx.up"]]
bl_sat_surv_volcano$inx.down <- mSet_blsatsurv[["analSet"]][["volcano"]][["inx.down"]]
bl_sat_surv_volcano$inx.p <- mSet_blsatsurv[["analSet"]][["volcano"]][["inx.p"]]

bl_sat_surv_volcano <- bl_sat_surv_volcano %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

p_volcano_bl_sat_surv <- ggplot(bl_sat_surv_volcano, aes(x = log_fc, y = log_p))+
  geom_point(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 30)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p-value"))+
  guides(color = "none") +
  # xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

p_volcano_bl_sat_surv

bl_sat_surv_volcano_metabolites <- bl_sat_surv_volcano %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()

### Adding PCA analysis
mSet_blsatsurv<-PCA.Anal(mSet_blsatsurv)

## Extraction of PCA component values
bl_sat_surv_PCA <- as.data.frame(mSet_blsatsurv[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_blsatsurv order for samples (used for PLSDA)
bl_sat_surv_PCA_sample_order <- bl_sat_surv_PCA %>%
  mutate(Sample = rownames(bl_sat_surv_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
bl_sat_surv_PCA$Sample <- rownames(bl_sat_surv_PCA)
bl_sat_surv_PCA <- bl_sat_surv_PCA %>%
  select(Sample, everything())

bl_sat_surv_PCA <- left_join(bl_sat_surv_PCA, bl_sat_surv_sample_id, by = "Sample")

bl_sat_surv_PCA <- bl_sat_surv_PCA %>%
  select(Sample, SAT_survival, everything())

bl_sat_surv_PCA$SAT_survival <- as.character(bl_sat_surv_PCA$SAT_survival)

## Visualization of PCA comp1 vs comp2
p_pca_bl_sat_surv <- bl_sat_surv_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = SAT_survival)) +
  geom_jitter() +  
  geom_text_repel(aes(label = Sample), max.overlaps = 30)+
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  #  coord_fixed(ratio = 1)+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

p_pca_bl_sat_surv


## Perform PLSDA
mSet_blsatsurv<-PLSR.Anal(mSet_blsatsurv, reg=TRUE)

## Extraction of PLSDA component values
bl_sat_surv_PLSDA <- as.matrix.data.frame(bl_sat_surv_PLSDA <- mSet_blsatsurv[["analSet"]][["plsr"]][["scores"]])
bl_sat_surv_PLSDA <- as.data.frame(bl_sat_surv_PLSDA)
# Sample information lost in mSet_blsatsurv upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
bl_sat_surv_PLSDA <- bl_sat_surv_PLSDA %>%
  mutate(Sample = bl_sat_surv_PCA_sample_order)

bl_sat_surv_PLSDA <- left_join(bl_sat_surv_PLSDA, bl_sat_surv_sample_id, by = "Sample")

bl_sat_surv_PLSDA$SAT_survival <- as.character(bl_sat_surv_PLSDA$SAT_survival)

## Visualization of PLSDA comp1 vs comp2
p_plsda_bl_sat_surv <- bl_sat_surv_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = SAT_survival)) +
  geom_jitter() +
  geom_text_repel(aes(label = Sample), max.overlaps = 30)+
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("black", "darkmagenta"))+
  scale_fill_manual(values = c("black", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  #guides(color = "none")+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

p_plsda_bl_sat_surv

bl_sat_surv_PLSDA_VIP <- as.data.frame(mSet_blsatsurv[["analSet"]][["plsr"]][["vip.mat"]])
bl_sat_surv_PLSDA_VIP <- bl_sat_surv_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

p_plasda_vip_bl_sat_surv <- bl_sat_surv_PLSDA_VIP %>%
  arrange(desc(`Comp. 1`)) %>%
  filter(`Comp. 1` > 2) %>%
  ggplot(aes(x = reorder(metabolite,`Comp. 1`), y = `Comp. 1`, fill = `Comp. 1`)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "PLS-DA Comp.1 VIP Score")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill = "none")+
  coord_flip() +
  scale_fill_gradient(low = "lightgrey", high = "darkmagenta")+
  theme_classic()

p_plasda_vip_bl_sat_surv

bl_sat_surv_plsda_metabolites <- bl_sat_surv_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) %>%
  select(metabolite)|>
  unlist() %>%
  as.vector()

## Combining volcano and plsda metabolites

bl_sat_surv_metabolites <- union(bl_sat_surv_plsda_metabolites, bl_sat_surv_volcano_metabolites)


###
### Identification of metabolites from TAT differences for survival and feature selection ----
###
## 1. Analysis of day 0 differences to extract metabolites

## 1.1 Loading of data for normalization and analyses
mSet_bltatsurv<-InitDataObjects("pktable", "stat", FALSE)
mSet_bltatsurv<-Read.TextData(mSet_bltatsurv, "Input_files/bl_cohort_TAT_survival.csv", "rowu", "disc");
mSet_bltatsurv<-SanityCheckData(mSet_bltatsurv)
mSet_bltatsurv<-ReplaceMin(mSet_bltatsurv);
mSet_bltatsurv<-SanityCheckData(mSet_bltatsurv)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_bltatsurv<-FilterVariable(mSet_bltatsurv, "median", 0, "F")
mSet_bltatsurv<-PreparePrenormData(mSet_bltatsurv)

## Normalization by sum and data scaling based on auto-scaling
mSet_bltatsurv<-Normalization(mSet_bltatsurv, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

# Extraction and saving of the normalized data into a new tibble
bl_tat_surv <- as.data.frame(mSet_bltatsurv[["dataSet"]][["norm"]])

# Load original table and cbind sample name and time point label
bl_tat_surv_original <- read.csv("Input_files/bl_cohort_TAT_survival.csv", na = "NA")

## Left join original data and normalized to link vat groups
bl_tat_surv$Sample <- row.names(bl_tat_surv)
bl_tat_surv <- bl_tat_surv %>% select(Sample, everything()) %>%
  arrange(Sample)

bl_tat_surv_sample_id <- bl_tat_surv_original %>% select(Sample, TAT_survival)

bl_tat_surv_norm <- left_join(bl_tat_surv, bl_tat_surv_sample_id, by = "Sample")
bl_tat_surv_norm <- bl_tat_surv_norm %>%
  select(Sample, TAT_survival, everything())

bl_tat_surv_norm$TAT_survival <- as.character(bl_tat_surv_norm$TAT_survival)

## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 0.05, and FC > 1.5
mSet_bltatsurv<-Volcano.Anal(mSet_bltatsurv, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

bl_tat_surv_volcano <- as.data.frame(mSet_bltatsurv[["analSet"]][["volcano"]][["fc.log"]])
bl_tat_surv_volcano$metabolite <- rownames(bl_tat_surv_volcano)
bl_tat_surv_volcano$log_p <- mSet_bltatsurv[["analSet"]][["volcano"]][["p.log"]]
bl_tat_surv_volcano$log_fc <- mSet_bltatsurv[["analSet"]][["volcano"]][["fc.log"]]
bl_tat_surv_volcano$inx.up <- mSet_bltatsurv[["analSet"]][["volcano"]][["inx.up"]]
bl_tat_surv_volcano$inx.down <- mSet_bltatsurv[["analSet"]][["volcano"]][["inx.down"]]
bl_tat_surv_volcano$inx.p <- mSet_bltatsurv[["analSet"]][["volcano"]][["inx.p"]]

bl_tat_surv_volcano <- bl_tat_surv_volcano %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

p_volcano_bl_tat_surv <- ggplot(bl_tat_surv_volcano, aes(x = log_fc, y = log_p))+
  geom_point(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 30)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p-value"))+
  guides(color = "none") +
  # xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

p_volcano_bl_tat_surv

bl_tat_surv_volcano_metabolites <- bl_tat_surv_volcano %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()

### Adding PCA analysis
mSet_bltatsurv<-PCA.Anal(mSet_bltatsurv)

## Extraction of PCA component values
bl_tat_surv_PCA <- as.data.frame(mSet_bltatsurv[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_bltatsurv order for samples (used for PLSDA)
bl_tat_surv_PCA_sample_order <- bl_tat_surv_PCA %>%
  mutate(Sample = rownames(bl_tat_surv_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
bl_tat_surv_PCA$Sample <- rownames(bl_tat_surv_PCA)
bl_tat_surv_PCA <- bl_tat_surv_PCA %>%
  select(Sample, everything())

bl_tat_surv_PCA <- left_join(bl_tat_surv_PCA, bl_tat_surv_sample_id, by = "Sample")

bl_tat_surv_PCA <- bl_tat_surv_PCA %>%
  select(Sample, TAT_survival, everything())

bl_tat_surv_PCA$TAT_survival <- as.character(bl_tat_surv_PCA$TAT_survival)

## Visualization of PCA comp1 vs comp2
p_pca_bl_tat_surv <- bl_tat_surv_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = TAT_survival)) +
  geom_jitter() +  
  geom_text_repel(aes(label = Sample), max.overlaps = 30)+
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  #  coord_fixed(ratio = 1)+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

p_pca_bl_tat_surv


## Perform PLSDA
mSet_bltatsurv<-PLSR.Anal(mSet_bltatsurv, reg=TRUE)

## Extraction of PLSDA component values
bl_tat_surv_PLSDA <- as.matrix.data.frame(bl_tat_surv_PLSDA <- mSet_bltatsurv[["analSet"]][["plsr"]][["scores"]])
bl_tat_surv_PLSDA <- as.data.frame(bl_tat_surv_PLSDA)
# Sample information lost in mSet_bltatsurv upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
bl_tat_surv_PLSDA <- bl_tat_surv_PLSDA %>%
  mutate(Sample = bl_tat_surv_PCA_sample_order)

bl_tat_surv_PLSDA <- left_join(bl_tat_surv_PLSDA, bl_tat_surv_sample_id, by = "Sample")

bl_tat_surv_PLSDA$TAT_survival <- as.character(bl_tat_surv_PLSDA$TAT_survival)

## Visualization of PLSDA comp1 vs comp2
p_plsda_bl_tat_surv <- bl_tat_surv_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = TAT_survival)) +
  geom_jitter() +
  geom_text_repel(aes(label = Sample), max.overlaps = 30)+
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("black", "darkmagenta"))+
  scale_fill_manual(values = c("black", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  #guides(color = "none")+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

p_plsda_bl_tat_surv

bl_tat_surv_PLSDA_VIP <- as.data.frame(mSet_bltatsurv[["analSet"]][["plsr"]][["vip.mat"]])
bl_tat_surv_PLSDA_VIP <- bl_tat_surv_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

p_plasda_vip_bl_tat_surv <- bl_tat_surv_PLSDA_VIP %>%
  arrange(desc(`Comp. 1`)) %>%
  filter(`Comp. 1` > 2) %>%
  ggplot(aes(x = reorder(metabolite,`Comp. 1`), y = `Comp. 1`, fill = `Comp. 1`)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "PLS-DA Comp.1 VIP Score")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill = "none")+
  coord_flip() +
  scale_fill_gradient(low = "lightgrey", high = "darkmagenta")+
  theme_classic()

p_plasda_vip_bl_tat_surv

bl_tat_surv_plsda_metabolites <- bl_tat_surv_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) %>%
  select(metabolite)|>
  unlist() %>%
  as.vector()

## Combining volcano and plsda metabolites

bl_tat_surv_metabolites <- union(bl_tat_surv_plsda_metabolites, bl_tat_surv_volcano_metabolites)


###
### Survival: Correlation of filtered survival metabolites with BC levels (Figure 1F, S4) ----
###

## Comparison of VAT metabolites with VAT measurements
bl_vat_surv_norm_bc <- left_join(bl_vat_surv_norm, all_master, by="Sample")
str(bl_vat_surv_norm_bc)

# Create an empty data frame to store correlation results
bl_vat_surv_corr <- data.frame(
  Metabolite = character(),
  Pearson_r  = numeric(),
  p_value    = numeric(),
  stringsAsFactors = FALSE
)

# Create a list to hold all ggplots
plot_list_bl_vat_surv_corr <- list()

# Loop over each metabolite column
for (met in bl_vat_surv_metabolites) {
  # Perform a Pearson correlation test against VAT
  corr_test <- cor.test(
    bl_vat_surv_norm_bc[[met]],
    bl_vat_surv_norm_bc[["VAT"]],
    method = "pearson"
  )
  
  # Extract the correlation coefficient (r) and p-value
  r_val <- corr_test$estimate
  p_val <- corr_test$p.value
  
  # Append the numeric results to 'results' data frame
  bl_vat_surv_corr <- rbind(
    bl_vat_surv_corr,
    data.frame(
      Metabolite = met,
      Pearson_r  = r_val,
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
  )
  
  # Create a scatter plot with regression line
  # Use aes(x = VAT, y = .data[[met]]) instead of aes_string()
  p <- ggplot(bl_vat_surv_norm_bc, aes(x = VAT, y = .data[[met]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
   # ggtitle(paste("Correlation of", met, "with VAT")) +
    labs(
      subtitle = paste0(
        "Pearson r = ", round(r_val, 3), 
        ", p-value = ", signif(p_val, 3)
      )
    )
  
  # Store each plot in the list, naming it by metabolite
  plot_list_bl_vat_surv_corr[[met]] <- p
}

# Create the arranged grid
grid_combined <- arrangeGrob(grobs = plot_list_bl_vat_surv_corr, ncol = 4)

# Save as SVG
ggsave("Figures_Manuscript/met_corr_VAT.svg", plot = grid_combined, width = 20, height = 24)

## SAT correlations
# Create an empty data frame to store correlation results
bl_sat_surv_corr <- data.frame(
  Metabolite = character(),
  Pearson_r  = numeric(),
  p_value    = numeric(),
  stringsAsFactors = FALSE
)

# Create a list to hold all ggplots
plot_list_bl_sat_surv_corr <- list()

# Loop over each metabolite column
for (met in bl_sat_surv_metabolites) {
  # Perform a Pearson correlation test against VAT
  corr_test <- cor.test(
    bl_vat_surv_norm_bc[[met]],
    bl_vat_surv_norm_bc[["SAT"]],
    method = "pearson"
  )
  
  # Extract the correlation coefficient (r) and p-value
  r_val <- corr_test$estimate
  p_val <- corr_test$p.value
  
  # Append the numeric results to 'results' data frame
  bl_sat_surv_corr <- rbind(
    bl_sat_surv_corr,
    data.frame(
      Metabolite = met,
      Pearson_r  = r_val,
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
  )
  
  # Create a scatter plot with regression line
  # Use aes(x = VAT, y = .data[[met]]) instead of aes_string()
  p <- ggplot(bl_vat_surv_norm_bc, aes(x = VAT, y = .data[[met]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    # ggtitle(paste("Correlation of", met, "with VAT")) +
    labs(
      subtitle = paste0(
        "Pearson r = ", round(r_val, 3), 
        ", p-value = ", signif(p_val, 3)
      )
    )
  
  # Store each plot in the list, naming it by metabolite
  plot_list_bl_sat_surv_corr[[met]] <- p
}

# Create the arranged grid
grid_combined <- arrangeGrob(grobs = plot_list_bl_sat_surv_corr, ncol = 4)

# Save as SVG
ggsave("Figures_Manuscript/met_corr_SAT.svg", plot = grid_combined, width = 20, height = 24)


## TAT correlations
# Create an empty data frame to store correlation results
bl_tat_surv_corr <- data.frame(
  Metabolite = character(),
  Pearson_r  = numeric(),
  p_value    = numeric(),
  stringsAsFactors = FALSE
)

# Create a list to hold all ggplots
plot_list_bl_tat_surv_corr <- list()

# Loop over each metabolite column
for (met in bl_tat_surv_metabolites) {
  # Perform a Pearson correlation test against VAT
  corr_test <- cor.test(
    bl_vat_surv_norm_bc[[met]],
    bl_vat_surv_norm_bc[["TAT"]],
    method = "pearson"
  )
  
  # Extract the correlation coefficient (r) and p-value
  r_val <- corr_test$estimate
  p_val <- corr_test$p.value
  
  # Append the numeric results to 'results' data frame
  bl_tat_surv_corr <- rbind(
    bl_tat_surv_corr,
    data.frame(
      Metabolite = met,
      Pearson_r  = r_val,
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
  )
  
  # Create a scatter plot with regression line
  # Use aes(x = TAT, y = .data[[met]]) instead of aes_string()
  p <- ggplot(bl_vat_surv_norm_bc, aes(x = TAT, y = .data[[met]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    # ggtitle(paste("Correlation of", met, "with VAT")) +
    labs(
      subtitle = paste0(
        "Pearson r = ", round(r_val, 3), 
        ", p-value = ", signif(p_val, 3)
      )
    )
  
  # Store each plot in the list, naming it by metabolite
  plot_list_bl_tat_surv_corr[[met]] <- p
}

# Create the arranged grid
grid_combined <- arrangeGrob(grobs = plot_list_bl_tat_surv_corr, ncol = 4)

# Save as SVG
ggsave("Figures_Manuscript/met_corr_TAT.svg", plot = grid_combined, width = 20, height = 24)

## Filtering for metabolites which are associated with BCs with a Pearson of at least 0.2 or -0.2

bl_vat_surv_corr_filter <- bl_vat_surv_corr |>
  filter(Pearson_r <= -0.2 | Pearson_r >= 0.2)|>
  pull(Metabolite) |>
  unname()
bl_sat_surv_corr_filter <- bl_sat_surv_corr |>
  filter(Pearson_r <= -0.2 | Pearson_r >= 0.2)|>
  pull(Metabolite) |>
  unname()
bl_tat_surv_corr_filter <- bl_tat_surv_corr |>
  filter(Pearson_r <= -0.2 | Pearson_r >= 0.2)|>
  pull(Metabolite) |>
  unname()

###
### Survival: Common and different metabolites to the VAT adipose depot ----
###

# Common to all three
surv_common_all <- sort(Reduce(intersect, list(
  bl_vat_surv_corr_filter,
  bl_sat_surv_corr_filter,
  bl_tat_surv_corr_filter
)))

# Unique to each
surv_only_in_tat <- sort(setdiff(bl_tat_surv_corr_filter, union(bl_vat_surv_corr_filter, bl_sat_surv_corr_filter)))
surv_only_in_vat <- sort(setdiff(bl_vat_surv_corr_filter, union(bl_tat_surv_corr_filter, bl_sat_surv_corr_filter)))
surv_only_in_sat <- sort(setdiff(bl_sat_surv_corr_filter, union(bl_tat_surv_corr_filter, bl_vat_surv_corr_filter)))

# Shared between any two but not all three
surv_tat_vat_shared <- sort(setdiff(intersect(bl_tat_surv_corr_filter, bl_vat_surv_corr_filter), surv_common_all))
surv_tat_sat_shared <- sort(setdiff(intersect(bl_tat_surv_corr_filter, bl_sat_surv_corr_filter), surv_common_all))
surv_vat_sat_shared <- sort(setdiff(intersect(bl_vat_surv_corr_filter, bl_sat_surv_corr_filter), surv_common_all))

# Organize output
list(
  common_to_all_three = surv_common_all,
  only_in_tat = surv_only_in_tat,
  only_in_vat = surv_only_in_vat,
  only_in_sat = surv_only_in_sat,
  shared_tat_vat = surv_tat_vat_shared,
  shared_tat_sat = surv_tat_sat_shared,
  shared_vat_sat = surv_vat_sat_shared
)

bl_all_surv_corr_filter <- sort(unique(c(
  bl_vat_surv_corr_filter,
  bl_sat_surv_corr_filter,
  bl_tat_surv_corr_filter
)))

###
### Associations with survival using COX models(Figure 1C, S2) ----
###

str(bl_vat_surv_norm_bc)
## Calculation of Cox model for day 0 metabolite levels
bl_all_cox_adjusted <- data.frame(marker = character(),
                                HR = numeric(),
                                lower95 = numeric(),
                                higher95 = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (i in bl_all_surv_corr_filter) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(PFS_days, PFS_event)  ~ `",i,"` + Costim + STLV + Geschlecht")
  model_cox <- coxph(as.formula(formula_str), data = bl_vat_surv_norm_bc)
  
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
  bl_all_cox_adjusted  <- rbind(bl_all_cox_adjusted, marker_results_cox)
}

# bl_all_cox_adjusted <- bl_all_cox_adjusted |>
#   mutate(FDR = p.adjust(p_value, method = "fdr"))


str(bl_all_cox_adjusted)

bl_all_cox_adjusted <- bl_all_cox_adjusted %>%
  mutate(group = sapply(marker, get_group_new)) 

bl_all_cox_adjusted <- bl_all_cox_adjusted |> filter(!(marker %in% c("PI-(38:07)", "PI-(40:03)","PI-(40:08)","PI-(40:09)")))

bl_all_cox_adjusted <- bl_all_cox_adjusted |>
  arrange(group, desc(HR)) |> 
  mutate(marker = factor(marker, levels = unique(marker)))

p_surv_group <- bl_all_cox_adjusted |>
  #filter(FDR <= 0.1)|>
  ggplot()+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = higher95), alpha = 0.8)+
  geom_point(aes(x = marker, y = HR, color = group),
             size = 3, shape = 19, alpha = 0.8)+
  coord_flip()+ #ylim = c(-1, 4)
  labs(y = "HR (95%CI)", x = "")+
  scale_x_discrete(limits = rev)+
  labs(color = "Group")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

print(p_surv_group)

ggsave("Figures_Manuscript/meta_surv_group.svg", plot = p_surv_group, width = 6, height = 9)


# Now create the columns needed for meta-analysis
meta_survival <- bl_all_cox_adjusted %>%
  mutate(
    logHR    = log(HR),
    logLower = log(lower95),
    logUpper = log(higher95),
    # approximate standard error for log(HR)
    SE_logHR = (logUpper - logLower) / (2 * 1.96)
  )


# Double-check that meta_df
meta_survival <- meta_survival |> filter(!(marker %in% c("PI-(38:07)", "PI-(40:03)","PI-(40:08)","PI-(40:09)")))

meta_results_survival <- meta_survival %>%
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

meta_surv_significant <- c("Acylcarnitine", "Phosphatidylethanolamine", "Diacylglycerol", "Sphingomyelin", "Plasmalogen" )   # groups to be red


p_meta_surv_group <- meta_results_survival |>
  mutate(col_group = if_else(group %in% meta_surv_significant, "highlight", "other")) %>% 
  filter(group != "Other/Unclassified")|>
  #filter(FDR <= 0.1)|>
  ggplot()+
  geom_linerange(aes(x = reorder(group, combined_HR), ymin = lower95, 
                     ymax = upper95), size=0.8, alpha = 0.8)+
  geom_point(aes(x = reorder(group, combined_HR), y = combined_HR, size = n_mets, color=col_group),
             shape = 19, alpha = 0.8)+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #ylim = c(-1, 4)
  labs(y = "Combined HR (95%CI)", x = "")+
  scale_color_manual(values = c("highlight" = "#CC0000", "other" = "black"), guide = "none") +
  scale_size_continuous(breaks = c(1,5,10), range = c(2,7))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  guides(color = "none")+
  labs(size="Metabolites")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

print(p_meta_surv_group) #(Figure 1C)
ggsave("Figures_Manuscript/meta_surv.svg", plot = p_meta_surv_group, width = 6, height = 4.5)



###
### Survival: Distribution of lipids based on adipose tissue compartment (Fig. 1E) ----
###

# Distribution of ACs
get_ACs <- function(x) {
  grep("^AC-", x, value = TRUE)
}

surv_ac_common_all <- get_ACs(surv_common_all)
surv_ac_only_tat <- get_ACs(surv_only_in_tat)
surv_ac_only_vat <- get_ACs(surv_only_in_vat)
surv_ac_only_sat <- get_ACs(surv_only_in_sat)
surv_ac_tat_vat <- get_ACs(surv_tat_vat_shared)
surv_ac_tat_sat <- get_ACs(surv_tat_sat_shared)
surv_ac_vat_sat <- get_ACs(surv_vat_sat_shared)

surv_ac_df <- data.frame(
  group = c(
    rep("common_all", length(surv_ac_common_all)),
    rep("only_tat", length(surv_ac_only_tat)),
    rep("only_vat", length(surv_ac_only_vat)),
    rep("only_sat", length(surv_ac_only_sat)),
    rep("tat_vat", length(surv_ac_tat_vat)),
    rep("tat_sat", length(surv_ac_tat_sat)),
    rep("vat_sat", length(surv_ac_vat_sat))
  ),
  metabolite = c(
    surv_ac_common_all,
    surv_ac_only_tat,
    surv_ac_only_vat,
    surv_ac_only_sat,
    surv_ac_tat_vat,
    surv_ac_tat_sat,
    surv_ac_vat_sat
  )
)

surv_ac_summary <- surv_ac_df %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) |>
  mutate(met_group = "AC")

surv_ac_vector <- surv_ac_df %>%
  select(metabolite)|>
  unlist()|>
  as.vector()


# Distribution of Phosphatidylethanolamine
get_PEAs <- function(x) {
  grep("^PEA-", x, value = TRUE)
}

surv_pea_common_all <- get_PEAs(surv_common_all)
surv_pea_only_tat <- get_PEAs(surv_only_in_tat)
surv_pea_only_vat <- get_PEAs(surv_only_in_vat)
surv_pea_only_sat <- get_PEAs(surv_only_in_sat)
surv_pea_tat_vat <- get_PEAs(surv_tat_vat_shared)
surv_pea_tat_sat <- get_PEAs(surv_tat_sat_shared)
surv_pea_vat_sat <- get_PEAs(surv_vat_sat_shared)

surv_pea_df <- data.frame(
  group = c(
    rep("common_all", length(surv_pea_common_all)),
    rep("only_tat", length(surv_pea_only_tat)),
    rep("only_vat", length(surv_pea_only_vat)),
    rep("only_sat", length(surv_pea_only_sat)),
    rep("tat_vat", length(surv_pea_tat_vat)),
    rep("tat_sat", length(surv_pea_tat_sat)),
    rep("vat_sat", length(surv_pea_vat_sat))
  ),
  metabolite = c(
    surv_pea_common_all,
    surv_pea_only_tat,
    surv_pea_only_vat,
    surv_pea_only_sat,
    surv_pea_tat_vat,
    surv_pea_tat_sat,
    surv_pea_vat_sat
  )
)

surv_pea_summary <- surv_pea_df %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) |>
  mutate(met_group = "PEA")

surv_pea_vector <- surv_pea_df %>%
  select(metabolite)|>
  unlist()|>
  as.vector()

# Distribution of Sphingomyeline
get_SMs <- function(x) {
  grep("^SM-", x, value = TRUE)
}

surv_sm_common_all <- get_SMs(surv_common_all)
surv_sm_only_tat <- get_SMs(surv_only_in_tat)
surv_sm_only_vat <- get_SMs(surv_only_in_vat)
surv_sm_only_sat <- get_SMs(surv_only_in_sat)
surv_sm_tat_vat <- get_SMs(surv_tat_vat_shared)
surv_sm_tat_sat <- get_SMs(surv_tat_sat_shared)
surv_sm_vat_sat <- get_SMs(surv_vat_sat_shared)

surv_sm_df <- data.frame(
  group = c(
    rep("common_all", length(surv_sm_common_all)),
    rep("only_tat", length(surv_sm_only_tat)),
    rep("only_vat", length(surv_sm_only_vat)),
    rep("only_sat", length(surv_sm_only_sat)),
    rep("tat_vat", length(surv_sm_tat_vat)),
    rep("tat_sat", length(surv_sm_tat_sat)),
    rep("vat_sat", length(surv_sm_vat_sat))
  ),
  metabolite = c(
    surv_sm_common_all,
    surv_sm_only_tat,
    surv_sm_only_vat,
    surv_sm_only_sat,
    surv_sm_tat_vat,
    surv_sm_tat_sat,
    surv_sm_vat_sat
  )
)

surv_sm_summary <- surv_sm_df %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) |>
  mutate(met_group = "SM")

surv_sm_vector <- surv_sm_df %>%
  select(metabolite)|>
  unlist()|>
  as.vector()

# Distribution of Sphingomyelins
get_Plass <- function(x) {
  grep("^Plas", x, value = TRUE)
}

surv_plas_common_all <- get_Plass(surv_common_all)
surv_plas_only_tat <- get_Plass(surv_only_in_tat)
surv_plas_only_vat <- get_Plass(surv_only_in_vat)
surv_plas_only_sat <- get_Plass(surv_only_in_sat)
surv_plas_tat_vat <- get_Plass(surv_tat_vat_shared)
surv_plas_tat_sat <- get_Plass(surv_tat_sat_shared)
surv_plas_vat_sat <- get_Plass(surv_vat_sat_shared)

surv_plas_df <- data.frame(
  group = c(
    rep("common_all", length(surv_plas_common_all)),
    rep("only_tat", length(surv_plas_only_tat)),
    rep("only_vat", length(surv_plas_only_vat)),
    rep("only_sat", length(surv_plas_only_sat)),
    rep("tat_vat", length(surv_plas_tat_vat)),
    rep("tat_sat", length(surv_plas_tat_sat)),
    rep("vat_sat", length(surv_plas_vat_sat))
  ),
  metabolite = c(
    surv_plas_common_all,
    surv_plas_only_tat,
    surv_plas_only_vat,
    surv_plas_only_sat,
    surv_plas_tat_vat,
    surv_plas_tat_sat,
    surv_plas_vat_sat
  )
)

surv_plas_summary <- surv_plas_df %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) |>
  mutate(met_group = "Plas")

surv_plas_vector <- surv_plas_df %>%
  select(metabolite)|>
  unlist()|>
  as.vector()


## Combining data frames and create a bar plot
surv_at_distr <- rbind(surv_ac_summary, surv_plas_summary, surv_sm_summary, surv_pea_summary) %>%
  group_by(met_group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(met_group = factor(met_group, levels = c("AC", "SM", "Plas",  "PEA")))

# New grouping variable
surv_at_distr<- surv_at_distr |>
  mutate(group_new = case_when(
    group == "only_tat" ~ "AT-shared",
    group == "only_vat" ~ "VAT-correlated",
    group == "only_sat" ~ "SAT-correlated",
    group == "tat_vat" ~ "VAT-enriched",
    group == "tat_sat" ~ "SAT-enriched",
    group == "vat_sat" ~ "AT-shared",
    group == "common_all" ~ "AT-shared",
  ))

group_colors_AT <- c(
  `AT-shared` = "#2196F3",
  `VAT-correlated` = "#4CAF50",
  `VAT-enriched` = "#009688",
  `SAT-correlated` = "#673AB7",
  `SAT-enriched` = "#9C27B0"
)

# Plot from the correctly prepared data
ggplot(surv_at_distr, aes(x = met_group, y = percentage, fill = group)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Metabolites [%]", fill = "Group") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

surv_at_distr <- surv_at_distr |>
  mutate(group_new = factor(group_new, levels = c("AT-shared" ,
                                                  "VAT-correlated",
                                                  "VAT-enriched",
                                                  "SAT-correlated",
                                                  "SAT-enriched" )))

p_surv_distr <- ggplot(surv_at_distr, aes(x = met_group, y = percentage, fill = group_new)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Metabolites [%]", fill = "Source") +
  scale_fill_manual(values = group_colors_AT) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


print(p_surv_distr)
ggsave("Figures_Manuscript/p_surv_distr.svg", plot = p_surv_distr, width = 5, height = 3)

###
### Kaplan-Meier plots for representative metabolites for each group (Fig. 1D, S3) ---- 
###
bl_vat_surv_norm_bc$PFS_days
bl_vat_surv_norm_bc$PFS_event
bl_vat_surv_norm_bc$`Plas-(d18:1/18:01)`
bl_vat_surv_norm_bc$EASIX

coxph(Surv(PFS_days, PFS_event) ~ `PEA-(34:03)`, data=bl_vat_surv_norm_bc)

pfs_PEA3403 <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`PEA-(34:03)` > median(`PEA-(34:03)`), "PEA_high", "PEA_low"), data=bl_vat_surv_norm_bc)
summary(pfs_PEA3403)
summary(pfs_PEA3403, times = 365)
p_km_pfs_PEA3403 <- ggsurvplot(pfs_PEA3403,
           ylab = "Estimated PFS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
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
           legend.title = "PEA-(34:03)",
           legend.labs= c("High", "Low"),
           palette = c("black","darkgrey")
)
ggsave(filename = "Figures_Manuscript/km/p_km_pfs_PEA3403.svg", plot = p_km_pfs_PEA3403$plot,
       width = 3, height = 2.5)

os_PEA3403 <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`PEA-(34:03)` > median(`PEA-(34:03)`), "PEA_high", "PEA_low"), data=bl_vat_surv_norm_bc)
summary(os_PEA3403)
summary(os_PEA3403, times = 365)
p_km_os_PEA3403 <- ggsurvplot(os_PEA3403,
           ylab = "Estimated OS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
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
           legend.title = "PEA-(34:03)",
           legend.labs= c("High", "Low"),
           palette = c("black","darkgrey")
)
ggsave(filename = "Figures_Manuscript/km/p_km_os_PEA3403.svg", plot = p_km_os_PEA3403$plot,
       width = 3, height = 2.5)

pfs_SM181 <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`SM-(d18:1/18:02)` > median(`SM-(d18:1/18:02)` ), "SMhigh", "SMlow"), data=bl_vat_surv_norm_bc)
summary(pfs_SM181)
summary(pfs_SM181, times = 365)
p_km_pfs_SM181 <- ggsurvplot(pfs_SM181,
           ylab = "Estimated PFS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
           censor.size = 5,
           pval = T,
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
           legend.title = "SM-(d18:1/18:02)",
           legend.labs= c("High", "Low"),
           palette = c("black","darkgrey")
)
ggsave(filename = "Figures_Manuscript/km/p_km_pfs_SM181.svg", plot = p_km_pfs_SM181$plot,
       width = 3, height = 2.5)

os_SM181 <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`SM-(d18:1/18:02)` > median(`SM-(d18:1/18:02)` ), "SMhigh", "SMlow"), data=bl_vat_surv_norm_bc)
summary(os_SM181)
summary(os_SM181, times = 365)
p_km_os_SM181 <- ggsurvplot(os_SM181,
           ylab = "Estimated OS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
           censor.size = 5,
           pval = T,
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
           legend.title = "SM-(d18:1/18:02)",
           legend.labs= c("High", "Low"),
           palette = c("black","darkgrey")
)
ggsave(filename = "Figures_Manuscript/km/p_km_os_SM181.svg", plot = p_km_os_SM181$plot,
       width = 3, height = 2.5)

pfs_Plas181 <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`PlasEA-(38:06)` > median(`PlasEA-(38:06)` ), "Plas_high", "Plas_low"), data=bl_vat_surv_norm_bc)
summary(pfs_Plas181)
summary(pfs_Plas181, times = 365)
p_km_pfs_Plas181 <- ggsurvplot(pfs_Plas181,
           ylab = "Estimated PFS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
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
           legend.title = "PlasEA-(38:06)",
           legend.labs= c("High", "Low"),
           palette = c("black","darkgrey")
)
ggsave(filename = "Figures_Manuscript/km/p_km_pfs_Plas181.svg", plot = p_km_pfs_Plas181$plot,
       width = 3, height = 2.5)

os_Plas181 <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`PlasEA-(38:06)` > median(`PlasEA-(38:06)` ), "Plas_high", "Plas_low"), data=bl_vat_surv_norm_bc)
summary(os_Plas181)
summary(os_Plas181, times = 365)
p_km_os_Plas181 <- ggsurvplot(os_Plas181,
           ylab = "Estimated OS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
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
           legend.title = "PlasEA-(38:06)",
           legend.labs= c("High", "Low"),
           palette = c("black","darkgrey")
)
ggsave(filename = "Figures_Manuscript/km/p_km_os_Plas181.svg", plot = p_km_os_Plas181$plot,
       width = 3, height = 2.5)

pfs_ac221 <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`AC-(22:1)` > median(`AC-(22:1)` ), "Plas_high", "Plas_low"), data=bl_vat_surv_norm_bc)
summary(pfs_ac221)
summary(pfs_ac221, times = 365)
p_km_pfs_ac221 <- ggsurvplot(pfs_ac221,
           ylab = "Estimated PFS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
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
           legend.title = "AC-(22:1)",
           legend.labs= c("High", "Low"),
           palette = c("black","darkgrey")
)
ggsave(filename = "Figures_Manuscript/km/p_km_pfs_ac221.svg", plot = p_km_pfs_ac221$plot,
       width = 3, height = 2.5)


os_ac221 <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`AC-(22:1)` > median(`AC-(22:1)` ), "Plas_high", "Plas_low"), data=bl_vat_surv_norm_bc)
summary(os_ac221)
summary(os_ac221, times = 365)
p_km_os_ac221 <- ggsurvplot(os_ac221,
           ylab = "Estimated OS",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
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
           legend.title = "AC-(22:1)",
           legend.labs= c("High", "Low"),
           palette = c("black","darkgrey")
)

ggsave(filename = "Figures_Manuscript/km/p_km_os_ac221.svg", plot = p_km_os_ac221$plot,
       width = 3, height = 2.5)


pfs_dag <- survfit(Surv(PFS_days/30.44, PFS_event) ~ ifelse(`DAG-(36:04)` > median(`DAG-(36:04)` ), "DAG_high", "DAG_low"), data=bl_vat_surv_norm_bc)
summary(pfs_dag)
summary(pfs_dag, times = 365)
p_km_pfs_dag <- ggsurvplot(pfs_dag,
                             ylab = "Estimated PFS",
                             xlab = "Months after CAR-T infusion",
                             break.time.by = 3,
                             xlim = c(0,26),
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
                             legend.title = "DAG-(36:04)",
                             legend.labs= c("High", "Low"),
                             palette = c("black","darkgrey")
)
p_km_pfs_dag
ggsave(filename = "Figures_Manuscript/km/p_km_pfs_dag.svg", plot = p_km_pfs_dag$plot,
       width = 3, height = 2.5)

os_dag <- survfit(Surv(OS_days/30.44, OS_event) ~ ifelse(`DAG-(36:04)` > median(`DAG-(36:04)` ), "DAG_high", "DAG_low"), data=bl_vat_surv_norm_bc)
summary(os_dag)
summary(os_dag, times = 365)
p_km_os_dag <- ggsurvplot(os_dag,
                            ylab = "Estimated OS",
                            xlab = "Months after CAR-T infusion",
                            break.time.by = 3,
                            xlim = c(0,26),
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
                            legend.title = "DAG-(36:04)",
                            legend.labs= c("High", "Low"),
                            palette = c("black","darkgrey")
)

ggsave(filename = "Figures_Manuscript/km/p_km_os_dag.svg", plot = p_km_os_dag$plot,
       width = 3, height = 2.5)

# ### Function to save plots ----
# Retrieve the plots using the pattern
pattern_km <- "p_km_[a-zA-Z0-9]+_[a-zA-Z0-9]"
plots_km <- mget(ls(pattern = pattern_km))

pattern_kml <- "p_kml_[a-zA-Z0-9]+_[a-zA-Z0-9]" #large
plots_kml <- mget(ls(pattern = pattern_kml))

# Define the directory to save the plots
output_dir <- "Figures_Manuscript/km"
dir.create(output_dir, showWarnings = FALSE)


for (name in names(plots_km)) {
  p_obj <- plots_km[[name]]
  combined_plot <- ggarrange(
    p_obj$plot
  )
  
  # Define filename/path
  subset_file_name <- paste0(name, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  message("Saving plot to: ", subset_file_path)
  
  # Save the combined figure
  ggsave(
    filename = subset_file_path,
    plot     = combined_plot,
    width    = 3,
    height   = 2.5
  )
}

for (name in names(plots_kml)) {
  p_obj <- plots_kml[[name]]
  combined_plot <- ggarrange(
    p_obj$plot,
    p_obj$table,
    nrow = 2,
    heights = c(3, 0.7) # controls relative space
  )
  
  # Define filename/path
  subset_file_name <- paste0(name, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  message("Saving plot to: ", subset_file_path)
  
  # Save the combined figure
  ggsave(
    filename = subset_file_path,
    plot     = combined_plot,
    width    = 4,
    height   = 3.5
  )
}


