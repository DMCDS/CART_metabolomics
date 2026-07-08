## Revisions

# cutoff table
# Define cutoffs
vat_crs_cutoff <- 162
sat_crs_cutoff <- 166   # replace with Youden-derived SAT cutoff if different
tat_crs_cutoff <- 310

# Prepare adiposity grouping dataset
adiposity_df <- all_master %>%
  filter(Entity != "MM") %>%
  mutate(
    VAT_survival_group = ifelse(VAT > median(VAT, na.rm = TRUE), "High", "Low"),
    SAT_survival_group = ifelse(SAT > median(SAT, na.rm = TRUE), "High", "Low"),
    TAT_survival_group = ifelse(TAT > median(TAT, na.rm = TRUE), "High", "Low"),
    
    VAT_CRS_group = ifelse(VAT >= vat_crs_cutoff, "High", "Low"),
    SAT_CRS_group = ifelse(SAT >= sat_crs_cutoff, "High", "Low"),
    TAT_CRS_group = ifelse(TAT >= tat_crs_cutoff, "High", "Low")
  )

median_iqr <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  med <- median(x)
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  
  paste0(
    round(med, 1),
    " (",
    round(q1, 1),
    ", ",
    round(q3, 1),
    ")"
  )
}


make_adiposity_summary <- function(data, group_var, group_label) {
  
  data %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      Grouping_variable = group_label,
      Group = first(.data[[group_var]]),
      n = n(),
      VAT_median_IQR = median_iqr(VAT),
      SAT_median_IQR = median_iqr(SAT),
      TAT_median_IQR = median_iqr(TAT),
      BMI_median_IQR = median_iqr(BMI),
      .groups = "drop"
    ) %>%
    select(
      Grouping_variable,
      Group,
      n,
      VAT_median_IQR,
      SAT_median_IQR,
      TAT_median_IQR,
      BMI_median_IQR
    )
}

supp_adiposity_table <- bind_rows(
  make_adiposity_summary(adiposity_df, "VAT_survival_group", "VAT survival group"),
  make_adiposity_summary(adiposity_df, "SAT_survival_group", "SAT survival group"),
  make_adiposity_summary(adiposity_df, "TAT_survival_group", "TAT survival group"),
  make_adiposity_summary(adiposity_df, "VAT_CRS_group", "VAT CRS group"),
  make_adiposity_summary(adiposity_df, "SAT_CRS_group", "SAT CRS group"),
  make_adiposity_summary(adiposity_df, "TAT_CRS_group", "TAT CRS group")
) %>%
  mutate(
    Group = factor(Group, levels = c("Low", "High"))
  ) %>%
  arrange(Grouping_variable, Group)

supp_adiposity_table