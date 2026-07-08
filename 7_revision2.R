
library(dplyr)
library(purrr)
library(tibble)
library(pROC)

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
