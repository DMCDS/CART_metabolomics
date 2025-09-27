library(readxl)
library(openxlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(purrr)
library(plotly)
library(stringr)
library(ComplexHeatmap)
library(ggpubr)
library(gridExtra)
library(hrbrthemes)
library(ggprism)
library(stringr)
library(MetaboAnalystR)
library(pls)
library(reshape2)
library(igraph)

library(RJSONIO)
library(pROC)

library(ggvenn)
library(ggVennDiagram)
library(RColorBrewer)
library(cutpointr)
library(survival)
library(survminer)
library(tidycmprsk)
library(ggsurvfit)
library(randomForest)
library(datasets)
library(caret)
library(metafor)
library(ggpattern)

#install.packages("Rtsne")
#library(Rtsne)
#install.packages("mdatools")
#library(mdatools)
#library(rapportools)
#install.packages("rapport")
# install.packages("remotes")
# library(remotes)
# remotes::install_github("xia-lab/MetaboAnalystR")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#library(preprocessCore)
# library(fuzzyjoin)

### Functions

## Building groups
get_group_new <- function(x) {
  x_upper <- toupper(x)          # unify case
  
  ## --------------- existing groups --------------- ##
  # 1. Acyl-carnitines (incl. free carnitine)
  if (grepl("^AC-\\(", x_upper) || x_upper == "CARNITINE") {
    return("Acylcarnitine")
    
    # 2. Sphingomyelin
  } else if (grepl("^SM-\\(", x_upper)) {
    return("Sphingomyelin")
    
    # 3. Phosphatidylcholine
  } else if (grepl("^PC-\\(", x_upper)) {
    return("Phosphatidylcholine")
    
    # 4. Phosphatidylethanolamine
  } else if (grepl("^PEA-\\(", x_upper)) {
    return("Phosphatidylethanolamine")
    
    # 5. Phosphatidylinositol
  } else if (grepl("^PI-\\(", x_upper)) {
    return("Phosphatidylinositol")
    
    # 6. Phosphatidylserine
  } else if (grepl("^PS-\\(", x_upper)) {
    return("Phosphatidylserine")
    
    # 7. Phosphatidic acid
  } else if (grepl("^PA-\\(", x_upper)) {
    return("Phosphatidic acid")
    
    # 8. Lyso-PAs / LPEAs / LPS / LPC / LPI -------------
  } else if (grepl("^LPA-\\(", x_upper)) {
    return("Lysophosphatidic acid")
  } else if (grepl("^LPEA-\\(", x_upper)) {
    return("Lysophosphatidylethanolamine")
  } else if (grepl("^LPS-\\(", x_upper)) {
    return("Lysophosphatidylserine")
  } else if (grepl("^LPC-\\(", x_upper)) {
    return("Lysophosphatidylcholine")
  } else if (grepl("^LPI-\\(", x_upper)) {
    return("Lysophosphatidylinositol")
  } else if (grepl("^LYSOPAF-", x_upper)) {
    return("Lysophospholipid / PAF")
    
    # 9. Plasmalogens
  } else if (grepl("^PLASC-\\(", x_upper) || grepl("^PLASEA-\\(", x_upper)) {
    return("Plasmalogen")
    
    # 10. Neutral glycerolipids
  } else if (grepl("^TAG-", x_upper)) {
    return("Triacylglycerol")
  } else if (grepl("^DAG-\\(", x_upper)) {
    return("Diacylglycerol")
    
    # 11. Fatty / carboxylic acids
  } else if (grepl("^CA-\\(", x_upper) || grepl("^FA-\\(", x_upper)) {
    return("Fatty acid / Carboxylic acid")
    
    # 12. Amino-acid & N-derivatives
  } else if (x_upper %in% c("CREATINE","CREATININE","CYSTINE","CYSTATHIONINE",
                            "AMINOADIPIC ACID","GLUARG","GLYCINE","ACETYLGLYCINE",
                            "DIMETHYLARGININE","HYPOTAURINE","N-AC-MET")) {
    return("Amino acid & derivatives")
    
    # 13. TCA / central-carbon intermediates
  } else if (x_upper %in% c("A-KETOGLUTARATE","ACONITATE","PYRUVATE")) {
    return("TCA cycle")
    
    # 14. Ketone bodies
  } else if (x_upper == "ACETOACETATE") {
    return("Ketone body")
    
    # 15. Organic acids (keep existing label)
  } else if (x_upper %in% c("3-HYDROXYMETHYLGLUTARATE")) {
    return("Organic acid")
    
    # 16. Vitamin derivatives
  } else if (x_upper %in% c("DEHYDROASCORBATE","PYRIDOXAL","MNAM")) {
    return("Vitamin derivative")
    
    # 17. Carbohydrates / sugars
  } else if (x_upper %in% c("ALDOHEXOSE","GLUCOSAMINE","FUCOSE","DRIBOSE-P")) {
    return("Carbohydrate / Sugar")
    
    # 18. Phospho intermediates
  } else if (x_upper %in% c("CDP-ETHANOLAMINE","PHOSPHOETHANOLAMINE")) {
    return("Phospho intermediate")
    
    # 19. Catch-all: any new sterol-lipid classes (CE-, CER-, THCA) or anything else
  } else {
    return("Other/Unclassified")
  }
}
