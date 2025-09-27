### Regradeting analyses for CRS
###
### Identification of metabolites from VAT survival analysis ----
### 

## 1. Analysis of day 0 differences to extract metabolites

## 1.1 Loading of data for normalization and analyses
mSet_t2vatcrs<-InitDataObjects("pktable", "stat", FALSE)
mSet_t2vatcrs<-Read.TextData(mSet_t2vatcrs, "Input_files/t2_cohort_VAT_crs.csv", "rowu", "disc");
mSet_t2vatcrs<-SanityCheckData(mSet_t2vatcrs)
mSet_t2vatcrs<-ReplaceMin(mSet_t2vatcrs);
mSet_t2vatcrs<-SanityCheckData(mSet_t2vatcrs)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_t2vatcrs<-FilterVariable(mSet_t2vatcrs, "median", 0, "F")
mSet_t2vatcrs<-PreparePrenormData(mSet_t2vatcrs)

## Normalization by sum and data scaling based on auto-scaling
mSet_t2vatcrs<-Normalization(mSet_t2vatcrs, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

# Extraction and saving of the normalized data into a new tibble
t2_vat_crs <- as.data.frame(mSet_t2vatcrs[["dataSet"]][["norm"]])

# Load original table and cbind sample name and time point label
t2_vat_crs_original <- read.csv("Input_files/t2_cohort_VAT_crs.csv", na = "NA")

## Left join original data and normalized to link vat groups
t2_vat_crs$Sample <- row.names(t2_vat_crs)
t2_vat_crs <- t2_vat_crs %>% select(Sample, everything()) %>%
  arrange(Sample)

t2_vat_crs_sample_id <- t2_vat_crs_original %>% select(Sample, VAT_CRS)

t2_vat_crs_norm <- left_join(t2_vat_crs, t2_vat_crs_sample_id, by = "Sample")
t2_vat_crs_norm <- t2_vat_crs_norm %>%
  select(Sample, VAT_CRS, everything())

t2_vat_crs_norm$VAT_CRS <- as.character(t2_vat_crs_norm$VAT_CRS)

## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 0.05, and FC > 1.5
mSet_t2vatcrs<-Volcano.Anal(mSet_t2vatcrs, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

t2_vat_crs_volcano <- as.data.frame(mSet_t2vatcrs[["analSet"]][["volcano"]][["fc.log"]])
t2_vat_crs_volcano$metabolite <- rownames(t2_vat_crs_volcano)
t2_vat_crs_volcano$log_p <- mSet_t2vatcrs[["analSet"]][["volcano"]][["p.log"]]
t2_vat_crs_volcano$log_fc <- mSet_t2vatcrs[["analSet"]][["volcano"]][["fc.log"]]
t2_vat_crs_volcano$inx.up <- mSet_t2vatcrs[["analSet"]][["volcano"]][["inx.up"]]
t2_vat_crs_volcano$inx.down <- mSet_t2vatcrs[["analSet"]][["volcano"]][["inx.down"]]
t2_vat_crs_volcano$inx.p <- mSet_t2vatcrs[["analSet"]][["volcano"]][["inx.p"]]

t2_vat_crs_volcano <- t2_vat_crs_volcano %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

p_volcano_t2_vat_crs <- ggplot(t2_vat_crs_volcano, aes(x = log_fc, y = log_p))+
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

p_volcano_t2_vat_crs

t2_vat_crs_volcano_metabolites <- t2_vat_crs_volcano %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()

### Adding PCA analysis
mSet_t2vatcrs<-PCA.Anal(mSet_t2vatcrs)

## Extraction of PCA component values
t2_vat_crs_PCA <- as.data.frame(mSet_t2vatcrs[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_t2vatcrs order for samples (used for PLSDA)
t2_vat_crs_PCA_sample_order <- t2_vat_crs_PCA %>%
  mutate(Sample = rownames(t2_vat_crs_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
t2_vat_crs_PCA$Sample <- rownames(t2_vat_crs_PCA)
t2_vat_crs_PCA <- t2_vat_crs_PCA %>%
  select(Sample, everything())

t2_vat_crs_PCA <- left_join(t2_vat_crs_PCA, t2_vat_crs_sample_id, by = "Sample")

t2_vat_crs_PCA <- t2_vat_crs_PCA %>%
  select(Sample, VAT_CRS, everything())

t2_vat_crs_PCA$VAT_CRS <- as.character(t2_vat_crs_PCA$VAT_CRS)

## Visualization of PCA comp1 vs comp2
p_pca_t2_vat_crs <- t2_vat_crs_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = VAT_CRS)) +
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

p_pca_t2_vat_crs


## Perform PLSDA
mSet_t2vatcrs<-PLSR.Anal(mSet_t2vatcrs, reg=TRUE)

## Extraction of PLSDA component values
t2_vat_crs_PLSDA <- as.matrix.data.frame(t2_vat_crs_PLSDA <- mSet_t2vatcrs[["analSet"]][["plsr"]][["scores"]])
t2_vat_crs_PLSDA <- as.data.frame(t2_vat_crs_PLSDA)
# Sample information lost in mSet_t2vatcrs upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
t2_vat_crs_PLSDA <- t2_vat_crs_PLSDA %>%
  mutate(Sample = t2_vat_crs_PCA_sample_order)

t2_vat_crs_PLSDA <- left_join(t2_vat_crs_PLSDA, t2_vat_crs_sample_id, by = "Sample")

t2_vat_crs_PLSDA$VAT_CRS <- as.character(t2_vat_crs_PLSDA$VAT_CRS)

## Visualization of PLSDA comp1 vs comp2
p_plsda_t2_vat_crs <- t2_vat_crs_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = VAT_CRS)) +
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

p_plsda_t2_vat_crs

t2_vat_crs_PLSDA_VIP <- as.data.frame(mSet_t2vatcrs[["analSet"]][["plsr"]][["vip.mat"]])
t2_vat_crs_PLSDA_VIP <- t2_vat_crs_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

p_plasda_vip_t2_vat_crs <- t2_vat_crs_PLSDA_VIP %>%
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

p_plasda_vip_t2_vat_crs

t2_vat_crs_plsda_metabolites <- t2_vat_crs_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) %>%
  select(metabolite)|>
  unlist() %>%
  as.vector()

## Combining volcano and plsda metabolites

t2_vat_crs_metabolites <- sort(union(t2_vat_crs_plsda_metabolites, t2_vat_crs_volcano_metabolites))

###
### Identification of metabolites from VAT survival analysis ----
### 

## 1. Analysis of day 0 differences to extract metabolites

## 1.1 Loading of data for normalization and analyses
mSet_t2satcrs<-InitDataObjects("pktable", "stat", FALSE)
mSet_t2satcrs<-Read.TextData(mSet_t2satcrs, "Input_files/t2_cohort_SAT_crs.csv", "rowu", "disc");
mSet_t2satcrs<-SanityCheckData(mSet_t2satcrs)
mSet_t2satcrs<-ReplaceMin(mSet_t2satcrs);
mSet_t2satcrs<-SanityCheckData(mSet_t2satcrs)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_t2satcrs<-FilterVariable(mSet_t2satcrs, "median", 0, "F")
mSet_t2satcrs<-PreparePrenormData(mSet_t2satcrs)

## Normalization by sum and data scaling based on auto-scaling
mSet_t2satcrs<-Normalization(mSet_t2satcrs, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

# Extraction and saving of the normalized data into a new tibble
t2_sat_crs <- as.data.frame(mSet_t2satcrs[["dataSet"]][["norm"]])

# Load original table and cbind sample name and time point label
t2_sat_crs_original <- read.csv("Input_files/t2_cohort_SAT_crs.csv", na = "NA")

## Left join original data and normalized to link vat groups
t2_sat_crs$Sample <- row.names(t2_sat_crs)
t2_sat_crs <- t2_sat_crs %>% select(Sample, everything()) %>%
  arrange(Sample)

t2_sat_crs_sample_id <- t2_sat_crs_original %>% select(Sample, SAT_CRS)

t2_sat_crs_norm <- left_join(t2_sat_crs, t2_sat_crs_sample_id, by = "Sample")
t2_sat_crs_norm <- t2_sat_crs_norm %>%
  select(Sample, SAT_CRS, everything())

t2_sat_crs_norm$SAT_CRS <- as.character(t2_sat_crs_norm$SAT_CRS)

## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 0.05, and FC > 1.5
mSet_t2satcrs<-Volcano.Anal(mSet_t2satcrs, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

t2_sat_crs_volcano <- as.data.frame(mSet_t2satcrs[["analSet"]][["volcano"]][["fc.log"]])
t2_sat_crs_volcano$metabolite <- rownames(t2_sat_crs_volcano)
t2_sat_crs_volcano$log_p <- mSet_t2satcrs[["analSet"]][["volcano"]][["p.log"]]
t2_sat_crs_volcano$log_fc <- mSet_t2satcrs[["analSet"]][["volcano"]][["fc.log"]]
t2_sat_crs_volcano$inx.up <- mSet_t2satcrs[["analSet"]][["volcano"]][["inx.up"]]
t2_sat_crs_volcano$inx.down <- mSet_t2satcrs[["analSet"]][["volcano"]][["inx.down"]]
t2_sat_crs_volcano$inx.p <- mSet_t2satcrs[["analSet"]][["volcano"]][["inx.p"]]

t2_sat_crs_volcano <- t2_sat_crs_volcano %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

p_volcano_t2_sat_crs <- ggplot(t2_sat_crs_volcano, aes(x = log_fc, y = log_p))+
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

p_volcano_t2_sat_crs

t2_sat_crs_volcano_metabolites <- t2_sat_crs_volcano %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()

### Adding PCA analysis
mSet_t2satcrs<-PCA.Anal(mSet_t2satcrs)

## Extraction of PCA component values
t2_sat_crs_PCA <- as.data.frame(mSet_t2satcrs[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_t2satcrs order for samples (used for PLSDA)
t2_sat_crs_PCA_sample_order <- t2_sat_crs_PCA %>%
  mutate(Sample = rownames(t2_sat_crs_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
t2_sat_crs_PCA$Sample <- rownames(t2_sat_crs_PCA)
t2_sat_crs_PCA <- t2_sat_crs_PCA %>%
  select(Sample, everything())

t2_sat_crs_PCA <- left_join(t2_sat_crs_PCA, t2_sat_crs_sample_id, by = "Sample")

t2_sat_crs_PCA <- t2_sat_crs_PCA %>%
  select(Sample, SAT_CRS, everything())

t2_sat_crs_PCA$SAT_CRS <- as.character(t2_sat_crs_PCA$SAT_CRS)

## Visualization of PCA comp1 vs comp2
p_pca_t2_sat_crs <- t2_sat_crs_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = SAT_CRS)) +
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

p_pca_t2_sat_crs


## Perform PLSDA
mSet_t2satcrs<-PLSR.Anal(mSet_t2satcrs, reg=TRUE)

## Extraction of PLSDA component values
t2_sat_crs_PLSDA <- as.matrix.data.frame(t2_sat_crs_PLSDA <- mSet_t2satcrs[["analSet"]][["plsr"]][["scores"]])
t2_sat_crs_PLSDA <- as.data.frame(t2_sat_crs_PLSDA)
# Sample information lost in mSet_t2satcrs upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
t2_sat_crs_PLSDA <- t2_sat_crs_PLSDA %>%
  mutate(Sample = t2_sat_crs_PCA_sample_order)

t2_sat_crs_PLSDA <- left_join(t2_sat_crs_PLSDA, t2_sat_crs_sample_id, by = "Sample")

t2_sat_crs_PLSDA$SAT_CRS <- as.character(t2_sat_crs_PLSDA$SAT_CRS)

## Visualization of PLSDA comp1 vs comp2
p_plsda_t2_sat_crs <- t2_sat_crs_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = SAT_CRS)) +
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

p_plsda_t2_sat_crs

t2_sat_crs_PLSDA_VIP <- as.data.frame(mSet_t2satcrs[["analSet"]][["plsr"]][["vip.mat"]])
t2_sat_crs_PLSDA_VIP <- t2_sat_crs_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

p_plasda_vip_t2_sat_crs <- t2_sat_crs_PLSDA_VIP %>%
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

p_plasda_vip_t2_sat_crs

t2_sat_crs_plsda_metabolites <- t2_sat_crs_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) %>%
  select(metabolite)|>
  unlist() %>%
  as.vector()

## Combining volcano and plsda metabolites

t2_sat_crs_metabolites <- sort(union(t2_sat_crs_plsda_metabolites, t2_sat_crs_volcano_metabolites))
###


### Identification of metabolites from VAT survival analysis ----
### 

## 1. Analysis of day 0 differences to extract metabolites

## 1.1 Loading of data for normalization and analyses
mSet_t2tatcrs<-InitDataObjects("pktable", "stat", FALSE)
mSet_t2tatcrs<-Read.TextData(mSet_t2tatcrs, "Input_files/t2_cohort_TAT_crs.csv", "rowu", "disc");
mSet_t2tatcrs<-SanityCheckData(mSet_t2tatcrs)
mSet_t2tatcrs<-ReplaceMin(mSet_t2tatcrs);
mSet_t2tatcrs<-SanityCheckData(mSet_t2tatcrs)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_t2tatcrs<-FilterVariable(mSet_t2tatcrs, "median", 0, "F")
mSet_t2tatcrs<-PreparePrenormData(mSet_t2tatcrs)

## Normalization by sum and data scaling based on auto-scaling
mSet_t2tatcrs<-Normalization(mSet_t2tatcrs, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

# Extraction and saving of the normalized data into a new tibble
t2_tat_crs <- as.data.frame(mSet_t2tatcrs[["dataSet"]][["norm"]])

# Load original table and cbind sample name and time point label
t2_tat_crs_original <- read.csv("Input_files/t2_cohort_TAT_crs.csv", na = "NA")

## Left join original data and normalized to link vat groups
t2_tat_crs$Sample <- row.names(t2_tat_crs)
t2_tat_crs <- t2_tat_crs %>% select(Sample, everything()) %>%
  arrange(Sample)

t2_tat_crs_sample_id <- t2_tat_crs_original %>% select(Sample, TAT_CRS)

t2_tat_crs_norm <- left_join(t2_tat_crs, t2_tat_crs_sample_id, by = "Sample")
t2_tat_crs_norm <- t2_tat_crs_norm %>%
  select(Sample, TAT_CRS, everything())

t2_tat_crs_norm$TAT_CRS <- as.character(t2_tat_crs_norm$TAT_CRS)

## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 0.05, and FC > 1.5
mSet_t2tatcrs<-Volcano.Anal(mSet_t2tatcrs, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

t2_tat_crs_volcano <- as.data.frame(mSet_t2tatcrs[["analSet"]][["volcano"]][["fc.log"]])
t2_tat_crs_volcano$metabolite <- rownames(t2_tat_crs_volcano)
t2_tat_crs_volcano$log_p <- mSet_t2tatcrs[["analSet"]][["volcano"]][["p.log"]]
t2_tat_crs_volcano$log_fc <- mSet_t2tatcrs[["analSet"]][["volcano"]][["fc.log"]]
t2_tat_crs_volcano$inx.up <- mSet_t2tatcrs[["analSet"]][["volcano"]][["inx.up"]]
t2_tat_crs_volcano$inx.down <- mSet_t2tatcrs[["analSet"]][["volcano"]][["inx.down"]]
t2_tat_crs_volcano$inx.p <- mSet_t2tatcrs[["analSet"]][["volcano"]][["inx.p"]]

t2_tat_crs_volcano <- t2_tat_crs_volcano %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

p_volcano_t2_tat_crs <- ggplot(t2_tat_crs_volcano, aes(x = log_fc, y = log_p))+
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

p_volcano_t2_tat_crs

t2_tat_crs_volcano_metabolites <- t2_tat_crs_volcano %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()

### Adding PCA analysis
mSet_t2tatcrs<-PCA.Anal(mSet_t2tatcrs)

## Extraction of PCA component values
t2_tat_crs_PCA <- as.data.frame(mSet_t2tatcrs[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_t2tatcrs order for samples (used for PLSDA)
t2_tat_crs_PCA_sample_order <- t2_tat_crs_PCA %>%
  mutate(Sample = rownames(t2_tat_crs_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
t2_tat_crs_PCA$Sample <- rownames(t2_tat_crs_PCA)
t2_tat_crs_PCA <- t2_tat_crs_PCA %>%
  select(Sample, everything())

t2_tat_crs_PCA <- left_join(t2_tat_crs_PCA, t2_tat_crs_sample_id, by = "Sample")

t2_tat_crs_PCA <- t2_tat_crs_PCA %>%
  select(Sample, TAT_CRS, everything())

t2_tat_crs_PCA$TAT_CRS <- as.character(t2_tat_crs_PCA$TAT_CRS)

## Visualization of PCA comp1 vs comp2
p_pca_t2_tat_crs <- t2_tat_crs_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = TAT_CRS)) +
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

p_pca_t2_tat_crs


## Perform PLSDA
mSet_t2tatcrs<-PLSR.Anal(mSet_t2tatcrs, reg=TRUE)

## Extraction of PLSDA component values
t2_tat_crs_PLSDA <- as.matrix.data.frame(t2_tat_crs_PLSDA <- mSet_t2tatcrs[["analSet"]][["plsr"]][["scores"]])
t2_tat_crs_PLSDA <- as.data.frame(t2_tat_crs_PLSDA)
# Sample information lost in mSet_t2tatcrs upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
t2_tat_crs_PLSDA <- t2_tat_crs_PLSDA %>%
  mutate(Sample = t2_tat_crs_PCA_sample_order)

t2_tat_crs_PLSDA <- left_join(t2_tat_crs_PLSDA, t2_tat_crs_sample_id, by = "Sample")

t2_tat_crs_PLSDA$TAT_CRS <- as.character(t2_tat_crs_PLSDA$TAT_CRS)

## Visualization of PLSDA comp1 vs comp2
p_plsda_t2_tat_crs <- t2_tat_crs_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = TAT_CRS)) +
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

p_plsda_t2_tat_crs

t2_tat_crs_PLSDA_VIP <- as.data.frame(mSet_t2tatcrs[["analSet"]][["plsr"]][["vip.mat"]])
t2_tat_crs_PLSDA_VIP <- t2_tat_crs_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

p_plasda_vip_t2_tat_crs <- t2_tat_crs_PLSDA_VIP %>%
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

p_plasda_vip_t2_tat_crs

t2_tat_crs_plsda_metabolites <- t2_tat_crs_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) %>%
  select(metabolite)|>
  unlist() %>%
  as.vector()

## Combining volcano and plsda metabolites

t2_tat_crs_metabolites <- sort(union(t2_tat_crs_plsda_metabolites, t2_tat_crs_volcano_metabolites))
###


### Comparison of filtered survival metabolites with BC levels ----
###

## Comparison of VAT metabolites with VAT measurements
t2_vat_crs_norm_bc <- left_join(t2_vat_crs_norm, all_master, by="Sample")
str(t2_vat_crs_norm_bc)

# Create an empty data frame to store correlation results
t2_vat_crs_corr <- data.frame(
  Metabolite = character(),
  Pearson_r  = numeric(),
  p_value    = numeric(),
  stringsAsFactors = FALSE
)

# Create a list to hold all ggplots
plot_list_t2_vat_crs_corr <- list()

# Loop over each metabolite column
for (met in t2_vat_crs_metabolites) {
  # Perform a Pearson correlation test against VAT
  corr_test <- cor.test(
    t2_vat_crs_norm_bc[[met]],
    t2_vat_crs_norm_bc[["VAT"]],
    method = "graderson"
  )
  
  # Extract the correlation coefficient (r) and p-value
  r_val <- corr_test$estimate
  p_val <- corr_test$p.value
  
  # Append the numeric results to 'results' data frame
  t2_vat_crs_corr <- rbind(
    t2_vat_crs_corr,
    data.frame(
      Metabolite = met,
      Pearson_r  = r_val,
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
  )
  
  # Create a scatter plot with regression line
  # Use aes(x = VAT, y = .data[[met]]) instead of aes_string()
  p <- ggplot(t2_vat_crs_norm_bc, aes(x = VAT, y = .data[[met]])) +
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
  plot_list_t2_vat_crs_corr[[met]] <- p
}

# If you want to see all the plots at once, you can do:
# Create the arranged grid
grid_combined <- arrangeGrob(grobs = plot_list_t2_vat_crs_corr, ncol = 4)

# Save as SVG
ggsave("Figures_Manuscript/met_corr_VAT_crs.svg", plot = grid_combined, width = 20, height = 36)

## SAT correlations
# Create an empty data frame to store correlation results
t2_sat_crs_corr <- data.frame(
  Metabolite = character(),
  Pearson_r  = numeric(),
  p_value    = numeric(),
  stringsAsFactors = FALSE
)

# Create a list to hold all ggplots
plot_list_t2_sat_crs_corr <- list()

# Loop over each metabolite column
for (met in t2_sat_crs_metabolites) {
  # Perform a Pearson correlation test against VAT
  corr_test <- cor.test(
    t2_vat_crs_norm_bc[[met]],
    t2_vat_crs_norm_bc[["VAT"]],
    method = "graderson"
  )
  
  # Extract the correlation coefficient (r) and p-value
  r_val <- corr_test$estimate
  p_val <- corr_test$p.value
  
  # Append the numeric results to 'results' data frame
  t2_sat_crs_corr <- rbind(
    t2_sat_crs_corr,
    data.frame(
      Metabolite = met,
      Pearson_r  = r_val,
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
  )
  
  # Create a scatter plot with regression line
  # Use aes(x = VAT, y = .data[[met]]) instead of aes_string()
  p <- ggplot(t2_vat_crs_norm_bc, aes(x = VAT, y = .data[[met]])) +
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
  plot_list_t2_sat_crs_corr[[met]] <- p
}

# If you want to see all the plots at once, you can do:
# Create the arranged grid
grid_combined <- arrangeGrob(grobs = plot_list_t2_sat_crs_corr, ncol = 4)

# Save as SVG
ggsave("Figures_Manuscript/met_corr_SAT_crs.svg", plot = grid_combined, width = 20, height = 36)


## TAT correlations
# Create an empty data frame to store correlation results
t2_tat_crs_corr <- data.frame(
  Metabolite = character(),
  Pearson_r  = numeric(),
  p_value    = numeric(),
  stringsAsFactors = FALSE
)

# Create a list to hold all ggplots
plot_list_t2_tat_crs_corr <- list()

# Loop over each metabolite column
for (met in t2_tat_crs_metabolites) {
  # Perform a Pearson correlation test against VAT
  corr_test <- cor.test(
    t2_vat_crs_norm_bc[[met]],
    t2_vat_crs_norm_bc[["VAT"]],
    method = "graderson"
  )
  
  # Extract the correlation coefficient (r) and p-value
  r_val <- corr_test$estimate
  p_val <- corr_test$p.value
  
  # Append the numeric results to 'results' data frame
  t2_tat_crs_corr <- rbind(
    t2_tat_crs_corr,
    data.frame(
      Metabolite = met,
      Pearson_r  = r_val,
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
  )
  
  # Create a scatter plot with regression line
  # Use aes(x = VAT, y = .data[[met]]) instead of aes_string()
  p <- ggplot(t2_vat_crs_norm_bc, aes(x = VAT, y = .data[[met]])) +
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
  plot_list_t2_tat_crs_corr[[met]] <- p
}

# If you want to see all the plots at once, you can do:
# Create the arranged grid
grid_combined <- arrangeGrob(grobs = plot_list_t2_tat_crs_corr, ncol = 4)

# Save as SVG
ggsave("Figures_Manuscript/met_corr_TAT_crs.svg", plot = grid_combined, width = 20, height = 36)


##
## Filtering for metabolites which are associated with BCs with a Pearson of at least 0.2 or -0.2

t2_vat_crs_corr_filter <- t2_vat_crs_corr |>
  filter(Pearson_r <= -0.2 | Pearson_r >= 0.2)|>
  pull(Metabolite) |>
  unname()
t2_sat_crs_corr_filter <- t2_sat_crs_corr |>
  filter(Pearson_r <= -0.2 | Pearson_r >= 0.2)|>
  pull(Metabolite) |>
  unname()
t2_tat_crs_corr_filter <- t2_tat_crs_corr |>
  filter(Pearson_r <= -0.2 | Pearson_r >= 0.2)|>
  pull(Metabolite) |>
  unname()

### Common and different metabolites to the VAT adipose depot ----
###

# Common to all three
crs_common_all <- sort(Reduce(intersect, list(
  t2_vat_crs_corr_filter,
  t2_sat_crs_corr_filter,
  t2_tat_crs_corr_filter
)))

# Unique to each
crs_only_in_tat <- sort(setdiff(t2_tat_crs_corr_filter, union(t2_vat_crs_corr_filter, t2_sat_crs_corr_filter)))
crs_only_in_vat <- sort(setdiff(t2_vat_crs_corr_filter, union(t2_tat_crs_corr_filter, t2_sat_crs_corr_filter)))
crs_only_in_sat <- sort(setdiff(t2_sat_crs_corr_filter, union(t2_tat_crs_corr_filter, t2_vat_crs_corr_filter)))

# Shared between any two but not all three
crs_tat_vat_shared <- sort(setdiff(intersect(t2_tat_crs_corr_filter, t2_vat_crs_corr_filter), crs_common_all))
crs_tat_sat_shared <- sort(setdiff(intersect(t2_tat_crs_corr_filter, t2_sat_crs_corr_filter), crs_common_all))
crs_vat_sat_shared <- sort(setdiff(intersect(t2_vat_crs_corr_filter, t2_sat_crs_corr_filter), crs_common_all))

# Organize output
list(
  common_to_all_three = crs_common_all,
  only_in_tat = crs_only_in_tat,
  only_in_vat = crs_only_in_vat,
  only_in_sat = crs_only_in_sat,
  shared_tat_vat = crs_tat_vat_shared,
  shared_tat_sat = crs_tat_sat_shared,
  shared_vat_sat = crs_vat_sat_shared
)

t2_all_crs_corr_filter <- sort(unique(c(
  t2_vat_crs_corr_filter,
  t2_sat_crs_corr_filter,
  t2_tat_crs_corr_filter
)))

###
### Associations with survival ----
###


## Performing a multivariable log regression analysis for CRS grade >= 2
str(t2_vat_crs_norm_bc)
t2_vat_crs_norm_bc$CRS_high

t2_all_glm_adjusted <- data.frame(marker = character(),
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
  model <- glm(formula_str, data = t2_vat_crs_norm_bc, family = binomial)
  
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
  t2_all_glm_adjusted <- rbind(t2_all_glm_adjusted, marker_results)
}

t2_all_glm_adjusted <- t2_all_glm_adjusted|>
  mutate(OR = exp(coefficient), FDR = p.adjust(p_value, method = "BH"))

t2_all_glm_adjusted$group <- sapply(as.character(t2_all_glm_adjusted$marker), get_group_new)

# Double-check that meta_df has the columns you expect:
t2_all_glm_adjusted <- t2_all_glm_adjusted |> filter(!(marker %in% c("PI-(38:07)", "PI-(40:03)","PI-(40:08)","PI-(40:09)")))

t2_all_glm_adjusted <- t2_all_glm_adjusted |>
  arrange(group, desc(OR)) |>  # You can change `desc(HR)` to another column if needed
  mutate(marker = factor(marker, levels = unique(marker)))


p_crs_group <- t2_all_glm_adjusted |>
  #filter(FDR <= 0.1)|>
  ggplot()+
  geom_point(aes(x = marker, y = OR, color = group),
             size = 3, shape = 19, alpha = 0.6)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = upper95))+
  coord_flip()+ 
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  labs(color = "Group")+
  scale_x_discrete(limits = rev)+
  labs(color = "Group")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

print(p_crs_group)
ggsave("Figures_Manuscript/meta_crs_group.svg", plot = p_crs_group, width =9, height = 9)


## Meta-analysis of group effects

# Now create the columns needed for meta-analysis
meta_crs <- t2_all_glm_adjusted %>%
  mutate(
    logOR    = log(OR),
    logLower = log(lower95),
    logUpper = log(upper95),
    # approximate standard error for log(HR)
    SE_logOR = (logUpper - logLower) / (2 * 1.96)
  )

## !! Multiple values are the same for different metabolites - need to be removed
str(meta_crs)


meta_results_crs <- meta_crs %>%
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

meta_crs_significant <- c("Acylcarnitine", "Phosphatidylethanolamine", 
                          "Carbohydrate / Sugar", "Lysophosphatidylcholine", 
                          "Lysophospholipid / PAF" )

p_meta_crs <- meta_results_crs |>
  mutate(col_group = if_else(group %in% meta_crs_significant, "highlight", "other")) %>% 
  filter(group != "Other/Unclassified")|>
  ggplot()+
  geom_linerange(aes(x = reorder(group, combined_OR), ymin = lower95, 
                     ymax = upper95), size=0.8, alpha = 0.8)+
  geom_point(aes(x = reorder(group, combined_OR), y = combined_OR, size = n_mets, color=col_group),
             shape = 19, alpha = 0.8)+
  coord_flip()+ 
  labs(y = "Combined OR (95%CI)", x = "")+
  scale_color_manual(values = c("highlight" = "#CC0000", "other" = "black"), guide = "none") +
  scale_size_continuous(breaks = c(1,5,10), range = c(2,7))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  guides(color = "none")+
  labs(size="Metabolites")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )
  
print(p_meta_crs)
ggsave("Figures_Manuscript/meta_crs.svg", plot = p_meta_crs, width = 6, height = 4.5)


###
### Test distribution analysis
###

group_colors <- c(
  only_vat = "#4CAF50",
  only_sat = "#9C27B0",
  only_tat = "#2196F3",
  tat_vat = "#009688",
  tat_sat = "#673AB7",
  vat_sat = "#8BC34A",
  common_all = "#9E9E9E"
)


# Distribution of ACs
get_ACs <- function(x) {
  grep("^AC-", x, value = TRUE)
}

crs_ac_common_all <- get_ACs(crs_common_all)
crs_ac_only_tat <- get_ACs(crs_only_in_tat)
crs_ac_only_vat <- get_ACs(crs_only_in_vat)
crs_ac_only_sat <- get_ACs(crs_only_in_sat)
crs_ac_tat_vat <- get_ACs(crs_tat_vat_shared)
crs_ac_tat_sat <- get_ACs(crs_tat_sat_shared)
crs_ac_vat_sat <- get_ACs(crs_vat_sat_shared)

crs_ac_df <- data.frame(
  group = c(
    rep("common_all", length(crs_ac_common_all)),
    rep("only_tat", length(crs_ac_only_tat)),
    rep("only_vat", length(crs_ac_only_vat)),
    rep("only_sat", length(crs_ac_only_sat)),
    rep("tat_vat", length(crs_ac_tat_vat)),
    rep("tat_sat", length(crs_ac_tat_sat)),
    rep("vat_sat", length(crs_ac_vat_sat))
  ),
  metabolite = c(
    crs_ac_common_all,
    crs_ac_only_tat,
    crs_ac_only_vat,
    crs_ac_only_sat,
    crs_ac_tat_vat,
    crs_ac_tat_sat,
    crs_ac_vat_sat
  )
)

crs_ac_summary <- crs_ac_df %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count))|>
  mutate(met_group = "AC")

crs_ac_vector <- crs_ac_df %>%
  select(metabolite)|>
  unlist()|>
  as.vector()

ggplot(crs_ac_summary, aes(x = "", y = count, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Acetylcarnitines (ACs) by Adipose Tissue Site") +
  scale_fill_manual(values = group_colors) +
  scale_pattern_manual(values = c(
    only_vat = "none",
    only_sat = "none",
    only_tat = "none",
    tat_vat = "stripe",
    tat_sat = "crosshatch",
    vat_sat = "stripe",
    common_all = "circle"
  ))

# Distribution of PEAs

crs_grade_common_all <- get_PEAs(crs_common_all)
crs_pea_only_tat <- get_PEAs(crs_only_in_tat)
crs_pea_only_vat <- get_PEAs(crs_only_in_vat)
crs_pea_only_sat <- get_PEAs(crs_only_in_sat)
crs_pea_tat_vat <- get_PEAs(crs_tat_vat_shared)
crs_pea_tat_sat <- get_PEAs(crs_tat_sat_shared)
crs_pea_vat_sat <- get_PEAs(crs_vat_sat_shared)

crs_pea_df <- data.frame(
  group = c(
    rep("common_all", length(crs_pea_common_all)),
    rep("only_tat", length(crs_pea_only_tat)),
    rep("only_vat", length(crs_pea_only_vat)),
    rep("only_sat", length(crs_pea_only_sat)),
    rep("tat_vat", length(crs_pea_tat_vat)),
    rep("tat_sat", length(crs_pea_tat_sat)),
    rep("vat_sat", length(crs_pea_vat_sat))
  ),
  metabolite = c(
    crs_pea_common_all,
    crs_pea_only_tat,
    crs_pea_only_vat,
    crs_pea_only_sat,
    crs_pea_tat_vat,
    crs_pea_tat_sat,
    crs_pea_vat_sat
  )
)

crs_pea_summary <- crs_pea_df %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count))|>
  mutate(met_group = "PEA")

crs_pea_vector <- crs_pea_df %>%
  select(metabolite)|>
  unlist()|>
  as.vector()

ggplot(crs_pea_summary, aes(x = "", y = count, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Phosphatidylethanolamine (PEAs) by Adipose Tissue Site") +
  scale_fill_manual(values = group_colors) +
  scale_pattern_manual(values = c(
    only_vat = "none",
    only_sat = "none",
    only_tat = "none",
    tat_vat = "stripe",
    tat_sat = "crosshatch",
    vat_sat = "stripe",
    common_all = "circle"
  ))


# Distribution of LPCs
get_LPCs <- function(x) {
  grep("^LPC-", x, value = TRUE)
}

crs_lpc_common_all <- get_LPCs(crs_common_all)
crs_lpc_only_tat <- get_LPCs(crs_only_in_tat)
crs_lpc_only_vat <- get_LPCs(crs_only_in_vat)
crs_lpc_only_sat <- get_LPCs(crs_only_in_sat)
crs_lpc_tat_vat <- get_LPCs(crs_tat_vat_shared)
crs_lpc_tat_sat <- get_LPCs(crs_tat_sat_shared)
crs_lpc_vat_sat <- get_LPCs(crs_vat_sat_shared)

crs_lpc_df <- data.frame(
  group = c(
    rep("common_all", length(crs_lpc_common_all)),
    rep("only_tat", length(crs_lpc_only_tat)),
    rep("only_vat", length(crs_lpc_only_vat)),
    rep("only_sat", length(crs_lpc_only_sat)),
    rep("tat_vat", length(crs_lpc_tat_vat)),
    rep("tat_sat", length(crs_lpc_tat_sat)),
    rep("vat_sat", length(crs_lpc_vat_sat))
  ),
  metabolite = c(
    crs_lpc_common_all,
    crs_lpc_only_tat,
    crs_lpc_only_vat,
    crs_lpc_only_sat,
    crs_lpc_tat_vat,
    crs_lpc_tat_sat,
    crs_lpc_vat_sat
  )
)

crs_lpc_summary <- crs_lpc_df %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count))|>
  mutate(met_group = "LPC")

crs_lpc_vector <- crs_lpc_df|>
  select(metabolite)|>
  unlist()|>
  as.vector()

ggplot(crs_lpc_summary, aes(x = "", y = count, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Lysophosphatidylcholine (LPCs) by Adipose Tissue Site") +
  scale_fill_manual(values = group_colors) +
  scale_pattern_manual(values = c(
    only_vat = "none",
    only_sat = "none",
    only_tat = "none",
    tat_vat = "stripe",
    tat_sat = "crosshatch",
    vat_sat = "stripe",
    common_all = "circle"
  ))


# Distribution of LPCs
get_lysoPAF <- function(x) {
  grep("^lysoPAF-", x, value = TRUE)
}

crs_lysoPAF_common_all <- get_lysoPAF(crs_common_all)
crs_lysoPAF_only_tat <- get_lysoPAF(crs_only_in_tat)
crs_lysoPAF_only_vat <- get_lysoPAF(crs_only_in_vat)
crs_lysoPAF_only_sat <- get_lysoPAF(crs_only_in_sat)
crs_lysoPAF_tat_vat <- get_lysoPAF(crs_tat_vat_shared)
crs_lysoPAF_tat_sat <- get_lysoPAF(crs_tat_sat_shared)
crs_lysoPAF_vat_sat <- get_lysoPAF(crs_vat_sat_shared)

crs_lysoPAF_df <- data.frame(
  group = c(
    rep("common_all", length(crs_lysoPAF_common_all)),
    rep("only_tat", length(crs_lysoPAF_only_tat)),
    rep("only_vat", length(crs_lysoPAF_only_vat)),
    rep("only_sat", length(crs_lysoPAF_only_sat)),
    rep("tat_vat", length(crs_lysoPAF_tat_vat)),
    rep("tat_sat", length(crs_lysoPAF_tat_sat)),
    rep("vat_sat", length(crs_lysoPAF_vat_sat))
  ),
  metabolite = c(
    crs_lysoPAF_common_all,
    crs_lysoPAF_only_tat,
    crs_lysoPAF_only_vat,
    crs_lysoPAF_only_sat,
    crs_lysoPAF_tat_vat,
    crs_lysoPAF_tat_sat,
    crs_lysoPAF_vat_sat
  )
)

crs_lysoPAF_summary <- crs_lysoPAF_df %>%
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) |>
  mutate(met_group = "lysoPAF")

crs_lysoPAF_vector <- crs_lysoPAF_df|>
  select(metabolite)|>
  unlist()|>
  as.vector()

ggplot(crs_lysoPAF_summary, aes(x = "", y = count, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Lysophosphatidylcholine (LPCs) by Adipose Tissue Site") +
  scale_fill_manual(values = group_colors) +
  scale_pattern_manual(values = c(
    only_vat = "none",
    only_sat = "none",
    only_tat = "none",
    tat_vat = "stripe",
    tat_sat = "crosshatch",
    vat_sat = "stripe",
    common_all = "circle"
  ))

## Combining data frames and create a bar plot
crs_at_distr <- rbind(crs_ac_summary, crs_pea_summary, crs_lpc_summary, crs_lysoPAF_summary) %>%
  group_by(met_group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(met_group = factor(met_group, levels = c("AC", "PEA", "LPC", "lysoPAF")))

# New grouping variable
crs_at_distr<- crs_at_distr |>
  mutate(group_new = case_when(
    group == "only_tat" ~ "AT-shared",
    group == "only_vat" ~ "VAT-derived",
    group == "only_sat" ~ "SAT-derived",
    group == "tat_vat" ~ "VAT-enriched",
    group == "tat_sat" ~ "SAT-enriched",
    group == "vat_sat" ~ "AT-shared",
    group == "common_all" ~ "AT-shared",
  ))


# Plot from the correctly prepared data
ggplot(crs_at_distr, aes(x = met_group, y = percentage, fill = group)) +
  geom_bar(stat = "identity") +
  labs(x = "Metabolite Group", y = "Percentage", fill = "Group") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

crs_at_distr <- crs_at_distr |>
  mutate(group_new = factor(group_new, levels = c("AT-shared" ,
                                                  "VAT-derived",
                                                  "VAT-enriched",
                                                  "SAT-derived",
                                                  "SAT-enriched" )))

p_crs_distr <- ggplot(crs_at_distr, aes(x = met_group, y = percentage, fill = group_new)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Metabolites [%]", fill = "Source") +
  scale_fill_manual(values = group_colors_AT) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


print(p_crs_distr)
ggsave("Figures_Manuscript/p_crs_distr.svg", plot = p_crs_distr, width = 5, height = 3)

### Analysis of CRS distribution and CRS development, based on metabolite examples

#LPC
p_crs_lpc <- t2_vat_crs_norm_bc %>%
  mutate(LPC_group = ifelse(`LPC-(20:01)` >= median(`LPC-(20:01)`), "LPC_high", "LPC_low")) %>%
  ggplot(aes(x = LPC_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  scale_x_discrete(limits = c("LPC_low", "LPC_high"), labels = c("Low", "High"))+
  scale_fill_manual(values=c("#FFF5E1","#FFDAB9", "#FFC04C" ,"#FF8C00" ))+
  scale_y_continuous(labels = scales::percent_format(scale = 100),expand = c(0,0))+
  labs(x="LPC", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

ggsave(filename = "Figures_Manuscript/p_crs_lpc.svg", plot = p_crs_lpc,
       width = 3, height = 2.5)

t2_vat_crs_norm_bc %>%
  mutate(LPC_group = ifelse(`LPC-(20:01)` >= median(`LPC-(20:01)`), "LPC_high", "LPC_low")) %>%
  with(table(LPC_group, Maximaler.CRS.Grad)) |>
  fisher.test(.)

ggplot(t2_vat_crs_norm_bc, aes(x=`LPC-(20:01)`, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

t2_vat_crs_norm_bc$crs
t2_vat_crs_norm_bc
crs_ci_lpc <- survfit2(Surv(CRS_onset_high, as.numeric(`CRS_high`)) ~ ifelse(`LPC-(20:01)` >= median(`LPC-(20:01)`), "LPC_high", "LPC_low"), data = t2_vat_crs_norm_bc)
p_crs_ci_lpc <- ggsurvplot(crs_ci_lpc,
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
            legend.title = c("LPC-(20:01)"),
            palette = c("#FF8C00", "darkgrey")
           )

ggsave(filename = "Figures_Manuscript/p_crs_ci_lpc.svg", plot = p_crs_ci_lpc$plot,
       width = 3, height = 2.5)

#lysoPAF
p_crs_lysoPAF <- t2_vat_crs_norm_bc %>%
  mutate(lysoPAF_group = ifelse(`lysoPAF-(18:1)` >= median(`lysoPAF-(18:1)`), "lysoPAF_high", "lysoPAF_low")) %>%
  ggplot(aes(x = lysoPAF_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  scale_x_discrete(limits = c("lysoPAF_low", "lysoPAF_high"), labels = c("Low", "High"))+
  scale_fill_manual(values=c("#FFF5E1","#FFDAB9", "#FFC04C" ,"#FF8C00" ))+
  scale_y_continuous(labels = scales::percent_format(scale = 100),expand = c(0,0))+
  labs(x="Lyso PAF", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

ggsave(filename = "Figures_Manuscript/p_crs_lysoPAF.svg", plot = p_crs_lysoPAF,
       width = 3, height = 2.5)

t2_vat_crs_norm_bc %>%
  mutate(lysoPAF_group = ifelse(`lysoPAF-(18:1)` >= median(`lysoPAF-(18:1)`), "lysoPAF_high", "lysoPAF_low")) %>%
  with(table(lysoPAF_group, Maximaler.CRS.Grad)) |>
  fisher.test(.)

ggplot(t2_vat_crs_norm_bc, aes(x=`lysoPAF-(18:1)`, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

ggplot(t2_vat_crs_norm_bc, aes(x=`lysoPAF-(18:1)`, y=`LPC-(20:01)`))+
  geom_point()+
  geom_smooth(method = "lm")

cor.test(t2_vat_crs_norm_bc$`lysoPAF-(18:1)`, t2_vat_crs_norm_bc$`LPC-(20:01)`, method = "pearson")

crs_ci_lysoPAF <- survfit2(Surv(CRS_onset_high, as.numeric(`CRS_high`)) ~ ifelse(`lysoPAF-(18:1)` >= median(`lysoPAF-(18:1)`), "lysoPAF_high", "lysoPAF_low"), data = t2_vat_crs_norm_bc)
p_crs_ci_lysoPAF <- ggsurvplot(crs_ci_lysoPAF,
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
           legend.title = c("lysoPAF-(18:1)"),
           palette = c("#FF8C00", "darkgrey")
)
ggsave(filename = "Figures_Manuscript/p_crs_ci_lysoPAF.svg", plot = p_crs_ci_lysoPAF$plot,
       width = 3, height = 2.5)
#AC
p_crs_ac <- t2_vat_crs_norm_bc %>%
  mutate(AC_group = ifelse(`AC-(10:0)` >= median(`AC-(10:0)`), "AC_high", "AC_low")) %>%
  ggplot(aes(x = AC_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  scale_x_discrete(limits = c("AC_low", "AC_high"), labels = c("Low", "High"))+
  scale_fill_manual(values=c("#FFF5E1","#FFDAB9", "#FFC04C" ,"#FF8C00" ))+
  scale_y_continuous(labels = scales::percent_format(scale = 100),expand = c(0,0))+
  labs(x="AC", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

ggsave(filename = "Figures_Manuscript/p_crs_ac.svg", plot = p_crs_ac,
       width = 3, height = 2.5)

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


cor.test(t2_vat_crs_norm_bc$`AC-(10:0)`, t2_vat_crs_norm_bc$STLV, method = "pearson")
cor.test(t2_vat_crs_norm_bc$`AC-(10:0)`, t2_vat_crs_norm_bc$VAT, method = "pearson")

 crs_ci_ac <- survfit2(Surv(CRS_onset_high, as.numeric(`CRS_high`)) ~ ifelse(`AC-(10:0)` >= median(`AC-(10:0)`), "AC_high", "AC_low"), data = t2_vat_crs_norm_bc)
 p_crs_ci_ac <- ggsurvplot(crs_ci_ac,
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
           legend.title = c("AC-(10:0)"),
           palette = c("#FF8C00", "darkgrey")
)
ggsave(filename = "Figures_Manuscript/p_crs_ci_ac.svg", plot = p_crs_ci_ac$plot,
       width = 3, height = 2.5)

#PEA
p_crs_pea <- t2_vat_crs_norm_bc %>%
  mutate(PEA_group = ifelse(`PEA-(38:06)` >= median(`PEA-(38:06)`), "PEA_high", "PEA_low")) %>%
  ggplot(aes(x = PEA_group, fill = factor(Maximaler.CRS.Grad))) +
  geom_bar(position = "fill") +
  scale_x_discrete(limits = c("PEA_low", "PEA_high"), labels = c("Low", "High"))+
  scale_fill_manual(values=c("#FFF5E1","#FFDAB9", "#FFC04C" ,"#FF8C00" ))+
  scale_y_continuous(labels = scales::percent_format(scale = 100),expand = c(0,0))+
  labs(x="PEA", y="Number of patients [%]", fill="CRS grade")+
  theme_classic()

ggsave(filename = "Figures_Manuscript/p_crs_pea.svg", plot = p_crs_pea,
       width = 3, height = 2.5)


t2_vat_crs_norm_bc %>%
  mutate(PEA_group = ifelse(`PEA-(38:06)` >= median(`PEA-(38:06)`), "PEA_high", "PEA_low")) %>%
  with(table(PEA_group, Maximaler.CRS.Grad)) |>
  fisher.test(.)

ggplot(t2_vat_crs_norm_bc, aes(x=`PEA-(38:06)`, y=STLV))+
  geom_point()+
  geom_smooth(method = "lm")

ggplot(t2_vat_crs_norm_bc, aes(x=`PEA-(38:06)`, y=VAT))+
  geom_point()+
  geom_smooth(method = "lm")


cor.test(t2_vat_crs_norm_bc$`PEA-(38:06)`, t2_vat_crs_norm_bc$STLV, method = "pearson")
cor.test(t2_vat_crs_norm_bc$`PEA-(38:06)`, t2_vat_crs_norm_bc$VAT, method = "pearson")

cor.test(t2_vat_crs_norm_bc$`PEA-(38:06)`, t2_vat_crs_norm_bc$`AC-(10:0)`, method = "pearson")


cor.test(t2_vat_crs_norm_bc$STLV, t2_vat_crs_norm_bc$SAT, method = "pearson")

ggplot(t2_vat_crs_norm_bc, aes(x=`PEA-(38:06)`, y=`AC-(10:0)`))+
  geom_point()+
  geom_smooth(method = "lm")

crs_ci_pea <- survfit2(Surv(CRS_onset_high, as.numeric(`CRS_high`)) ~ ifelse(`PEA-(38:06)` >= median(`PEA-(38:06)`), "PEA_high", "PEA_low"), data = t2_vat_crs_norm_bc)
p_crs_ci_pea <- ggsurvplot(crs_ci_pea,
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
           legend.title = c("PEA-(38:06)"),
           palette = c("#FF8C00", "darkgrey")
)
ggsave(filename = "Figures_Manuscript/p_crs_ci_pea.svg", plot = p_crs_ci_pea$plot,
       width = 3, height = 2.5)



## Loading Werner HMDB dictionary
hmdb_dict <- read.xlsx("Input_files/Werner_HMDB_translation.xlsx")


###
### Clinical analysis for CRS
###

p_crs_training_pie <- all_master |> 
  filter(cohort == "training") |> 
  ggplot(aes(x = "", fill = factor(Maximaler.CRS.Grad, levels = c(0, 1, 2, 3)))) +
  geom_bar(width = 1) +
  coord_polar(theta = "y", start = 0, direction = -1) +  # start at 12 o'clock and go clockwise
  scale_fill_manual(values = c("#FFF5E1", "#FFDAB9", "#FFC04C", "#FF8C00"), name = "CRS grade") +
  theme_void()  # clean theme for pie chart
  #theme(legend.title = element_blank())

ggsave("Figures_Manuscript/crs_training_pie.svg", plot = p_crs_training_pie, width = 3, height = 3)

all_master |> 
  filter(cohort == "training") |> 
  group_by(Maximaler.CRS.Grad) |>
  count() |>
  ungroup()|>
  mutate(perc = n/sum(n))


p_crs_training_vat <- all_master |> 
  filter(cohort == "training") |> 
  ggplot(aes(x = factor(Maximaler.CRS.Grad, levels = c(0, 1, 2, 3)), y = VAT)) +
  geom_jitter(aes(fill = factor(Maximaler.CRS.Grad), color = factor(Maximaler.CRS.Grad)), width = 0.2, size = 3, alpha = 0.9) +
  geom_boxplot(aes(fill =  factor(Maximaler.CRS.Grad)), outlier.shape = NA, alpha = 0.7) +
  geom_pwc(ref.group = "0", method = "wilcox.test", label.size = 4, bracket.nudge.y = -0.06, step.increase = 0.12
    #tip.length = 0.01,
    ) +
  scale_fill_manual(values = c("#FFF5E1", "#FFDAB9", "#FFC04C", "#FF8C00")) +
  scale_color_manual(values = c("#FFF5E1", "#FFDAB9", "#FFC04C", "#FF8C00")) +
  labs(x = "CRS grade", y = "VAT [cm²]") +
  guides(color = "none", fill = "none")+
  ylim(0,400)+
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )
p_crs_training_vat
ggsave("Figures_Manuscript/crs_training_vat.svg", plot = p_crs_training_vat, width = 3, height = 3)


crs_ci_vat <- survfit2(Surv(CRS_onset_high, as.numeric(`CRS_high`)) ~ VAT_CRS.y, data = t2_vat_crs_norm_bc)
p_crs_ci_vat <- ggsurvplot(crs_ci_vat,
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
p_crs_ci_vat
ggsave(filename = "Figures_Manuscript/p_crs_ci_vat.svg", plot = p_crs_ci_vat$plot,
       width = 3, height = 2.5)


# all_master |> 
#   filter(cohort == "training") |> 
#   ggplot(aes(x = factor(Maximaler.CRS.Grad, levels = c(0, 1, 2, 3)), y = SAT)) +
#   geom_jitter(aes(fill = factor(Maximaler.CRS.Grad), color = factor(Maximaler.CRS.Grad)), width = 0.2, size = 3, alpha = 0.9) +
#   geom_boxplot(aes(fill =  factor(Maximaler.CRS.Grad)), outlier.shape = NA, alpha = 0.7) +
#   geom_pwc(ref.group = "0", method = "wilcox.test",
#            #tip.length = 0.01,
#            label.size = 4) +
#   scale_fill_manual(values = c("#FFF5E1", "#FFDAB9", "#FFC04C", "#FF8C00")) +
#   scale_color_manual(values = c("#FFF5E1", "#FFDAB9", "#FFC04C", "#FF8C00")) +
#   labs(x = "CRS grade", y = "VAT [cm²]") +
#   guides(color = "none", fill = "none")+
#   theme_classic() +
#   theme(
#     axis.text = element_text(size = 9),
#     axis.title = element_text(size = 9),
#     legend.title = element_text(size = 9),
#     legend.text = element_text(size = 9)
#   )



