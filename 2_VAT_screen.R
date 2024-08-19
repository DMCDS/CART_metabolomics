## New Script for final vat analysis
## 1. Analysis of day 0 differences to extract metabolites ----

## 1.1 Loading of data for normalization and analyses
mSet_vat0<-InitDataObjects("pktable", "stat", FALSE)
mSet_vat0<-Read.TextData(mSet_vat0, "Input_files/CART_VAT_Tag0.csv", "rowu", "disc");
mSet_vat0<-SanityCheckData(mSet_vat0)
mSet_vat0<-ReplaceMin(mSet_vat0);
mSet_vat0<-SanityCheckData(mSet_vat0)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_vat0<-FilterVariable(mSet_vat0, "median", 0, "F")
mSet_vat0<-PreparePrenormData(mSet_vat0)

## Normalization by sum and data scaling based on auto-scaling
mSet_vat0<-Normalization(mSet_vat0, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)


## Extraction and saving of the normalized data into a new tibble
cart_vat_0 <- as.data.frame(mSet_vat0[["dataSet"]][["norm"]])

## Load original table and cbind sample name and time point label
cart_vat_0_original <- read_xlsx("Input_files/CART_VAT_Tag0.xlsx", na = "NA")

## Left join original data and normalized to link vat groups
cart_vat_0$Sample <- row.names(cart_vat_0)
cart_vat_0 <- cart_vat_0 %>% select(Sample, everything()) %>%
  arrange(Sample)

Sample_vat_0 <- cart_vat_0_original %>% select(Sample, VAT_high)

cart_vat_0_norm <- left_join(cart_vat_0, Sample_vat_0, by = "Sample")
cart_vat_0_norm <- cart_vat_0_norm %>%
  select(Sample, VAT_high, everything())

cart_vat_0_norm$VAT_high <- as.character(cart_vat_0_norm$VAT_high)

## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 0.05, and FC > 1.5
mSet_vat0<-Volcano.Anal(mSet_vat0, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

cart_volcano_vat_0_raw <- as.data.frame(mSet_vat0[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_vat_0_raw$metabolite <- rownames(cart_volcano_vat_0_raw)
cart_volcano_vat_0_raw$log_p <- mSet_vat0[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_vat_0_raw$log_fc <- mSet_vat0[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_vat_0_raw$inx.up <- mSet_vat0[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_vat_0_raw$inx.down <- mSet_vat0[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_vat_0_raw$inx.p <- mSet_vat0[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_vat_0_raw <- cart_volcano_vat_0_raw %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

plot_volc_vat_0_raw <- ggplot(cart_volcano_vat_0_raw, aes(x = log_fc, y = log_p))+
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

plot_volc_vat_0_raw

## FDR values with p-threshold of 0.1, and FC > 1.5, comparison vat 0 vs vat 1
mSet_vat0<-Volcano.Anal(mSet_vat0, FALSE, 1.2, 1, F, 0.1, F, "fdr")

cart_volcano_vat_0_fdr <- as.data.frame(mSet_vat0[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_vat_0_fdr$metabolite <- rownames(cart_volcano_vat_0_fdr)
cart_volcano_vat_0_fdr$log_p <- mSet_vat0[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_vat_0_fdr$log_fc <- mSet_vat0[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_vat_0_fdr$inx.up <- mSet_vat0[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_vat_0_fdr$inx.down <- mSet_vat0[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_vat_0_fdr$inx.p <- mSet_vat0[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_vat_0_fdr <- cart_volcano_vat_0_fdr %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

plot_volc_vat_0_fdr <- ggplot(cart_volcano_vat_0_fdr, aes(x = log_fc, y = log_p))+
  geom_jitter(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 50)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))+
  guides(color = "none") +
  # xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_volc_vat_0_fdr

## Extraction of siginficant metabolties into a new data frame for later 
cart_vat_0_sig_metabolites <- cart_volcano_vat_0_raw %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()

# plot_vat_0_sig_metabolite <- cart_vat_0_norm %>%
#   pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
#   filter(metabolite %in% cart_vat_0_sig_metabolites) %>%
#   #  filter(metabolite == "Alanine") %>%
#   group_by(VAT_high) %>%
#   ggplot(aes(x = VAT_high, y = level, color = VAT_high, fill = VAT_high)) +
#   geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA)+
#   geom_pwc(stat = "pwc", method = "t.test", label = "p.signif",
#            bracket.nudge.y = -0.08) +
#   geom_jitter(alpha = 0.3)+
#   scale_x_discrete(label = c("vat 0-1", "vat > 1")) +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
#   scale_color_manual(values = c("grey", "darkmagenta")) +
#   scale_fill_manual(values = c("grey", "darkmagenta")) +
#   guides(color = "none", fill = "none") +
#   xlab("") +
#   ylab("Concentration [norm.]") +
#   facet_wrap(vars(metabolite), scales = "free")+
#   theme_classic()+
#   theme(strip.background = element_blank())
# 
# #plot_vat_0_sig_metabolite
# 
# plot_vat_0_sig_metabolite_adj <- ggadjust_pvalue(plot_vat_0_sig_metabolite,
#                                                    p.adjust.method = "BH", 
#                                                    label = "p.adj.signif")
# 
# plot_vat_0_sig_metabolite_adj

### 1.3 Performing PLSDA to get VIP > 1.5/2

mSet_vat0<-PCA.Anal(mSet_vat0)

## Extraction of PCA component values
cart_vat_0_PCA <- as.data.frame(mSet_vat0[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_vat0 order for samples (used for PLSDA)
Sample_order_PCA_vat_0 <- cart_vat_0_PCA %>%
  mutate(Sample = rownames(cart_vat_0_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
cart_vat_0_PCA$Sample <- rownames(cart_vat_0_PCA)
cart_vat_0_PCA <- cart_vat_0_PCA %>%
  select(Sample, everything())

cart_vat_0_PCA <- left_join(cart_vat_0_PCA, Sample_vat_0, by = "Sample")

cart_vat_0_PCA <- cart_vat_0_PCA %>%
  select(Sample, VAT_high, everything())

cart_vat_0_PCA$VAT_high <- as.character(cart_vat_0_PCA$VAT_high)

## Visualization of PCA comp1 vs comp2
plot_vat_0_PCA <- cart_vat_0_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = VAT_high, fill = VAT_high)) +
  geom_jitter() +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  #  coord_fixed(ratio = 1)+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

plot_vat_0_PCA

## Perform PLSDA
mSet_vat0<-PLSR.Anal(mSet_vat0, reg=TRUE)

## Extraction of PLSDA component values
cart_vat_0_PLSDA <- as.matrix.data.frame(PLSDA_vat_0 <- mSet_vat0[["analSet"]][["plsr"]][["scores"]])
cart_vat_0_PLSDA <- as.data.frame(cart_vat_0_PLSDA)
# Sample information lost in mSet_vat0 upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
cart_vat_0_PLSDA <- cart_vat_0_PLSDA %>%
  mutate(Sample = Sample_order_PCA_vat_0)

cart_vat_0_PLSDA <- left_join(cart_vat_0_PLSDA, Sample_vat_0, by = "Sample")

cart_vat_0_PLSDA$VAT_high <- as.character(cart_vat_0_PLSDA$VAT_high)

## Visualization of PLSDA comp1 vs comp2
plot_vat_0_PLSDA <- cart_vat_0_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = VAT_high, fill = VAT_high)) +
  geom_jitter() +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

plot_vat_0_PLSDA

## PLSDA VIP Scores
cart_vat_0_PLSDA_VIP <- as.data.frame(mSet_vat0[["analSet"]][["plsr"]][["vip.mat"]])
cart_vat_0_PLSDA_VIP <- cart_vat_0_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

plot_vat_0_PLSDA_VIP <- cart_vat_0_PLSDA_VIP %>%
  arrange(desc(`Comp. 1`)) %>%
  filter(`Comp. 1` > 1.5) %>%
  ggplot(aes(x = reorder(metabolite,`Comp. 1`), y = `Comp. 1`, fill = `Comp. 1`)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "PLS-DA Comp.1 VIP Score")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill = "none")+
  coord_flip() +
  scale_fill_gradient(low = "lightgrey", high = "darkmagenta")+
  theme_classic()

plot_vat_0_PLSDA_VIP

vat_0_VIP_2 <- cart_vat_0_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) |>
  select(metabolite) |>
  unlist()|>
  as.vector()

# 1.4 Performing MEBA analysis to extract 
mSet_vat_dyn<-InitDataObjects("pktable", "mf", FALSE)
mSet_vat_dyn<-SetDesignType(mSet_vat_dyn, "time")
mSet_vat_dyn<-Read.TextDataTs(mSet_vat_dyn, "Input_files/CART_VAT_Time.csv", "rowmf");
mSet_vat_dyn<-ReadMetaData(mSet_vat_dyn, "Input_files/CART_VAT_Time_Metafile.csv");
mSet_vat_dyn<-SanityCheckData(mSet_vat_dyn)
mSet_vat_dyn<-ReplaceMin(mSet_vat_dyn);
mSet_vat_dyn<-SanityCheckMeta(mSet_vat_dyn, 1)
mSet_vat_dyn<-SetDataTypeOfMeta(mSet_vat_dyn);
mSet_vat_dyn<-SanityCheckData(mSet_vat_dyn)
#mSet_vat_dyn<-FilterVariable(mSet_vat_dyn, "median", 0, "F", 25, F)
mSet_vat_dyn<-PreparePrenormData(mSet_vat_dyn)
mSet_vat_dyn<-Normalization(mSet_vat_dyn, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

## Extracting normalized values from mSet_vat_dyn object
## Extraction and saving of the normalized data into a new tibble
cart_vat_timeseries <- as.data.frame(mSet_vat_dyn[["dataSet"]][["norm"]])

cart_vat_timeseries <- cart_vat_timeseries %>%
  mutate(Sample = rownames(cart_vat_timeseries)) %>%
  select(-Subject, -Phenotype, -Time) %>%
  select(Sample, everything())

cart_vat_timeseries_meta <- as.data.frame(mSet_vat_dyn[["dataSet"]][["orig.meta.info"]])

cart_vat_timeseries_meta <- cart_vat_timeseries_meta %>%
  mutate(Sample = rownames(cart_vat_timeseries_meta)) %>%
  select(Sample, everything())

cart_vat_timeseries <- left_join(cart_vat_timeseries, cart_vat_timeseries_meta, by = "Sample")

cart_vat_timeseries <- cart_vat_timeseries %>%
  select(Sample, Phenotype, Time, Subject, everything())


#### Multivariate Empirical Bayes Analysis of Variance (MEBA) for Time Series
meta.vec.mb <- c("Phenotype", "Time")
mSet_vat_dyn<-performMB(mSet_vat_dyn, 10)

CART_vat_timeadj_hotelling <- as.data.frame(mSet_vat_dyn[["analSet"]][["MB"]][["stats"]])
CART_vat_timeadj_hotelling <- CART_vat_timeadj_hotelling %>%
  mutate(Sample = rownames(CART_vat_timeadj_hotelling)) %>%
  select(Sample, everything())

CART_vat_timeadj_hotelling  <- CART_vat_timeadj_hotelling %>% 
  filter(Sample != "Phenotype" & Sample != "Subject" & Sample != "Time")

## Calculation of F-statistic and p-value from T-square
## number of dependent variables p = 2
## number of items per group n1 = 32, n2 = 15
## F = n1+n2-p-1 / p(n1+n2-2) * T2 <=> F = 63/128 * T2
## df1 = p, df2 = n1+n2-p-1

CART_vat_timeadj_hotelling <- CART_vat_timeadj_hotelling %>%
  mutate(F_Statistic = 44/90 * CART_vat_timeadj_hotelling$`Hotelling-T2`) %>%
  mutate(p_value = pf(F_Statistic, 2, 44, lower.tail = F))


## Adjust for multiple testing
CART_vat_timeadj_hotelling <- CART_vat_timeadj_hotelling %>%
  mutate(fdr = p.adjust(p_value, method = "fdr"),
         BH = p.adjust(p_value, method = "BH"),
         BF = p.adjust(p_value, method = "bonferroni"))


## Plotting of significant metabolites in MEBA
MEBA_sign_vat_0.1 <- CART_vat_timeadj_hotelling %>%
  filter(fdr < 0.1) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

plot_vat_MEBA_0.1 <- cart_vat_timeseries %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% MEBA_sign_vat_0.1) %>%
  group_by(Time, Phenotype, metabolite) %>%
  summarize(mean_values = mean(level), SD = sd(level), CI = confint(lm(level ~ 1))) %>%
  mutate(CI_lower = CI[,1], CI_upper = CI[,2]) %>%
  ggplot()+
  geom_ribbon(aes(x = as.double(Time), ymin = CI_lower, ymax = CI_upper, fill = Phenotype), alpha = 0.1) +
  geom_point(aes(x = as.factor(Time), y = mean_values, color = Phenotype), alpha = 0.4, size = 2) +
  geom_line(aes(x = as.factor(Time), y = mean_values, group = Phenotype, color = Phenotype), linewidth = 1, alpha = 0.5)+
  #geom_errorbar(aes(ymin = mean_values - SD, ymax = mean_values + SD), width = 0.2)+
  scale_x_discrete(label = c("Day 0", "Day 3", "Day 14"), expand = c(0,0.05)) +
  labs(x = "", y = "Abundance [norm.]")+
  scale_color_manual(values = c("darkgrey", "darkmagenta"),
                     name = "vat",
                     labels = c("vat 0-1", "vat > 1"))+
  scale_fill_manual(values = c("darkgrey", "darkmagenta"),
                    name = "vat",
                    labels = c("vat 0-1", "vat > 1"))+
  guides(fill = "none") +
  # xlab("Time after CAR-T transfusion [days]") +
  # ylab("Concentration [norm.]") +
  facet_wrap(vars(metabolite), scales = "free") +
  # geom_text(data = MEBA_sign_vat_fdr, aes(label= fdr), x = 1, y = 1)+
  theme_classic() +
  theme(strip.background = element_blank())

#1 vor dem ersten Jahr Progress -> in R 1 kein Progress vor dem ersten Jahr -> 1 = responder

plot_vat_MEBA_0.1

# 1.5 Creating one data frame with all singificantly changed metabolites
venn_vat_0 <- list(cart_vat_0_sig_metabolites,
                   vat_0_VIP_2,
                   MEBA_sign_vat_0.1)

p_venn_vat_0 <- ggVennDiagram(venn_vat_0, label_alpha = 0,
                              category.names = c("DEM", "VIP>2", "MEBA"),
                              label = "count") +
  scale_fill_distiller(palette = "Blues") +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none")

p_venn_vat_0

vat_0_metabolites <- c(cart_vat_0_sig_metabolites, vat_0_VIP_2)
vat_0_metabolites <- unique(vat_0_metabolites)

vat_features_0 <- cart_vat_0_norm |>
  pivot_longer(cols = ('Alanine':'Cer-(24:01)'), names_to = "metabolite", values_to = "level")|>
  filter(metabolite %in% vat_0_metabolites) |>
  pivot_wider(names_from = metabolite, values_from = level) |>
  mutate(time = 0) |>
  select(Sample, time, VAT_high, everything())

#Repeat Analysis for dazy 3/5
mSet_vat3<-InitDataObjects("pktable", "stat", FALSE)
mSet_vat3<-Read.TextData(mSet_vat3, "Input_files/CART_VAT_Tag3.csv", "rowu", "disc");
mSet_vat3<-SanityCheckData(mSet_vat3)
mSet_vat3<-ReplaceMin(mSet_vat3);
mSet_vat3<-SanityCheckData(mSet_vat3)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_vat3<-FilterVariable(mSet_vat3, "median", 3, "F")
mSet_vat3<-PreparePrenormData(mSet_vat3)

## Normalization by sum and data scaling based on auto-scaling
mSet_vat3<-Normalization(mSet_vat3, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=23)


## Extraction and saving of the normalized data into a new tibble
cart_vat_3 <- as.data.frame(mSet_vat3[["dataSet"]][["norm"]])

## Load original table and cbind sample name and time point label
cart_vat_3_original <- read_xlsx("Input_files/CART_VAT_Tag3.xlsx", na = "NA")

## Left join original data and normalized to link vat groups
cart_vat_3$Sample <- row.names(cart_vat_3)
cart_vat_3 <- cart_vat_3 %>% select(Sample, everything()) %>%
  arrange(Sample)

Sample_vat_3 <- cart_vat_3_original %>% select(Sample, VAT_high)

cart_vat_3_norm <- left_join(cart_vat_3, Sample_vat_3, by = "Sample")
cart_vat_3_norm <- cart_vat_3_norm %>%
  select(Sample, VAT_high, everything())

cart_vat_3_norm$VAT_high <- as.character(cart_vat_3_norm$VAT_high)

## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 3.35, and FC > 1.5
mSet_vat3<-Volcano.Anal(mSet_vat3, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

cart_volcano_vat_3_raw <- as.data.frame(mSet_vat3[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_vat_3_raw$metabolite <- rownames(cart_volcano_vat_3_raw)
cart_volcano_vat_3_raw$log_p <- mSet_vat3[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_vat_3_raw$log_fc <- mSet_vat3[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_vat_3_raw$inx.up <- mSet_vat3[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_vat_3_raw$inx.down <- mSet_vat3[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_vat_3_raw$inx.p <- mSet_vat3[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_vat_3_raw <- cart_volcano_vat_3_raw %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

plot_volc_vat_3_raw <- ggplot(cart_volcano_vat_3_raw, aes(x = log_fc, y = log_p))+
  geom_point(aes(color = gene_type), alpha = 3.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3", "gray53", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 33)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[13]*"p-value"))+
  guides(color = "none") +
  # xlim(-13, 13)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_volc_vat_3_raw

## FDR values with p-threshold of 3.1, and FC > 1.5, comparison vat 3 vs vat 1
mSet_vat3<-Volcano.Anal(mSet_vat3, FALSE, 1.2, 1, F, 0.1, F, "fdr")

cart_volcano_vat_3_fdr <- as.data.frame(mSet_vat3[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_vat_3_fdr$metabolite <- rownames(cart_volcano_vat_3_fdr)
cart_volcano_vat_3_fdr$log_p <- mSet_vat3[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_vat_3_fdr$log_fc <- mSet_vat3[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_vat_3_fdr$inx.up <- mSet_vat3[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_vat_3_fdr$inx.down <- mSet_vat3[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_vat_3_fdr$inx.p <- mSet_vat3[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_vat_3_fdr <- cart_volcano_vat_3_fdr %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

plot_volc_vat_3_fdr <- ggplot(cart_volcano_vat_3_fdr, aes(x = log_fc, y = log_p))+
  geom_jitter(aes(color = gene_type), alpha = 3.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3","gray53")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 53)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[13]*"FDR"))+
  guides(color = "none") +
  # xlim(-13, 13)
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_volc_vat_3_fdr

## Extraction of siginficant metabolties into a new data frame for later 
cart_vat_3_sig_metabolites <- cart_volcano_vat_3_raw %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()

# plot_vat_3_sig_metabolite <- cart_vat_3_norm %>%
#   pivot_longer(cols = Alanine:`Cer-(24:31)`, names_to = "metabolite", values_to = "level") %>%
#   filter(metabolite %in% cart_vat_3_sig_metabolites) %>%
#   #  filter(metabolite == "Alanine") %>%
#   group_by(VAT_high) %>%
#   ggplot(aes(x = VAT_high, y = level, color = VAT_high, fill = VAT_high)) +
#   geom_boxplot(width = 3.3, alpha = 3.5, outlier.shape = NA)+
#   geom_pwc(stat = "pwc", method = "t.test", label = "p.signif",
#            bracket.nudge.y = -3.38) +
#   geom_jitter(alpha = 3.3)+
#   scale_x_discrete(label = c("vat 3-1", "vat > 1")) +
#   scale_y_continuous(expand = expansion(mult = c(3.35, 3.15)))+
#   scale_color_manual(values = c("grey", "darkmagenta")) +
#   scale_fill_manual(values = c("grey", "darkmagenta")) +
#   guides(color = "none", fill = "none") +
#   xlab("") +
#   ylab("Concentration [norm.]") +
#   facet_wrap(vars(metabolite), scales = "free")+
#   theme_classic()+
#   theme(strip.background = element_blank())
# 
# #plot_vat_3_sig_metabolite
# 
# plot_vat_3_sig_metabolite_adj <- ggadjust_pvalue(plot_vat_3_sig_metabolite,
#                                                    p.adjust.method = "BH", 
#                                                    label = "p.adj.signif")
# 
# plot_vat_3_sig_metabolite_adj

### 1.3 Performing PLSDA to get VIP > 1.5/2

mSet_vat3<-PCA.Anal(mSet_vat3)

## Extraction of PCA component values
cart_vat_3_PCA <- as.data.frame(mSet_vat3[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_vat3 order for samples (used for PLSDA)
Sample_order_PCA_vat_3 <- cart_vat_3_PCA %>%
  mutate(Sample = rownames(cart_vat_3_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
cart_vat_3_PCA$Sample <- rownames(cart_vat_3_PCA)
cart_vat_3_PCA <- cart_vat_3_PCA %>%
  select(Sample, everything())

cart_vat_3_PCA <- left_join(cart_vat_3_PCA, Sample_vat_3, by = "Sample")

cart_vat_3_PCA <- cart_vat_3_PCA %>%
  select(Sample, VAT_high, everything())

cart_vat_3_PCA$VAT_high <- as.character(cart_vat_3_PCA$VAT_high)

## Visualization of PCA comp1 vs comp2
plot_vat_3_PCA <- cart_vat_3_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = VAT_high, fill = VAT_high)) +
  geom_jitter() +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  #  coord_fixed(ratio = 1)+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

plot_vat_3_PCA

## Perform PLSDA
mSet_vat3<-PLSR.Anal(mSet_vat3, reg=TRUE)

## Extraction of PLSDA component values
cart_vat_3_PLSDA <- as.matrix.data.frame(PLSDA_vat_3 <- mSet_vat3[["analSet"]][["plsr"]][["scores"]])
cart_vat_3_PLSDA <- as.data.frame(cart_vat_3_PLSDA)
# Sample information lost in mSet_vat3 upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
cart_vat_3_PLSDA <- cart_vat_3_PLSDA %>%
  mutate(Sample = Sample_order_PCA_vat_3)

cart_vat_3_PLSDA <- left_join(cart_vat_3_PLSDA, Sample_vat_3, by = "Sample")

cart_vat_3_PLSDA$VAT_high <- as.character(cart_vat_3_PLSDA$VAT_high)

## Visualization of PLSDA comp1 vs comp2
plot_vat_3_PLSDA <- cart_vat_3_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = VAT_high, fill = VAT_high)) +
  geom_jitter() +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

plot_vat_3_PLSDA

## PLSDA VIP Scores
cart_vat_3_PLSDA_VIP <- as.data.frame(mSet_vat3[["analSet"]][["plsr"]][["vip.mat"]])
cart_vat_3_PLSDA_VIP <- cart_vat_3_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

plot_vat_3_PLSDA_VIP <- cart_vat_3_PLSDA_VIP %>%
  arrange(desc(`Comp. 1`)) %>%
  filter(`Comp. 1` > 1.5) %>%
  ggplot(aes(x = reorder(metabolite,`Comp. 1`), y = `Comp. 1`, fill = `Comp. 1`)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "PLS-DA Comp.1 VIP Score")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill = "none")+
  coord_flip() +
  scale_fill_gradient(low = "lightgrey", high = "darkmagenta")+
  theme_classic()

vat_3_VIP_2 <- cart_vat_3_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) |>
  select(metabolite) |>
  unlist()|>
  as.vector()

# 1.5 Creating one data frame with all singificantly changed metabolites
venn_vat_3 <- list(cart_vat_3_sig_metabolites,
                   vat_3_VIP_2,
                   MEBA_sign_vat_0.1)

p_venn_vat_3 <- ggVennDiagram(venn_vat_3, label_alpha = 0,
                              category.names = c("DEM", "VIP>2", "MEBA"),
                              label = "count") +
  scale_fill_distiller(palette = "Greys") +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none") 

p_venn_vat_3

#Test all days in one Venn

venn_vat_total <- list(cart_vat_3_sig_metabolites,
                       cart_vat_0_sig_metabolites,
                       vat_0_VIP_2,
                       vat_3_VIP_2,
                       MEBA_sign_vat_0.1)

p_venn_vat_total <- ggVennDiagram(venn_vat_total, label_alpha = 0,
                                  category.names = c("DEM_0", "DEM_3", "VIP>2_0", "VIP>2_3", "MEBA"),
                                  label = "count") +
  # scale_fill_distiller(palette = "Greys") +
  #scale_color_manual(values = c("black", "black", "black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none") 

p_venn_vat_total


vat_3_metabolites <- c(cart_vat_3_sig_metabolites, vat_3_VIP_2)
vat_3_metabolites <- unique(vat_3_metabolites)

vat_features_3 <- cart_vat_3_norm |>
  pivot_longer(cols = ('Alanine':'Cer-(24:01)'), names_to = "metabolite", values_to = "level")|>
  filter(metabolite %in% vat_3_metabolites) |>
  pivot_wider(names_from = metabolite, values_from = level) |>
  mutate(time = 3) |>
  select(Sample, time, VAT_high, everything())




## 1.1 Loading of data for normalization and analyses
mSet_vat14<-InitDataObjects("pktable", "stat", FALSE)
mSet_vat14<-Read.TextData(mSet_vat14, "Input_files/CART_VAT_Tag14.csv", "rowu", "disc");
mSet_vat14<-SanityCheckData(mSet_vat14)
mSet_vat14<-ReplaceMin(mSet_vat14);
mSet_vat14<-SanityCheckData(mSet_vat14)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet_vat14<-FilterVariable(mSet_vat14, "median", 0, "F")
mSet_vat14<-PreparePrenormData(mSet_vat14)

## Normalization by sum and data scaling based on auto-scaling
mSet_vat14<-Normalization(mSet_vat14, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)


## Extraction and saving of the normalized data into a new tibble
cart_vat_14 <- as.data.frame(mSet_vat14[["dataSet"]][["norm"]])

## Load original table and cbind sample name and time point label
cart_vat_14_original <- read_xlsx("Input_files/CART_VAT_Tag14.xlsx", na = "NA")

## Left join original data and normalized to link CRS groups
cart_vat_14$Sample <- row.names(cart_vat_14)
cart_vat_14 <- cart_vat_14 %>% select(Sample, everything()) %>%
  arrange(Sample)

Sample_vat_14 <- cart_vat_14_original %>% select(Sample, VAT_high)

cart_vat_14_norm <- left_join(cart_vat_14, Sample_vat_14, by = "Sample")
cart_vat_14_norm <- cart_vat_14_norm %>%
  select(Sample, VAT_high, everything())

cart_vat_14_norm$VAT_high <- as.character(cart_vat_14_norm$VAT_high)


## 1.2 Fold-change analysis and extraction of significantly changed metabolites with p-threshold of 0.05, and FC > 1.5
mSet_vat14<-Volcano.Anal(mSet_vat14, FALSE, 1.2, 1, F, 0.05, FALSE, "raw")

cart_volcano_vat_14_raw <- as.data.frame(mSet_vat14[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_vat_14_raw$metabolite <- rownames(cart_volcano_vat_14_raw)
cart_volcano_vat_14_raw$log_p <- mSet_vat14[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_vat_14_raw$log_fc <- mSet_vat14[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_vat_14_raw$inx.up <- mSet_vat14[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_vat_14_raw$inx.down <- mSet_vat14[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_vat_14_raw$inx.p <- mSet_vat14[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_vat_14_raw <- cart_volcano_vat_14_raw %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

plot_volc_vat_14_raw <- ggplot(cart_volcano_vat_14_raw, aes(x = log_fc, y = log_p))+
  geom_point(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 20)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p-value"))+
  guides(color = "none") +
  # xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_volc_vat_14_raw

## FDR values with p-threshold of 0.1, and FC > 1.5, comparison CRS 0 vs CRS 1
mSet_vat14<-Volcano.Anal(mSet_vat14, FALSE, 1.2, 1, F, 0.1, F, "fdr")

cart_volcano_vat_14_fdr <- as.data.frame(mSet_vat14[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_vat_14_fdr$metabolite <- rownames(cart_volcano_vat_14_fdr)
cart_volcano_vat_14_fdr$log_p <- mSet_vat14[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_vat_14_fdr$log_fc <- mSet_vat14[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_vat_14_fdr$inx.up <- mSet_vat14[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_vat_14_fdr$inx.down <- mSet_vat14[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_vat_14_fdr$inx.p <- mSet_vat14[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_vat_14_fdr <- cart_volcano_vat_14_fdr %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

plot_volc_vat_14_fdr <- ggplot(cart_volcano_vat_14_fdr, aes(x = log_fc, y = log_p))+
  geom_jitter(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3","gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 50)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))+
  guides(color = "none") +
  # xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_volc_vat_14_fdr

## Plotting the significant metabolites 
cart_vat_14_sig_metabolites <- cart_volcano_vat_14_raw %>%
  filter(inx.p == T) %>%
  filter(inx.down == T | inx.up == T) %>%
  select(metabolite) %>%
  unlist() %>%
  as.vector()


### 1.3 Performing PLSDA to get VIP > 1.5/2

mSet_vat14<-PCA.Anal(mSet_vat14)

## Extraction of PCA component values
cart_vat_14_PCA <- as.data.frame(mSet_vat14[["analSet"]][["pca"]][["x"]])

## Extraction of mSet_vat14 order for samples (used for PLSDA)
Sample_order_PCA_vat_14 <- cart_vat_14_PCA %>%
  mutate(Sample = rownames(cart_vat_14_PCA)) %>%
  select(Sample) %>%
  unlist() %>%
  as.vector()

## Building PCA data frame
cart_vat_14_PCA$Sample <- rownames(cart_vat_14_PCA)
cart_vat_14_PCA <- cart_vat_14_PCA %>%
  select(Sample, everything())

cart_vat_14_PCA <- left_join(cart_vat_14_PCA, Sample_vat_14, by = "Sample")

cart_vat_14_PCA <- cart_vat_14_PCA %>%
  select(Sample, VAT_high, everything())

cart_vat_14_PCA$VAT_high <- as.character(cart_vat_14_PCA$VAT_high)

## Visualization of PCA comp1 vs comp2
plot_vat_14_PCA <- cart_vat_14_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = VAT_high, fill = VAT_high)) +
  geom_jitter() +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  #  coord_fixed(ratio = 1)+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

plot_vat_14_PCA

## Perform PLSDA
mSet_vat14<-PLSR.Anal(mSet_vat14, reg=TRUE)

## Extraction of PLSDA component values
cart_vat_14_PLSDA <- as.matrix.data.frame(PLSDA_vat_14 <- mSet_vat14[["analSet"]][["plsr"]][["scores"]])
cart_vat_14_PLSDA <- as.data.frame(cart_vat_14_PLSDA)
# Sample information lost in mSet_vat14 upon PLSDA analysis
# Sample order from PCA extracted and saved in Sample_order_PCA
cart_vat_14_PLSDA <- cart_vat_14_PLSDA %>%
  mutate(Sample = Sample_order_PCA_vat_14)

cart_vat_14_PLSDA <- left_join(cart_vat_14_PLSDA, Sample_vat_14, by = "Sample")

cart_vat_14_PLSDA$VAT_high <- as.character(cart_vat_14_PLSDA$VAT_high)

## Visualization of PLSDA comp1 vs comp2
plot_vat_14_PLSDA <- cart_vat_14_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = VAT_high, fill = VAT_high)) +
  geom_jitter() +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.8) +
  scale_color_manual(values = c("grey", "darkmagenta"))+
  scale_fill_manual(values = c("grey", "darkmagenta"), 
                    labels = c("vat 0-1", "vat > 1"))+
  guides(color = "none")+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

plot_vat_14_PLSDA

## PLSDA VIP Scores
cart_vat_14_PLSDA_VIP <- as.data.frame(mSet_vat14[["analSet"]][["plsr"]][["vip.mat"]])
cart_vat_14_PLSDA_VIP <- cart_vat_14_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

plot_vat_14_PLSDA_VIP <- cart_vat_14_PLSDA_VIP %>%
  arrange(desc(`Comp. 1`)) %>%
  filter(`Comp. 1` > 1.5) %>%
  ggplot(aes(x = reorder(metabolite,`Comp. 1`), y = `Comp. 1`, fill = `Comp. 1`)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "PLS-DA Comp.1 VIP Score")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill = "none")+
  coord_flip() +
  scale_fill_gradient(low = "lightgrey", high = "darkmagenta")+
  theme_classic()

vat_14_VIP_2 <- cart_vat_14_PLSDA_VIP %>%
  filter(`Comp. 1` > 1.5) |>
  select(metabolite) |>
  unlist()|>
  as.vector()

# 1.5 Creating one data frame with all singificantly changed metabolites
venn_vat_14 <- list(cart_vat_14_sig_metabolites,
                    vat_14_VIP_2,
                    MEBA_sign_vat_0.1)

p_venn_vat_14 <- ggVennDiagram(venn_vat_14, label_alpha = 0,
                               category.names = c("DEM", "VIP>2", "MEBA"),
                               label = "count") +
  scale_fill_distiller(palette = "Greys") +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none") 

p_venn_vat_14

#Test all days in one Venn

venn_vat_total <- list(cart_vat_14_sig_metabolites,
                       cart_vat_0_sig_metabolites,
                       vat_0_VIP_2,
                       vat_14_VIP_2,
                       MEBA_sign_vat_0.1)

p_venn_vat_total <- ggVennDiagram(venn_vat_total, label_alpha = 0,
                                  category.names = c("DEM_0", "DEM_3", "VIP>2_0", "VIP>2_3", "MEBA"),
                                  label = "count") +
  # scale_fill_distiller(palette = "Greys") +
  #scale_color_manual(values = c("black", "black", "black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none") 

p_venn_vat_total


vat_14_metabolites <- c(cart_vat_14_sig_metabolites, vat_14_VIP_2)
vat_14_metabolites <- unique(vat_14_metabolites)

vat_features_14 <- cart_vat_14_norm |>
  pivot_longer(cols = ('Alanine':'Cer-(24:01)'), names_to = "metabolite", values_to = "level")|>
  filter(metabolite %in% vat_14_metabolites) |>
  pivot_wider(names_from = metabolite, values_from = level) |>
  mutate(time = 14) |>
  select(Sample, time, VAT_high, everything())



# plot_vat_14_sig_metabolite <- cart_vat_14_norm %>%
#   pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
#   filter(metabolite %in% cart_vat_14_sig_metabolites) %>%
#   #  filter(metabolite == "Alanine") %>%
#   group_by(VAT_high) %>%
#   ggplot(aes(x = VAT_high, y = level, color = VAT_high, fill = VAT_high)) +
#   geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA)+
#   geom_pwc(stat = "pwc", method = "t.test", label = "p.signif",
#            bracket.nudge.y = -0.08) +
#   geom_jitter(alpha = 0.3)+
#   scale_x_discrete(label = c("NR", "R")) +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
#   scale_color_manual(values = c("grey", "darkgreen")) +
#   scale_fill_manual(values = c("grey", "darkgreen")) +
#   guides(color = "none", fill = "none") +
#   xlab("") +
#   ylab("Concentration [norm.]") +
#   facet_wrap(vars(metabolite), scales = "free")+
#   theme_classic()+
#   theme(strip.background = element_blank())

# #plot_vat_14_sig_metabolite
# 
# plot_vat_14_sig_metabolite_adj <- ggadjust_pvalue(plot_vat_14_sig_metabolite,
#                                                   p.adjust.method = "BH", 
#                                                   label = "p.adj.signif")

#plot_vat_14_sig_metabolite_adj

#### Correlation analysis between metabolite levels and clinical data ----
## Loading meta data and analysis of correlation between significantly altered metabolites and bodycomp / inflammatory markers
meta_master <- read_xlsx("Input_files/CART_meta_master_neu.xlsx", na = "NA")
str(meta_master)

meta_master$Ferritin <- as.numeric(meta_master$Ferritin)
meta_master$Neutro <- as.numeric(meta_master$Neutro)
meta_master$IPI <- as.character(meta_master$IPI)
meta_master$ECOG <- as.character(meta_master$ECOG)
meta_master$`Sample ID Number` <- as.character(meta_master$`Sample ID Number`)

##combining metabolites of selected features with clinical data
vat_features_0_clinic <- vat_features_0 |>
  separate(Sample, into = c("Sample", "Letter", "Number"), sep = "_")

vat_features_0_clinic <- left_join(vat_features_0_clinic, meta_master, by = c("Number" = "Sample ID Number"))

str(vat_features_0_clinic)


vat_features_0_clinic$CRS_high <- as.double(vat_features_0_clinic$CRS_high)


vat_features_3_clinic <- vat_features_3 |>
  separate(Sample, into = c("Sample", "Letter", "Number"), sep = "_")

vat_features_3_clinic <- left_join(vat_features_3_clinic, meta_master, by = c("Number" = "Sample ID Number"))
str(vat_features_3_clinic)

vat_features_3_clinic$CRS_high <- as.double(vat_features_3_clinic$CRS_high)

vat_features_14_clinic <- vat_features_14 |>
  separate(Sample, into = c("Sample", "Letter", "Number"), sep = "_")

vat_features_14_clinic <- left_join(vat_features_14_clinic, meta_master, by = c("Number" = "Sample ID Number"))
str(vat_features_14_clinic)

vat_features_14_clinic$CRS_high <- as.double(vat_features_14_clinic$CRS_high)

# clinic_corr <- cart_icans_0_correlation |>
#   select(CRP, Albumin, Ferritin, Leuko, Neutro, Lympho, VAT, SAT, BMI, SMI,
#          vat, WtHR, TAT, PMI, STLV)|>
#   colnames()|>
#   as.vector()
# 

## Log Regression MVA ----

## MULTIVARIABLE log regression analysis focusing on icans day 0 metabolites and lab values
#Metadata can be loaded from cart_crs0_correlation table

#Defining an optimal baseline model for icans development

summary(glm(CRS_high ~ Construct, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ CRP, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ Albumin, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ Geschlecht, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ STLV, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ Ferritin, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ LDH, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ Costim, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ EASIX, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ mEASIX, vat_features_0_clinic, family =binomial))
summary(glm(as.factor(CRS_high) ~ PLT, vat_features_0_clinic, family =binomial))


crs_baseline_model <- glm(CRS_high ~ Costim + CRP  + Geschlecht + STLV + LDH + EASIX, vat_features_0_clinic, family =binomial,
                          na.action = na.omit)
summary(crs_baseline_model)
step(crs_baseline_model, direction = "both") 

#Optimal baseline model for CRS does not contain any covariates

# Testing metabolites against CAR-T model

crs_vat_0_glm_MVA <- data.frame(marker = character(),
                                coefficient = numeric(),
                                std_error = numeric(),
                                p_value = numeric(),
                                lower95 = numeric(),
                                upper95 = numeric(),
                                stringsAsFactors = FALSE)


for(i in vat_0_metabolites) {
  # Fit logistic regression model
  formula_str <- paste0("as.factor(CRS_high) ~ `",i,"` + EASIX")
  
  #print(formula_str)
  model <- glm(formula_str, data = vat_features_0_clinic, family = binomial)
  
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
  crs_vat_0_glm_MVA <- rbind(crs_vat_0_glm_MVA, marker_results)
}

crs_vat_0_glm_MVA <- crs_vat_0_glm_MVA|>
  mutate(OR = exp(coefficient), FDR = p.adjust(p_value, method = "BH"))

p_crs_vat_0_glm <- crs_vat_0_glm_MVA |>
  ggplot()+
  geom_point(aes(x = reorder(marker, coefficient) , y = OR),
             size = 4, shape = 19, color = "darkred", alpha = 0.7)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = upper95), color = "darkred")+
  geom_text(aes(x = marker, y = -0.7, 
                label = paste("q =", round(p_value,3))))+
  coord_flip(ylim = c(-1, 6))+
  labs(y = "Odds Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

p_crs_vat_0_glm

## Log Regression Day 3

crs_vat_3_glm_MVA <- data.frame(marker = character(),
                                coefficient = numeric(),
                                std_error = numeric(),
                                p_value = numeric(),
                                lower95 = numeric(),
                                upper95 = numeric(),
                                stringsAsFactors = FALSE)


for(i in vat_3_metabolites) {
  # Fit logistic regression model
  formula_str <- paste0("as.factor(CRS_high) ~ `",i,"` + EASIX")
  
  #print(formula_str)
  model <- glm(formula_str, data = vat_features_3_clinic, family = binomial)
  
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
  crs_vat_3_glm_MVA <- rbind(crs_vat_3_glm_MVA, marker_results)
}

crs_vat_3_glm_MVA <- crs_vat_3_glm_MVA|>
  mutate(OR = exp(coefficient), FDR = p.adjust(p_value, method = "BH"))

p_crs_vat_3_glm <- crs_vat_3_glm_MVA |>
  ggplot()+
  geom_point(aes(x = reorder(marker, coefficient) , y = OR),
             size = 4, shape = 19, color = "darkred", alpha = 3.7)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = upper95), color = "darkred")+
  geom_text(aes(x = marker, y = -3.7, 
                label = paste("q =", round(p_value,3))))+
  coord_flip()+
  labs(y = "Odds Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

p_crs_vat_3_glm


## Combined results of day 0 and day 3

crs_vat_0_glm_MVA <- crs_vat_0_glm_MVA |>
  mutate(time = 0) |>
  mutate(FDR = p.adjust(p_value, method = "BH"),
         time = as.factor(time))

crs_vat_3_glm_MVA <- crs_vat_3_glm_MVA |>
  mutate(time =3)|>
  mutate(FDR = p.adjust(p_value, method = "BH"),
         time = as.factor(time))

crs_vat_glm_combined <- rbind(crs_vat_0_glm_MVA, crs_vat_3_glm_MVA)

crs_vat_glm_combined_0.05 <- crs_vat_glm_combined |>
  filter(FDR < 0.1)

p_crs_vat_glm <- crs_vat_glm_combined_0.05  |>
  ggplot()+
  geom_jitter(aes(x = reorder(marker, coefficient) , y = OR, color = time),
              size = 4, shape = 19, width = 0.2, alpha = 0.6)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = upper95))+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = paste0("p=", round(p_value,3))))+
  coord_flip(ylim = c(0, 6))+
  scale_color_manual(values = c("darkblue"))+
  labs(y = "Odds Ratio (95%CI)", x = "", color = "Timepoint")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

p_crs_vat_glm

## Calculation of Cox regression models for crs
## Setting the baseline model

tox_master <- meta_master |>
  select(ID_CARTBC, ICANS_onset, ICANS_high, CRS_onset, CRS_high)|>
  mutate(start = 0)|>
  rename(id = ID_CARTBC)

crs_cox_bm <- meta_master |>
  select(ID_CARTBC, CRS_onset, CRS_high,
         STLV,Ferritin, LDH, Geschlecht, Costim, CRP, Albumin, EASIX, mEASIX, PLT)|>
  mutate(start = 0)|>
  rename(id = ID_CARTBC)

str(crs_cox_bm)

crs_cox_bm$CRS_high <- as.logical(crs_cox_bm$CRS_high)

cox_crs_base_model <- coxph(Surv(CRS_onset, CRS_high)~ STLV + CRP + Costim + LDH + Geschlecht , crs_cox_bm)

summary(cox_crs_base_model)

step(coxph(Surv(CRS_onset, CRS_high)~ STLV + CRP + LDH + Costim + EASIX + PLT +
             Geschlecht, na.omit(crs_cox_bm)), direction = "both")

#optimales Modell ohne jegliche Variablen
summary(coxph(Surv(CRS_onset, CRS_high)~ STLV, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ CRP, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ Albumin, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ Geschlecht, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ LDH, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ Ferritin, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ Costim, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ EASIX, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ mEASIX, crs_cox_bm))
summary(coxph(Surv(CRS_onset, CRS_high)~ PLT, crs_cox_bm))

## Calculation of Cox model for day 0 metabolite levels
crs_vat_0_cox_MVA <- data.frame(marker = character(),
                                HR = numeric(),
                                lower95 = numeric(),
                                higher95 = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (i in vat_0_metabolites) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(CRS_onset, CRS_high)  ~ `",i,"` + EASIX")
  model_cox <- coxph(as.formula(formula_str), data = vat_features_0_clinic)
  
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
  crs_vat_0_cox_MVA  <- rbind(crs_vat_0_cox_MVA, marker_results_cox)
}

crs_vat_0_cox_MVA |>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkred", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95,
                    ymax = higher95), width = 0.1,
                color = "darkred")+
  geom_text(aes(x = marker, y = 0.1,
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip(ylim=c(0, 3))+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

crs_vat_3_cox_MVA <- data.frame(marker = character(),
                                HR = numeric(),
                                lower95 = numeric(),
                                higher95 = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)


for (i in vat_3_metabolites) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(CRS_onset, CRS_high)  ~ `",i,"` + EASIX")
  model_cox <- coxph(as.formula(formula_str), data = vat_features_3_clinic)
  
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
  crs_vat_3_cox_MVA  <- rbind(crs_vat_3_cox_MVA, marker_results_cox)
}

crs_vat_3_cox_MVA |>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkred", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95,
                    ymax = higher95), width = 0.1,
                color = "darkred")+
  geom_text(aes(x = marker, y = 0.1,
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip(ylim=c(0, 3))+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

## Combining Hazard Ratio tables from day 0 and day 3

crs_vat_0_cox_MVA <- crs_vat_0_cox_MVA |>
  mutate(time = 0) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))

crs_vat_3_cox_MVA <- crs_vat_3_cox_MVA |>
  mutate(time =3) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))

crs_vat_cox_combined <- rbind(crs_vat_0_cox_MVA, crs_vat_3_cox_MVA)

crs_vat_cox_combined_0.05 <- crs_vat_cox_combined |>
  filter(FDR < 0.1)|>
  mutate(time = as.factor(time))

p_crs_vat_cox <- crs_vat_cox_combined_0.05  |>
  ggplot()+
  geom_jitter(aes(x = reorder(marker, HR) , y = HR, color = time),
              size = 4, shape = 19, width = 0.2, alpha = 0.6)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = higher95))+
  # geom_text(aes(x = marker, y = -0.7,
  #               label = paste0("p=", round(p_value,3))))+
  coord_flip(ylim = c(0, 4))+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  labs(color = "Timepoint")+
  scale_color_manual(values = c("darkblue"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

p_crs_vat_cox

## Cumulative incidence for top metabolites as visualization

# crs_vat_cox_0 <- crs_vat_cox_combined_0.05 |>
#   filter(time == 0) |>
#   select(marker) |>
#   unlist() |>
#   as.vector()
# 
# p_ci_crs_vat_0 <- list()
# 
# for (metabolite in crs_vat_cox_0) {
#   # Directly use the formula in survfit2
#   fit <- survfit2(Surv(CRS_onset, as.numeric(`CRS_high`)) ~ ifelse(get(metabolite, vat_features_0_clinic) > mean(get(metabolite, vat_features_0_clinic)), 1, 0), data = vat_features_0_clinic)
#   
#   # Generate the plot
#   p <- ggsurvplot(fit,
#                   fun = "event",
#                   ylim= c(0,1), xlim=c(0, 19), break.x.by = 3, ylab="CRS CI", xlab="Days after CAR-T infusion",
#                   pval= TRUE, pval.coord = c(1, 0.9), pval.size = 3,
#                   size = 1.15,
#                   axes.offset = FALSE,
#                   risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
#                   tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
#                   conf.int = FALSE,
#                   ggtheme = theme_classic2(10),
#                   font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
#                   fontsize=3,
#                   legend.labs = c("Low", "High"),
#                   legend.title = c(paste(metabolite)),
#                   palette = c("lightgrey", "darkgrey"))
#   # Add the plot to the list
#   p_ci_crs_vat_0[[metabolite]] <- p
# }
# 
# grobs_ci_crs_vat_0 <- lapply(p_ci_crs_vat_0, function(x) ggplotGrob(x$plot))
# do.call(grid.arrange, c(grobs_ci_crs_vat_0, ncol = 5))


crs_vat_cox_3 <- crs_vat_cox_combined_0.05 |>
  filter(time == 3) |>
  select(marker) |>
  unlist() |>
  as.vector()

p_ci_crs_vat_3 <- list()

for (metabolite in crs_vat_cox_3) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(CRS_onset, as.numeric(`CRS_high`)) ~ ifelse(get(metabolite, vat_features_3_clinic) > mean(get(metabolite, vat_features_3_clinic)), 1, 0), data = vat_features_3_clinic)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  fun = "event",
                  ylim= c(0,1), xlim=c(0, 19), break.x.by = 3, ylab="CRS CI", xlab="Days after CAR-T infusion",
                  pval= TRUE, pval.coord = c(1, 0.9), pval.size = 3,
                  size = 1.15,
                  axes.offset = FALSE,
                  risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                  tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                  conf.int = FALSE,
                  ggtheme = theme_classic2(10),
                  font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                  fontsize=3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste(metabolite)),
                  palette = c("lightblue", "darkblue"))
  # Add the plot to the list
  p_ci_crs_vat_3[[metabolite]] <- p
}

grobs_ci_crs_vat_3 <- lapply(p_ci_crs_vat_3, function(x) ggplotGrob(x$plot))
do.call(grid.arrange, c(grobs_ci_crs_vat_3, ncol = 5))

p_ci_crs_vat_3$`PEA-(40:07)`

p_ci_crs_vat_3$`LPC-(20:03)`

#### PFS and OS analyses based on metabolites different in vat ----
survival_master <- meta_master |>
  select(`Sample ID Number`, PFS_days, PFS_event, OS_days, OS_event)|>
  mutate(start = 0)|>
  rename(id = `Sample ID Number`)

str(survival_master)


#### Calculation of Cox proportional hazard model with time varying covariates based on MEBA
## Table with significant metabolites in form> id, time, value
cox_surv_vat_dyn_metabolites <- cart_vat_timeseries %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% MEBA_sign_vat_0.1)|>
  pivot_wider(names_from = metabolite, values_from = level)|>
  arrange(Subject)

cox_surv_vat_dyn_metabolites <- cox_surv_vat_dyn_metabolites |>
  mutate(Time = case_when(
    Time == 0 ~ 0,
    Time == 1 ~ 5,
    Time == 2 ~ 14
  ))


cox_surv_vat_dyn_metabolites <- cox_surv_vat_dyn_metabolites|>
  separate(Sample, into = c("Sample", "suffix", "id"), sep = "_", remove = FALSE) |>
  select(-Sample, -suffix, -Phenotype, -Subject)

survival_base_model <- meta_master |>
  select(`Sample ID Number`, PFS_days, PFS_event, OS_days, OS_event,
         STLV, CRP, LDH, Geschlecht, Costim, TAT)|>
  mutate(start = 0)|>
  rename(id = `Sample ID Number`)

#survival_base_model$`STLV [ml]` <- as.numeric(survival_base_model$`STLV [ml]`)

survival_base_model

#cox_pfs_base_model 

coxph(Surv(PFS_days, PFS_event)~ STLV , survival_base_model)
coxph(Surv(PFS_days, PFS_event)~ CRP, survival_base_model) #signifikant
coxph(Surv(PFS_days, PFS_event)~ log(LDH), survival_base_model) #signifikant
coxph(Surv(PFS_days, PFS_event)~ Geschlecht, survival_base_model)
coxph(Surv(PFS_days, PFS_event)~ Costim, survival_base_model)
coxph(Surv(PFS_days, PFS_event)~ TAT, survival_base_model)

#step(cox_pfs_base_model, direction = "both") #nur LDH bleibt

resp_ids <- cox_surv_vat_dyn_metabolites |>
  select(id)|>
  unique() |>
  unlist() |>
  as.vector()

survival_resp <- survival_master|>
  filter(id %in% resp_ids)

resp_base_model <- survival_base_model|>
  filter(id %in% resp_ids)

vat_pfs_bm_cox <- tmerge(resp_base_model, resp_base_model, id=id, 
                         endpt=event(survival_resp$PFS_days,survival_resp$PFS_event))

vat_os_bm_cox <- tmerge(resp_base_model, resp_base_model, id=id, 
                        endpt=event(survival_resp$OS_days,survival_resp$OS_event))

code_lines <- sapply(MEBA_sign_vat_0.1, function(metabolite) {
  paste0("`", metabolite, "`", "=tdc(Time,`", metabolite, "`)")
})

code <- paste("resp_pfs_bm_cox <- tmerge(resp_pfs_bm_cox, cox_resp_metabolites, id=id, ", 
              paste(code_lines, collapse=", "), 
              ")", 
              sep="")

writeLines(code)

vat_pfs_bm_cox <- tmerge(vat_pfs_bm_cox, cox_surv_vat_dyn_metabolites, id=id, `AC-(06:0)`=tdc(Time,`AC-(06:0)`), `AC-(02:0)`=tdc(Time,`AC-(02:0)`), `CE-(16:00) NH4`=tdc(Time,`CE-(16:00) NH4`), `CE-(18:01) NH4`=tdc(Time,`CE-(18:01) NH4`), `CE-(20:03) NH4`=tdc(Time,`CE-(20:03) NH4`), `AC-(06:1)`=tdc(Time,`AC-(06:1)`), `CA-(14:00)`=tdc(Time,`CA-(14:00)`), `AC-(08:1)`=tdc(Time,`AC-(08:1)`), `CE-(20:04) NH4`=tdc(Time,`CE-(20:04) NH4`), `CE-(18:02) NH4`=tdc(Time,`CE-(18:02) NH4`), `CA-(24:00)`=tdc(Time,`CA-(24:00)`), `AC-(13:0)`=tdc(Time,`AC-(13:0)`), `CE-(22:05) NH4`=tdc(Time,`CE-(22:05) NH4`), `AC-(14:1)`=tdc(Time,`AC-(14:1)`), `Xanthine`=tdc(Time,`Xanthine`), `Cholesterol`=tdc(Time,`Cholesterol`), `PlasC-(34:01)`=tdc(Time,`PlasC-(34:01)`), `CA-(22:00)`=tdc(Time,`CA-(22:00)`), `LPC-(24:03)`=tdc(Time,`LPC-(24:03)`), `Taurine`=tdc(Time,`Taurine`), `CE-(18:03) NH4`=tdc(Time,`CE-(18:03) NH4`), `PA-(38:04)`=tdc(Time,`PA-(38:04)`), `LPC-(20:05)`=tdc(Time,`LPC-(20:05)`), `PC-(36:03)`=tdc(Time,`PC-(36:03)`), `AC-(16:1)`=tdc(Time,`AC-(16:1)`), `CA-(20:00)`=tdc(Time,`CA-(20:00)`), `AC-(08:0)`=tdc(Time,`AC-(08:0)`), `AC-(05:0)`=tdc(Time,`AC-(05:0)`), `LPC-(14:00)`=tdc(Time,`LPC-(14:00)`), `TAG-(56:06) NH4`=tdc(Time,`TAG-(56:06) NH4`), `Creatine`=tdc(Time,`Creatine`), `TAG-(54:05) NH4`=tdc(Time,`TAG-(54:05) NH4`), `Creatinine`=tdc(Time,`Creatinine`), `CA-(16:01)`=tdc(Time,`CA-(16:01)`))

vat_os_bm_cox <- tmerge(vat_os_bm_cox, cox_surv_vat_dyn_metabolites, id=id, `AC-(06:0)`=tdc(Time,`AC-(06:0)`), `AC-(02:0)`=tdc(Time,`AC-(02:0)`), `CE-(16:00) NH4`=tdc(Time,`CE-(16:00) NH4`), `CE-(18:01) NH4`=tdc(Time,`CE-(18:01) NH4`), `CE-(20:03) NH4`=tdc(Time,`CE-(20:03) NH4`), `AC-(06:1)`=tdc(Time,`AC-(06:1)`), `CA-(14:00)`=tdc(Time,`CA-(14:00)`), `AC-(08:1)`=tdc(Time,`AC-(08:1)`), `CE-(20:04) NH4`=tdc(Time,`CE-(20:04) NH4`), `CE-(18:02) NH4`=tdc(Time,`CE-(18:02) NH4`), `CA-(24:00)`=tdc(Time,`CA-(24:00)`), `AC-(13:0)`=tdc(Time,`AC-(13:0)`), `CE-(22:05) NH4`=tdc(Time,`CE-(22:05) NH4`), `AC-(14:1)`=tdc(Time,`AC-(14:1)`), `Xanthine`=tdc(Time,`Xanthine`), `Cholesterol`=tdc(Time,`Cholesterol`), `PlasC-(34:01)`=tdc(Time,`PlasC-(34:01)`), `CA-(22:00)`=tdc(Time,`CA-(22:00)`), `LPC-(24:03)`=tdc(Time,`LPC-(24:03)`), `Taurine`=tdc(Time,`Taurine`), `CE-(18:03) NH4`=tdc(Time,`CE-(18:03) NH4`), `PA-(38:04)`=tdc(Time,`PA-(38:04)`), `LPC-(20:05)`=tdc(Time,`LPC-(20:05)`), `PC-(36:03)`=tdc(Time,`PC-(36:03)`), `AC-(16:1)`=tdc(Time,`AC-(16:1)`), `CA-(20:00)`=tdc(Time,`CA-(20:00)`), `AC-(08:0)`=tdc(Time,`AC-(08:0)`), `AC-(05:0)`=tdc(Time,`AC-(05:0)`), `LPC-(14:00)`=tdc(Time,`LPC-(14:00)`), `TAG-(56:06) NH4`=tdc(Time,`TAG-(56:06) NH4`), `Creatine`=tdc(Time,`Creatine`), `TAG-(54:05) NH4`=tdc(Time,`TAG-(54:05) NH4`), `Creatinine`=tdc(Time,`Creatinine`), `CA-(16:01)`=tdc(Time,`CA-(16:01)`))

# resp_pfs_bm_cox <- tmerge(vat_pfs_bm_cox, cox_surv_vat_dyn_metabolites, id=id, 
#                           MEBA_sign_vat_0.1=tdc(Time,`MEBA_sign_vat_0.1`))

# resp_os_bm_cox <- tmerge(resp_os_bm_cox, cox_resp_metabolites, id=id, 
#                          `PlasC-(38:05)`=tdc(Time,`PlasC-(38:05)`),
#                          `CE-(22:06) NH4`=tdc(Time,`CE-(22:06) NH4`),
#                          `LPEA-(22:04)`=tdc(Time,`LPEA-(22:04)`),
#                          `LPI-(16:00)`=tdc(Time,`LPI-(16:00)`),
#                          `LPI-(16:01)`=tdc(Time,`LPI-(16:01)`),
#                          `PC-(38:06)`=tdc(Time,`PC-(38:06)`),
#                          `PC-(38:05)`=tdc(Time,`PC-(38:05)`),
#                          `PI-(34:01)`=tdc(Time,`PI-(34:01)`),
#                          `PI-(36:03)`=tdc(Time,`PI-(36:03)`),
#                          `PI-(38:05)`=tdc(Time,`PI-(38:05)`),
#                          `PlasC-(38:05)`=tdc(Time,`PlasC-(38:05)`),
#                          `PlasC-(40:06)`=tdc(Time,`PlasC-(40:06)`),
#                          `Cer-(16:00)`=tdc(Time,`Cer-(16:00)`),
#                          `Cer-(18:00)`=tdc(Time,`Cer-(18:00)`),
#                          Pyruvate=tdc(Time,Pyruvate),
#                          GluAsn=tdc(Time,GluAsn),
#                          Phenylalanine=tdc(Time,Phenylalanine))



###Multivariate Cox model with survival base model

vat_dyn_pfs_cox <- data.frame(marker = character(),
                              HR = numeric(),
                              lower95 = numeric(),
                              higher95 = numeric(),
                              p_value = numeric(),
                              stringsAsFactors = FALSE)

for (i in MEBA_sign_vat_0.1) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(tstart, tstop, endpt)  ~ `",i,"`+ log(LDH) + CRP")
  model_cox <- coxph(as.formula(formula_str), data = vat_pfs_bm_cox)
  
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
  vat_dyn_pfs_cox  <- rbind(vat_dyn_pfs_cox , marker_results_cox)
}

vat_dyn_pfs_cox  |>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkgreen", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95, 
                    ymax = higher95), width = 0.1,
                color = "darkgreen")+
  geom_text(aes(x = marker, y = 3.5, 
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )


vat_dyn_os_cox  <- data.frame(marker = character(),
                              HR = numeric(),
                              lower95 = numeric(),
                              higher95 = numeric(),
                              p_value = numeric(),
                              stringsAsFactors = FALSE)

for (i in MEBA_sign_vat_0.1) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(tstart, tstop, endpt)  ~ `",i,"`+ log(LDH) + CRP")
  model_cox <- coxph(as.formula(formula_str), data = vat_os_bm_cox)
  
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
  vat_dyn_os_cox <- rbind(vat_dyn_os_cox, marker_results_cox)
}

vat_dyn_os_cox |>
  # filter(marker != "LPI-(16:01)")|>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkgreen", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95, 
                    ymax = higher95), width = 0.1,
                color = "darkgreen")+
  geom_text(aes(x = marker, y = 6.5, 
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

### Analysis of different time points for PFS
vat_0_pfs_cox_MVA <- data.frame(marker = character(),
                                HR = numeric(),
                                lower95 = numeric(),
                                higher95 = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (i in vat_0_metabolites) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(PFS_days, PFS_event)  ~ `",i,"`+ log(LDH) + CRP")
  model_cox <- coxph(as.formula(formula_str), data = vat_features_0_clinic)
  
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
  vat_0_pfs_cox_MVA  <- rbind(vat_0_pfs_cox_MVA, marker_results_cox)
}

vat_0_pfs_cox_MVA |>
  # filter(marker != "LPI-(16:01)")|>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkgreen", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95, 
                    ymax = higher95), width = 0.1,
                color = "darkgreen")+
  geom_text(aes(x = marker, y = 6.5, 
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

vat_3_pfs_cox_MVA <- data.frame(marker = character(),
                                HR = numeric(),
                                lower95 = numeric(),
                                higher95 = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (i in vat_3_metabolites) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(PFS_days, PFS_event)  ~ `",i,"`+ log(LDH) + CRP")
  model_cox <- coxph(as.formula(formula_str), data = vat_features_3_clinic)
  
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
  vat_3_pfs_cox_MVA  <- rbind(vat_3_pfs_cox_MVA, marker_results_cox)
}

vat_3_pfs_cox_MVA |>
  # filter(marker != "LPI-(16:01)")|>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkgreen", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95, 
                    ymax = higher95), width = 0.1,
                color = "darkgreen")+
  geom_text(aes(x = marker, y = 6.5, 
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )


vat_14_pfs_cox_MVA <- data.frame(marker = character(),
                                 HR = numeric(),
                                 lower95 = numeric(),
                                 higher95 = numeric(),
                                 p_value = numeric(),
                                 stringsAsFactors = FALSE)

for (i in vat_14_metabolites) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(PFS_days, PFS_event)  ~ `",i,"`+ log(LDH) + CRP")
  model_cox <- coxph(as.formula(formula_str), data = vat_features_14_clinic)
  
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
  vat_14_pfs_cox_MVA  <- rbind(vat_14_pfs_cox_MVA, marker_results_cox)
}

vat_14_pfs_cox_MVA |>
  # filter(marker != "LPI-(16:01)")|>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkgreen", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95, 
                    ymax = higher95), width = 0.1,
                color = "darkgreen")+
  geom_text(aes(x = marker, y = 6.5, 
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )


### Analysis of different time points for OS
vat_0_os_cox_MVA <- data.frame(marker = character(),
                               HR = numeric(),
                               lower95 = numeric(),
                               higher95 = numeric(),
                               p_value = numeric(),
                               stringsAsFactors = FALSE)

for (i in vat_0_metabolites) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(OS_days, OS_event)  ~ `",i,"`+ log(LDH) + CRP")
  model_cox <- coxph(as.formula(formula_str), data = vat_features_0_clinic)
  
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
  vat_0_os_cox_MVA  <- rbind(vat_0_os_cox_MVA, marker_results_cox)
}

vat_0_os_cox_MVA |>
  # filter(marker != "LPI-(16:01)")|>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkgreen", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95, 
                    ymax = higher95), width = 0.1,
                color = "darkgreen")+
  geom_text(aes(x = marker, y = 6.5, 
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

vat_3_os_cox_MVA <- data.frame(marker = character(),
                               HR = numeric(),
                               lower95 = numeric(),
                               higher95 = numeric(),
                               p_value = numeric(),
                               stringsAsFactors = FALSE)

for (i in vat_3_metabolites) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(OS_days, OS_event)  ~ `",i,"`+ log(LDH) + CRP")
  model_cox <- coxph(as.formula(formula_str), data = vat_features_3_clinic)
  
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
  vat_3_os_cox_MVA  <- rbind(vat_3_os_cox_MVA, marker_results_cox)
}

vat_3_os_cox_MVA |>
  # filter(marker != "LPI-(16:01)")|>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkgreen", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95, 
                    ymax = higher95), width = 0.1,
                color = "darkgreen")+
  geom_text(aes(x = marker, y = 6.5, 
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )


vat_14_os_cox_MVA <- data.frame(marker = character(),
                                HR = numeric(),
                                lower95 = numeric(),
                                higher95 = numeric(),
                                p_value = numeric(),
                                stringsAsFactors = FALSE)

for (i in vat_14_metabolites) {
  # Fit Cox proportional hazards model
  formula_str <- paste0("Surv(OS_days, OS_event)  ~ `",i,"`+ log(LDH) + CRP")
  model_cox <- coxph(as.formula(formula_str), data = vat_features_14_clinic)
  
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
  vat_14_os_cox_MVA  <- rbind(vat_14_os_cox_MVA, marker_results_cox)
}

vat_14_os_cox_MVA |>
  # filter(marker != "LPI-(16:01)")|>
  ggplot()+
  geom_point(aes(x = reorder(marker, HR) , y = HR),
             size = 4, shape = 19, color = "darkgreen", alpha = 0.7)+
  geom_errorbar(aes(x = marker, ymin = lower95, 
                    ymax = higher95), width = 0.1,
                color = "darkgreen")+
  geom_text(aes(x = marker, y = 6.5, 
                label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

## Combining PFS Hazard Ratio tables from different days and dynamic evaluation

vat_0_pfs_cox_MVA <- vat_0_pfs_cox_MVA |>
  mutate(time = 0) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))

vat_3_pfs_cox_MVA <- vat_3_pfs_cox_MVA |>
  mutate(time =3) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))

vat_14_pfs_cox_MVA <- vat_14_pfs_cox_MVA |>
  mutate(time =14)|>
  mutate(FDR = p.adjust(p_value, method = "BH"))

vat_dyn_pfs_cox <- vat_dyn_pfs_cox |>
  mutate(time ="dyn")|>
  mutate(FDR = p.adjust(p_value, method = "BH"))

vat_pfs_cox_combined <- rbind(vat_0_pfs_cox_MVA, vat_3_pfs_cox_MVA, vat_14_pfs_cox_MVA, vat_dyn_pfs_cox)

vat_pfs_cox_combined_0.05 <- vat_pfs_cox_combined |>
  filter(p_value < 0.05)|>
  mutate(time = as.factor(time))

#p_pfs_cox_combined

p_vat_pfs_cox_comibned_0.05 <- vat_pfs_cox_combined_0.05  |>
  ggplot()+
  geom_jitter(aes(x = reorder(marker, HR) , y = HR, color = factor(time, levels = c("0", "3", "14", "dyn"))),
              size = 4, shape = 19, alpha = 0.6, width = 0.2)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = higher95))+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #ylim = c(-1, 4)
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  labs(color = "Timepoint")+
  scale_color_manual(values = c("darkgrey", "darkblue", "orange", "pink"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

p_vat_pfs_cox_comibned_0.05
## Combining PFS Hazard Ratio tables from different days and dynamic evaluation

vat_0_os_cox_MVA <- vat_0_os_cox_MVA |>
  mutate(time = 0) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))

vat_3_os_cox_MVA <- vat_3_os_cox_MVA |>
  mutate(time =3) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))

vat_14_os_cox_MVA <- vat_14_os_cox_MVA |>
  mutate(time =14) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))

vat_dyn_os_cox <- vat_dyn_os_cox |>
  mutate(time ="dyn")|>
  mutate(FDR = p.adjust(p_value, method = "BH"))

vat_os_cox_combined <- rbind(vat_0_os_cox_MVA, vat_3_os_cox_MVA, vat_14_os_cox_MVA, vat_dyn_os_cox)

vat_os_cox_combined_0.05 <- vat_os_cox_combined |>
  filter(FDR < 0.1)|>
  mutate(time = as.factor(time))

#p_os_cox_combined

p_vat_os_cox_combined_0.05 <- vat_os_cox_combined_0.05  |>
  filter(marker != 'LPI-(16:01)')|>
  ggplot()+
  geom_jitter(aes(x = reorder(marker, HR) , y = HR, color = factor(time, levels = c("0", "3", "14", "dyn"))),
              size = 4, shape = 19, alpha = 0.6, width = 0.2, height = 0.05)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = higher95))+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #ylim = c(-1, 4)
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  labs(color = "Timepoint")+
  scale_color_manual(values = c("darkgrey", "darkblue", "pink"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

p_vat_os_cox_combined_0.05

## Visualization of top metabolites in Kaplan-Meier survival curves ----
## OS

vat_os_cox_3 <- vat_os_cox_combined_0.05 |>
  filter(time == 3) |>
  select(marker)|>
  unlist()|>
  as.vector()

p_km_vat_os_cox_3 <- list()

for (metabolite in vat_os_cox_3) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(OS_days, OS_event) ~ ifelse(get(metabolite, vat_features_3_clinic) > mean(get(metabolite, vat_features_3_clinic)), 1, 0), data = vat_features_3_clinic)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  xlim=c(0, 1110), break.x.by = 180, ylab="OS probability", xlab="Days after CAR-T infusion",
                  pval= TRUE, pval.coord = c(10, 0.08), pval.size = 3,
                  size = 1.15,
                  axes.offset = FALSE,
                  risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                  tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                  conf.int = FALSE,
                  ggtheme = theme_classic2(10),
                  font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                  fontsize=3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste(metabolite)),
                  palette = c("lightblue", "darkblue"))
  
  
  # Add the plot to the list
  p_km_vat_os_cox_3[[metabolite]] <- p
}

grobs_vat_os_cox_3 <- lapply(p_km_vat_os_cox_3, function(x) ggplotGrob(x$plot))
p_grobs_vat_os_cox_3 <- do.call(grid.arrange, c(grobs_vat_os_cox_3, ncol = 3))

p_km_vat_os_cox_3$`LPC-(16:01)`
p_km_vat_os_cox_3$`LPI-(22:05)`

vat_os_cox_14 <- vat_os_cox_combined_0.05 |>
  filter(time == 14) |>
  select(marker)|>
  unlist()|>
  as.vector()

# p_km_vat_os_cox_14 <- list()
# 
# for (metabolite in vat_os_cox_14) {
#   # Directly use the formula in survfit2
#   fit <- survfit2(Surv(OS_days, OS_event) ~ ifelse(get(metabolite, vat_features_14_clinic) > mean(get(metabolite, vat_features_14_clinic)), 1, 0), data = vat_features_14_clinic)
#   
#   # Generate the plot
#   p <- ggsurvplot(fit,
#                   xlim=c(0, 1110), break.x.by = 180, ylab="OS probability", xlab="Days after CAR-T infusion",
#                   pval= TRUE, pval.coord = c(10, 0.08), pval.size = 3,
#                   size = 1.15,
#                   axes.offset = FALSE,
#                   risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
#                   tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
#                   conf.int = FALSE,
#                   ggtheme = theme_classic2(10),
#                   font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
#                   fontsize=3,
#                   legend.labs = c("Low", "High"),
#                   legend.title = c(paste(metabolite)),
#                   palette = c("orange", "orange3"))
#   
#   
#   # Add the plot to the list
#   p_km_vat_os_cox_14[[metabolite]] <- p
# }
# 
# grobs_vat_os_cox_14 <- lapply(p_km_vat_os_cox_14, function(x) ggplotGrob(x$plot))
# p_grobs_vat_os_cox_14 <- do.call(grid.arrange, c(grobs_vat_os_cox_14, ncol = 5))

## Visualization of top metabolites in Kaplan-Meier survival curves
## PFS analysis
# Function to create a Kaplan-Meier plot for a given metabolite
# Initialize an empty list to store plots
pfs_cox_0 <- vat_pfs_cox_combined_0.05 |>
  filter(time == 0) |>
  select(marker)|>
  unlist()|>
  as.vector()

p_km_pfs_cox_0 <- list()

for (metabolite in pfs_cox_0) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(PFS_days, PFS_event) ~ ifelse(get(metabolite, vat_features_0_clinic) > mean(get(metabolite, vat_features_0_clinic)), 1, 0), data = vat_features_0_clinic)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  xlim=c(0, 1110), break.x.by = 180, ylab="PFS probability", xlab="Days after CAR-T infusion",
                  pval= TRUE, pval.coord = c(10, 0.08), pval.size = 3,
                  size = 1.15,
                  axes.offset = FALSE,
                  risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                  tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                  conf.int = FALSE,
                  ggtheme = theme_classic2(10),
                  font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                  fontsize=3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste(metabolite)),
                  palette = c("lightgrey", "darkgrey"))
  
  
  # Add the plot to the list
  p_km_pfs_cox_0[[metabolite]] <- p
}

grobs_pfs_cox_0 <- lapply(p_km_pfs_cox_0, function(x) ggplotGrob(x$plot))
p_grobs_pfs_cox_0 <- do.call(grid.arrange, c(grobs_pfs_cox_0))

pfs_cox_3 <- vat_pfs_cox_combined_0.05 |>
  filter(time == 3) |>
  select(marker)|>
  unlist()|>
  as.vector()

p_km_pfs_cox_3 <- list()

for (metabolite in pfs_cox_3) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(PFS_days, PFS_event) ~ ifelse(get(metabolite, vat_features_3_clinic) > mean(get(metabolite, vat_features_3_clinic)), 1, 0), data = vat_features_3_clinic)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  xlim=c(0, 1110), break.x.by = 180, ylab="PFS probability", xlab="Days after CAR-T infusion",
                  pval= TRUE, pval.coord = c(10, 0.08), pval.size = 3,
                  size = 1.15,
                  axes.offset = FALSE,
                  risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                  tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                  conf.int = FALSE,
                  ggtheme = theme_classic2(10),
                  font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                  fontsize=3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste(metabolite)),
                  palette = c("lightblue", "darkblue"))
  
  
  # Add the plot to the list
  p_km_pfs_cox_3[[metabolite]] <- p
}

grobs_pfs_cox_3 <- lapply(p_km_pfs_cox_3, function(x) ggplotGrob(x$plot))
p_grobs_pfs_cox_3 <- do.call(grid.arrange, c(grobs_pfs_cox_3, ncol = 5))

pfs_cox_14 <- vat_pfs_cox_combined_0.05 |>
  filter(time == 14) |>
  select(marker)|>
  unlist()|>
  as.vector()

p_km_pfs_cox_14 <- list()

for (metabolite in pfs_cox_14) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(PFS_days, PFS_event) ~ ifelse(get(metabolite, vat_features_14_clinic) > mean(get(metabolite, vat_features_14_clinic)), 1, 0), data = vat_features_14_clinic)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  xlim=c(0, 1110), break.x.by = 180, ylab="PFS probability", xlab="Days after CAR-T infusion",
                  pval= TRUE, pval.coord = c(10, 0.08), pval.size = 3,
                  size = 1.15,
                  axes.offset = FALSE,
                  risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                  tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                  conf.int = FALSE,
                  ggtheme = theme_classic2(10),
                  font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                  fontsize=3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste(metabolite)),
                  palette = c("orange", "orange3"))
  
  
  # Add the plot to the list
  p_km_pfs_cox_14[[metabolite]] <- p
}

grobs_pfs_cox_14 <- lapply(p_km_pfs_cox_14, function(x) ggplotGrob(x$plot))
p_grobs_pfs_cox_14 <- do.call(grid.arrange, c(grobs_pfs_cox_14, ncol = 4))

p_km_pfs_cox_14$`CE-(22:06) NH4`



### MSEA analysis ----
### Enrichment analysis of significant metabolites from log and cox regression models

#Preparation of data sets for MSEA
metabolites_HMDB <- read_xlsx("Input_files/Werner_HMDB_translation.xlsx", na = "NA")

metabolites_HMDB <- metabolites_HMDB %>%
  rename(metabolite = `...1`)

#crs_cox beinhaltet alles von crs_glm
vat_crs_cox_HMDB <- crs_vat_cox_combined_0.05 |>
  select(marker, HR)

vat_crs_cox_HMDB <- left_join(vat_crs_cox_HMDB, metabolites_HMDB, by = c("marker" = "metabolite"))

vat_crs_cox_b1_HMDB <- vat_crs_cox_HMDB |>
  filter(HR<1)|>
  select(HMDB)|>
  unlist()|>
  as.vector()

vat_crs_cox_o1_HMDB <- vat_crs_cox_HMDB |>
  filter(HR>1)|>
  select(HMDB)|>
  unlist()|>
  as.vector()

##pfs does not include all of os
vat_pfs_cox_HMDB <- vat_pfs_cox_combined_0.05 |>
  select(marker, HR)

vat_os_cox_HMDB <- vat_os_cox_combined_0.05 |>
  select(marker, HR)

vat_survival_cox_HMDB <- rbind(vat_pfs_cox_HMDB, vat_os_cox_HMDB)

vat_survival_cox_HMDB <- vat_survival_cox_HMDB |>
  distinct(marker, .keep_all = T)

vat_survival_cox_HMDB <- left_join(vat_survival_cox_HMDB, metabolites_HMDB, by = c("marker" = "metabolite"))

vat_survival_cox_b1_HMDB <- vat_survival_cox_HMDB |>
  filter(HR<1)|>
  select(HMDB)|>
  unlist()|>
  as.vector()

vat_survival_cox_o1_HMDB <- vat_survival_cox_HMDB |>
  filter(HR>1)|>
  select(HMDB)|>
  unlist()|>
  as.vector()

## Enrichment analysis of CRS
#Enrichment of HR<1
mset_vat_crs_cox_b1_msea<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("HMDB0010379","HMDB0010393","HMDB0010397","HMDB0010392","HMDB0010387","HMDB0010406","HMDB0011514","HMDB0011512")
mset_vat_crs_cox_b1_msea<-Setup.MapData(mset_vat_crs_cox_b1_msea, cmpd.vec);
mset_vat_crs_cox_b1_msea<-CrossReferencing(mset_vat_crs_cox_b1_msea, "hmdb", lipid = T);
mset_vat_crs_cox_b1_msea<-CreateMappingResultTable(mset_vat_crs_cox_b1_msea)
mset_vat_crs_cox_b1_msea<-SetMetabolomeFilter(mset_vat_crs_cox_b1_msea, F);
mset_vat_crs_cox_b1_msea<-SetCurrentMsetLib(mset_vat_crs_cox_b1_msea, "sub_class", 2);
mset_vat_crs_cox_b1_msea<-CalculateHyperScore(mset_vat_crs_cox_b1_msea)

vat_crs_cox_b1_msea <- as.data.frame(mset_vat_crs_cox_b1_msea[["analSet"]][["ora.mat"]])
vat_crs_cox_b1_msea <- vat_crs_cox_b1_msea %>%
  mutate(pathway = rownames(vat_crs_cox_b1_msea), ratio = hits/expected) %>%
  rename("Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(vat_crs_cox_b1_msea)

vat_crs_cox_b1_msea$"Raw p" <- colnames("Raw_p")


#MSEA of HR>1 metabolites
mset_vat_crs_cox_o1_msea<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-vat_crs_cox_o1_HMDB
mset_vat_crs_cox_o1_msea<-Setup.MapData(mset_vat_crs_cox_o1_msea, cmpd.vec);
mset_vat_crs_cox_o1_msea<-CrossReferencing(mset_vat_crs_cox_o1_msea, "hmdb", lipid = T);
mset_vat_crs_cox_o1_msea<-CreateMappingResultTable(mset_vat_crs_cox_o1_msea)
mset_vat_crs_cox_o1_msea<-SetMetabolomeFilter(mset_vat_crs_cox_o1_msea, F);
mset_vat_crs_cox_o1_msea<-SetCurrentMsetLib(mset_vat_crs_cox_o1_msea, "sub_class", 2);
mset_vat_crs_cox_o1_msea<-CalculateHyperScore(mset_vat_crs_cox_o1_msea)

vat_crs_cox_o1_msea <- as.data.frame(mset_vat_crs_cox_o1_msea[["analSet"]][["ora.mat"]])
vat_crs_cox_o1_msea <- vat_crs_cox_o1_msea %>%
  mutate(pathway = rownames(vat_crs_cox_o1_msea), ratio = hits/expected) %>%
  rename("Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(vat_crs_cox_o1_msea)

vat_crs_cox_o1_msea$"Raw p" <- colnames("Raw_p")

vat_crs_cox_b1_msea <- vat_crs_cox_b1_msea |>
  mutate(direction = "HR<1")
vat_crs_cox_o1_msea <- vat_crs_cox_o1_msea |>
  mutate(direction = "HR>1")

vat_crs_cox_msea_combined <- rbind(vat_crs_cox_b1_msea, vat_crs_cox_o1_msea)

p_vat_crs_cox_msea_combined <- vat_crs_cox_msea_combined |>
  filter(Raw_p < 0.01) |>
  # filter(hits >= 1)|>
  ggplot()+
  geom_point(aes(x=direction, y=pathway, size=ratio, color=Raw_p))+
  scale_color_gradient(low = "darkorange", high = "darkblue") +
  #scale_size(limits = c(1,100), breaks = c(1,10,100))+
  #scale_x_discrete(limits = c("123", "121", "321"), labels = c("Increase \n pattern", "In-Decrease \n pattern", "Decrease \n pattern"))+
  scale_y_discrete(limits = rev(levels(vat_crs_cox_msea_combined))) +
  labs(y="", x="", size = "Enrichment \n score", color = "P-value")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a black frame
  ) 

p_vat_crs_cox_msea_combined

## Enrichment analysis of PFS
#Enrichment of HR<1
mset_vat_survival_cox_b1_msea<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-vat_survival_cox_b1_HMDB
mset_vat_survival_cox_b1_msea<-Setup.MapData(mset_vat_survival_cox_b1_msea, cmpd.vec);
mset_vat_survival_cox_b1_msea<-CrossReferencing(mset_vat_survival_cox_b1_msea, "hmdb", lipid = T);
mset_vat_survival_cox_b1_msea<-CreateMappingResultTable(mset_vat_survival_cox_b1_msea)
mset_vat_survival_cox_b1_msea<-SetMetabolomeFilter(mset_vat_survival_cox_b1_msea, F);
mset_vat_survival_cox_b1_msea<-SetCurrentMsetLib(mset_vat_survival_cox_b1_msea, "sub_class", 2);
mset_vat_survival_cox_b1_msea<-CalculateHyperScore(mset_vat_survival_cox_b1_msea)

vat_survival_cox_b1_msea <- as.data.frame(mset_vat_survival_cox_b1_msea[["analSet"]][["ora.mat"]])
vat_survival_cox_b1_msea <- vat_survival_cox_b1_msea %>%
  mutate(pathway = rownames(vat_survival_cox_b1_msea), ratio = hits/expected) %>%
  rename("Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(vat_survival_cox_b1_msea)

vat_survival_cox_b1_msea$"Raw p" <- colnames("Raw_p")


#MSEA of HR>1 metabolites
mset_vat_survival_cox_o1_msea<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-vat_survival_cox_o1_HMDB
mset_vat_survival_cox_o1_msea<-Setup.MapData(mset_vat_survival_cox_o1_msea, cmpd.vec);
mset_vat_survival_cox_o1_msea<-CrossReferencing(mset_vat_survival_cox_o1_msea, "hmdb", lipid = T);
mset_vat_survival_cox_o1_msea<-CreateMappingResultTable(mset_vat_survival_cox_o1_msea)
mset_vat_survival_cox_o1_msea<-SetMetabolomeFilter(mset_vat_survival_cox_o1_msea, F);
mset_vat_survival_cox_o1_msea<-SetCurrentMsetLib(mset_vat_survival_cox_o1_msea, "sub_class", 2);
mset_vat_survival_cox_o1_msea<-CalculateHyperScore(mset_vat_survival_cox_o1_msea)

vat_survival_cox_o1_msea <- as.data.frame(mset_vat_survival_cox_o1_msea[["analSet"]][["ora.mat"]])
vat_survival_cox_o1_msea <- vat_survival_cox_o1_msea %>%
  mutate(pathway = rownames(vat_survival_cox_o1_msea), ratio = hits/expected) %>%
  rename("Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(vat_survival_cox_o1_msea)

vat_survival_cox_o1_msea$"Raw p" <- colnames("Raw_p")

vat_survival_cox_b1_msea <- vat_survival_cox_b1_msea |>
  mutate(direction = "HR<1")
vat_survival_cox_o1_msea <- vat_survival_cox_o1_msea |>
  mutate(direction = "HR>1")

vat_survival_cox_msea_combined <- rbind(vat_survival_cox_b1_msea, vat_survival_cox_o1_msea)

p_vat_survival_cox_msea_combined <- vat_survival_cox_msea_combined |>
  # filter(FDR < 0.2) |>
  # filter(hits >= 1)|>
  ggplot()+
  geom_point(aes(x=direction, y=pathway, size=ratio, color=Raw_p))+
  scale_color_gradient(low = "darkorange", high = "darkblue") +
  #scale_size(limits = c(1,100), breaks = c(1,10,100))+
  #scale_x_discrete(limits = c("123", "121", "321"), labels = c("Increase \n pattern", "In-Decrease \n pattern", "Decrease \n pattern"))+
  scale_y_discrete(limits = rev(levels(vat_survival_cox_msea_combined))) +
  labs(y="", x="", size = "Enrichment \n score", color = "P-value")+
  theme_classic()+  # Chopfse a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a black frame
  ) 

p_vat_survival_cox_msea_combined



### Correlation analysis with Olink data from V. Blumenberg ----
## Loading of Olink data and linking with baseline data

olink <- read_xlsx("Input_files/olink.xlsx", na = "NA")

vat_0_olink <- left_join(cart_vat_0_norm, olink)

vat_0_olink_master <- vat_0_olink|>
  pivot_longer(cols = c(ADA:TNFRSF21), names_to = "cytokines", values_to = "level")|>
  pivot_longer(cols = c(Alanine:'Cer-(24:01)'), names_to = "metabolite", values_to = "abundance")

cytokine <- vat_0_olink|>
  select(ADA:TNFRSF21) |>
  colnames()|>
  unlist()|>
  as.vector()

crs_vat_markers <- crs_vat_cox_combined_0.05 |>
  select(marker) |>
  unlist() |>
  as.vector()

vat0_crs_olink_corr <- data.frame()
for (i in crs_vat_markers) {
  for (k in cytokine) {
    # Extract the relevant columns from the data frame
    subset_data <- vat_0_olink[, c(k, i)]
    
    # Remove NA values
    subset_data <- na.omit(subset_data)
    
    # Calculate Spearman correlation
    correlation_result <- cor.test(subset_data[[k]], subset_data[[i]], method = "pearson")
    
    # Create a data frame with the results
    result_row <- data.frame(
      Metabolite = i,
      Cytokine = k,
      Correlation = correlation_result$estimate,
      P_Value = correlation_result$p.value
    )
    
    # Append the results to the main data frame
    vat0_crs_olink_corr <- rbind(vat0_crs_olink_corr, result_row)
  }
}

vat0_crs_olink_corr 

vat0_crs_olink_corr  |>
  ggplot(aes(x=Metabolite, y=Cytokine, color=Correlation))+
  geom_point(size = 5)+
  geom_text(aes(label = ifelse(P_Value < 0.05, "*", "")), size = 6, color = "black",
            nudge_x = 0.015, nudge_y = -0.06)+
  scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue" )+
  scale_y_discrete(limits = rev(levels(vat0_crs_olink_corr$Cytokine)))+
  labs(x="",y="")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black frame
    axis.text.x = element_text(angle = 45, hjust=0.9))


## Correlation between day 3 metabolites and 
vat_3_olink <- left_join(cart_vat_3_norm, olink)

vat_0_olink_master <- vat_0_olink|>
  pivot_longer(cols = c(ADA:TNFRSF21), names_to = "cytokines", values_to = "level")|>
  pivot_longer(cols = c(Alanine:'Cer-(24:01)'), names_to = "metabolite", values_to = "abundance")

cytokine <- vat_0_olink|>
  select(ADA:TNFRSF21) |>
  colnames()|>
  unlist()|>
  as.vector()

survival_vat_markers <- vat_survival_cox_HMDB |>
  select(marker) |>
  unlist() |>
  as.vector()

vat0_survival_olink_corr <- data.frame()
for (i in survival_vat_markers) {
  for (k in cytokine) {
    # Extract the relevant columns from the data frame
    subset_data <- vat_0_olink[, c(k, i)]
    
    # Remove NA values
    subset_data <- na.omit(subset_data)
    
    # Calculate Spearman correlation
    correlation_result <- cor.test(subset_data[[k]], subset_data[[i]], method = "pearson")
    
    # Create a data frame with the results
    result_row <- data.frame(
      Metabolite = i,
      Cytokine = k,
      Correlation = correlation_result$estimate,
      P_Value = correlation_result$p.value
    )
    
    # Append the results to the main data frame
    vat0_survival_olink_corr <- rbind(vat0_survival_olink_corr, result_row)
  }
}

vat0_survival_olink_corr 

vat0_survival_olink_corr   |>
  ggplot(aes(x=Metabolite, y=Cytokine, color=Correlation))+
  geom_point(size = 5)+
  geom_text(aes(label = ifelse(P_Value < 0.05, "*", "")), size = 6, color = "black",
            nudge_x = 0.015, nudge_y = -0.06)+
  scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue" )+
  scale_y_discrete(limits = rev(levels(vat0_survival_olink_corr$Cytokine)))+
  labs(x="",y="")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black frame
    axis.text.x = element_text(angle = 45, hjust=0.9))





### Analysis to which time pattern the selected features belong ----
## Left_join of CRS relevant metabolites and time pattern with keeping all on the left side

vat_crs_all <- rbind(vat_crs_cox_HMDB, vat_survival_cox_HMDB)

vat_crs_all <- vat_crs_all |>
  distinct(marker)

vat_crs_all <- left_join(vat_crs_all, cart_time_pattern, by = join_by("marker" == "metabolite"))

vat_crs_all$direction <- as.character(vat_crs_all$direction)

vat_crs_all |>
  count(direction) |>
  ggplot(aes(x = "", y = n, fill = direction)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(fill = "Category")


### Export of figures ----
# Save plots

save_plot <- function(plot_list, prefix = "plot", height = 7, width = 7) {
  for (plot_name in names(plot_list)) {
    filename <- paste0(prefix, plot_name, ".svg")
    svg(filename = filename, height = height, width = width)
    print(plot_list[[plot_name]])
    dev.off()
  }
}

plot_pattern_PCA <- "plot_[a-zA-Z]+_[0-9]+_PCA"
plot_pattern_PLSDA <- "plot_[a-zA-Z]+_[0-9]+_PLSDA"
plot_pattern_PLSDA_VIP <- "plot_[a-zA-Z]+_[0-9]+_PLSDA_VIP"
#plot_pattern_sig_metabolite <- "plot_[a-zA-Z]+_sig_metabolite"
#plot_pattern_sig_metabolite_adj <- "plot_[a-zA-Z]+_sig_metabolite_adj"
plot_pattern_volc <- "plot_volc_vat_[0-9]+_[a-zA-Z]"
plot_pattern_CI <- "p_CI_[a-zA-Z]+_[a-zA-Z]+_[a-zA-Z]"
plot_pattern_glm <- "p_[a-zA-Z]+_[a-zA-Z]+_glm"
plot_pattern_cox <- "p_[a-zA-Z]+_[a-zA-Z]+_cox"
plot_pattern_venn <- "p_venn_[a-zA-Z]+_[0-9]"
plot_pattern_corr <- "p_crs_[a-zA-Z]+_corr_[a-zA-Z]"
plot_pattern_msea <- "p_crs_msea_combined"

plots_PCA <- mget(ls(pattern = plot_pattern_PCA))
plots_PLSDA <- mget(ls(pattern = plot_pattern_PLSDA))
plots_PLSDA_VIP <- mget(ls(pattern = plot_pattern_PLSDA_VIP))
#plots_sig_metabolite <- mget(ls(pattern = plot_pattern_sig_metabolite))
#plots_sig_metabolite_adj <- mget(ls(pattern = plot_pattern_sig_metabolite_adj))
plots_volc <- mget(ls(pattern = plot_pattern_volc))
plots_CI <- mget(ls(pattern = plot_pattern_CI))
plots_glm <- mget(ls(pattern = plot_pattern_glm))
plots_cox <- mget(ls(pattern = plot_pattern_cox))
plots_venn <- mget(ls(pattern = plot_pattern_venn))
plots_corr <- mget(ls(pattern = plot_pattern_corr))
plots_msea <- mget(ls(pattern = plot_pattern_msea))

# plots_sig_metabolite <- list(plot_ABX_sig_metabolite = plot_ABX_sig_metabolite,
#                              plot_Dys_sig_metabolite = plot_Dys_sig_metabolite,
#                              plot_Exp_sig_metabolite = plot_Exp_sig_metabolite)

save_plot(plots_PCA)
save_plot(plots_PLSDA)
save_plot(plots_volc)
save_plot(plots_PLSDA_VIP, height = 8, width = 5)
save_plot(plots_CI)
save_plot(plots_glm, height = 12, width = 6)
save_plot(plots_cox, height = 12, width = 6)
save_plot(plots_venn)
save_plot(plots_corr, height = 3, width = 12)
save_plot(plots_msea, height = 4.5, width = 5)
#save_plot(plots_sig_metabolite, height = 10, width = 15)
#save_plot(plots_sig_metabolite_adj, height = 10, width = 15)

'darkgreen' '' 'red4' 'red2'