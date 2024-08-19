##### Time point analysis of paired samples ----
### Analysis workflow copied from MetaboAnalyst website
### Original 153 samples with 525 features (qc before based on Lamivudin - Werner)


### Loading, filtering and normalization of data ----

rm(mSet)
## Loading/sanity check
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "Input_files/230418_CART_Metabolomics_Data_qc_paired_samples.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet<-FilterVariable(mSet, "median", 0, "F")
mSet<-PreparePrenormData(mSet)

## Normalization by sum and data scaling based on auto-scaling
mSet<-Normalization(mSet, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
# mSet<-PlotNormSummary(mSet, "norm_1_", "png", 72, width=NA)
# mSet<-PlotSampleNormSummary(mSet, "snorm_1_", "png", 72, width=NA)

## Extraction and saving of the normalized data into a new tibble
cart_time_norm <- as_tibble(mSet[["dataSet"]][["norm"]])

## Load original table and cbind sample name and time point label
cart_original <- read_xlsx("Input_files/230418_CART_Metabolomics_Data_qc_paired_samples.xlsx", na = "NA")
cart_time_norm <- cbind(cart_original[,1:2], cart_time_norm)

# ###### ASH Abstract Calculations
# 
# cart_time_norm %>%
#    pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
#    filter(metabolite == "Alanine")
# 
# alanine_ttest <- cart_time_norm %>%
#   pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
#   filter(metabolite == "Alanine", Timepoint != 2)
# 
# alanine_ttest_result <- t.test(level ~ Timepoint, data = alanine_ttest)
# 
# print(alanine_ttest_result)
# 
# dRiboseP_ttest <- cart_time_norm %>%
#   pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
#   filter(metabolite == "dRibose-P", Timepoint != 2)
# 
# dRiboseP_ttest_result <- t.test(level ~ Timepoint, data = dRiboseP_ttest)
# 
# print(dRiboseP_ttest_result)

### Calculation of ANOVA for detection of significantly different metabolites, FDR cutoff 0.05 ----
# 
## Adding ANOVA to object
mSet<-ANOVA.Anal(mSet, F, 0.05, FALSE)
# mSet<-PlotANOVA(mSet, "aov_0_", "png", 72, width=NA)
# 
# 
# cart_time_norm %>%
#   pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
#   filter(metabolite %in% cart_time_123_sig$metabolite) %>%
#   #  filter(metabolite == "Alanine") %>%
#   group_by(Timepoint) %>%
#   ggplot() +
#   geom_violin(aes(x = as.factor(Timepoint), y = level))+
#   geom_jitter(aes(x = as.factor(Timepoint), y = level, alpha= 0.5, color = as.factor(Timepoint)))+
#   guides(alpha = "none", color = "none")+
#   scale_color_manual(values = c("grey", "darkblue", "orange")) +
#   scale_fill_manual(values = c("grey", "darkblue", "orange")) +
#   scale_x_discrete(label = c("0", "3-5", "14")) +
#   xlab("Time after CAR-T transfusion [days]") +
#   ylab("Concentration [norm.]") +
#   facet_wrap(vars(metabolite), scales = "free")+
#   theme_classic()

## Extraction of metabolites and p-values in new tibble
cart_time_ANOVA <- tibble(mSet[["dataSet"]][["prenorm.feat.nms"]])
colnames(cart_time_ANOVA)[1] <- "metabolite"
cart_time_ANOVA$p_value <- mSet[["analSet"]][["aov"]][["p.value"]]

## Sanity check: uniform p-value distribution?
plot_pvalue_distribution <- cart_time_ANOVA %>%
  ggplot() +
  geom_histogram(aes(x=p_value))+
  theme_minimal()

plot_pvalue_distribution

### Pattern analysis with Spearman rank correlation ----

## Adding 1-2-3/3-2-1 pattern to object
mSet<-Match.Pattern(mSet, "pearson", "1-2-3")
#mSet<-PlotCorr(mSet, "ptn_1_", "feature", "png", 72, width=NA)

## Extracting 1-2-3 pattern information and saving top 123 and top 321 correlations
cart_time_123 <- as.data.frame(mSet[["analSet"]][["corr"]][["cor.mat"]])
cart_time_123$metabolite <- rownames(cart_time_123)
cart_time_123 <- cart_time_123 %>% select(metabolite, everything())

cart_time_123_sig <- cart_time_123 %>%
  arrange(FDR) %>%
  filter(FDR < 0.1, correlation > 0.3)

plot_123_sig_correlation <- cart_time_123_sig %>%
  ggplot() +
  geom_point(aes(x = reorder(metabolite, correlation), y = correlation), 
             stat = "identity", size = 3) +
  labs(x = "Metabolite", y = "Pearson Correlation", title = "Increase Pattern")+
  coord_flip() +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

plot_123_sig_correlation

plot_123_all_grid <- cart_time_norm %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% cart_time_123_sig$metabolite) %>%
  #  filter(metabolite == "Alanine") %>%
  group_by(Timepoint) %>%
  ggplot() +
  geom_violin(aes(x = as.factor(Timepoint), y = level))+
  geom_jitter(aes(x = as.factor(Timepoint), y = level, alpha= 0.5, color = as.factor(Timepoint)))+
  guides(alpha = "none", color = "none")+
  scale_color_manual(values = c("grey", "darkblue", "orange")) +
  scale_fill_manual(values = c("grey", "darkblue", "orange")) +
  scale_x_discrete(label = c("0", "3-5", "14")) +
  xlab("Time after CAR-T transfusion [days]") +
  ylab("Concentration [norm.]") +
  facet_wrap(vars(metabolite), scales = "free")+
  theme_classic()

#plot_123_all_grid

cart_time_321_sig <- cart_time_123 %>%
  arrange(FDR) %>%
  filter(FDR < 0.1, correlation < -0.3)

plot_321_sig_correlation <- cart_time_321_sig %>%
  ggplot() +
  geom_point(aes(x = reorder(metabolite, correlation), y = correlation), 
             stat = "identity", size = 3) +
  labs(x = "Metabolite", y = "Pearson Correlation", title = "Decrease Pattern")+
  coord_flip() +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

#plot_321_sig_correlation

plot_321_all_grid <- cart_time_norm %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% cart_time_321_sig$metabolite) %>%
  #  filter(metabolite == "Alanine") %>%
  group_by(Timepoint) %>%
  ggplot() +
  geom_violin(aes(x = as.factor(Timepoint), y = level))+
  geom_jitter(aes(x = as.factor(Timepoint), y = level, alpha= 0.5, color = as.factor(Timepoint)))+
  guides(alpha = "none", color = "none")+
  scale_color_manual(values = c("grey", "darkblue", "orange")) +
  scale_fill_manual(values = c("grey", "darkblue", "orange")) +
  scale_x_discrete(label = c("0", "3-5", "14")) +
  xlab("Time after CAR-T transfusion [days]") +
  ylab("Concentration [norm.]") +
  facet_wrap(vars(metabolite), scales = "free")+
  theme_classic()

#plot_321_all_grid

## Overwriting 1-2-1/2-1-2 pattern to object
mSet<-Match.Pattern(mSet, "pearson", "1-2-1")
# mSet<-PlotCorr(mSet, "ptn_1_", "feature", "png", 72, width=NA)

## Extracting 1-2-1/2-1-2 pattern information and saving top 121 and top 211 correlations
cart_time_121 <- as.data.frame(mSet[["analSet"]][["corr"]][["cor.mat"]])
cart_time_121$metabolite <- rownames(cart_time_121)
cart_time_121 <- cart_time_121 %>% select(metabolite, everything())

cart_time_121_sig <- cart_time_121 %>%
  arrange(FDR) %>%
  filter(FDR < 0.1, correlation > 0.3)

plot_121_sig_correlation <- cart_time_121_sig %>%
  ggplot() +
  geom_point(aes(x = reorder(metabolite, correlation), y = correlation), 
             stat = "identity", size = 3) +
  labs(x = "Metabolite", y = "Pearson Correlation", title = "Increase-Decrease Pattern")+
  coord_flip() +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

#plot_121_sig_correlation

plot_121_all_grid <- cart_time_norm %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% cart_time_121_sig$metabolite) %>%
  #  filter(metabolite == "Alanine") %>%
  group_by(Timepoint) %>%
  ggplot() +
  geom_violin(aes(x = as.factor(Timepoint), y = level))+
  geom_jitter(aes(x = as.factor(Timepoint), y = level, alpha= 0.5, color = as.factor(Timepoint)))+
  guides(alpha = "none", color = "none")+
  scale_color_manual(values = c("grey", "darkblue", "orange")) +
  scale_fill_manual(values = c("grey", "darkblue", "orange")) +
  scale_x_discrete(label = c("0", "3-5", "14")) +
  xlab("Time after CAR-T transfusion [days]") +
  ylab("Concentration [norm.]") +
  facet_wrap(vars(metabolite), scales = "free")+
  theme_classic()

#plot_121_all_grid

cart_time_212_sig <- cart_time_121 %>%
  arrange(FDR) %>%
  filter(FDR < 0.1, correlation < -0.3)

plot_212_sig_correlation <- cart_time_212_sig %>%
  ggplot() +
  geom_point(aes(x = reorder(metabolite, correlation), y = correlation), 
             stat = "identity", size = 3) +
  labs(x = "Metabolite", y = "Pearson Correlation", title = "Decrease-Increase Pattern")+
  coord_flip() +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

#plot_212_sig_correlation

plot_212_all_grid <- cart_time_norm %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% cart_time_212_sig$metabolite) %>%
  #  filter(metabolite == "Alanine") %>%
  group_by(Timepoint) %>%
  ggplot() +
  geom_violin(aes(x = as.factor(Timepoint), y = level))+
  geom_jitter(aes(x = as.factor(Timepoint), y = level, alpha= 0.5, color = as.factor(Timepoint)))+
  guides(alpha = "none", color = "none")+
  scale_color_manual(values = c("grey", "darkblue", "orange")) +
  scale_fill_manual(values = c("grey", "darkblue", "orange")) +
  scale_x_discrete(label = c("0", "3-5", "14")) +
  xlab("Time after CAR-T transfusion [days]") +
  ylab("Concentration [norm.]") +
  facet_wrap(vars(metabolite), scales = "free")+
  theme_classic()

#plot_212_all_grid

## Top significantly patterns shown together

## Problem: Some patterns fall into two categories, so they have to be divided based
## on the higher correlation coefficient. 

cart_time_123_sig <- cart_time_123_sig %>%
  mutate(direction = 123)
cart_time_321_sig <- cart_time_321_sig %>%
  mutate(direction = 321)
cart_time_121_sig <- cart_time_121_sig %>%
  mutate(direction = 121)
cart_time_212_sig <- cart_time_212_sig %>%
  mutate(direction = 212)

cart_time_pattern <- rbind(cart_time_123_sig, cart_time_321_sig, cart_time_121_sig,
                           cart_time_212_sig)

# cart_time_pattern <- cart_time_pattern %>%
#   arrange(desc(abs(correlation))) %>%
#   mutate(duplicate = duplicated(metabolite)) %>%
#   subset(duplicate == FALSE) %>%
#   arrange(direction)

# Removing duplicated metabolites
cart_time_pattern <- cart_time_pattern %>%
  mutate(duplicate = duplicated(metabolite))

cart_time_pattern %>%
  mutate(duplicate = duplicated(metabolite))|>
  filter(duplicate == T)|>
  select(metabolite)|>
  as.vector()

duplicated_metabolites <- cart_time_pattern %>%
  mutate(duplicate = duplicated(metabolite))|>
  filter(duplicate == T)|>
  select(metabolite)|>
  as.vector()

cart_time_pattern <- cart_time_pattern %>%
  group_by(metabolite) %>%
  filter(correlation == max(correlation)) %>%
  ungroup()

# Subsetting unique metabolites following direction for enrichment analysis
cart_time_123_sig_unique <- cart_time_pattern %>%
  filter(direction == 123) %>%
  select(metabolite) %>%
  as.data.frame()

cart_time_321_sig_unique <- cart_time_pattern %>%
  filter(direction == 321) %>%
  select(metabolite) %>%
  as.data.frame()

cart_time_121_sig_unique <- cart_time_pattern %>%
  filter(direction == 121) %>%
  select(metabolite) %>%
  as.data.frame()

cart_time_212_sig_unique <- cart_time_pattern %>%
  filter(direction == 212) %>%
  select(metabolite) %>%
  as.data.frame()

# Visualization of top 4 patterns

top4_123 <- cart_time_pattern %>%
  filter(direction == "123") %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:4) %>%
  select(metabolite) %>%
  as.vector()

top4_321 <- cart_time_pattern %>%
  filter(direction == "321") %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:4) %>%
  select(metabolite) %>%
  as.vector()

top4_121 <- cart_time_pattern %>%
  filter(direction == "121") %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:4) %>%
  select(metabolite) %>%
  as.vector()

top4_212 <- cart_time_pattern %>%
  filter(direction == "212") %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:4) %>%
  select(metabolite) %>%
  as.vector()


plot_top123 <- cart_time_norm %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% top4_123$metabolite) %>%
  group_by(Timepoint) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = as.factor(Timepoint), fill = as.factor(Timepoint)))+
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_pwc(stat = "pwc", method = "t.test", label = "p.signif", 
           bracket.nudge.y = -0.25, step.increase = 0.1,
           tip.length = 0.01) +
  geom_jitter(alpha = 0.3) +
  scale_x_discrete(label = c("0", "3-5", "14")) +
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_color_manual(values = c("grey", "darkblue", "orange")) +
  scale_fill_manual(values = c("grey", "darkblue", "orange")) +
  guides(color = "none", fill = "none") +
  xlab("Time after CAR-T transfusion [days]") +
  ylab("Concentration [norm.]") +
  facet_grid(cols = vars(metabolite), scales = "free") +
  theme_classic()

plot_top123

plot_top123_adj <- ggadjust_pvalue(plot_top123, p.adjust.method = "BH", label = "p.adj.signif")

plot_top123_adj

plot_top321 <- cart_time_norm %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% top4_321$metabolite) %>%
  group_by(Timepoint) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = as.factor(Timepoint), fill = as.factor(Timepoint)))+
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_pwc(stat = "pwc", method = "t.test", label = "p.signif", 
           bracket.nudge.y = -0.25, step.increase = 0.1,
           tip.length = 0.01) +
  geom_jitter(alpha = 0.3, ) +
  scale_x_discrete(label = c("0", "3-5", "14")) +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_color_manual(values = c("grey", "darkblue", "orange")) +
  scale_fill_manual(values = c("grey", "darkblue", "orange")) +
  guides(color = "none", fill = "none") +
  xlab("Time after CAR-T transfusion [days]") +
  ylab("Concentration [norm.]") +
  facet_grid(cols = vars(metabolite), scales = "free") +
  theme_classic()

plot_top321

plot_top321_adj <- ggadjust_pvalue(plot_top321, p.adjust.method = "BH", label = "p.adj.signif")

plot_top321_adj

plot_top121 <- cart_time_norm %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% top4_121$metabolite) %>%
  group_by(Timepoint) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = as.factor(Timepoint), fill = as.factor(Timepoint)))+
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_pwc(stat = "pwc", method = "t.test", label = "p.signif", 
           bracket.nudge.y = -0.25, step.increase = 0.1,
           tip.length = 0.01) +
  geom_jitter(alpha = 0.3, ) +
  scale_x_discrete(label = c("0", "3-5", "14")) +
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_color_manual(values = c("grey", "darkblue", "orange")) +
  scale_fill_manual(values = c("grey", "darkblue", "orange")) +
  guides(color = "none", fill = "none") +
  xlab("Time after CAR-T transfusion [days]") +
  ylab("Concentration [norm.]") +
  facet_grid(cols = vars(metabolite), scales = "free") +
  theme_classic()

plot_top121_adj <- ggadjust_pvalue(plot_top121, p.adjust.method = "BH", label = "p.adj.signif")

plot_top121_adj

plot_top212 <- cart_time_norm %>%
  pivot_longer(cols = Alanine:`Cer-(24:01)`, names_to = "metabolite", values_to = "level") %>%
  filter(metabolite %in% top4_212$metabolite) %>%
  group_by(Timepoint) %>%
  ggplot(aes(x = as.factor(Timepoint), y = level, color = as.factor(Timepoint), fill = as.factor(Timepoint)))+
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_pwc(stat = "pwc", method = "t.test", label = "p.signif", 
           bracket.nudge.y = -0.25, step.increase = 0.1,
           tip.length = 0.01) +
  geom_jitter(alpha = 0.3, ) +
  scale_x_discrete(label = c("0", "3-5", "14")) +
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_color_manual(values = c("grey", "darkblue", "orange")) +
  scale_fill_manual(values = c("grey", "darkblue", "orange")) +
  guides(color = "none", fill = "none") +
  xlab("Time after CAR-T transfusion [days]") +
  ylab("Concentration [norm.]") +
  facet_grid(cols = vars(metabolite), scales = "free") +
  theme_classic()

plot_top212_adj <- ggadjust_pvalue(plot_top212, p.adjust.method = "BH", label = "p.adj.signif")

plot_top212_adj

#grid.arrange(plot_top123, plot_top321, plot_top121, plot_top212,
#            ncol = 1, nrow = 4)

#grid.arrange(plot_top123_adj, plot_top321_adj, plot_top121_adj, plot_top212_adj,
#             ncol = 1, nrow = 4)


### Principal component analysis (PCA) ----
mSet<-PCA.Anal(mSet)

## Extraction of PCA component values
cart_time_PCA <- as.data.frame(mSet[["analSet"]][["pca"]][["x"]])
cart_time_PCA <- cbind(cart_time_norm[,1:2], cart_time_PCA)
cart_time_PCA <- cart_time_PCA %>%
  mutate(timepoint = as.character(Timepoint)) %>%
  select(-Timepoint) %>%
  select(Sample, timepoint, everything())

## Visualization of PCA comp1 vs comp2
plot_time_PCA <- cart_time_PCA %>%
  ggplot(aes(x = PC1, y = PC2, color = timepoint, fill = timepoint)) +
  geom_jitter() +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.9) +
  scale_color_manual(values = c("grey", "darkblue", "orange"))+
  scale_fill_manual(values = c("grey", "darkblue", "orange"), 
                    name = "Days",
                    labels = c("0", "3-5", "14"))+
  guides(color = "none")+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

plot_time_PCA

## 3D Visualization of Comp 1-3
#plot_ly(data = cart_time_PCA, x = ~PC1, y = ~PC2, z = ~PC3, color = ~timepoint, type = "scatter3d", 
#        mode = "markers", alpha = 0.7)

### Partial Least Square Discriminant Analysis (PLS-DA) ----
mSet<-PLSR.Anal(mSet, reg=TRUE)

## Extraction of PLSDA component values
cart_time_PLSDA <- as.matrix.data.frame(PLSDA_timepoints <- mSet[["analSet"]][["plsr"]][["scores"]])
cart_time_PLSDA <- as.data.frame(cart_time_PLSDA)
cart_time_PLSDA <- cbind(cart_time_norm[,1:2], cart_time_PLSDA)
cart_time_PLSDA <- cart_time_PLSDA %>%
  mutate(timepoint = as.character(Timepoint)) %>%
  select(-Timepoint) %>%
  select(Sample, timepoint, everything())

## 2D Visualization
plot_time_PLSDA <- cart_time_PLSDA %>%
  ggplot(aes(x = V1, y = V2, color = timepoint, fill = timepoint)) +
  geom_jitter() +
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.9) +
  ylim(-20, 15)+
  scale_color_manual(values = c("grey", "darkblue", "orange"))+
  scale_fill_manual(values = c("grey", "darkblue", "orange"), 
                    name = "Days",
                    labels = c("0", "3-5", "14"))+
  guides(color = "none")+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1),
        aspect.ratio = 1)

plot_time_PLSDA


## 3D Visualization of Comp 1-3
#plot_ly(data = cart_time_PLSDA, x = ~V1, y = ~V2, z = ~V3, color = ~timepoint, type = "scatter3d", 
#        mode = "markers", alpha = 0.7)

## PLSDA VIP Scores
cart_time_PLSDA_VIP <- as.data.frame(mSet[["analSet"]][["plsr"]][["vip.mat"]])
cart_time_PLSDA_VIP <- cart_time_PLSDA_VIP %>%
  tibble::rownames_to_column(var = "metabolite")

plot_PLSDA_VIP <- cart_time_PLSDA_VIP %>%
  arrange(desc(`Comp. 1`)) %>%
  filter(`Comp. 1` > 2) %>%
  ggplot(aes(x = reorder(metabolite,`Comp. 1`), y = `Comp. 1`, fill = `Comp. 1`)) +
  geom_bar(stat = "identity") +
  guides(fill = "none")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = "", y = "PLS-DA Comp. 1 VIP Score") +
  coord_flip() +
  scale_fill_gradient(low = "lightgrey", high = "darkblue")+
  theme_classic()

plot_PLSDA_VIP

##### Comparison of time points 1 and 2 ----

### Selecting groups in mSet
mSet<-GetGroupNames(mSet, "")
feature.nm.vec <- c("")
smpl.nm.vec <- c("")
grp.nm.vec <- c("1","2")
mSet<-UpdateData(mSet)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet<-FilterVariable(mSet, "median", 0, "F")
mSet<-PreparePrenormData(mSet)

## Normalization by sum and data scaling based on auto-scaling
mSet<-Normalization(mSet, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

## Volcano Plot / T tests
mSet<-Volcano.Anal(mSet, FALSE, 2, 1, F, 0.05, TRUE, "fdr")
#mSet<-PlotVolcano(mSet, "volcano_0_",1, 0, "png", 72, width=NA)

cart_volcano_2vs1 <- as.data.frame(mSet[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_2vs1$metabolite <- rownames(cart_volcano_2vs1)
cart_volcano_2vs1$log_p <- mSet[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_2vs1$log_fc <- mSet[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_2vs1$inx.up <- mSet[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_2vs1$inx.down <- mSet[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_2vs1$inx.p <- mSet[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_2vs1 <- cart_volcano_2vs1 %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

plot_volc_2vs1 <- ggplot(cart_volcano_2vs1, aes(x = log_fc, y = log_p))+
  geom_point(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 20)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))+
  guides(color = "none")+
  #  xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_volc_2vs1

### Comparison of time points 2 and 3 ----

### Selecting groups in mSet
mSet<-GetGroupNames(mSet, "")
feature.nm.vec <- c("")
smpl.nm.vec <- c("")
grp.nm.vec <- c("2","3")
mSet<-UpdateData(mSet)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet<-FilterVariable(mSet, "median", 0, "F")
mSet<-PreparePrenormData(mSet)

## Normalization by sum and data scaling based on auto-scaling
mSet<-Normalization(mSet, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

## Volcano Plot / T tests
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 1, F, 0.05, TRUE, "fdr")
#mSet<-PlotVolcano(mSet, "volcano_0_",1, 0, "png", 72, width=NA)

cart_volcano_3vs2 <- as.data.frame(mSet[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_3vs2$metabolite <- rownames(cart_volcano_2vs1)
cart_volcano_3vs2$log_p <- mSet[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_3vs2$log_fc <- mSet[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_3vs2$inx.up <- mSet[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_3vs2$inx.down <- mSet[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_3vs2$inx.p <- mSet[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_3vs2 <- cart_volcano_3vs2 %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))


plot_volc_3vs2 <- ggplot(cart_volcano_3vs2, aes(x = log_fc, y = log_p))+
  geom_point(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 10)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))+
  guides(color = "none")+
  # xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_volc_3vs2

### Comparison of time points 1 and 3 ----
### Selecting groups in mSet
mSet<-GetGroupNames(mSet, "")
feature.nm.vec <- c("")
smpl.nm.vec <- c("")
grp.nm.vec <- c("1","3")
mSet<-UpdateData(mSet)

## Filter based on median, 15% of features filtered to delete variables with very low values (features = 446)
#mSet<-FilterVariable(mSet, "median", 0, "F")
mSet<-PreparePrenormData(mSet)

## Normalization by sum and data scaling based on auto-scaling
mSet<-Normalization(mSet, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

## Volcano Plot / T tests
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 1, F, 0.05, TRUE, "fdr")
#mSet<-PlotVolcano(mSet, "volcano_0_",1, 0, "png", 72, width=NA)

cart_volcano_3vs1 <- as.data.frame(mSet[["analSet"]][["volcano"]][["fc.log"]])
cart_volcano_3vs1$metabolite <- rownames(cart_volcano_2vs1)
cart_volcano_3vs1$log_p <- mSet[["analSet"]][["volcano"]][["p.log"]]
cart_volcano_3vs1$log_fc <- mSet[["analSet"]][["volcano"]][["fc.log"]]
cart_volcano_3vs1$inx.up <- mSet[["analSet"]][["volcano"]][["inx.up"]]
cart_volcano_3vs1$inx.down <- mSet[["analSet"]][["volcano"]][["inx.down"]]
cart_volcano_3vs1$inx.p <- mSet[["analSet"]][["volcano"]][["inx.p"]]

cart_volcano_3vs1 <- cart_volcano_3vs1 %>%
  mutate(gene_type = case_when(inx.down == T & inx.p == T ~ "down",
                               inx.up == T & inx.p == T ~ "up",
                               inx.p == F ~ "ns"))

plot_volc_3vs1 <- ggplot(cart_volcano_3vs1, aes(x = log_fc, y = log_p))+
  geom_point(aes(color = gene_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(aes(label = ifelse(gene_type != "ns", metabolite, ""), color = gene_type),
                  max.overlaps = 10)+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))+
  guides(color = "none")+
  # xlim(-10, 10)+
  #  coord_fixed()+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


plot_volc_3vs1


#### Enrichment Analysis ----

### Loading of HDMB file
metabolites_HMDB <- read_xlsx("Input_files/Werner_HMDB_translation.xlsx", na = "NA")

metabolites_HMDB <- metabolites_HMDB %>%
  rename(metabolite = `...1`)

### Create a list with significant metabolites from time analysis
cart_time_ANOVA_sig <- as.data.frame(mSet[["analSet"]][["aov"]][["sig.fdr"]])

cart_time_ANOVA_sig<- cart_time_ANOVA_sig %>%
  mutate(metabolite = rownames(cart_time_ANOVA_sig))%>%
  rename(FDR = `mSet[["analSet"]][["aov"]][["sig.fdr"]]`)

cart_time_ANOVA_sig <- cart_time_ANOVA_sig %>%
  filter(FDR < 0.1)

cart_time_ANOVA_sig_ID <- left_join(cart_time_ANOVA_sig, 
                                    metabolites_HMDB, by = "metabolite")

## Vector with HMDB IDs from time analysis without NAs for enrichment analysis
ANOVA_HMDB <-  as.vector(na.omit(cart_time_ANOVA_sig_ID$HMDB))

### Create list with significant 123 pattern metabolites for enrichment
cart_time_123_sig_unique <- left_join(cart_time_123_sig_unique, 
                                      metabolites_HMDB, by = "metabolite")

HMDB_123 <- as.vector(na.omit(cart_time_123_sig_unique$HMDB))

### Create list with significant 321 pattern metabolites for enrichment
cart_time_321_sig_unique <- left_join(cart_time_321_sig_unique, 
                                      metabolites_HMDB, by = "metabolite")

HMDB_321 <- as.vector(na.omit(cart_time_321_sig_unique$HMDB))

### Create list with significant 121 pattern metabolites for enrichment
cart_time_121_sig_unique <- left_join(cart_time_121_sig_unique, 
                                      metabolites_HMDB, by = "metabolite")

HMDB_121 <- as.vector(na.omit(cart_time_121_sig_unique$HMDB))

### Create list with significant 212 pattern metabolites for enrichment
cart_time_212_sig_unique <- left_join(cart_time_212_sig_unique, 
                                      metabolites_HMDB, by = "metabolite")

HMDB_212 <- as.vector(na.omit(cart_time_212_sig_unique$HMDB))

rm(mSet)

mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-ANOVA_HMDB
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)

cart_time_enrichment <- as.data.frame(mSet[["analSet"]][["ora.mat"]])
cart_time_enrichment <- cart_time_enrichment %>%
  mutate(pathway = rownames(cart_time_enrichment), ratio = hits/expected) %>%
  rename("Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(cart_time_enrichment)

cart_time_enrichment$"Raw p" <- colnames("Raw_p")

plot_time_enrichment <- cart_time_enrichment %>%
  filter(FDR < 0.1) %>%
  ggplot(aes(x = reorder(pathway, -FDR), y = FDR, color = as.integer(reorder(pathway, -FDR)))) +
  geom_segment(aes(x = reorder(pathway, -FDR), xend = reorder(pathway, FDR), y = 0.1, yend = FDR),
               alpha = 0.5, linewidth = 1) +
  geom_point(aes(size = ratio))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)), trans = "reverse")+
  #scale_y_reverse()+
  scale_color_gradientn(name = 'category', colours = c('lightgrey', 'darkblue')) +
  labs(x = "", size = "Enrichment Score")+
  coord_flip()+
  guides(color = "none") +
  theme_classic()+
  theme(legend.position = "bottom")

plot_time_enrichment

### Enrichment 123 pattern
rm(mSet)

mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-HMDB_123
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)

cart_time_enrichment_123 <- data.frame()
cart_time_enrichment_123 <- read_xlsx("Input_files/HMDB_123_metabolites.xlsx", na = "NA")
cart_time_enrichment_123 <- cart_time_enrichment_123 %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(cart_time_enrichment_123)

cart_time_enrichment_123$"Raw p" <- colnames("Raw_p")

plot_time_enrichment_123 <- cart_time_enrichment_123 %>%
  filter(FDR < 0.1) %>%
  ggplot(aes(x = reorder(pathway, -FDR), y = FDR, color = as.integer(reorder(pathway, -FDR)))) +
  geom_segment(aes(x = reorder(pathway, -FDR), xend = reorder(pathway, FDR), y = 0.1, yend = FDR),
               alpha = 0.5, linewidth = 1) +
  geom_point(aes(size = ratio))+
  # scale_y_continuous(expand = expansion(mult = c(0, .1)), trans = "reverse")+
  #scale_y_reverse()+
  scale_color_gradientn(name = 'category', colours = c('lightgrey', 'darkblue')) +
  labs(x = "", size = "Enrichment Score")+
  coord_flip()+
  guides(color = "none") +
  theme_classic()+
  theme(legend.position = "bottom")

plot_time_enrichment_123

### Enrichment 321 pattern
rm(mSet)

mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-HMDB_321
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)

cart_time_enrichment_321 <- read_xlsx("Input_files/HMDB_321_metabolites.xlsx", na = "NA")
cart_time_enrichment_321 <- cart_time_enrichment_321 %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(cart_time_enrichment_321)

cart_time_enrichment_321$"Raw p" <- colnames("Raw_p")

plot_time_enrichment_321 <- cart_time_enrichment_321 %>%
  filter(FDR < 0.1) %>%
  ggplot(aes(x = reorder(pathway, -FDR), y = FDR, color = as.integer(reorder(pathway, -FDR)))) +
  geom_segment(aes(x = reorder(pathway, -FDR), xend = reorder(pathway, FDR), y = 0.1, yend = FDR),
               alpha = 0.5, linewidth = 1) +
  geom_point(aes(size = ratio))+
  # scale_y_continuous(expand = expansion(mult = c(0, .1)), trans = "reverse")+
  #scale_y_reverse()+
  scale_color_gradientn(name = 'category', colours = c('lightgrey', 'darkblue')) +
  labs(x = "", size = "Enrichment Score")+
  coord_flip()+
  guides(color = "none") +
  theme_classic()+
  theme(legend.position = "bottom")

plot_time_enrichment_321

### Enrichment 121 pattern
rm(mSet)

mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-HMDB_121
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)

cart_time_enrichment_121 <- read_xlsx("Input_files/HMDB_121_metabolites.xlsx", na = "NA")
cart_time_enrichment_121 <- cart_time_enrichment_121 %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(cart_time_enrichment_121)

cart_time_enrichment_121$"Raw p" <- colnames("Raw_p")

plot_time_enrichment_121 <- cart_time_enrichment_121 %>%
  filter(FDR < 0.1) %>%
  ggplot(aes(x = reorder(pathway, -FDR), y = FDR, color = as.integer(reorder(pathway, -FDR)))) +
  geom_segment(aes(x = reorder(pathway, -FDR), xend = reorder(pathway, FDR), y = 0.1, yend = FDR),
               alpha = 0.5, linewidth = 1) +
  geom_point(aes(size = ratio))+
  # scale_y_continuous(expand = expansion(mult = c(0, .1)), trans = "reverse")+
  #scale_y_reverse()+
  scale_color_gradientn(name = 'category', colours = c('lightgrey', 'darkblue')) +
  labs(x = "", size = "Enrichment Score")+
  coord_flip()+
  guides(color = "none") +
  theme_classic()+
  theme(legend.position = "bottom")

plot_time_enrichment_121

### Enrichment 212 pattern
rm(mSet)

mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-HMDB_212
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)

cart_time_enrichment_212 <- read_xlsx("Input_files/HMDB_212_metabolites.xlsx", na = "NA")
cart_time_enrichment_212 <- cart_time_enrichment_212 %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(cart_time_enrichment_212)

cart_time_enrichment_212$"Raw p" <- colnames("Raw_p")

plot_time_enrichment_212 <- cart_time_enrichment_212 %>%
  filter(FDR < 0.1) %>%
  ggplot(aes(x = reorder(pathway, -FDR), y = FDR, color = as.integer(reorder(pathway, -FDR)))) +
  geom_segment(aes(x = reorder(pathway, -FDR), xend = reorder(pathway, FDR), y = 0.1, yend = FDR),
               alpha = 0.5, linewidth = 1) +
  geom_point(aes(size = ratio))+
  # scale_y_continuous(expand = expansion(mult = c(0, .1)), trans = "reverse")+
  #scale_y_reverse()+
  scale_color_gradientn(name = 'category', colours = c('lightgrey', 'darkblue')) +
  labs(x = "", size = "Enrichment Score")+
  coord_flip()+
  guides(color = "none") +
  theme_classic()+
  theme(legend.position = "bottom")

plot_time_enrichment_212

### Combined visualization of Enrichment pathways
cart_time_enrichment_123 <- cart_time_enrichment_123 |>
  mutate(direction = "123")
cart_time_enrichment_321 <- cart_time_enrichment_321 |>
  mutate(direction = "321")
cart_time_enrichment_121 <- cart_time_enrichment_121 |>
  mutate(direction = "121")
cart_time_enrichment_212 <- cart_time_enrichment_212 |>
  mutate(direction = "212")

cart_time_enrichment_dir <- rbind(cart_time_enrichment_123, cart_time_enrichment_121,
                                  cart_time_enrichment_321, cart_time_enrichment_212)

cart_time_enrichment_dir |>
  filter(FDR < 0.2) |>
  filter(hits >= 2)|>
  ggplot()+
  geom_point(aes(x=direction, y=pathway, size=ratio, color=Raw_p))+
  scale_color_gradient(low = "darkorange", high = "darkblue") +
  #scale_size(limits = c(1,100), breaks = c(1,10,100))+
  scale_x_discrete(limits = c("123", "121", "321"), labels = c("Increase \n pattern", "In-Decrease \n pattern", "Decrease \n pattern"))+
  scale_y_discrete(limits = rev(levels(cart_time_enrichment_dir$pathway))) +
  labs(y="", x="", size = "Enrichment \n score", color = "P-value")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a black frame
  ) 

## Separate enrichment analysis for lipids


cart_time_enrichment_123_lipid <- read_xlsx("Input_files/HMDB_123_lipids.xlsx", na = "NA")
cart_time_enrichment_123_lipid <- cart_time_enrichment_123_lipid %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

cart_time_enrichment_121_lipid <- read_xlsx("Input_files/HMDB_121_lipids.xlsx", na = "NA")
cart_time_enrichment_121_lipid <- cart_time_enrichment_121_lipid %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

cart_time_enrichment_321_lipid <- read_xlsx("Input_files/HMDB_321_lipids.xlsx", na = "NA")
cart_time_enrichment_321_lipid <- cart_time_enrichment_321_lipid %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

cart_time_enrichment_123_lipid <- cart_time_enrichment_123_lipid |>
  mutate(direction = "123")
cart_time_enrichment_321_lipid <- cart_time_enrichment_321_lipid |>
  mutate(direction = "321")
cart_time_enrichment_121_lipid <- cart_time_enrichment_121_lipid |>
  mutate(direction = "121")

cart_time_enrichment_dir_lipid <- rbind(cart_time_enrichment_123_lipid, cart_time_enrichment_121_lipid,
                                        cart_time_enrichment_321_lipid)

cart_time_enrichment_dir_lipid |>
  filter(FDR < 0.2) |>
  filter(hits >= 2)|>
  ggplot()+
  geom_point(aes(x=direction, y=pathway, size=ratio, color=Raw_p))+
  scale_color_gradient(low = "darkorange", high = "darkblue") +
  #scale_size(limits = c(1,100), breaks = c(1,10,100))+
  scale_x_discrete(limits = c("123", "121", "321"), labels = c("Increase \n pattern", "In-Decrease \n pattern", "Decrease \n pattern"))+
  scale_y_discrete(limits = rev(levels(cart_time_enrichment_dir$pathway))) +
  labs(y="", x="", size = "Enrichment \n score", color = "P-value")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a black frame
  ) 
