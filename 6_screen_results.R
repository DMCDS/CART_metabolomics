### Collection of all CRS and survival data from single analyses
## CRS GLM
crs_vat_glm_combined <- crs_vat_glm_combined |>
  mutate(bc = "vat")

crs_sat_glm_combined <- crs_sat_glm_combined |>
  mutate(bc = "bc")

crs_pmi_glm_combined <- crs_pmi_glm_combined |>
  mutate(bc = "pmi")

crs_smi_glm_combined <- crs_smi_glm_combined |>
  mutate(bc = "smi")

crs_bc_glm_combined <- rbind(crs_vat_glm_combined,
                             crs_sat_glm_combined,
                             crs_pmi_glm_combined,
                             crs_smi_glm_combined)

crs_bc_glm_combined |>
  filter(FDR < 0.1) |>
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

## CRS COX
crs_vat_cox_combined <- crs_vat_cox_combined |>
  mutate(bc = "vat")

crs_sat_cox_combined <- crs_sat_cox_combined |>
  mutate(bc = "sat")

crs_pmi_cox_combined <- crs_pmi_cox_combined |>
  mutate(bc = "pmi")

crs_smi_cox_combined <- crs_smi_cox_combined |>
  mutate(bc = "smi")

crs_bc_cox_combined <- rbind(crs_vat_cox_combined,
                             crs_sat_cox_combined,
                             crs_pmi_cox_combined,
                             crs_smi_cox_combined)

crs_bc_cox_combined |>
  filter(FDR <0.1) |>
  ggplot()+
  geom_jitter(aes(x = reorder(marker, HR) , y = HR, color = bc),
              size = 4, shape = 19, width = 0.2, alpha = 0.6)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = higher95))+
  # geom_text(aes(x = marker, y = -0.7,
  #               label = paste0("p=", round(p_value,3))))+
  coord_flip(ylim = c(0, 4))+
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  labs(color = "Timepoint")+
  #scale_color_manual(values = c("darkgrey", "darkblue", "orange", "pink"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

## PFS COX
vat_pfs_cox_combined <- vat_pfs_cox_combined |>
  mutate(bc = "vat")

sat_pfs_cox_combined <- sat_pfs_cox_combined |>
  mutate(bc = "sat")

pmi_pfs_cox_combined <- pmi_pfs_cox_combined |>
  mutate(bc = "pmi")

smi_pfs_cox_combined <- smi_pfs_cox_combined |>
  mutate(bc = "smi")

bc_pfs_cox_combined <- rbind(vat_pfs_cox_combined,
                             sat_pfs_cox_combined,
                             pmi_pfs_cox_combined,
                             smi_pfs_cox_combined)
bc_pfs_cox_combined |>
  filter(FDR <= 0.1)|>
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

bc_pfs_cox_combined |>
  filter(FDR <= 0.1)|>
  ggplot()+
  geom_jitter(aes(x = reorder(marker, HR) , y = HR, color = bc),
              size = 4, shape = 19, alpha = 0.6, width = 0.2)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = higher95))+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #ylim = c(-1, 4)
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  labs(color = "Timepoint")+
  #scale_color_manual(values = c("darkgrey", "darkblue", "orange", "pink"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

## OS COX
vat_os_cox_combined <- vat_os_cox_combined |>
  mutate(bc = "vat")

sat_os_cox_combined <- sat_os_cox_combined |>
  mutate(bc = "sat")

pmi_os_cox_combined <- pmi_os_cox_combined |>
  mutate(bc = "pmi")

smi_os_cox_combined <- smi_os_cox_combined |>
  mutate(bc = "smi")

bc_os_cox_combined <- rbind(vat_os_cox_combined,
                            sat_os_cox_combined,
                            pmi_os_cox_combined,
                            smi_os_cox_combined)
bc_os_cox_combined |>
  filter(FDR <= 0.1)|>
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

bc_os_cox_combined |>
  filter(FDR <= 0.1)|>
  ggplot()+
  geom_jitter(aes(x = reorder(marker, HR) , y = HR, color = bc),
              size = 4, shape = 19, alpha = 0.6, width = 0.2)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = higher95))+
  # geom_text(aes(x = marker, y = -0.7, 
  #               label = ifelse(p_value < 0.001, paste("p < 0.001"), paste("p =", round(p_value, 3)))))+
  coord_flip()+ #ylim = c(-1, 4)
  labs(y = "Hazard Ratio (95%CI)", x = "")+
  labs(color = "Timepoint")+
  #scale_color_manual(values = c("darkgrey", "darkblue", "orange", "pink"))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")  # Adjust the size and style as needed
  )

## Overlap between CRS

vat_crs_0.1 <- crs_vat_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

sat_crs_0.1 <- crs_sat_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

pmi_crs_0.1 <- crs_pmi_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

smi_crs_0.1 <- crs_smi_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

venn_crs_0.1 <- list(vat_crs_0.1,
                     sat_crs_0.1,
                     pmi_crs_0.1,
                     smi_crs_0.1)

ggVennDiagram(venn_crs_0.1, label_alpha = 0,
              category.names = c("VAT", "SAT", "PMI", "SMI"),
              label = "count") +
  scale_fill_distiller(palette = "Blues") +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none")



## Overlap between body composition parameters for OS
vat_os_0.1 <- vat_os_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

sat_os_0.1 <- sat_os_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

pmi_os_0.1 <- pmi_os_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

smi_os_0.1 <- smi_os_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

venn_os_0.1 <- list(vat_os_0.1,
                    sat_os_0.1,
                    pmi_os_0.1,
                    smi_os_0.1)

ggVennDiagram(venn_os_0.1, label_alpha = 0,
              category.names = c("VAT", "SAT", "PMI", "SMI"),
              label = "count") +
  scale_fill_distiller(palette = "Blues") +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none")

intersect(vat_os_0.1, smi_os_0.1)
intersect(sat_os_0.1, smi_os_0.1)

## Overlap between metabolites included in CRS COX analysis and Survival COX analysis
bc_crs_cox_0.1 <- crs_bc_cox_combined |>
  filter(FDR <0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

bc_os_0.1 <- bc_os_cox_combined |>
  filter(FDR <0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

bc_pfs_0.1 <- bc_pfs_cox_combined |>
  filter(FDR <0.1)|>
  select(marker) |>
  unique() |>
  unlist() |>
  as.vector()

venn_crs_os <- list(bc_crs_cox_0.1,
                    bc_os_0.1)

ggVennDiagram(venn_crs_os, label_alpha = 0,
              category.names = c("CRS", "OS"),
              label = "count") +
  scale_fill_distiller(palette = "Blues") +
  scale_color_manual(values = c("black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none")

intersect(bc_crs_cox_0.1, bc_os_0.1)

##Overlap between OS and PFS metabolites

venn_os_pfs <- list(bc_os_0.1,
                    bc_pfs_0.1)

ggVennDiagram(venn_os_pfs, label_alpha = 0,
              category.names = c("OS", "PFS"),
              label = "count") +
  scale_fill_distiller(palette = "Blues") +
  scale_color_manual(values = c("black", "black"))+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill = "none")

intersect(bc_crs_cox_0.1, bc_os_0.1)

## Kaplan Meier curves

crs_bc_3_metabolites <- crs_bc_cox_combined |>
  filter(FDR <0.1) |>
  select(marker) |>
  unlist() |>
  as.vector() |>
  sort() |>
  unique()

bc_pfs_3_metabolites <- bc_pfs_cox_combined |>
  filter(time == 3)|>
  filter(FDR <0.1) |>
  select(marker) |>
  unlist() |>
  as.vector() |>
  sort() |>
  unique()

bc_pfs_14_metabolites <- bc_pfs_cox_combined |>
  filter(time == 14)|>
  filter(FDR <0.1) |>
  select(marker) |>
  unlist() |>
  as.vector() |>
  sort() |>
  unique()

bc_os_0_metabolites <- bc_os_cox_combined |>
  filter(time == 0)|>
  filter(FDR <0.1) |>
  select(marker) |>
  unlist() |>
  as.vector() |>
  sort() |>
  unique()

bc_os_3_metabolites <- bc_os_cox_combined |>
  filter(time == 3)|>
  filter(FDR <0.1) |>
  select(marker) |>
  unlist() |>
  as.vector() |>
  sort() |>
  unique()

bc_os_14_metabolites <- bc_os_cox_combined |>
  filter(time == 14)|>
  filter(FDR <0.1) |>
  select(marker) |>
  unlist() |>
  as.vector() |>
  sort() |>
  unique()


crs_bc_features_3 <- cart_pmi_3_norm |>
  pivot_longer(cols = ('Alanine':'Cer-(24:01)'), names_to = "metabolite", values_to = "level")|>
  filter(metabolite %in% crs_bc_3_metabolites) |>
  pivot_wider(names_from = metabolite, values_from = level) |>
  mutate(time = 3) |>
  select(Sample, time, PMI_high, everything())

##combining metabolites of selected features with clinical data
crs_bc_features_3_clinic <- crs_bc_features_3 |>
  separate(Sample, into = c("Sample", "Letter", "Number"), sep = "_")

crs_bc_features_3_clinic <- left_join(crs_bc_features_3_clinic, meta_master, by = c("Number" = "Sample ID Number"))

str(crs_bc_features_3_clinic)


crs_bc_features_3_clinic$CRS_high <- as.double(crs_bc_features_3_clinic$CRS_high)

bc_features_0 <- cart_pmi_0_norm |>
  pivot_longer(cols = ('Alanine':'Cer-(24:01)'), names_to = "metabolite", values_to = "level")|>
  pivot_wider(names_from = metabolite, values_from = level) |>
  mutate(time = 0) |>
  select(Sample, time, PMI_high, everything())

bc_features_0_clinic <- bc_features_0  |>
  separate(Sample, into = c("Sample", "Letter", "Number"), sep = "_")

bc_features_0_clinic<- left_join(bc_features_0_clinic, meta_master, by = c("Number" = "Sample ID Number"))

str(bc_features_0_clinic)

bc_features_3 <- cart_pmi_3_norm |>
  pivot_longer(cols = ('Alanine':'Cer-(24:01)'), names_to = "metabolite", values_to = "level")|>
  pivot_wider(names_from = metabolite, values_from = level) |>
  mutate(time = 3) |>
  select(Sample, time, PMI_high, everything())

bc_features_3_clinic <- bc_features_3  |>
  separate(Sample, into = c("Sample", "Letter", "Number"), sep = "_")

bc_features_3_clinic<- left_join(bc_features_3_clinic, meta_master, by = c("Number" = "Sample ID Number"))

str(bc_features_3_clinic)

bc_features_14 <- cart_pmi_14_norm |>
  pivot_longer(cols = ('Alanine':'Cer-(24:01)'), names_to = "metabolite", values_to = "level")|>
  pivot_wider(names_from = metabolite, values_from = level) |>
  mutate(time = 14) |>
  select(Sample, time, PMI_high, everything())

bc_features_14_clinic <- bc_features_14  |>
  separate(Sample, into = c("Sample", "Letter", "Number"), sep = "_")

bc_features_14_clinic<- left_join(bc_features_14_clinic, meta_master, by = c("Number" = "Sample ID Number"))

str(bc_features_14_clinic)


p_ci_crs_bc_3 <- list()

for (metabolite in crs_bc_3_metabolites) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(CRS_onset, as.numeric(`CRS_high`)) ~ ifelse(get(metabolite, crs_bc_features_3_clinic) > mean(get(metabolite, crs_bc_features_3_clinic)), 1, 0), data = crs_bc_features_3_clinic)
  
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
  p_ci_crs_bc_3[[metabolite]] <- p
}

grobs_ci_crs_bc_3 <- lapply(p_ci_crs_bc_3, function(x) ggplotGrob(x$plot))
do.call(grid.arrange, c(grobs_ci_crs_bc_3, ncol = 5))

p_ci_crs_bc_3$`PEA-(40:07)`

p_ci_crs_bc_3$`AC-(10:0)`

p_ci_crs_bc_3$`LPC-(20:03)`


p_km_bc_pfs_cox_3 <- list()

for (metabolite in bc_pfs_3_metabolites) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(PFS_days, PFS_event) ~ ifelse(get(metabolite, bc_features_3_clinic) > mean(get(metabolite, bc_features_3_clinic)), 1, 0), data = bc_features_3_clinic)
  
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
  p_km_bc_pfs_cox_3[[metabolite]] <- p
}

grobs_bc_pfs_cox_3 <- lapply(p_km_bc_pfs_cox_3, function(x) ggplotGrob(x$plot))
p_grobs_bc_pfs_cox_3 <- do.call(grid.arrange, c(grobs_bc_pfs_cox_3))

p_km_bc_pfs_cox_14 <- list()

for (metabolite in bc_pfs_14_metabolites) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(PFS_days, PFS_event) ~ ifelse(get(metabolite, bc_features_14_clinic) > mean(get(metabolite, bc_features_14_clinic)), 1, 0), data = bc_features_14_clinic)
  
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
  p_km_bc_pfs_cox_14[[metabolite]] <- p
}

grobs_bc_pfs_cox_14 <- lapply(p_km_bc_pfs_cox_14, function(x) ggplotGrob(x$plot))
p_grobs_bc_pfs_cox_14 <- do.call(grid.arrange, c(grobs_bc_pfs_cox_14))


p_km_bc_os_cox_0 <- list()

for (metabolite in bc_os_0_metabolites) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(OS_days, OS_event) ~ ifelse(get(metabolite, bc_features_0_clinic) > mean(get(metabolite, bc_features_0_clinic)), 1, 0), data = bc_features_0_clinic)
  
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
                  palette = c("lightgrey", "darkgrey"))
  
  
  # Add the plot to the list
  p_km_bc_os_cox_0[[metabolite]] <- p
}

grobs_bc_os_cox_0 <- lapply(p_km_bc_os_cox_0, function(x) ggplotGrob(x$plot))
p_grobs_bc_os_cox_0 <- do.call(grid.arrange, c(grobs_bc_os_cox_0))



p_km_bc_os_cox_3 <- list()

for (metabolite in bc_os_3_metabolites) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(OS_days, OS_event) ~ ifelse(get(metabolite, bc_features_3_clinic) > mean(get(metabolite, bc_features_3_clinic)), 1, 0), data = bc_features_3_clinic)
  
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
  p_km_bc_os_cox_3[[metabolite]] <- p
}

grobs_bc_os_cox_3 <- lapply(p_km_bc_os_cox_3, function(x) ggplotGrob(x$plot))
p_grobs_bc_os_cox_3 <- do.call(grid.arrange, c(grobs_bc_os_cox_3))

p_km_bc_os_cox_14 <- list()

for (metabolite in bc_os_14_metabolites) {
  # Directly use the formula in survfit2
  fit <- survfit2(Surv(OS_days, OS_event) ~ ifelse(get(metabolite, bc_features_14_clinic) > mean(get(metabolite, bc_features_14_clinic)), 1, 0), data = bc_features_14_clinic)
  
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
                  palette = c("orange", "orange3"))
  
  
  # Add the plot to the list
  p_km_bc_os_cox_14[[metabolite]] <- p
}

grobs_bc_os_cox_14 <- lapply(p_km_bc_os_cox_14, function(x) ggplotGrob(x$plot))
p_grobs_bc_os_cox_14 <- do.call(grid.arrange, c(grobs_bc_os_cox_14))


p_km_bc_os_cox_0$`PlasC-(36:00)`
p_km_bc_os_cox_3$`FA-(18:00)`
p_km_bc_os_cox_3$`LPI-(22:05)`
p_km_bc_os_cox_14$Lactate


### MSEA analysis ----
### Enrichment analysis of significant metabolites from log and cox regression models

#Preparation of data sets for MSEA
metabolites_HMDB <- read_xlsx("C:/Users/dcs54/Desktop/Projekt_CART_Metabolomics/MetaboAnalyst_Input_Data/Finale_Inpute_Dateien/Werner_HMDB_translation.xlsx", na = "NA")

metabolites_HMDB <- metabolites_HMDB %>%
  rename(metabolite = `...1`)

#crs_cox beinhaltet alles von crs_glm
bc_crs_cox_HMDB <- crs_bc_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker, HR)

bc_crs_cox_HMDB <- left_join(bc_crs_cox_HMDB, metabolites_HMDB, by = c("marker" = "metabolite"))

bc_crs_cox_b1_HMDB <- bc_crs_cox_HMDB |>
  filter(HR<1)|>
  select(HMDB)|>
  unlist()|>
  as.vector() |>
  unique()

bc_crs_cox_o1_HMDB <- bc_crs_cox_HMDB |>
  filter(HR>1)|>
  select(HMDB)|>
  unlist()|>
  as.vector() |>
  unique()

##pfs does not include all of os
bc_pfs_cox_HMDB <- bc_pfs_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker, HR)

bc_os_cox_HMDB <- bc_os_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker, HR)

bc_survival_cox_HMDB <- rbind(bc_pfs_cox_HMDB, bc_os_cox_HMDB)

bc_survival_cox_HMDB <- bc_survival_cox_HMDB |>
  distinct(marker, .keep_all = T)

bc_survival_cox_HMDB <- left_join(bc_survival_cox_HMDB, metabolites_HMDB, by = c("marker" = "metabolite"))

bc_survival_cox_b1_HMDB <- bc_survival_cox_HMDB |>
  filter(HR<1)|>
  select(HMDB)|>
  unlist()|>
  as.vector() |>
  unique()

bc_survival_cox_o1_HMDB <- bc_survival_cox_HMDB |>
  filter(HR>1)|>
  select(HMDB)|>
  unlist()|>
  as.vector() |>
  unique()

## Enrichment analysis of CRS
#Enrichment of HR<1
mset_bc_crs_cox_b1_msea<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-bc_crs_cox_b1_HMDB
mset_bc_crs_cox_b1_msea<-Setup.MapData(mset_bc_crs_cox_b1_msea, cmpd.vec);
mset_bc_crs_cox_b1_msea<-CrossReferencing(mset_bc_crs_cox_b1_msea, "hmdb", lipid = T);
mset_bc_crs_cox_b1_msea<-CreateMappingResultTable(mset_bc_crs_cox_b1_msea)
mset_bc_crs_cox_b1_msea<-SetMetabolomeFilter(mset_bc_crs_cox_b1_msea, F);
mset_bc_crs_cox_b1_msea<-SetCurrentMsetLib(mset_bc_crs_cox_b1_msea, "sub_class", 2);
mset_bc_crs_cox_b1_msea<-CalculateHyperScore(mset_bc_crs_cox_b1_msea)

bc_crs_cox_b1_msea <- read_xlsx("C:/Users/dcs54/Desktop/Projekt_CART_Metabolomics/MetaboAnalyst_Input_Data/Finale_Inpute_Dateien/HMDB_crs_b1.xlsx", na = "NA")
bc_crs_cox_b1_msea <- bc_crs_cox_b1_msea %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(bc_crs_cox_b1_msea)

bc_crs_cox_b1_msea$"Raw p" <- colnames("Raw_p")

# 
# mSet<-InitDataObjects("conc", "msetora", FALSE)
# cmpd.vec<-c("HMDB0008730","HMDB0008732","HMDB0010379","HMDB0010387","HMDB0010392","HMDB0010393","HMDB0010397","HMDB0010406","HMDB0011512","HMDB0011514")
# mSet<-Setup.MapData(mSet, cmpd.vec);
# mSet<-CrossReferencing(mSet, "hmdb", lipid = T);
# mSet<-CreateMappingResultTable(mSet)
# mSet<-SetMetabolomeFilter(mSet, F);
# mSet<-SetCurrentMsetLib(mSet, "sub_class", 2);
# mSet<-CalculateHyperScore(mSet)
# 
# mSet[["analSet"]][["ora.mat"]]


#MSEA of HR>1 metabolites
mset_bc_crs_cox_o1_msea<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-bc_crs_cox_o1_HMDB
mset_bc_crs_cox_o1_msea<-Setup.MapData(mset_bc_crs_cox_o1_msea, cmpd.vec);
mset_bc_crs_cox_o1_msea<-CrossReferencing(mset_bc_crs_cox_o1_msea, "hmdb", lipid = T);
mset_bc_crs_cox_o1_msea<-CreateMappingResultTable(mset_bc_crs_cox_o1_msea)
mset_bc_crs_cox_o1_msea<-SetMetabolomeFilter(mset_bc_crs_cox_o1_msea, F);
mset_bc_crs_cox_o1_msea<-SetCurrentMsetLib(mset_bc_crs_cox_o1_msea, "sub_class", 2);
mset_bc_crs_cox_o1_msea<-CalculateHyperScore(mset_bc_crs_cox_o1_msea)

bc_crs_cox_o1_msea <- read_xlsx("C:/Users/dcs54/Desktop/Projekt_CART_Metabolomics/MetaboAnalyst_Input_Data/Finale_Inpute_Dateien/HMDB_crs_o1.xlsx", na = "NA")
bc_crs_cox_o1_msea <- bc_crs_cox_o1_msea %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(bc_crs_cox_o1_msea)

bc_crs_cox_o1_msea$"Raw p" <- colnames("Raw_p")

bc_crs_cox_b1_msea <- bc_crs_cox_b1_msea |>
  mutate(direction = "HR<1")
bc_crs_cox_o1_msea <- bc_crs_cox_o1_msea |>
  mutate(direction = "HR>1")

bc_crs_cox_msea_combined <- rbind(bc_crs_cox_b1_msea, bc_crs_cox_o1_msea)

p_bc_crs_cox_msea_combined <- bc_crs_cox_msea_combined |>
  filter(FDR < 0.1) |>
  # filter(hits >= 1)|>
  ggplot()+
  geom_point(aes(x=direction, y=pathway, size=ratio, color=Raw_p))+
  scale_color_gradient(low = "darkorange", high = "darkblue") +
  #scale_size(limits = c(1,100), breaks = c(1,10,100))+
  #scale_x_discrete(limits = c("123", "121", "321"), labels = c("Increase \n pattern", "In-Decrease \n pattern", "Decrease \n pattern"))+
  scale_y_discrete(limits = rev(levels(bc_crs_cox_msea_combined))) +
  labs(y="", x="", size = "Enrichment \n score", color = "P-value")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a black frame
  ) 

p_bc_crs_cox_msea_combined

## Enrichment analysis of Survival
#Enrichment of HR<1
mset_bc_survival_cox_b1_msea<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-bc_survival_cox_b1_HMDB
mset_bc_survival_cox_b1_msea<-Setup.MapData(mset_bc_survival_cox_b1_msea, cmpd.vec);
mset_bc_survival_cox_b1_msea<-CrossReferencing(mset_bc_survival_cox_b1_msea, "hmdb");
mset_bc_survival_cox_b1_msea<-CreateMappingResultTable(mset_bc_survival_cox_b1_msea)
mset_bc_survival_cox_b1_msea<-SetMetabolomeFilter(mset_bc_survival_cox_b1_msea, F);
mset_bc_survival_cox_b1_msea<-SetCurrentMsetLib(mset_bc_survival_cox_b1_msea, "sub_class", 2);
mset_bc_survival_cox_b1_msea<-CalculateHyperScore(mset_bc_survival_cox_b1_msea)

bc_survival_cox_b1_msea <- read_xlsx("C:/Users/dcs54/Desktop/Projekt_CART_Metabolomics/MetaboAnalyst_Input_Data/Finale_Inpute_Dateien/HMDB_survival_b1.xlsx", na = "NA")
bc_survival_cox_b1_msea <- bc_survival_cox_b1_msea %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(bc_survival_cox_b1_msea)

bc_survival_cox_b1_msea$"Raw p" <- colnames("Raw_p")


#MSEA of HR>1 metabolites
mset_bc_survival_cox_o1_msea<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-bc_survival_cox_o1_HMDB
mset_bc_survival_cox_o1_msea<-Setup.MapData(mset_bc_survival_cox_o1_msea, cmpd.vec);
mset_bc_survival_cox_o1_msea<-CrossReferencing(mset_bc_survival_cox_o1_msea, "hmdb", lipid = T);
mset_bc_survival_cox_o1_msea<-CreateMappingResultTable(mset_bc_survival_cox_o1_msea)
mset_bc_survival_cox_o1_msea<-SetMetabolomeFilter(mset_bc_survival_cox_o1_msea, F);
mset_bc_survival_cox_o1_msea<-SetCurrentMsetLib(mset_bc_survival_cox_o1_msea, "sub_class", 2);
mset_bc_survival_cox_o1_msea<-CalculateHyperScore(mset_bc_survival_cox_o1_msea)

bc_survival_cox_o1_msea <- read_xlsx("C:/Users/dcs54/Desktop/Projekt_CART_Metabolomics/MetaboAnalyst_Input_Data/Finale_Inpute_Dateien/HMDB_crs_o1.xlsx", na = "NA")
bc_survival_cox_o1_msea <- bc_survival_cox_o1_msea %>%
  mutate(ratio = hits/expected) %>%
  rename("pathway" = "...1", "Raw_p" = "Raw p", "Holm_p" = "Holm p") %>%
  select(pathway, ratio, everything())

str(bc_survival_cox_o1_msea)

bc_survival_cox_o1_msea$"Raw p" <- colnames("Raw_p")

bc_survival_cox_b1_msea <- bc_survival_cox_b1_msea |>
  mutate(direction = "HR<1")
bc_survival_cox_o1_msea <- bc_survival_cox_o1_msea |>
  mutate(direction = "HR>1")

bc_survival_cox_msea_combined <- rbind(bc_survival_cox_b1_msea, bc_survival_cox_o1_msea)

p_bc_survival_cox_msea_combined <- bc_survival_cox_msea_combined |>
  filter(FDR < 0.1) |>
  # filter(hits >= 1)|>
  ggplot()+
  geom_point(aes(x=direction, y=pathway, size=ratio, color=Raw_p))+
  scale_color_gradient(low = "darkorange", high = "darkblue") +
  #scale_size(limits = c(1,100), breaks = c(1,10,100))+
  #scale_x_discrete(limits = c("123", "121", "321"), labels = c("Increase \n pattern", "In-Decrease \n pattern", "Decrease \n pattern"))+
  scale_y_discrete(limits = rev(levels(bc_survival_cox_msea_combined))) +
  labs(y="", x="", size = "Enrichment \n score", color = "P-value")+
  theme_classic()+  # Chopfse a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a black frame
  ) 

p_bc_survival_cox_msea_combined


### Correlation analysis with Olink data from V. Blumenberg ----
## Loading of Olink data and linking with baseline data

olink <- read_xlsx("C:/Users/dcs54/Desktop/Projekt_CART_Metabolomics/MetaboAnalyst_Input_Data/Finale_Inpute_Dateien/olink.xlsx", na = "NA")

bc_0_olink <- left_join(cart_sat_0_norm, olink)

bc_0_olink_master <- bc_0_olink|>
  pivot_longer(cols = c(ADA:TNFRSF21), names_to = "cytokines", values_to = "level")|>
  pivot_longer(cols = c(Alanine:'Cer-(24:01)'), names_to = "metabolite", values_to = "abundance")

cytokine <- bc_0_olink|>
  select(ADA:TNFRSF21) |>
  colnames()|>
  unlist()|>
  as.vector()

crs_bc_markers <- crs_bc_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unlist() |>
  as.vector()

bc_0_crs_olink_corr <- data.frame()
for (i in crs_bc_markers) {
  for (k in cytokine) {
    # Extract the relevant columns from the data frame
    subset_data <- bc_0_olink[, c(k, i)]
    
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
    bc_0_crs_olink_corr <- rbind(bc_0_crs_olink_corr, result_row)
  }
}

bc_0_crs_olink_corr <- bc_0_crs_olink_corr |>
  mutate(FDR = p.adjust(P_Value, method = "BH"))

bc_0_crs_olink_corr  |>
  ggplot(aes(x=Metabolite, y=Cytokine, color=Correlation))+
  geom_point(size = 5)+
  geom_text(aes(label = ifelse(P_Value < 0.05, "*", "")), size = 6, color = "black",
            nudge_x = 0.015, nudge_y = -0.06)+
  scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue" )+
  scale_y_discrete(limits = rev(levels(bc_0_crs_olink_corr$Cytokine)))+
  labs(x="",y="")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black frame
    axis.text.x = element_text(angle = 45, hjust=0.9))


#dummy column for joining tables
cart_sat_3_norm_join <- cart_sat_3_norm |>
  mutate(common_key = gsub("_B_", "_A_", Sample))

bc_3_olink <- left_join(cart_sat_3_norm_join, olink, by = join_by(common_key == Sample))

bc_3_crs_olink_corr <- data.frame()
for (i in crs_bc_markers) {
  for (k in cytokine) {
    # Extract the relevant columns from the data frame
    subset_data <- bc_3_olink[, c(k, i)]
    
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
    bc_3_crs_olink_corr <- rbind(bc_3_crs_olink_corr, result_row)
  }
}

bc_3_crs_olink_corr <- bc_3_crs_olink_corr |>
  mutate(FDR = p.adjust(P_Value, method = "BH"))

bc_3_crs_olink_corr  |>
  ggplot(aes(x=Metabolite, y=Cytokine, color=Correlation))+
  geom_point(size = 5)+
  geom_text(aes(label = ifelse(FDR < 0.1, "*", "")), size = 6, color = "black",
            nudge_x = 0.015, nudge_y = -0.06)+
  scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue" )+
  scale_y_discrete(limits = rev(levels(bc_3_crs_olink_corr$Cytokine)))+
  labs(x="",y="")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black frame
    axis.text.x = element_text(angle = 45, hjust=0.9))

## Correlation between day 3 metabolites and 

os_bc_markers <- bc_os_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unlist() |>
  unique() |>
  as.vector()

pfs_bc_markers <- bc_pfs_cox_combined |>
  filter(FDR < 0.1)|>
  select(marker) |>
  unlist() |>
  unique() |>
  as.vector()

survival_bc_markers <- c(os_bc_markers, pfs_bc_markers) |>
  unique()

bc_0_survival_olink_corr <- data.frame()

for (i in survival_bc_markers) {
  for (k in cytokine) {
    # Extract the relevant columns from the data frame
    subset_data <- bc_0_olink[, c(k, i)]
    
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
    bc_0_survival_olink_corr <- rbind(bc_0_survival_olink_corr, result_row)
  }
}

bc_0_survival_olink_corr <- bc_0_survival_olink_corr |>
  mutate(FDR = p.adjust(P_Value, method = "BH"))

bc_0_survival_olink_corr  |>
  ggplot(aes(x=Metabolite, y=Cytokine, color=Correlation))+
  geom_point(size = 5)+
  geom_text(aes(label = ifelse(P_Value < 0.05, "*", "")), size = 6, color = "black",
            nudge_x = 0.015, nudge_y = -0.06)+
  scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue" )+
  scale_y_discrete(limits = rev(levels(bc_0_survival_olink_corr$Cytokine)))+
  labs(x="",y="")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black frame
    axis.text.x = element_text(angle = 45, hjust=0.9))


bc_3_survival_olink_corr <- data.frame()

for (i in survival_bc_markers) {
  for (k in cytokine) {
    # Extract the relevant columns from the data frame
    subset_data <- bc_3_olink[, c(k, i)]
    
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
    bc_3_survival_olink_corr <- rbind(bc_3_survival_olink_corr, result_row)
  }
}

bc_3_survival_olink_corr <- bc_3_survival_olink_corr |>
  mutate(FDR = p.adjust(P_Value, method = "BH"))

bc_3_survival_olink_corr  |>
  ggplot(aes(x=Metabolite, y=Cytokine, color=Correlation))+
  geom_point(size = 5)+
  geom_text(aes(label = ifelse(FDR < 0.1, "*", "")), size = 6, color = "black",
            nudge_x = 0.015, nudge_y = -0.06)+
  scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue" )+
  scale_y_discrete(limits = rev(levels(bc_3_survival_olink_corr$Cytokine)))+
  labs(x="",y="")+
  theme_classic()+  # Choose a theme as a starting point
  theme(
    axis.ticks = element_blank(),  # Remove ticks from both axes
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black frame
    axis.text.x = element_text(angle = 45, hjust=0.9))

### Analysis to which time pattern the selected features belong ----
## Left_join of CRS relevant metabolites and time pattern with keeping all on the left side

sat_crs_all <- rbind(sat_crs_cox_HMDB, sat_survival_cox_HMDB)

sat_crs_all <- sat_crs_all |>
  distinct(marker)

crs_bc_cox_combined_marker <- crs_bc_cox_combined |>
  filter(FDR < 0.1) |>
  distinct(marker)

bc_crs_time_pattern <- left_join(crs_bc_cox_combined_marker, cart_time_pattern, by = join_by("marker" == "metabolite"))

bc_crs_time_pattern$direction <- as.character(bc_crs_time_pattern$direction)

bc_crs_time_pattern |>
  count(direction) |>
  ggplot(aes(x = "", y = n, fill = direction)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(fill = "Category")
