
validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  group_by(Entity)|>
  count()

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  group_by(Geschlecht)|>
  count()

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  summarise(
    median_time = median(Column1, na.rm = TRUE),                     
    iqr = IQR(Column1, na.rm = TRUE),                            
    mean_time = mean(Column1, na.rm = TRUE),,                        
    q25 = quantile(Column1, 0.25, na.rm = TRUE),                   
    q75 = quantile(Column1, 0.75, na.rm = TRUE),
    max = max(Column1, na.rm = T),
    min = min(Column1, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  summarise(
    median_time = median(BMI, na.rm = TRUE),                     
    iqr = IQR(BMI, na.rm = TRUE),                            
    mean_time = mean(BMI, na.rm = TRUE),,                        
    q25 = quantile(BMI, 0.25, na.rm = TRUE),                   
    q75 = quantile(BMI, 0.75, na.rm = TRUE),
    max = max(BMI, na.rm = T),
    min = min(BMI, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|> 
  summarise(
    median_time = median(Waist, na.rm = TRUE),                     
    iqr = IQR(Waist, na.rm = TRUE),                            
    mean_time = mean(Waist, na.rm = TRUE),,                        
    q25 = quantile(Waist, 0.25, na.rm = TRUE),                   
    q75 = quantile(Waist, 0.75, na.rm = TRUE),
    max = max(Waist, na.rm = T),
    min = min(Waist, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|> 
  summarise(
    median_time = median(SAT, na.rm = TRUE),                     
    iqr = IQR(SAT, na.rm = TRUE),                            
    mean_time = mean(SAT, na.rm = TRUE),,                        
    q25 = quantile(SAT, 0.25, na.rm = TRUE),                   
    q75 = quantile(SAT, 0.75, na.rm = TRUE),
    max = max(SAT, na.rm = T),
    min = min(SAT, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|> 
  summarise(
    median_time = median(VAT, na.rm = TRUE),                     
    iqr = IQR(VAT, na.rm = TRUE),                            
    mean_time = mean(VAT, na.rm = TRUE),,                        
    q25 = quantile(VAT, 0.25, na.rm = TRUE),                   
    q75 = quantile(VAT, 0.75, na.rm = TRUE),
    max = max(VAT, na.rm = T),
    min = min(VAT, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  summarise(
    median_time = median(STLV, na.rm = TRUE),                     
    iqr = IQR(STLV, na.rm = TRUE),                            
    mean_time = mean(STLV, na.rm = TRUE),,                        
    q25 = quantile(STLV, 0.25, na.rm = TRUE),                   
    q75 = quantile(STLV, 0.75, na.rm = TRUE),
    max = max(STLV, na.rm = T),
    min = min(STLV, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  summarise(
    median_time = median(Leuko, na.rm = TRUE),                     
    iqr = IQR(Leuko, na.rm = TRUE),                            
    mean_time = mean(Leuko, na.rm = TRUE),,                        
    q25 = quantile(Leuko, 0.25, na.rm = TRUE),                   
    q75 = quantile(Leuko, 0.75, na.rm = TRUE),
    max = max(Leuko, na.rm = T),
    min = min(Leuko, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  summarise(
    median_time = median(Ferritin, na.rm = TRUE),                     
    iqr = IQR(Ferritin, na.rm = TRUE),                            
    mean_time = mean(Ferritin, na.rm = TRUE),,                        
    q25 = quantile(Ferritin, 0.25, na.rm = TRUE),                   
    q75 = quantile(Ferritin, 0.75, na.rm = TRUE),
    max = max(Ferritin, na.rm = T),
    min = min(Ferritin, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  summarise(
    median_time = median(LDH, na.rm = TRUE),                     
    iqr = IQR(LDH, na.rm = TRUE),                            
    mean_time = mean(LDH, na.rm = TRUE),,                        
    q25 = quantile(LDH, 0.25, na.rm = TRUE),                   
    q75 = quantile(LDH, 0.75, na.rm = TRUE),
    max = max(LDH, na.rm = T),
    min = min(LDH, na.rm = T))

validation_master |>
  filter(Entity != "ALL")|>
  filter(Time == "Day 0")|>
  summarise(
    median_time = median(CRP, na.rm = TRUE),                     
    iqr = IQR(CRP, na.rm = TRUE),                            
    mean_time = mean(CRP, na.rm = TRUE),,                        
    q25 = quantile(CRP, 0.25, na.rm = TRUE),                   
    q75 = quantile(CRP, 0.75, na.rm = TRUE),
    max = max(CRP, na.rm = T),
    min = min(CRP, na.rm = T))




t2_vat_crs_norm_bc |>
  group_by(Construct)|>
  count()

t2_vat_crs_norm_bc |>
  group_by(Geschlecht)|>
  count()

t2_vat_crs_norm_bc |>
  group_by(Costim)|>
  count()

all_master|>
  filter(cohort == "training")|>
  summarise(
    median_time = median(Age, na.rm = TRUE),                     
    iqr = IQR(Age, na.rm = TRUE),                            
    mean_time = mean(Age, na.rm = TRUE),,                        
    q25 = quantile(Age, 0.25, na.rm = TRUE),                   
    q75 = quantile(Age, 0.75, na.rm = TRUE),
    max = max(Age, na.rm = T),
    min = min(Age, na.rm = T))


t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(Column1, na.rm = TRUE),                     
    iqr = IQR(Column1, na.rm = TRUE),                            
    mean_time = mean(Column1, na.rm = TRUE),,                        
    q25 = quantile(Column1, 0.25, na.rm = TRUE),                   
    q75 = quantile(Column1, 0.75, na.rm = TRUE),
    max = max(Column1, na.rm = T),
    min = min(Column1, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(BMI, na.rm = TRUE),                     
    iqr = IQR(BMI, na.rm = TRUE),                            
    mean_time = mean(BMI, na.rm = TRUE),,                        
    q25 = quantile(BMI, 0.25, na.rm = TRUE),                   
    q75 = quantile(BMI, 0.75, na.rm = TRUE),
    max = max(BMI, na.rm = T),
    min = min(BMI, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(Waist, na.rm = TRUE),                     
    iqr = IQR(Waist, na.rm = TRUE),                            
    mean_time = mean(Waist, na.rm = TRUE),,                        
    q25 = quantile(Waist, 0.25, na.rm = TRUE),                   
    q75 = quantile(Waist, 0.75, na.rm = TRUE),
    max = max(Waist, na.rm = T),
    min = min(Waist, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(SAT, na.rm = TRUE),                     
    iqr = IQR(SAT, na.rm = TRUE),                            
    mean_time = mean(SAT, na.rm = TRUE),,                        
    q25 = quantile(SAT, 0.25, na.rm = TRUE),                   
    q75 = quantile(SAT, 0.75, na.rm = TRUE),
    max = max(SAT, na.rm = T),
    min = min(SAT, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(VAT, na.rm = TRUE),                     
    iqr = IQR(VAT, na.rm = TRUE),                            
    mean_time = mean(VAT, na.rm = TRUE),,                        
    q25 = quantile(VAT, 0.25, na.rm = TRUE),                   
    q75 = quantile(VAT, 0.75, na.rm = TRUE),
    max = max(VAT, na.rm = T),
    min = min(VAT, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(STLV, na.rm = TRUE),                     
    iqr = IQR(STLV, na.rm = TRUE),                            
    mean_time = mean(STLV, na.rm = TRUE),,                        
    q25 = quantile(STLV, 0.25, na.rm = TRUE),                   
    q75 = quantile(STLV, 0.75, na.rm = TRUE),
    max = max(STLV, na.rm = T),
    min = min(STLV, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(Leuko, na.rm = TRUE),                     
    iqr = IQR(Leuko, na.rm = TRUE),                            
    mean_time = mean(Leuko, na.rm = TRUE),,                        
    q25 = quantile(Leuko, 0.25, na.rm = TRUE),                   
    q75 = quantile(Leuko, 0.75, na.rm = TRUE),
    max = max(Leuko, na.rm = T),
    min = min(Leuko, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(Ferritin, na.rm = TRUE),                     
    iqr = IQR(Ferritin, na.rm = TRUE),                            
    mean_time = mean(Ferritin, na.rm = TRUE),,                        
    q25 = quantile(Ferritin, 0.25, na.rm = TRUE),                   
    q75 = quantile(Ferritin, 0.75, na.rm = TRUE),
    max = max(Ferritin, na.rm = T),
    min = min(Ferritin, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(LDH, na.rm = TRUE),                     
    iqr = IQR(LDH, na.rm = TRUE),                            
    mean_time = mean(LDH, na.rm = TRUE),,                        
    q25 = quantile(LDH, 0.25, na.rm = TRUE),                   
    q75 = quantile(LDH, 0.75, na.rm = TRUE),
    max = max(LDH, na.rm = T),
    min = min(LDH, na.rm = T))

t2_vat_crs_norm_bc |>
  summarise(
    median_time = median(CRP, na.rm = TRUE),                     
    iqr = IQR(CRP, na.rm = TRUE),                            
    mean_time = mean(CRP, na.rm = TRUE),,                        
    q25 = quantile(CRP, 0.25, na.rm = TRUE),                   
    q75 = quantile(CRP, 0.75, na.rm = TRUE),
    max = max(CRP, na.rm = T),
    min = min(CRP, na.rm = T))



