#remove all files on environment
rm(list = ls())

#check working directory
getwd()
#set working directory
setwd("/Users/nafisahalias/Documents/Imperial/Research project/R/Data")
#check working directory
getwd()
SEM_EDS<-read.csv("SEM_EDS_Data.csv", header=TRUE)

# Load necessary packages
library(tidyverse)

# Convert all nd to 0
SEM_EDS[SEM_EDS == "nd"] <- 0

# Convert all mol columns to numeric
mol_cols <- c("Si_mol", "O_mol", "Al_mol", "Zn_mol", "P_mol")
SEM_EDS <- SEM_EDS %>%
  mutate(across(all_of(mol_cols), ~ as.numeric(.)))

#add Al/Si column
SEM_EDS <- SEM_EDS %>%
  mutate(
    AlSi_ratio = Al_mol / Si_mol,
    ZnSi_ratio = Zn_mol / Si_mol)

# Pivot to long format
long_data <- SEM_EDS %>%
  pivot_longer(cols = c(AlSi_ratio, ZnSi_ratio),
               names_to = "Metal",
               values_to = "Ratio") %>%
  mutate(
    Conc_label = ifelse(is.na(Concentration), "Control", paste0(Concentration, "M"))
  )

# Filter to keep only Species TP and TW
filtered_data <- long_data %>%
  filter(Species %in% c("T_weissflogii", "T_pseudonana"))

# Filter only valid Metal–Concentration combinations
filtered_data <- filtered_data %>%
  filter(
    (Metal == "AlSi_ratio" & Conc_label %in% c("0M", "0.1M", "1M")) |
      (Metal == "ZnSi_ratio" & Conc_label %in% c("0M", "8M", "24M"))
  )

summary_data <- filtered_data %>%
  group_by(Species, Metal, Conc_label) %>%
  summarise(
    mean_ratio = mean(Ratio, na.rm = TRUE),
    se_ratio = sd(Ratio, na.rm = TRUE) / sqrt(n()),
    sd_ratio = sd(Ratio, na.rm=TRUE),
    median_ratio = median(Ratio, na.rm = TRUE),
    .groups = "drop"
  )

# Reorder x_group levels
summary_data <- summary_data %>%
  mutate(
    x_group = paste0(Metal, ".", Conc_label)
  )
summary_data <- summary_data %>%
  mutate(
    x_group = factor(
      x_group,
      levels = c(
        "AlSi_ratio.0M",
        "AlSi_ratio.0.1M",
        "AlSi_ratio.1M",
        " ", #spacer
        "ZnSi_ratio.0M",
        "ZnSi_ratio.8M",
        "ZnSi_ratio.24M"
      )
    )
  )

gap_rows <- summary_data %>%
  distinct(Species) %>%
  mutate(
    x_group = factor(" ", levels = levels(summary_data$x_group)),
    mean_ratio = NA,
    se_ratio = NA,
    Metal = NA
  )

# Add gap rows to your data
summary_data_with_gap <- bind_rows(summary_data, gap_rows)
str(summary_data_with_gap$Species)

summary_data_with_gap$Species <- dplyr::recode(summary_data_with_gap$Species,
                                               "T_pseudonana" = "T. pseudonana",
                                               "T_weissflogii" = "T. weissflogii")
species_labels <- c(
  "T. pseudonana" = "italic('T. pseudonana')",
  "T. weissflogii" = "italic('T. weissflogii')"
)

# Create significance annotations dataframe
sig_data <- data.frame(
  Species = c("T. pseudonana", "T. pseudonana", "T. pseudonana", "T. pseudonana",
              "T. weissflogii", "T. weissflogii", "T. weissflogii", "T. weissflogii"),
  x_group = c("AlSi_ratio.0.1M", "AlSi_ratio.1M", "ZnSi_ratio.8M", "ZnSi_ratio.24M",
              "AlSi_ratio.0.1M", "AlSi_ratio.1M", "ZnSi_ratio.8M", "ZnSi_ratio.24M"),
  Metal = c("Al", "Al", "Zn", "Zn",
            "Al", "Al", "Zn", "Zn"),
  y_pos = c(0.128, 0.28, 0.21, 0.019,
            0.09, 0.0699, 0.0195, 0.0665),
  significance = c("***", "***", "***", "*",
                   "ns", "ns", "ns", "***")
)

metal_colors <- c("AlSi_ratio" = "#E69F00", "ZnSi_ratio" = "#0173B2") 

library(ggplot2)
ggplot(summary_data_with_gap, aes(x = x_group, y = mean_ratio, fill = Metal)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
  geom_errorbar(aes(ymin = mean_ratio - se_ratio, ymax = mean_ratio + se_ratio),
                width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(
    ~Species,
    scales = "free_x",
    nrow = 1,
    labeller = labeller(Species = as_labeller(species_labels, label_parsed))
  ) +
  scale_x_discrete(labels = c(
    "AlSi_ratio.0M" = "Control",
    "AlSi_ratio.0.1M" = "0.1 \u03bcM Al",
    "AlSi_ratio.1M" = "1.0 \u03bcM Al",
    "ZnSi_ratio.0M" = "Control",
    "ZnSi_ratio.8M" = "8.0 \u03bcM Zn",
    "ZnSi_ratio.24M" = "24 \u03bcM Zn"
  )) +
  scale_fill_manual(values = metal_colors, na.translate = FALSE,
                    labels = c("AlSi_ratio" = "Al/Si ratio", 
                               "ZnSi_ratio" = "Zn/Si ratio")) +
  labs(
    x = "Treatment",
    y = "Metal/silicon ratio in diatoms",
    fill = "Metal/silicon ratio"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  geom_text(data = sig_data, 
            aes(x = x_group, y = y_pos, label = significance),
            size = 4, vjust = -0.5,
            inherit.aes = FALSE)

#STATISTICAL MODELS
#load libraries
library(lme4)
library(lmtest)
library(sjPlot)
library(performance)
library(MuMIn)

#TW vs Al: df and lmm
# Create data frame for T_weissflogii in Al (including control)
TW_Al <- SEM_EDS %>%
  filter(Species == "T_weissflogii" & (Treatment == "Al" | Treatment == "control")) %>%
  select(Species, CellID, PointID, Treatment, Concentration, AlSi_ratio) %>%
  arrange(Treatment, Concentration, CellID, PointID)
str(TW_Al)
hist(TW_Al$AlSi_ratio)

# Square root transform
TW_Al$sqrt_AlSi_ratio <- sqrt(TW_Al$AlSi_ratio)

#check transformed distribution
hist(TW_Al$sqrt_AlSi_ratio)

# use sqrt-transformed response for model
library(lmerTest)
TW_Al$Concentration <- factor(TW_Al$Concentration)
TW_Al$Concentration <- as.numeric(as.character(TW_Al$Concentration))
model_TW_Al <- lmer(sqrt_AlSi_ratio ~ Concentration + (1|CellID), data = TW_Al)
summary(model_TW_Al)


#check residuals
plot(model_TW_Al) 
qqnorm(residuals(model_TW_Al))
qqline(residuals(model_TW_Al)) 
check_model(model_TW_Al)
check_normality(model_TW_Al) 
## NOTE: Non-normality of residuals detected (p =0.027). might be due to the outlier in histogram which is biologically plausible
## visual inspection of q-q plot is good. model will be used.

#likelihood convergence test
null_model_TW_Al <- lmer(log_AlSi_ratio ~ 1 + (1|CellID), data = TW_Al)
anova(null_model_TW_Al, model_TW_Al)

#plot linear model of TW vs Al
TW_Al$Concentration <- as.numeric(as.character(TW_Al$Concentration))
LM3<-ggplot(TW_Al, aes(x = Concentration, y = sqrt_AlSi_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", 
              se = TRUE,
              colour = "#E69F00",
              formula = y ~ x) +
  theme_minimal() +
  ylab(expression("√(Al/Si ratio) in"~italic("T. weissflogii"))) +
  xlab("Concentration of Al added (µM)") 

plot(LM3)

##TW vs Zn: df and lmm
# Create data frame for T_weissflogii in Zn (including control)
TW_Zn <- SEM_EDS %>%
  filter(Species == "T_weissflogii" & (Treatment == "Zn" | Treatment == "control")) %>%
  select(Species, CellID, PointID, Treatment, Concentration, ZnSi_ratio) %>%
  arrange(Treatment, Concentration, CellID, PointID)
hist(TW_Zn$ZnSi_ratio)

# sqrt transform
TW_Zn$sqrt_ZnSi_ratio <- sqrt(TW_Zn$ZnSi_ratio)

#check transformed distribution
hist(TW_Zn$sqrt_ZnSi_ratio)

# use log-transformed response for model
library(lmerTest)
TW_Zn$Concentration <- factor(TW_Zn$Concentration)
model_TW_Zn <- lmer(sqrt_ZnSi_ratio ~ Concentration + (1|CellID), data = TW_Zn)
summary(model_TW_Zn)

#check residuals
plot(model_TW_Zn) 
qqnorm(residuals(model_TW_Zn)) 
qqline(residuals(model_TW_Zn)) 

#likelihood convergence test
null_model_TW_Zn <- lmer(log_ZnSi_ratio ~ 1 + (1|CellID), data = TW_Zn)
anova(null_model_TW_Zn, model_TW_Zn)


#plot linear model for TW vs Zn
TW_Zn$Concentration <- as.numeric(as.character(TW_Zn$Concentration))
LM4<-ggplot(TW_Zn, aes(x = Concentration, y = sqrt_ZnSi_ratio)) +
  geom_point() +               
  geom_smooth(method = "lm",         
              se = TRUE,               
              color = "#0173B2",         
              formula = y ~ x) +     
  ylab(expression("√(Zn/Si ratio) in"~italic("T. weissflogii"))) +
  xlab("Concentration of zinc added (µM)") +
  theme_minimal()

plot(LM4)

#TP vs Al: df and lmm
# Create data frame for T_pseudonana in Al (including control)
TP_Al <- SEM_EDS %>%
  filter(Species == "T_pseudonana" & (Treatment == "Al" | Treatment == "control")) %>%
  select(Species, CellID, PointID, Treatment, Concentration, AlSi_ratio) %>%
  arrange(Treatment, Concentration, CellID, PointID)

str(TP_Al)
hist(TP_Al$AlSi_ratio)

#histogram is right-skewed check zero, negative values
summary(TP_Al$AlSi_ratio)
sum(TP_Al$AlSi_ratio <= 0, na.rm = TRUE)

# square-root transform
TP_Al$sqrt_AlSi_ratio <- sqrt(TP_Al$AlSi_ratio)

#check transformed distribution
hist(TP_Al$sqrt_AlSi_ratio)

summary(TP_Al$sqrt_AlSi_ratio)
var(TP_Al$sqrt_AlSi_ratio, na.rm=TRUE)
sd_value <- sd(TP_Al$sqrt_AlSi_ratio, na.rm = TRUE)
print(sd_value)

# use log-transformed response for model
library(lmerTest)
TP_Al$Concentration <- factor(TP_Al$Concentration)
model_TP_Al <- lmer(sqrt_AlSi_ratio ~ Concentration + (1|CellID), data = TP_Al)
summary(model_TP_Al)

#check residuals
plot(model_TP_Al)  
qqnorm(residuals(model_TP_Al)) 
qqline(residuals(model_TP_Al)) 


#likelihood convergence test
null_model_TP_Al <- lmer(sqrt_AlSi_ratio ~ 1 + (1|CellID), data = TP_Al)
anova(null_model_TP_Al, model_TP_Al)

#plotting linear model of TP vs Al
TP_Al$Concentration <- as.numeric(as.character(TP_Al$Concentration))
LM1<-ggplot(TP_Al, aes(x = Concentration, y = sqrt_AlSi_ratio)) +
  geom_point() +                      
  geom_smooth(method = "lm",          
              se = TRUE,                
              color = "#E69F00",
              formula = y ~ x) +    
  ylab(expression("√(Al/Si ratio) in"~italic("T. pseudonana"))) +
  xlab("Concentration of Al added (µM)") +
  theme_minimal()

plot(LM1)


# Create data frame for T_pseudonana in Zn (including control)
TP_Zn <- SEM_EDS %>%
  filter(Species == "T_pseudonana" & (Treatment == "Zn" | Treatment == "control")) %>%
  select(Species, CellID, PointID, Treatment, Concentration, ZnSi_ratio) %>%
  arrange(Treatment, Concentration, CellID, PointID)

# square-root transform
TP_Zn$sqrt_ZnSi_ratio <- sqrt(TP_Zn$ZnSi_ratio)

#check transformed distribution
hist(TP_Zn$sqrt_ZnSi_ratio)

summary(TP_Zn$sqrt_ZnSi_ratio)
var(TP_Zn$sqrt_ZnSi_ratio, na.rm=TRUE)
sd_value <- sd(TP_Zn$sqrt_ZnSi_ratio, na.rm = TRUE)
print(sd_value)

# use sqrt transformed response for model
library(lmerTest)
TP_Zn$Concentration <- factor(TP_Zn$Concentration)
model_TP_Zn <- lmer(sqrt_ZnSi_ratio ~ Concentration + (1|CellID), data = TP_Zn)
summary(model_TP_Zn)

#check residuals
plot(model_TP_Zn)  
qqnorm(residuals(model_TP_Zn)) 
qqline(residuals(model_TP_Zn)) 


#likelihood convergence test
null_model_TP_Al <- lmer(sqrt_AlSi_ratio ~ 1 + (1|CellID), data = TP_Al)
anova(null_model_TP_Al, model_TP_Al)

#plotting linear model of TP vs Al
TP_Zn$Concentration <- as.numeric(as.character(TP_Zn$Concentration))
LM2<-ggplot(TP_Zn, aes(x = Concentration, y = sqrt_ZnSi_ratio)) +
  geom_point() +                        
  geom_smooth(method = "lm",            
              se = TRUE,                
              color = "#0173B2",       
              formula = y ~ x) +     
  ylab(expression("√(Zn/Si ratio) in"~italic("T. pseudonana"))) +
  xlab("Concentration of Zn added (µM)") +
  theme_minimal()

plot(LM2)

# Create data frame for C_pulvinata
CP <- SEM_EDS %>%
  filter(Species == "C_pulvinata" & (Treatment == "P")) %>%
  select(Species, CellID, PointID, Treatment, Concentration, P_mol) %>%
  arrange(Treatment, Concentration, CellID, PointID)

#histogram is right-skewed check zero, negative values
summary(CP$P_mol)
str(CP$P_mol)
CP$P_mol <- as.numeric(CP$P_mol)

library(ggplot2)
library(dplyr)

# Summarize the data: mean and standard error of X_mol grouped by Concentration
CP_summary <- CP %>%
  group_by(Concentration) %>%
  summarise(
    mean_P = mean(P_mol, na.rm = TRUE),
    se_P = sd(P_mol, na.rm = TRUE) / sqrt(n()),
    sd_P = sd(P_mol, na.rm=TRUE),
    median_P = median(P_mol, na.rm = TRUE)
  )
write.csv(CP_summary, "CP_EDS_summary.csv", row.names = FALSE)
#check distribution
hist(CP$P_mol)
CP$sqrt_P_mol <- sqrt(CP$P_mol)
hist(CP$sqrt_P_mol)

# use sqrt transformed response for model
model_CP <- lmer(sqrt_P_mol ~ Concentration + (1|CellID), data = CP)
summary(model_CP)
CP$Concentration <- factor(CP$Concentration)
summary(model_CP)

#check residuals
plot(model_CP) 
qqnorm(residuals(model_CP)) 
qqline(residuals(model_CP)) 
check_model(model_CP)
check_normality(model_CP) 
## NOTE: Non-normality of residuals detected (p =0.014)
## might be due to the outlier in histogram which is biologically plausible

#likelihood convergence test
null_model_CP <- lmer(sqrt_P_mol ~ 1 + (1|CellID), data = CP)
anova(null_model_CP, model_CP)

# Plot CP bar graph with error bars
ggplot(CP_summary, aes(x = factor(Concentration), y = mean_P)) +
  geom_bar(stat = "identity", width = 0.4, fill = "#669900") +
  geom_errorbar(aes(ymin = mean_P - se_P, ymax = mean_P + se_P), width = 0.1) +
  ylab(expression("% mol P in"~italic("C. pulvinata"))) +
  xlab("Concentration of phosphate added (µM)") +
  theme_minimal() +
  geom_text(data = data.frame(
    concentration = factor(c("20", "2000", "20000"), levels = c("0", "20", "2000", "20000")),
    mean_C = c(0.125, 0.295, 0.35), 
    label = c("ns", "***", "***")   
  ),
  aes(x = concentration, y = mean_C, label = label),
  inherit.aes = FALSE,
  vjust = -0.5,
  size = 4)
#plot linear model for CP
CP$Concentration <- as.numeric(as.character(CP$Concentration))
ggplot(CP, aes(x = Concentration, y = P_mol)) +
  geom_point() +                       
  geom_smooth(method = "lm",           
              se = TRUE,              
              color = "#669900",      
              formula = y ~ x) +     
  ylab(expression("%P mol in"~italic("C. pulvinata"))) +
  xlab("Concentration of phosphate added (µM)") +
  theme_minimal()
