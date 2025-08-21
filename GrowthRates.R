#remove all files on environment
rm(list = ls())

#check working directory
getwd()
#set working directory
setwd("/Users/nafisahalias/Documents/Imperial/Research project/R/Data")
#check working directory
getwd()
dailyrfu<-read.csv("Daily_RFU.csv", header=TRUE)

#load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)

#---GROWTH RATE FOR TW vs Al---#
# Filter for TW and Al only
TWAl <- dailyrfu %>%
  filter(Species == "T_weissflogii" & Treatment == "Al")

# plotting of growth graph
# Log transform RFU values
TWAl <- TWAl %>%
  mutate(log_rfu = log(RFU))
#Average replicates by flask and day
flask_TWAl <- TWAl %>%
  group_by(Species, Treatment, Concentration, Day, FlaskID) %>%
  summarize(
    log_rfu = mean(log_rfu, na.rm = TRUE),
    rfu_raw = mean(RFU, na.rm = TRUE),
    n_readings = n(),
    .groups = "drop"
  )
flask_TWAl <- flask_TWAl %>%
  mutate(
    concentration_factor = factor(Concentration),
    flask_id = paste(Concentration, FlaskID, sep = "_")
  )

# Plot 1: Individual flask growth curves
ggplot(flask_TWAl, aes(x = Day, y = log_rfu)) +
  geom_line(aes(group = flask_id, color = concentration_factor), 
            alpha = 0.7, size = 0.8) +
  geom_point(aes(color = concentration_factor), 
             alpha = 0.8, size = 1.5) +
  labs(
    title = "Individual Flask Growth Curves",
    subtitle = "TW with additional Al",
    x = "Day",
    y = "Log(RFU)",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(type = "qual", palette = "Set1")

# Calculate log averages and SD by concentration and day
concentration_averages <- flask_TWAl %>%
  group_by(Concentration, Day) %>%
  summarize(
    mean_log_rfu = mean(log_rfu, na.rm = TRUE),
    sd_log_rfu = sd(log_rfu, na.rm = TRUE),
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))

# Plot 2: Average growth curves with sd error bars
p1 <-ggplot(concentration_averages, 
                aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("T. weissflogii")*" in different concentrations of Al"),
    x = "Day",
    y = "Natural log of RFU",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +  
  scale_color_manual(values = c("#000000", "#DDAA33", "#BB5566"),  
                     labels = c("0 μM Al", "0.1 μM Al", "1.0 μM Al"))
print(p1)

# Calculation of gradient per flask (growth rate)
TWAl_growth_rates <- TWAl %>%
  filter(Day <= 6, RFU > 0) %>% 
  group_by(Concentration, FlaskID) %>%
  filter(RFU > 0) %>% 
  mutate(logRFU = log(RFU)) %>%
  do({
    model <- lm(logRFU ~ Day, data = .)
    tidy(model) %>%
      filter(term == "Day") %>%
      select(estimate)  # slope = growth rate
  }) %>%
  rename(growth_rate = estimate) %>%
  ungroup()

print(TWAl_growth_rates)

TWAl_growth_summary <- TWAl_growth_rates %>%
  group_by(Concentration) %>%
  summarise(
    mean_growth_rate = mean(growth_rate),
    se_growth_rate = sd(growth_rate) / sqrt(n()),
    sd_growth_rate = sd(growth_rate),
    median_growth_rate = median(growth_rate),
    n = n()
  )

print(TWAl_growth_summary)

library(lme4)
TWAl_lm <- lm(growth_rate ~ Concentration, data = TWAl_growth_rates)
summary(TWAl_lm)
par(mfrow = c(2, 2))
plot(TWAl_lm)


#plotting of linear model
TWAl_lmplot <- ggplot(TWAl_growth_rates, aes(x = Concentration, y = growth_rate)) +
  geom_point() +                        
  geom_smooth(method = "lm",            
              se = TRUE,                
              color = "#E69F00",           
              formula = y ~ x) +     
  ylab(expression("Growth rates of"~italic("T. weissflogii"))) +
  xlab("Concentration of Al added (µM)") +
  theme_minimal()

#---GROWTH RATE FOR TW vs Zn---#
# Filter for TW and Al only
TWZn <- dailyrfu %>%
  filter(Species == "T_weissflogii" & Treatment == "Zn")
#Log transform RFU values
TWZn <- TWZn %>%
  mutate(log_rfu = log(RFU))
#Average replicates by flask and day
flask_TWZn <- TWZn %>%
  group_by(Species, Treatment, Concentration, Day, FlaskID) %>%
  summarize(
    log_rfu = mean(log_rfu, na.rm = TRUE),
    rfu_raw = mean(RFU, na.rm = TRUE),
    n_readings = n(),
    .groups = "drop"
  )
flask_TWZn <- flask_TWZn %>%
  mutate(
    concentration_factor = factor(Concentration),
    flask_id = paste(Concentration, FlaskID, sep = "_")
  )

# Plot 1: Individual flask growth curves
ggplot(flask_TWZn, aes(x = Day, y = log_rfu)) +
  geom_line(aes(group = flask_id, color = concentration_factor), 
            alpha = 0.7, size = 0.8) +
  geom_point(aes(color = concentration_factor), 
             alpha = 0.8, size = 1.5) +
  labs(
    title = "Individual Flask Growth Curves",
    subtitle = "TW with additional Zn",
    x = "Day",
    y = "Log(RFU)",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(type = "qual", palette = "Set1")

# Calculate log averages and SE by concentration and day
TWZn_concentration_averages <- flask_TWZn %>%
  group_by(Concentration, Day) %>%
  summarize(
    mean_log_rfu = mean(log_rfu, na.rm = TRUE),
    sd_log_rfu = sd(log_rfu, na.rm = TRUE),
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))

# Plot 2: Average growth curves with error bars
p2 <- ggplot(TWZn_concentration_averages, 
                    aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("T. weissflogii")*" in different concentrations of Zn"),
    x = "Day",
    y = "Natural log of RFU",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +  
  scale_color_manual(values = c("#000000", "#DDAA33", "#BB5566"),  
                     labels = c("0 μM Zn", "8.0 μM Zn", "24 μM Zn"))  

# Step 1: Calculate slopes per flask
TWZn_growth_rates <- TWZn %>%
  filter(Day <= 6, RFU > 0) %>% 
  group_by(Concentration, FlaskID) %>%
  filter(RFU > 0) %>%  
  mutate(logRFU = log(RFU)) %>%
  do({
    model <- lm(logRFU ~ Day, data = .)
    tidy(model) %>%
      filter(term == "Day") %>%
      select(estimate)  # slope = growth rate
  }) %>%
  rename(growth_rate = estimate) %>%
  ungroup()

print(TWZn_growth_rates)

TWZn_growth_summary <- TWZn_growth_rates %>%
  group_by(Concentration) %>%
  summarise(
    mean_growth_rate = mean(growth_rate),
    se_growth_rate = sd(growth_rate) / sqrt(n()),
    sd_growth_rate = sd(growth_rate),
    median_growth_rate = median(growth_rate),
    n = n()
  )

print(TWZn_growth_summary)

#linear model
TWZn_lm <- lm(growth_rate ~ Concentration, data = TWZn_growth_rates)
summary(TWZn_lm)
plot(TWZn_lm)

TWZn_lmplot <- ggplot(TWZn_growth_rates, aes(x = Concentration, y = growth_rate)) +
  geom_point() +                       
  geom_smooth(method = "lm",           
              se = TRUE,               
              color = "#0173B2",          
              formula = y ~ x) +     
  ylab(expression("Growth rates of"~italic("T. weissflogii"))) +
  xlab("Concentration of Zn added (µM)") +
  theme_minimal()

#---GROWTH RATE FOR TP vs Al---#
# Filter for TP and Al only
TPAl <- dailyrfu %>%
  filter(Species == "T_pseudonana" & Treatment == "Al")

#Log transform RFU values
TPAl <- TPAl %>%
  mutate(log_rfu = log(RFU))
#Average replicates by flask and day
flask_TPAl <- TPAl %>%
  group_by(Species, Treatment, Concentration, Day, FlaskID) %>%
  summarize(
    log_rfu = mean(log_rfu, na.rm = TRUE),
    rfu_raw = mean(RFU, na.rm = TRUE),
    n_readings = n(),
    .groups = "drop"
  )
flask_TPAl <- flask_TPAl %>%
  mutate(
    concentration_factor = factor(Concentration),
    flask_id = paste(Concentration, FlaskID, sep = "_")
  )

# Plot 1: Individual flask growth curves
ggplot(flask_TPAl, aes(x = Day, y = log_rfu)) +
  geom_line(aes(group = flask_id, color = concentration_factor), 
            alpha = 0.7, size = 0.8) +
  geom_point(aes(color = concentration_factor), 
             alpha = 0.8, size = 1.5) +
  labs(
    title = "Individual Flask Growth Curves",
    subtitle = "TP with additional Al",
    x = "Day",
    y = "Log(RFU)",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(type = "qual", palette = "Set1")

# Calculate log averages and SE by concentration and day
TPAl_concentration_averages <- flask_TPAl %>%
  group_by(Concentration, Day) %>%
  summarize(
    mean_log_rfu = mean(log_rfu, na.rm = TRUE),
    sd_log_rfu = sd(log_rfu, na.rm = TRUE),
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))

# Plot 2: Average growth curves with error bars
p3 <- ggplot(TPAl_concentration_averages, 
                    aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("T. pseudonana")*" in different concentrations of Al"),
    x = "Day",
    y = "Natural log of RFU",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +  
  scale_color_manual(values = c("#000000", "#DDAA33", "#BB5566"),  
                     labels = c("0 μM Al", "0.1 μM Al", "1.0 μM Al"))  


# Step 1: Calculate slopes per flask
TPAl_growth_rates <- TPAl %>%
  filter(Day >= 1 & Day <= 6, RFU > 0) %>% 
  group_by(Concentration, FlaskID) %>%
  filter(RFU > 0) %>%  # Remove zero or negative RFU before log
  mutate(logRFU = log(RFU)) %>%
  do({
    model <- lm(logRFU ~ Day, data = .)
    tidy(model) %>%
      filter(term == "Day") %>%
      select(estimate)  # slope = growth rate
  }) %>%
  rename(growth_rate = estimate) %>%
  ungroup()

print(TPAl_growth_rates)

TPAl_growth_summary <- TPAl_growth_rates %>%
  group_by(Concentration) %>%
  summarise(
    mean_growth_rate = mean(growth_rate),
    se_growth_rate = sd(growth_rate) / sqrt(n()),
    sd_growth_rate = sd(growth_rate),
    median_growth_rate = median(growth_rate),
    n = n()
  )

print(TPAl_growth_summary)

library(lme4)
TPAl_lm <- lm(growth_rate ~ Concentration, data = TPAl_growth_rates)
summary(TPAl_lm)
plot(TPAl_lm)

TPAl_lmplot <- ggplot(TPAl_lm, aes(x = Concentration, y = growth_rate)) +
  geom_point() +                        
  geom_smooth(method = "lm",          
              se = TRUE,                
              color = "#E69F00",     
              formula = y ~ x) +    
  ylab(expression("Growth rates of"~italic("T. pseudonana"))) +
  xlab("Concentration of Al added (µM)") +
  theme_minimal()

#---GROWTH RATE FOR TP vs Zn---#
# Filter for TP and Al only
TPZn <- dailyrfu %>%
  filter(Species == "T_pseudonana" & Treatment == "Zn")

#Log transform RFU values
TPZn <- TPZn %>%
  mutate(log_rfu = log(RFU))
#Average replicates by flask and day
flask_TPZn <- TPZn %>%
  group_by(Species, Treatment, Concentration, Day, FlaskID) %>%
  summarize(
    log_rfu = mean(log_rfu, na.rm = TRUE),
    rfu_raw = mean(RFU, na.rm = TRUE),
    n_readings = n(),
    .groups = "drop"
  )
flask_TPZn <- flask_TPZn %>%
  mutate(
    concentration_factor = factor(Concentration),
    flask_id = paste(Concentration, FlaskID, sep = "_")
  )

# Plot 1: Individual flask growth curves
TPZnplot1 <- ggplot(flask_TPZn, aes(x = Day, y = log_rfu)) +
  geom_line(aes(group = flask_id, color = concentration_factor), 
            alpha = 0.7, size = 0.8) +
  geom_point(aes(color = concentration_factor), 
             alpha = 0.8, size = 1.5) +
  labs(
    title = "Individual Flask Growth Curves",
    subtitle = "TP with additional Zn",
    x = "Day",
    y = "Log(RFU)",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(type = "qual", palette = "Set1")

print(TPZnplot1)

# Calculate log averages and SE by concentration and day
TPZn_concentration_averages <- flask_TPZn %>%
  group_by(Concentration, Day) %>%
  summarize(
    mean_log_rfu = mean(log_rfu, na.rm = TRUE),
    sd_log_rfu = sd(log_rfu, na.rm = TRUE),
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))

# Plot 2: Average growth curves with error bars
p4 <- ggplot(TPZn_concentration_averages, 
                    aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("T. pseudonana")*" in different concentrations of Zn"),
    x = "Day",
    y = "Natural log of RFU",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) + 
  scale_color_manual(values = c("#000000", "#DDAA33", "#BB5566"), 
                     labels = c("0 μM Zn", "8.0 μM Zn", "24 μM Zn"))  


# Step 1: Calculate slopes per flask
TPZn_growth_rates <- TPZn %>%
  filter(Day <= 6, RFU > 0) %>% 
  group_by(Concentration, FlaskID) %>%
  filter(RFU > 0) %>%  # Remove zero or negative RFU before log
  mutate(logRFU = log(RFU)) %>%
  do({
    model <- lm(logRFU ~ Day, data = .)
    tidy(model) %>%
      filter(term == "Day") %>%
      select(estimate)  # slope = growth rate
  }) %>%
  rename(growth_rate = estimate) %>%
  ungroup()

print(TPZn_growth_rates)

TPZn_growth_summary <- TPZn_growth_rates %>%
  group_by(Concentration) %>%
  summarise(
    mean_growth_rate = mean(growth_rate),
    se_growth_rate = sd(growth_rate) / sqrt(n()),
    sd_growth_rate = sd(growth_rate),
    median_growth_rate = median(growth_rate),
    n = n()
  )

print(TPZn_growth_summary)

#linear model
TPZn_lm <- lm(growth_rate ~ Concentration, data = TPZn_growth_rates)
summary(TPZn_lm)
plot(TPZn_lm)

TPZn_lmplot <- ggplot(TPZn_growth_rates, aes(x = Concentration, y = growth_rate)) +
  geom_point() +                       
  geom_smooth(method = "lm",            
              se = TRUE,                
              color = "#0173B2",           
              formula = y ~ x) +     
  ylab(expression("Growth rates of"~italic("T. pseudonana"))) +
  xlab("Concentration of Zn added (µM)") +
  theme_minimal()

#---GROWTH RATE FOR CP vs Phosphate---#
# Filter for CP and P only
CP <- dailyrfu %>%
  filter(Species == "C_pulvinata" & Treatment == "P")

#Log transform RFU values
CP <- CP %>%
  mutate(log_rfu = log(RFU))

#Average replicates by flask and day
flask_CP <- CP %>%
  group_by(Species, Treatment, Concentration, Day, FlaskID) %>%
  summarize(
    log_rfu = mean(log_rfu, na.rm = TRUE),
    rfu_raw = mean(RFU, na.rm = TRUE),
    n_readings = n(),
    .groups = "drop"
  )
flask_CP <- flask_CP %>%
  mutate(
    concentration_factor = factor(Concentration),
    flask_id = paste(Concentration, FlaskID, sep = "_")
  )

# Plot 1: Individual flask growth curves
CPplot1 <- ggplot(flask_CP, aes(x = Day, y = log_rfu)) +
  geom_line(aes(group = flask_id, color = concentration_factor), 
            alpha = 0.7, size = 0.8) +
  geom_point(aes(color = concentration_factor), 
             alpha = 0.8, size = 1.5) +
  labs(
    title = "Individual Flask Growth Curves",
    subtitle = "CP with different concentrations of P",
    x = "Day",
    y = "Log(RFU)",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(type = "qual", palette = "Set1")

print(CPplot1)

# Calculate log averages and SE by concentration and day
concentration_averages <- flask_CP %>%
  group_by(Concentration, Day) %>%
  summarize(
    mean_log_rfu = mean(log_rfu, na.rm = TRUE),
    sd_log_rfu = sd(log_rfu, na.rm = TRUE),
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))


# Plot 2: Average growth curves with error bars
CPplot2 <- ggplot(concentration_averages, 
                  aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("C. pulvinata")*" in different concentrations of phosphate"),
    x = "Day",
    y = "Natural log of RFU",
    color = "Concentration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +
  scale_color_manual(values = c("0" = "#000000", "20" = "#E69F00", "2000" = "#9467BD", "20000" = "#007A52"),
                     labels = c("0" = "0 μM phosphate", "20" = "20 μM phosphate", "2000" = "2000 μM phosphate", "20000" = "20000 μM phosphate"))
print(CPplot2)

# Step 1: Calculate slopes per flask
CP_growth_rates <- CP %>%
  filter(Day >= 2 & Day <= 12, RFU > 0) %>% 
  group_by(Concentration, FlaskID) %>%
  filter(RFU > 0) %>%  # Remove zero or negative RFU before log
  mutate(logRFU = log(RFU)) %>%
  do({
    model <- lm(logRFU ~ Day, data = .)
    tidy(model) %>%
      filter(term == "Day") %>%
      select(estimate)  # slope = growth rate
  }) %>%
  rename(growth_rate = estimate) %>%
  ungroup()

print(CP_growth_rates)

CP_growth_summary <- CP_growth_rates %>%
  group_by(Concentration) %>%
  summarise(
    mean_growth_rate = mean(growth_rate),
    se_growth_rate = sd(growth_rate) / sqrt(n()),
    sd_growth_rate = sd(growth_rate),
    median_growth_rate = median(growth_rate),
    n = n()
  )

print(CP_growth_summary)

#RESCALING
CP_growth_rates <- CP_growth_rates %>%
  mutate(Concentration_k = Concentration / 1000)
CP_growth_rates %>% count(FlaskID)


CP_lm <- lm(growth_rate ~ Concentration_k, data = CP_growth_rates)
summary(CP_lm)
plot(CP_lm)

ggplot(CP_growth_rates, aes(x = Concentration, y = growth_rate)) +
  geom_point() +                        
  geom_smooth(method = "lm",           
              se = TRUE,               
              color = "darkgreen",          
              formula = y ~ x) +     
  ylab(expression("Growth rates of"~italic("C.pulvinata"))) +
  xlab("Concentration of phosphate added (µM)") +
  theme_minimal()

