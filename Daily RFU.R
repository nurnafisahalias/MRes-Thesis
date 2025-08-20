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

# =======================================
#        TW vs Al
# =======================================
# Filter for TW and Al only
TWAl <- dailyrfu %>%
  filter(Species == "T_weissflogii" & Treatment == "Al")
#Log transform RFU values
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
plot1 <- ggplot(flask_TWAl, aes(x = Day, y = log_rfu)) +
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

print(plot1)

# Calculate log averages and SE by concentration and day
concentration_averages <- flask_TWAl %>%
  group_by(Concentration, Day) %>%
  summarize(
    mean_log_rfu = mean(log_rfu, na.rm = TRUE),
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))

# Plot 2: Average growth curves with error bars
plot2 <- ggplot(concentration_averages, 
                aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("T. weissflogii")*" in different added concentrations of Al"),
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
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +  # NEW: Integer x-axis breaks
  scale_color_manual(values = c("#004488", "#DDAA33", "#BB5566"),  # CHANGE 1: Single color for all lines
                     labels = c("0 μM Al added", "0.1 μM Al added", "1.0 μM Al added"))  # CHANGE 3: New legend labels
  
print(plot2)

# Linear Mixed Effects Model for Growth Rate Analysis
#load libraries
library(lme4)
library(lmerTest)
library(performance)
library(sjPlot)
library(MuMIn)
library(ggplot2)
library(dplyr)

# Histogram of response variable
hist(flask_TWAl$log_rfu, main = "Distribution of log(RFU)", 
     xlab = "log(RFU)", breaks = 20)

#simple linear model worked, try nlme package for linear mixed model
library(nlme)
#fitting TW vs Al into linear mixed effects model using nlme package
# Model 1: Basic model - same growth rate for all concentrations
# Assumes concentration has NO effect on growth
# Hypothesis: All concentrations have the same growth rate
TWAl_model1 <- lme(log_rfu ~ Day, 
              random = ~ 1 | FlaskID, 
              data = flask_TWAl,
              method = "REML")

# Model 2: Different intercepts by concentration  
# Assumes concentration affects STARTING POINT but not growth rate
# Hypothesis: Concentrations start at different levels, same growth rate
TWAl_model2 <- lme(log_rfu ~ Day + Concentration, 
              random = ~ 1 | FlaskID, 
              data = flask_TWAl,
              method = "REML")

# Model 3: Different growth rates by concentration (interaction)
# Assumes concentration affects BOTH starting point AND growth rate
# Hypothesis: Concentrations have different starting points AND growth rates
TWAl_model3 <- lme(log_rfu ~ Day * Concentration, 
              random = ~ 1 | FlaskID, 
              data = flask_TWAl,
              method = "REML")

# Check that all models converged
models <- list(TWAl_model1, TWAl_model2, TWAl_model3)
model_names <- c("Model 1", "Model 2", "Model 3")

for(i in 1:length(models)) {
  cat("\n", model_names[i], "convergence: ")
  # nlme doesn't have the same convergence warnings as lme4
  # If the model ran without error, it converged
  cat("SUCCESS\n")
}

# Compare nested models using likelihood ratio tests
#model 1 vs model 2
TWAl_LRT_1vs2 <- anova (TWAl_model1, TWAl_model2)
print(TWAl_LRT_1vs2)

#model 2 vs model 3
TWAl_LRT_2vs3 <- anova (TWAl_model2, TWAl_model3)
print(TWAl_LRT_2vs3)

#model 1 vs model 3
TWAl_LRT_1vs3 <- anova (TWAl_model1, TWAl_model3)
print(TWAl_LRT_1vs3)

# Compare AIC values 
aic_values <- c(AIC(TWAl_model1), AIC(TWAl_model2), AIC(TWAl_model3))
names(aic_values) <- model_names
print(aic_values)

# final model diagnostic (model 2)
# Residual plots
plot(TWAl_model2, main = paste("Residuals vs Fitted -", "Model 2"))

# Q-Q plots
qqnorm(TWAl_model2, ~ resid(.), main = paste("Residuals Q-Q -", "Model 2"))
qqnorm(TWAl_model2, ~ ranef(.), main = paste("Random Effects Q-Q -", "Model 2"))

# Check residuals by concentration
plot(TWAl_model2, resid(.) ~ fitted(.) | Concentration, 
     main = paste("Residuals by Concentration -", "Model 2"))

#VIF analysis
library(car)
# Extract the fixed effects design matrix from nlme model
TWAl_model2_lm <- lm(log_rfu ~ Day + Concentration, data = flask_TWAl)

# Calculate VIF
vif_values <- vif(TWAl_model2_lm)
print(vif_values)

summary(TWAl_model2)

# Extract variance components
variance_table <- VarCorr(TWAl_model2)
print(variance_table)

# Get the numbers we need
random_variance <- as.numeric(variance_table[1, "Variance"])    # Between flasks
residual_variance <- as.numeric(variance_table[2, "Variance"])  # Within flasks

# Calculate ICC (percentage due to flask differences)
icc <- random_variance / (random_variance + residual_variance)

print("\nResults:")
print(paste("Flask-to-flask variance:", round(random_variance, 4)))
print(paste("Within-flask variance:", round(residual_variance, 4)))
print(paste("ICC (flask effect):", round(icc, 3)))
print(paste("Percentage due to flask differences:", round(icc * 100, 1), "%"))
print(paste("Percentage due to within-flask variation:", round((1-icc) * 100, 1), "%"))

summary(TWAl_model2)$tTable

TWAl$Concentration <- as.numeric(as.character(TWAl$Concentration))
ggplot(TWAl, aes(x = Concentration, y = log_rfu)) +
  geom_point() +                        # scatter points
  geom_smooth(method = "lm",            # fit linear model
              se = TRUE,                # show confidence interval as shaded area
              color = "#E69F00",           # line color
              formula = y ~ x) +     # shaded CI color
  ylab(expression("Natural log of RFU of"~italic("T. weissflogii"))) +
  xlab("Concentration of Al added (µM)") +
  theme_minimal()

# =======================================
#        TW vs Zn
# =======================================
# Filter for TW and Zn only
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
TWZnplot1 <- ggplot(flask_TWZn, aes(x = Day, y = log_rfu)) +
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

print(TWZnplot1)

# Calculate log averages and SE by concentration and day
TWZn_concentration_averages <- flask_TWZn %>%
  group_by(Concentration, Day) %>%
  summarize(
    mean_log_rfu = mean(log_rfu, na.rm = TRUE),
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))

# Plot 2: Average growth curves with error bars
TWZnplot2 <- ggplot(TWZn_concentration_averages, 
                    aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("T. weissflogii")*" in different added concentrations of Zn"),
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
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +  # NEW: Integer x-axis breaks
  scale_color_manual(values = c("#004488", "#DDAA33", "#BB5566"),  # CHANGE 1: Single color for all lines
                     labels = c("0 μM Zn added", "8.0 μM Al added", "24 μM Al added"))  # CHANGE 3: New legend labels

print(TWZnplot2)

#Linear mixed effects model using nlme
# Histogram of response variable
hist(flask_TWZn$log_rfu, main = "Distribution of log(RFU)", 
     xlab = "log(RFU)", breaks = 20)

# Model 1: Basic model - same growth rate for all concentrations
# Assumes concentration has NO effect on growth
# Hypothesis: All concentrations have the same growth rate
TWZn_model1 <- lme(log_rfu ~ Day, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TWZn,
                   method = "REML")

# Model 2: Different intercepts by concentration  
# Assumes concentration affects STARTING POINT but not growth rate
# Hypothesis: Concentrations start at different levels, same growth rate
TWZn_model2 <- lme(log_rfu ~ Day + Concentration, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TWZn,
                   method = "REML")

# Model 3: Different growth rates by concentration (interaction)
# Assumes concentration affects BOTH starting point AND growth rate
# Hypothesis: Concentrations have different starting points AND growth rates
TWZn_model3 <- lme(log_rfu ~ Day * Concentration, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TWZn,
                   method = "REML")

# Check that all models converged
models <- list(TWZn_model1, TWZn_model2, TWZn_model3)
model_names <- c("Model 1", "Model 2", "Model 3")

for(i in 1:length(models)) {
  cat("\n", model_names[i], "convergence: ")
  # nlme doesn't have the same convergence warnings as lme4
  # If the model ran without error, it converged
  cat("SUCCESS\n")
}

# Compare nested models using likelihood ratio tests
#model 1 vs model 2
TWZn_LRT_1vs2 <- anova (TWZn_model1, TWZn_model2)
print(TWZn_LRT_1vs2)

#model 2 vs model 3
TWZn_LRT_2vs3 <- anova (TWZn_model2, TWZn_model3)
print(TWZn_LRT_2vs3)

#model 1 vs model 3
TWZn_LRT_1vs3 <- anova (TWZn_model1, TWZn_model3)
print(TWZn_LRT_1vs3)

# Compare AIC values 
aic_values <- c(AIC(TWZn_model1), AIC(TWZn_model2), AIC(TWZn_model3))
names(aic_values) <- model_names
print(aic_values)

# final model diagnostic (model 3)
# Residual plots
plot(TWZn_model3, main = paste("Residuals vs Fitted -", "Model 3"))

# Q-Q plots
qqnorm(TWZn_model3, ~ resid(.), main = paste("Residuals Q-Q -", "Model 3"))
qqnorm(TWZn_model3, ~ ranef(.), main = paste("Random Effects Q-Q -", "Model 3"))

# Check residuals by concentration
plot(TWZn_model3, resid(.) ~ fitted(.) | Concentration, 
     main = paste("Residuals by Concentration -", "Model 3"))

#VIF analysis
library(car)
# Extract the fixed effects design matrix from nlme model
TWZn_model3_lm <- lm(log_rfu ~ Day + Concentration, data = flask_TWZn)

# Calculate VIF
vif_values <- vif(TWZn_model3_lm)
print(vif_values)

summary(TWZn_model3)

# Extract variance components
variance_table <- VarCorr(TWZn_model3)
print(variance_table)

# Get the numbers we need
random_variance <- as.numeric(variance_table[1, "Variance"])    # Between flasks
residual_variance <- as.numeric(variance_table[2, "Variance"])  # Within flasks

# Calculate ICC (percentage due to flask differences)
icc <- random_variance / (random_variance + residual_variance)

print("\nResults:")
print(paste("Flask-to-flask variance:", round(random_variance, 4)))
print(paste("Within-flask variance:", round(residual_variance, 4)))
print(paste("ICC (flask effect):", round(icc, 3)))
print(paste("Percentage due to flask differences:", round(icc * 100, 1), "%"))
print(paste("Percentage due to within-flask variation:", round((1-icc) * 100, 1), "%"))


TWZn$Concentration <- as.numeric(as.character(TWZn$Concentration))
ggplot(TWZn, aes(x = Concentration, y = log_rfu)) +
  geom_point() +                        # scatter points
  geom_smooth(method = "lm",            # fit linear model
              se = TRUE,                # show confidence interval as shaded area
              color = "#0173B2",           # line color
              formula = y ~ x) +     # shaded CI color
  ylab(expression("Natural log of RFU of"~italic("T. weissflogii"))) +
  xlab("Concentration of Zn added (µM)") +
  theme_minimal()


# =======================================
#        TP vs Al
# =======================================
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
TPAlplot1 <- ggplot(flask_TPAl, aes(x = Day, y = log_rfu)) +
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

print(TPAlplot1)

# Calculate log averages and SE by concentration and day
TPAl_concentration_averages <- flask_TPAl %>%
  group_by(Concentration, Day) %>%
  summarize(
    mean_log_rfu = mean(log_rfu, na.rm = TRUE),
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))

# Plot 2: Average growth curves with error bars
TPAlplot2 <- ggplot(TPAl_concentration_averages, 
                aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("T. pseudonana")*" in different added concentrations of Al"),
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
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +  # NEW: Integer x-axis breaks
  scale_color_manual(values = c("#004488", "#DDAA33", "#BB5566"),  # CHANGE 1: Single color for all lines
                     labels = c("0 μM Al added", "0.1 μM Al added", "1.0 μM Al added"))  # CHANGE 3: New legend labels


print(TPAlplot2)

# Linear Mixed Effects Model for Growth Rate Analysis
library(nlme)
result <- tryCatch({
  nlme_model <- lme(log_rfu ~ Day * Concentration, 
                    random = ~ 1 | FlaskID, 
                    data = flask_TPAl)
  print("nlme SUCCESS!")
  summary(nlme_model)
}, error = function(e) {
  print(paste("nlme also failed:", e$message))
})

#fitting TW vs Al into linear mixed effects model using nlme package
# Model 1: Basic model - same growth rate for all concentrations
# Assumes concentration has NO effect on growth
# Hypothesis: All concentrations have the same growth rate
TPAl_model1 <- lme(log_rfu ~ Day, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TPAl,
                   method = "REML")

# Model 2: Different intercepts by concentration  
# Assumes concentration affects STARTING POINT but not growth rate
# Hypothesis: Concentrations start at different levels, same growth rate
TPAl_model2 <- lme(log_rfu ~ Day + Concentration, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TPAl,
                   method = "REML")

# Model 3: Different growth rates by concentration (interaction)
# Assumes concentration affects BOTH starting point AND growth rate
# Hypothesis: Concentrations have different starting points AND growth rates
TPAl_model3 <- lme(log_rfu ~ Day * Concentration, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TPAl,
                   method = "REML")

# Check that all models converged
models <- list(TPAl_model1, TPAl_model2, TPAl_model3)
model_names <- c("Model 1", "Model 2", "Model 3")

for(i in 1:length(models)) {
  cat("\n", model_names[i], "convergence: ")
  # nlme doesn't have the same convergence warnings as lme4
  # If the model ran without error, it converged
  cat("SUCCESS\n")
}

# Compare nested models using likelihood ratio tests
#model 1 vs model 2
TPAl_LRT_1vs2 <- anova (TPAl_model1, TPAl_model2)
print(TPAl_LRT_1vs2)

#model 2 vs model 3
TPAl_LRT_2vs3 <- anova (TPAl_model2, TPAl_model3)
print(TPAl_LRT_2vs3)

#model 1 vs model 3
TPAl_LRT_1vs3 <- anova (TPAl_model1, TPAl_model3)
print(TPAl_LRT_1vs3)

# Compare AIC values 
aic_values <- c(AIC(TPAl_model1), AIC(TPAl_model2), AIC(TPAl_model3))
names(aic_values) <- model_names
print(aic_values)

# final model diagnostic (model 2)
# Residual plots
plot(TPAl_model2, main = paste("Residuals vs Fitted -", "Model 2"))

# Residuals Q-Q plot
residuals_vals <- residuals(TPAl_model2)
qqnorm(residuals_vals, main = paste("Residuals Q-Q -", "Model 2"))
qqline(residuals_vals, col = "red", lwd = 2)

# Check residuals by concentration
plot(TPAl_model2, resid(.) ~ fitted(.) | Concentration, 
     main = paste("Residuals by Concentration -", "Model 2"))

#VIF analysis
library(car)
# Extract the fixed effects design matrix from nlme model
TPAl_model2_lm <- lm(log_rfu ~ Day + Concentration, data = flask_TPAl)

# Calculate VIF
vif_values <- vif(TPAl_model2_lm)
print(vif_values)

summary(TPAl_model2)

# Extract variance components
variance_table <- VarCorr(TPAl_model2)
print(variance_table)

# Get the numbers we need
random_variance <- as.numeric(variance_table[1, "Variance"])    # Between flasks
residual_variance <- as.numeric(variance_table[2, "Variance"])  # Within flasks

# Calculate ICC (percentage due to flask differences)
icc <- random_variance / (random_variance + residual_variance)

print("\nResults:")
print(paste("Flask-to-flask variance:", round(random_variance, 4)))
print(paste("Within-flask variance:", round(residual_variance, 4)))
print(paste("ICC (flask effect):", round(icc, 3)))
print(paste("Percentage due to flask differences:", round(icc * 100, 1), "%"))
print(paste("Percentage due to within-flask variation:", round((1-icc) * 100, 1), "%"))


TPAl$Concentration <- as.numeric(as.character(TPAl$Concentration))
ggplot(TPAl, aes(x = Concentration, y = log_rfu)) +
  geom_point() +                        # scatter points
  geom_smooth(method = "lm",            # fit linear model
              se = TRUE,                # show confidence interval as shaded area
              color = "#E69F00",           # line color
              formula = y ~ x) +     # shaded CI color
  ylab(expression("Natural log of RFU of"~italic("T. pseudonana"))) +
  xlab("Concentration of Al added (µM)") +
  theme_minimal()

# =======================================
#        TP vs Zn
# =======================================
# Filter for TW and Al only
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
    se_log_rfu = sd(log_rfu, na.rm = TRUE) / sqrt(n()),
    n_flasks = n(),
    .groups = "drop"
  ) %>%
  mutate(concentration_factor = factor(Concentration))

# Plot 2: Average growth curves with error bars
TPZnplot2 <- ggplot(TPZn_concentration_averages, 
                    aes(x = Day, y = mean_log_rfu,color = concentration_factor)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_rfu - se_log_rfu, 
                    ymax = mean_log_rfu + se_log_rfu),
                width = 0.2, size = 0.8) +
  labs(
    title = expression("Growth of "*italic("T. pseudonana")*" in different added concentrations of Zn"),
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
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +  # NEW: Integer x-axis breaks
  scale_color_manual(values = c("#004488", "#DDAA33", "#BB5566"),  # CHANGE 1: Single color for all lines
                     labels = c("0 μM Zn added", "8.0 μM Zn added", "24 μM Zn added"))  # CHANGE 3: New legend labels

print(TPZnplot2)

#Linear mixed effects model using nlme
# Histogram of response variable
hist(flask_TPZn$log_rfu, main = "Distribution of log(RFU)", 
     xlab = "log(RFU)", breaks = 20)

# Model 1: Basic model - same growth rate for all concentrations
# Assumes concentration has NO effect on growth
# Hypothesis: All concentrations have the same growth rate
TPZn_model1 <- lme(log_rfu ~ Day, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TPZn,
                   method = "REML")

# Model 2: Different intercepts by concentration  
# Assumes concentration affects STARTING POINT but not growth rate
# Hypothesis: Concentrations start at different levels, same growth rate
TPZn_model2 <- lme(log_rfu ~ Day + Concentration, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TPZn,
                   method = "REML")

# Model 3: Different growth rates by concentration (interaction)
# Assumes concentration affects BOTH starting point AND growth rate
# Hypothesis: Concentrations have different starting points AND growth rates
TPZn_model3 <- lme(log_rfu ~ Day * Concentration, 
                   random = ~ 1 | FlaskID, 
                   data = flask_TPZn,
                   method = "REML")

# Check that all models converged
models <- list(TPZn_model1, TPZn_model2, TPZn_model3)
model_names <- c("Model 1", "Model 2", "Model 3")

for(i in 1:length(models)) {
  cat("\n", model_names[i], "convergence: ")
  # nlme doesn't have the same convergence warnings as lme4
  # If the model ran without error, it converged
  cat("SUCCESS\n")
}

# Compare nested models using likelihood ratio tests
#model 1 vs model 2
TPZn_LRT_1vs2 <- anova (TPZn_model1, TPZn_model2)
print(TPZn_LRT_1vs2)

#model 2 vs model 3
TPZn_LRT_2vs3 <- anova (TPZn_model2, TPZn_model3)
print(TPZn_LRT_2vs3)

#model 1 vs model 3
TPZn_LRT_1vs3 <- anova (TPZn_model1, TPZn_model3)
print(TPZn_LRT_1vs3)

# Compare AIC values 
aic_values <- c(AIC(TPZn_model1), AIC(TPZn_model2), AIC(TPZn_model3))
names(aic_values) <- model_names
print(aic_values)

# final model diagnostic (model 1)
# Residual plots
plot(TPZn_model1, main = paste("Residuals vs Fitted -", "Model 1"))

# Residuals Q-Q plot
residuals_vals <- residuals(TPZn_model1)
qqnorm(residuals_vals, main = paste("Residuals Q-Q -", "Model 1"))
qqline(residuals_vals, col = "red", lwd = 2)

# Check residuals by concentration
plot(TWZn_model1, resid(.) ~ fitted(.) | Concentration, 
     main = paste("Residuals by Concentration -", "Model 1"))

#VIF analysis
library(car)
# Extract the fixed effects design matrix from nlme model
TPZn_model1_lm <- lm(log_rfu ~ Day + Concentration, data = flask_TPZn)

# Calculate VIF
vif_values <- vif(TPZn_model1_lm)
print(vif_values)

summary(TWAl_model2)

# Extract variance components
variance_table <- VarCorr(TWZn_model1)
print(variance_table)

# Get the numbers we need
random_variance <- as.numeric(variance_table[1, "Variance"])    # Between flasks
residual_variance <- as.numeric(variance_table[2, "Variance"])  # Within flasks

# Calculate ICC (percentage due to flask differences)
icc <- random_variance / (random_variance + residual_variance)

print("\nResults:")
print(paste("Flask-to-flask variance:", round(random_variance, 4)))
print(paste("Within-flask variance:", round(residual_variance, 4)))
print(paste("ICC (flask effect):", round(icc, 3)))
print(paste("Percentage due to flask differences:", round(icc * 100, 1), "%"))
print(paste("Percentage due to within-flask variation:", round((1-icc) * 100, 1), "%"))

TPZn$Concentration <- as.numeric(as.character(TPZn$Concentration))
ggplot(TPZn, aes(x = Concentration, y = log_rfu)) +
  geom_point() +                        # scatter points
  geom_smooth(method = "lm",            # fit linear model
              se = TRUE,                # show confidence interval as shaded area
              color = "#0173B2",           # line color
              formula = y ~ x) +     # shaded CI color
  ylab(expression("Natural log of RFU of"~italic("T. pseudonana"))) +
  xlab("Concentration of Zn added (µM)") +
  theme_minimal()

# =======================================
#        CP vs P
# =======================================
# Filter for CP only
CP <- dailyrfu %>%
  filter(Species == "C_pulvinata")
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
    title = expression("Growth of "*italic("C. pulvinata")*" in different added concentrations of P"),
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
  scale_x_continuous(breaks = seq(1, max(concentration_averages$Day), by = 1)) +  # NEW: Integer x-axis breaks
  scale_color_manual(values = c("#0173B2", "#DE8F05", "#029E73", "#CC78BC"),  # CHANGE 1: Single color for all lines
                     labels = c("0 μM P added", "20 μM P added", "2000 μM P added", "20000 μM P added"))  # CHANGE 3: New legend labels


print(CPplot2)

#Linear mixed effects model using nlme
# Histogram of response variable
hist(flask_CP$log_rfu, main = "Distribution of log(RFU)", 
     xlab = "log(RFU)", breaks = 20)

# Model 1: Basic model - same growth rate for all concentrations
# Assumes concentration has NO effect on growth
# Hypothesis: All concentrations have the same growth rate
CP_model1 <- lme(log_rfu ~ Day, 
                   random = ~ 1 | FlaskID, 
                   data = flask_CP,
                   method = "REML")

# Model 2: Different intercepts by concentration  
# Assumes concentration affects STARTING POINT but not growth rate
# Hypothesis: Concentrations start at different levels, same growth rate
CP_model2 <- lme(log_rfu ~ Day + Concentration, 
                   random = ~ 1 | FlaskID, 
                   data = flask_CP,
                   method = "REML")

# Model 3: Different growth rates by concentration (interaction)
# Assumes concentration affects BOTH starting point AND growth rate
# Hypothesis: Concentrations have different starting points AND growth rates
CP_model3 <- lme(log_rfu ~ Day * Concentration, 
                   random = ~ 1 | FlaskID, 
                   data = flask_CP,
                   method = "REML")

# Check that all models converged
models <- list(CP_model1, CP_model2, CP_model3)
model_names <- c("Model 1", "Model 2", "Model 3")

for(i in 1:length(models)) {
  cat("\n", model_names[i], "convergence: ")
  # nlme doesn't have the same convergence warnings as lme4
  # If the model ran without error, it converged
  cat("SUCCESS\n")
}

# Compare nested models using likelihood ratio tests
#model 1 vs model 2
CP_LRT_1vs2 <- anova (CP_model1, CP_model2)
print(CP_LRT_1vs2)

#model 2 vs model 3
CP_LRT_2vs3 <- anova (CP_model2, CP_model3)
print(CP_LRT_2vs3)

#model 1 vs model 3
CP_LRT_1vs3 <- anova (CP_model1, CP_model3)
print(CP_LRT_1vs3)

# Compare AIC values 
aic_values <- c(AIC(CP_model1), AIC(CP_model2), AIC(CP_model3))
names(aic_values) <- model_names
print(aic_values)

# final model diagnostic (model 3)
# Residual plots
plot(CP_model3, main = paste("Residuals vs Fitted -", "Model 3"))

# Residuals Q-Q plot
residuals_vals <- residuals(CP_model3)
qqnorm(residuals_vals, main = paste("Residuals Q-Q -", "Model 3"))
qqline(residuals_vals, col = "red", lwd = 2)

# Check residuals by concentration
plot(CP_model3, resid(.) ~ fitted(.) | Concentration, 
     main = paste("Residuals by Concentration -", "Model 3"))

#VIF analysis
library(car)
# Extract the fixed effects design matrix from nlme model
CP_model3_lm <- lm(log_rfu ~ Day + Concentration, data = flask_CP)

# Calculate VIF
vif_values <- vif(CP_model3_lm)
print(vif_values)

summary(CP_model3)

# Extract variance components
variance_table <- VarCorr(CP_model3)
print(variance_table)

# Get the numbers we need
random_variance <- as.numeric(variance_table[1, "Variance"])    # Between flasks
residual_variance <- as.numeric(variance_table[2, "Variance"])  # Within flasks

# Calculate ICC (percentage due to flask differences)
icc <- random_variance / (random_variance + residual_variance)

print("\nResults:")
print(paste("Flask-to-flask variance:", round(random_variance, 4)))
print(paste("Within-flask variance:", round(residual_variance, 4)))
print(paste("ICC (flask effect):", round(icc, 3)))
print(paste("Percentage due to flask differences:", round(icc * 100, 1), "%"))
print(paste("Percentage due to within-flask variation:", round((1-icc) * 100, 1), "%"))

CP$Concentration <- as.numeric(as.character(CP$Concentration))
ggplot(CP, aes(x = Concentration, y = log_rfu)) +
  geom_point() +                        # scatter points
  geom_smooth(method = "lm",            # fit linear model
              se = TRUE,                # show confidence interval as shaded area
              color = "darkgreen",           # line color
              formula = y ~ x) +     # shaded CI color
  ylab(expression("Natural log of RFU of"~italic("C. pulvinata"))) +
  xlab("Concentration of P added (µM)") +
  theme_minimal()
