# This script was compiled by Charles Knowlton in the NDC Lab for conducting preliminary analyses of READ Study 2
# Load libraries----
library(readr)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(tidyr)
library(reghelper)
library(ggeffects)
library(interactions)
library(patchwork)

# Initialize director folders----
derivatives_dir <- "/home/data/NDClab/analyses/read-study2-alpha/derivatives/"

# Load in data----
read_df <- read_csv(file.path(derivatives_dir, "read_long_s1_r1_30_01_2026_15_00_02.csv"))

# Aim 1: Mixed-effect model predicting Social Anxiety from Social Observation, Age, and ERN----
## Aggregate variables of interest----
df_summary <- read_df %>%
  group_by(sub, soc) %>%
  summarise(
    ERN = amplitude[acc == 0] - amplitude[acc == 1],
    spaic = first(spaic_scrdTotal_s1_r1_e1),
    age = first(age_m),
    spaip = first(spaip_scrdTotal_s1_r1_e1),
    acc = first(accuracy_score),
    staic = first(staic_scrdTotal_s1_r1_e1)
  )

## Convert soc to factor w/ 2 levels----
df_summary$soc <- factor(df_summary$soc, levels = c(0, 1), labels = c("alone", "social"))

## lm----
df_wide <- df_summary %>%
  pivot_wider(
    names_from = soc,
    values_from = ERN,
    names_prefix = "ERN_"
  )

## Omit NaNs----
participants_remove <- c() # Stored vector of participants with Nan for both conditions 
p_i <- 1 # Looping variable for index
for (i in 1:nrow(df_wide)){
  p_id = as.character(df_wide$sub[i])
  if ((is.na(df_wide$ERN_alone[i]) == TRUE) & (is.na(df_wide$ERN_social[i]) == FALSE)){
    print(paste("Participant's ERN alone (social is valid) is NAN, ID:", p_id))
  }
  if ((is.na(df_wide$ERN_social[i]) == TRUE) & (is.na(df_wide$ERN_alone[i]) == FALSE)){
    print(paste("Participant's ERN social (alone is valid) is NAN, ID:", p_id))
  }
  if ((is.na(df_wide$ERN_social[i]) == TRUE) & (is.na(df_wide$ERN_alone[i]) == TRUE)){
    print(paste("Participant has ERN NaN for both conditions, ID:", p_id))
    participants_remove[p_i] <- p_id
    p_i <- p_i + 1
  }
}

# Remove participants from long and wide data frames that contained no ERN data
df_wide <- df_wide %>% filter(!(sub %in% participants_remove))
df_summary <- df_summary %>% filter(!(sub %in% participants_remove))
# Remove participants with no SPAIP data
df_wide <- df_wide %>% filter(!is.na(spaip))
df_summary <- df_summary %>% filter(!is.na(spaip))
# Remove participants with no SPAIC data
df_wide <- df_wide %>% filter(!is.na(spaic))
df_summary <- df_summary %>% filter(!is.na(spaic))

## Testing removing participant > -3 SD of age mean
df_summary_age_subset <- df_summary |> subset(age > 150)

## Scaling variables
df_summary$age_s    <- scale(df_summary$age)
df_summary$spaic_s  <- scale(df_summary$spaic)
df_summary$spaip_s  <- scale(df_summary$spaip)
df_summary$acc_s    <- scale(df_summary$acc)

## Mixed-effect modeling----
model_lm_mixed_spaip <- lmer(ERN ~ age_s + soc + spaip_s + (1|sub), data = df_summary)
summary(model_lm_mixed_spaip)

model_lm_mixed_spaic <- lmer(ERN ~ age_s + soc + spaic_s + (1|sub), data = df_summary)
summary(model_lm_mixed_spaic)

## Test for interaction----
model_lm_mixed_spaic_int <- lmer(ERN ~ age_s * spaic_s * soc * acc_s + (1|sub), data = df_summary)
summary(model_lm_mixed_spaic_int)

model_lm_mixed_spaip_int <- lmer(ERN ~ age_s * spaip_s * soc * acc_s + (1|sub), data = df_summary)
summary(model_lm_mixed_spaip_int)

## Simple slopes to interpret interaction better

simple_slopes(model_lm_mixed_spaip_int)
simple_slopes(model_lm_mixed_spaic_int)

## Plotting interaction + orig data
# Getting predictions for the slopes
preds_spaic <- ggpredict(model_lm_mixed_spaic_int, terms = c("spaic_s", "age_s", "soc", "acc_s"))
preds_spaic$soc <- preds_spaic$facet
# Plotting
ggplot() +
  geom_point(data = df_summary, 
             aes(x = spaic_s, y = ERN, size = acc_s, color = acc_s),
             alpha = 0.3) +
  scale_
  geom_line(data = preds_spaic,
            aes(x = x, y = predicted, color = group),
            linewidth = 1) +
  geom_ribbon(data = preds_spaic, 
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), 
              alpha = 0.1) +
  facet_wrap(~soc) +
  labs(x = "Social Anxiety (SPAIC)", y = "Delta ERN", color = "Age (months)", fill = "Age (months)") +
  theme_minimal()

# Getting predictions
preds_spaip <- ggpredict(model_lm_mixed_spaip_int, terms = c("spaic_s", "age_s", "soc", "acc_s"))
preds_spaip$soc <- preds_spaip$facet
# Plotting
ggplot() +
  geom_point(data = df_summary, 
             aes(x = spaip_s, y = ERN_s),
             alpha = 0.3) +
  geom_line(data = preds_spaip, 
            aes(x = x, y = predicted, color = group),
            linewidth = 1) +
  geom_ribbon(data = preds_spaip, 
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), 
              alpha = 0.1) +
  facet_wrap(~soc) +
  labs(x = "Social Anxiety (SPAIP)", y = "Delta ERN", color = "Age (months)", fill = "Age (months)") +
  theme_minimal()

## Test for interaction in STAIC----
## Scaling variables
read_df$age_s    <- scale(read_df$age_m)
read_df$spaic_s  <- scale(read_df$spaic_scrdTotal_s1_r1_e1)
read_df$spaip_s  <- scale(read_df$spaip_scrdTotal_s1_r1_e1)
read_df$rcadsanxs_s <- scale(read_df$rcads_scrdAnx_s1_r1_e1)

modelspaip <- lmer(amplitude ~ acc + spaip_s * soc + age_s + (1|sub), data = read_df)
summary(modelspaip)

modelspaic <- lmer(amplitude ~ age_s * spaic_s * acc + soc + (1|sub), data = read_df)
summary(modelspaic)


## Testing new plotting method
# 1. Define the variables for clarity
# Replace "Level_A" and "Level_B" with your actual acc1 factor levels (e.g., 0 and 1)
# If acc1 is numeric, choose specific values like mean +/- SD
my_labels <- c("Low age", "Mean age", "High age")
# --- Plot for Acc Level 1 ---
p1spaic <- interact_plot(
  modelspaic,
  pred = spaic_s,
  modx = soc,
  mod2 = age_s,
  mod2.labels = my_labels,
  at = list(acc = 0),  # <--- ISOLATE ACC LEVEL HERE
  #plot.points = TRUE,            # Optional: show raw data
  #data = read_df %>% filter(acc == 0)
) +
  ggtitle("Error") +
  theme_minimal()

# --- Plot for Acc Level 2 ---
p2spaic <- interact_plot(
  modelspaic,
  pred = spaic_s,
  modx = soc,
  mod2 = age_s,
  mod2.labels = my_labels,
  at = list(acc = 1),  # <--- ISOLATE ACC LEVEL HERE
  #plot.points = TRUE,
  #data = read_df %>% filter(acc == 1)
) +
  ggtitle("Correct") +
  theme_minimal()

# --- Combine them ---
# This creates a side-by-side comparison
# p1 + p2

combined_plot_spaic <- p1spaic + p2spaic + 
  plot_layout(guides = "collect") # Optional: merges legends to save even more space

p1spaip <- interact_plot(
  modelspaip,
  pred = spaip_s,
  modx = soc,
  mod2 = age_s,
  mod2.labels = my_labels,
  at = list(acc = 0),  # <--- ISOLATE ACC LEVEL HERE
  #plot.points = TRUE,            # Optional: show raw data
) +
  ggtitle("Error") +
  theme_minimal() +
  coord_cartesian(ylim = range(read_df$amplitude, na.rm = TRUE))

# --- Plot for Acc Level 2 ---
p2spaip <- interact_plot(
  modelspaip,
  pred = spaip_s,
  modx = soc,
  mod2 = age_s,
  mod2.labels = my_labels,
  at = list(acc = 1),  # <--- ISOLATE ACC LEVEL HERE
  #plot.points = TRUE,
) +
  ggtitle("Correct") +
  theme_minimal() +
  coord_cartesian(ylim = range(read_df$amplitude, na.rm = TRUE))

# --- Combine them ---
# This creates a side-by-side comparison
# p1 + p2

combined_plot_spaip <- p1spaip + p2spaip + 
  plot_layout(guides = "collect") # Optional: merges legends to save even more space

# RCADS as predictor ----

modelrcads <- lmer(amplitude ~ age_s * rcadsanxs_s * acc * soc + (1|sub), data = read_df)
summary(modelrcads)


## Testing new plotting method
# 1. Define the variables for clarity
# Replace "Level_A" and "Level_B" with your actual acc1 factor levels (e.g., 0 and 1)
# If acc1 is numeric, choose specific values like mean +/- SD
my_labels <- c("Low age", "Mean age", "High age")
# --- Plot for Acc Level 1 ---
p1spaic <- interact_plot(
  modelspaic,
  pred = spaic_s,
  modx = soc,
  mod2 = age_s,
  mod2.labels = my_labels,
  at = list(acc = 0),  # <--- ISOLATE ACC LEVEL HERE
  #plot.points = TRUE,            # Optional: show raw data
  #data = read_df %>% filter(acc == 0)
) +
  ggtitle("Error") +
  theme_minimal()

# --- Plot for Acc Level 2 ---
p2spaic <- interact_plot(
  modelspaic,
  pred = spaic_s,
  modx = soc,
  mod2 = age_s,
  mod2.labels = my_labels,
  at = list(acc = 1),  # <--- ISOLATE ACC LEVEL HERE
  #plot.points = TRUE,
  #data = read_df %>% filter(acc == 1)
) +
  ggtitle("Correct") +
  theme_minimal()

# --- Combine them ---
# This creates a side-by-side comparison
# p1 + p2

combined_plot_spaic <- p1spaic + p2spaic + 
  plot_layout(guides = "collect") # Optional: merges legends to save even more space
