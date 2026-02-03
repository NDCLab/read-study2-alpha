# This script was compiled by Charles Knowlton in the NDC Lab for conducting preliminary analyses of READ Study 2
# Load libraries----
library(readr)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(tidyr)
library(reghelper)
library(sjPlot)

# Initialize director folders----
derivatives_dir <- "/home/data/NDClab/analyses/read-study2-alpha/derivatives/"

# Load in data----
read_df <- read_csv(file.path(derivatives_dir, "read_long_s1_r1_30_01_2026_15_00_02.csv"))

# Separate by social and nonsocial
read_soc <- subset(read_df, soc == 1)
read_nonsoc <- subset(read_df, soc == 0)

# Summary statistics of reaction time per condition----
summary(read_soc$rt)
summary(read_nonsoc$rt)

#t.test(read_soc$rt, read_nonsoc$rt, paired = TRUE)

# Aim 1: Mixed-effect model predicting Social Anxiety from Social Observation, Age, and ERN----
## Aggregate variables of interest----
df_summary <- read_df %>%
  group_by(sub, soc) %>%
  summarise(
    ERN = amplitude[acc == 0] - amplitude[acc == 1],
    spaic = first(spaic_scrdTotal_s1_r1_e1),
    age = first(age_m),
    spaip = first(spaip_scrdTotal_s1_r1_e1)
  )

## Remove participants with just one condition----
df_summary <- df_summary %>%
  group_by(sub) %>%
  filter(n() == 2) %>%
  ungroup()

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
df_wide <- df_wide %>% filter(!is.na(ERN_alone))
df_wide <- df_wide %>% filter(!is.na(ERN_social))
df_wide <- df_wide %>% filter(!is.na(spaip))
df_summary <- df_summary %>% filter(!is.na(ERN))
df_summary <- df_summary %>% filter(!is.na(ERN))
df_summary <- df_summary %>% filter(!is.na(spaip))

## Scale Predictors----
df_summary$age <- scale(df_summary$age, center = TRUE, scale = TRUE)
df_summary$spaip <- scale(df_summary$spaip, center = TRUE, scale = TRUE)
df_summary$spaic <- scale(df_summary$spaic, center = TRUE, scale = TRUE)
df_wide$age <- scale(df_wide$age, center = TRUE, scale = TRUE)
df_wide$spaip <- scale(df_wide$spaip, center = TRUE, scale = TRUE)
df_wide$spaic <- scale(df_wide$spaic, center = TRUE, scale = TRUE)

## Mixed-effect modeling----
model_lm_mixed_spaip <- lmer(ERN ~ age + soc + spaip + (1|sub), data = df_summary_age_subset)
summary(model_lm_mixed_spaip)

model_lm_mixed_spaic <- lmer(ERN ~ age + soc + spaic + (1|sub), data = df_summary_age_subset)
summary(model_lm_mixed_spaic)

## Test for interaction----
model_lm_mixed_spaic_int <- lmer(ERN ~ age + spaic + soc + (1|sub), data = df_summary_age_subset)
summary(model_lm_mixed_spaic_int)

model_lm_mixed_spaip_int <- lmer(ERN ~ age + spaip + soc + (1|sub), data = df_summary_age_subset)
summary(model_lm_mixed_spaip_int)

## Simple slopes to interpret interaction better

simple_slopes(model_lm_mixed_spaip_int)
simple_slopes(model_lm_mixed_spaic_int)

## Plotting this interaction
plot_model(model_lm_mixed_spaic_int, type = "pred", terms = c("age", "spaic", "soc"),
           title = "Predicting ERN from Child Reported Social Anxiety, Age, and Condition",
           axis.title = c("Age", "Delta ERN"),
           legend.title = "Social Anxiety")
plot_model(model_lm_mixed_spaip_int, type = "pred", terms = c("age", "spaip", "soc"),
           title = "Predicting ERN from Parent Reported Social Anxiety, Age, and Condition",
           axis.title = c("Age", "Delta ERN"),
           legend.title = "Social Anxiety")
