# This script was compiled by Charles Knowlton in the NDC Lab for conducting preliminary analyses of READ Study 2
# Load libraries----
library(readr)
library(ggplot2)
library(glmmTMB)
library(lme4)
library(dplyr)
library(tidyr)

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

t.test(read_soc$rt, read_nonsoc$rt, paired = TRUE)

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

## Predict SPAIC from ERN diff and age----
model_lm_c <- lm(spaic ~ ERN_diff + age, data = df_wide)
summary(model_lm_c)

## Predict SPAIP from ERN diff and age----
model_lm_p <- lm(spaip ~ ERN_diff + age, data = df_wide)
summary(model_lm_p)


## Testing difference score----

df_wide <- df_wide %>%
  mutate(
    ERN_diff = ERN_social - ERN_alone
  )
### Predict SPAIC ERN (difference score) and age in each condition----
model <- lm(spaic ~ ERN_diff + age, data = df_wide)
summary(model)
### Predict SPAIP from ERN (difference score) and ERN (difference) and age in each condition----
model <- lm(spaip ~ ERN_diff + age, data = df_wide)
summary(model)

### Predict FNE----

model <- lm(spaic ~ ERN_diff + age, data = df_wide)
summary(model)

model <- lm(spaip ~ ERN_diff + age, data = df_wide)
summary(model)

