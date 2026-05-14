# This script was created by Charles Knowlton in the NDC Lab for conducting preliminary analyses of READ Study 2 flanker EEG + behavioral data
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

# Initialize directory folders and load data----
derivatives_dir <- "/home/data/NDClab/analyses/read-study2-alpha/derivatives/csv/s1_r1"
read_df <- read_csv(file.path(derivatives_dir, "read_long_s1_r1_14_04_2026_13_45_42.csv"))

#enter participants that failed deception check (i.e. they answered "very much" to did you believe that no one was watching you at all?)
#failed_deception <- c("3300118", "3300026", "3300071", "3300131", "3300123", "3300067", "3300052", "3300013")

#read_df <- read_df |> filter(!(sub %in% failed_deception))

# Prepare data----
## Aggregate variables of interest----
df_summary_long <- read_df |>
  group_by(sub, soc) |>
  summarise(
    deltaERN = first(ERN_min_CRN),
    thetapower = first(power_early_diff),
    spaic = first(spaic_scrdTotal_s1_r1_e1),
    age = first(age_m),
    spaip = first(spaip_scrdTotal_s1_r1_e1),
    sex = first(sex),
    #acc_score = first(accuracy_score),
    #invalid_rt_percent = first(invalid_rt_percent),
    #rt = first(rt),
    #rt_con = first(rt_con),
    #skipped_percent = first(skipped_percent),
    #acc_con = first(acc_con),
    #rt_incon = first(rt_incon),
    #acc_incon = first(acc_incon)
  )

#df_summary_long <- df_summary_long |> subset(
#  age > mean(age) - (3 * sd(age)) & age < mean(age) + (3 * sd(age))
#)

## Omit NaNs----
# NaN values is the data set were assigned during data aggregation, therefore they reflect either:
# 1.) Missing data
# 2.) Outlier based NaN assignment (+- 3 SD )
# participants_remove <- c() # Stored vector of participants with Nan for both conditions 
# p_i <- 1 # Looping variable for index
# for (i in 1:nrow(df_summary_wide)){
#   p_id = as.character(df_summary_wide$sub[i])
#   if ((is.na(df_summary_wide$ERN_alone[i]) == TRUE) & (is.na(df_summary_wide$ERN_social[i]) == FALSE)){
#     print(paste("Participant's ERN alone (social is valid) is NAN, ID:", p_id))
#   }
#   if ((is.na(df_summary_wide$ERN_social[i]) == TRUE) & (is.na(df_summary_wide$ERN_alone[i]) == FALSE)){
#     print(paste("Participant's ERN social (alone is valid) is NAN, ID:", p_id))
#   }
#   if ((is.na(df_summary_wide$ERN_social[i]) == TRUE) & (is.na(df_summary_wide$ERN_alone[i]) == TRUE)){
#     print(paste("Participant has ERN NaN for both conditions, ID:", p_id))
#     participants_remove[p_i] <- p_id
#     p_i <- p_i + 1
#   }
# }

# Remove participants from long and wide data frames that contained no ERN data
#df_wide <- df_wide %>% filter(!(sub %in% participants_remove))
#df_summary <- df_summary %>% filter(!(sub %in% participants_remove))
# Remove participants with no SPAIP data
#df_wide <- df_wide %>% filter(!is.na(spaip))
#df_summary <- df_summary %>% filter(!is.na(spaip))
# Remove participants with no SPAIC data
#df_wide <- df_wide %>% filter(!is.na(spaic))
#df_summary <- df_summary %>% filter(!is.na(spaic))

## Remove any data that is not recorded in both conditions
df_summary_long <- na.omit(df_summary_long) #note: data was compiled to assign NA to any variable that is +3 or -3 sd from mean as identifying outliers early
df_summary_long <- df_summary_long |>
  group_by(sub) |>
  filter(n() == 2) |>
  ungroup()

## Scaling variables----
df_summary_long$age_s <- scale(df_summary_long$age)
df_summary_long$spaic_s <- scale(df_summary_long$spaic)
df_summary_long$spaip_s <- scale(df_summary_long$spaip)
df_summary_long$deltaERN_s <- scale(df_summary_long$deltaERN)
df_summary_long$thetapower_s <- scale(df_summary_long$thetapower)

## Apply sum contrasts ----
df_summary_long$soc <- factor(df_summary_long$soc, 
                              levels = c(0, 1), 
                              labels = c("alone", "social"))

df_summary_long$sex <- factor(df_summary_long$sex, 
                              levels = c(1, 2), 
                              labels = c("male", "female"))

contrasts(df_summary_long$soc) <- rev(contr.sum(2))
contrasts(df_summary_long$sex) <- rev(contr.sum(2))

# Aim #1: Social Anxiety x Social Observation x Age predicting ERN model----
### Mixed effect modeling---- 
model_lm_spaic_int_lmer <- lmer(deltaERN_s ~ spaic_s * age_s * soc * sex + (1|sub), data = df_summary_long)
summary(model_lm_spaic_int_lmer)

model_lm_spaip_int_lmer <- lmer(deltaERN_s ~ age_s * spaip_s * soc * sex + (1|sub), data = df_summary_long)
summary(model_lm_spaip_int_lmer)

#### Interaction plots----

# First, check levels and contrasts to verify they are being generated correctly
print("Factor levels for soc:")
print(levels(df_summary_long$soc))
print("Contrasts:")
print(contrasts(df_summary_long$soc))

#Plot only significant relationships
plot_lm_spaic_int_lmer <- interact_plot(
  model_lm_spaic_int_lmer, 
  pred = spaic_s, # X - axis
  modx = age_s, # Interaction term
  modx.values = "mean-plus-minus", # Method of tercile subset
  plot.points = TRUE, 
  interval = TRUE, # Confidence bands
  main.title = "SPAIC",
  x.label = "SPAIC",
  y.label = "Delta ERN",
  legend.main = "Age"
)

p_large <- plot_lm_spaic_int_lmer + 
  theme(
    text = element_text(size = 16, family = "sans"),           
    axis.title = element_text(size = 18, family = "sans"),     
    axis.text = element_text(size = 14, family = "sans"),      
    plot.title = element_text(size = 20, family = "sans"),     
    legend.text = element_text(size = 14, family = "sans"),    
    legend.title = element_text(size = 16, family = "sans")    
  )

ggsave("model_lm_spaip_age_int_lmer.png", 
       plot = p_large, 
       width = 10, 
       height = 7, 
       units = "in",
       dpi = 600)

plot_lm_spaic_int_lmer <- interact_plot(
  model_lm_spaic_int_lmer, 
  pred = spaic_s, # X - axis
  modx = sex, # Interaction term
  plot.points = TRUE, 
  interval = TRUE, # Confidence bands
  main.title = "SPAIC",
  x.label = "SPAIC",
  y.label = "Delta ERN",
  mod2.labels = c("Male", "Female"),
  legend.main = "Sex"
)

p_large <- plot_lm_spaic_int_lmer + 
  theme(
    text = element_text(size = 16, family = "sans"),           
    axis.title = element_text(size = 18, family = "sans"),     
    axis.text = element_text(size = 14, family = "sans"),      
    plot.title = element_text(size = 20, family = "sans"),     
    legend.text = element_text(size = 14, family = "sans"),    
    legend.title = element_text(size = 16, family = "sans")    
  )

ggsave("model_lm_spaic_int_sex_lmer.png", 
       plot = p_large, 
       width = 10, 
       height = 7, 
       units = "in",
       dpi = 600)

plot_lm_spaip_int_lmer <- interact_plot(
  model_lm_spaip_int_lmer, 
  pred = spaip_s, # X - axis
  modx = age_s, # Interaction term
  modx.values = "mean-plus-minus", # Method of tercile subset
  mod2 = soc,
  plot.points = TRUE, 
  interval = TRUE, # Confidence bands
  main.title = "SPAICP",
  x.label = "SPAICP",
  y.label = "Delta ERN",
  mod2.labels = c("Alone", "Social"),
  legend.main = "Age"
)

p_large <- plot_lm_spaip_int_lmer + 
  theme(
    text = element_text(size = 16, family = "sans"),           
    axis.title = element_text(size = 18, family = "sans"),     
    axis.text = element_text(size = 14, family = "sans"),      
    plot.title = element_text(size = 20, family = "sans"),     
    legend.text = element_text(size = 14, family = "sans"),    
    legend.title = element_text(size = 16, family = "sans")    
  )

ggsave("model_lm_spaip_int_lmer.png", 
       plot = p_large, 
       width = 10, 
       height = 7, 
       units = "in",
       dpi = 600)

# Aim #2: Social Anxiety x Social Observation x Age predicting Midfrontal Theta power model----
model_lm_spaictheta_int_lmer <- lmer(thetapower_s ~ spaic_s * soc * age_s * sex + (1|sub), data = df_summary_long)
summary(model_lm_spaictheta_int_lmer)

model_lm_spaiptheta_int_lmer <- lmer(thetapower_s ~ sex + (1|sub), data = df_summary_long)
summary(model_lm_spaiptheta_int_lmer)

interact_plot(
  model_lm_spaictheta_int_lmer,
  pred = spaic_s,
  modx = sex,
  mod2 = soc,
  plot.points = TRUE, 
  interval = TRUE,
  main.title = "SPAIC",
  x.label = "SPAIC (scaled)",
  y.label = "Power (scaled)",
  legend.main = "Sex",
  mod2.labels = c("Alone", "Social")
)

interact_plot(
  model_lm_spaictheta_int_lmer,
  pred = spaic_s,
  modx = age_s,
  mod2 = soc,
  modx.values = "mean-plus-minus",
  plot.points = TRUE, 
  interval = TRUE,
  main.title = "SPAIC",
  x.label = "SPAIC (scaled)",
  y.label = "Power (scaled)",
  legend.main = "Age",
  mod2.labels = c("Alone", "Social")
)

## final plots

# ERN
model_lm_spaic_int_lmer <- lmer(deltaERN_s ~ spaic_s * age_s * soc * sex + (1|sub), data = df_summary_long)
summary(model_lm_spaic_int_lmer)

model_lm_spaip_int_lmer <- lmer(deltaERN_s ~ age_s * spaip_s * soc + (1|sub), data = df_summary_long)
summary(model_lm_spaip_int_lmer)

p_age <- interact_plot(
  model_lm_spaic_int_lmer, 
  pred = spaic_s,
  modx = age_s,
  modx.values = "mean-plus-minus",
  plot.points = TRUE, 
  interval = TRUE,
  main.title = "SPAIC",
  x.label = "SPAIC",
  y.label = "Delta ERN"
)

p_large <- p_age + 
  theme(
    text = element_text(size = 16, family = "sans"),           
    axis.title = element_text(size = 18, family = "sans"),     
    axis.text = element_text(size = 14, family = "sans"),      
    plot.title = element_text(size = 20, family = "sans"),     
    legend.text = element_text(size = 14, family = "sans"),    
    legend.title = element_text(size = 16, family = "sans")    
  )

ggsave("model_lm_deltaern_spaic_int_lmer_age.png", 
       plot = p_large, 
       width = 7, 
       height = 7, 
       units = "in",
       dpi = 600)

p_sex <- interact_plot(
  model_lm_spaic_int_lmer, 
  pred = spaic_s,
  modx = sex,
  plot.points = TRUE, 
  interval = TRUE,
  vary.lty = FALSE,
  main.title = "SPAIC",
  x.label = "SPAIC",
  y.label = "Delta ERN"
)

p_large <- p_sex + coord_cartesian(ylim = yr) + 
  theme(
    text = element_text(size = 16, family = "sans"),           
    axis.title = element_text(size = 18, family = "sans"),     
    axis.text = element_text(size = 14, family = "sans"),      
    plot.title = element_text(size = 20, family = "sans"),     
    legend.text = element_text(size = 14, family = "sans"),    
    legend.title = element_text(size = 16, family = "sans")    
  )

ggsave("model_lm_deltaern_spaic_int_lmer_sex.png", 
       plot = p_large, 
       width = 7, 
       height = 7, 
       units = "in",
       dpi = 600)


yr <- ggplot_build(p_age)$layout$panel_scales_y[[1]]$range$range

# theta
model_lm_spaictheta_int_lmer <- lmer(thetapower_s ~ spaic_s * soc + (1|sub), data = df_summary_long)
summary(model_lm_spaictheta_int_lmer)

model_lm_spaiptheta_int_lmer <- lmer(thetapower_s ~ spaip_s * soc + (1|sub), data = df_summary_long)
summary(model_lm_spaiptheta_int_lmer)

plot_lm_spaic_int_lmer <- interact_plot(
  model_lm_spaictheta_int_lmer,
  pred = spaic_s,
  modx = sex,
  mod2 = soc,
  plot.points = TRUE, 
  interval = TRUE,
  vary.lty = FALSE,
  main.title = "SPAIC",
  x.label = "SPAIC (scaled)",
  y.label = "Power (scaled)",
  legend.main = "Sex",
  mod2.labels = c("Alone", "Social")
)

p_large <- plot_lm_spaic_int_lmer + 
  coord_cartesian(ylim = age_ylim) +
  theme(
    text = element_text(size = 16, family = "sans"),           
    axis.title = element_text(size = 18, family = "sans"),     
    axis.text = element_text(size = 14, family = "sans"),      
    plot.title = element_text(size = 20, family = "sans"),     
    legend.text = element_text(size = 14, family = "sans"),    
    legend.title = element_text(size = 16, family = "sans")    
  )

ggsave("model_lm_spaic_int_lmer_sex.png", 
       plot = p_large, 
       width = 10, 
       height = 7, 
       units = "in",
       dpi = 600)

plot_lm_spaic_int_lmer2 <- interact_plot(
  model_lm_spaictheta_int_lmer,
  pred = spaic_s,
  modx = age_s,
  mod2 = soc,
  modx.values = "mean-plus-minus",
  plot.points = TRUE, 
  interval = TRUE,
  main.title = "SPAIC",
  x.label = "SPAIC (scaled)",
  y.label = "Power (scaled)",
  legend.main = "Age",
  mod2.labels = c("Alone", "Social")
)

p_large <- plot_lm_spaic_int_lmer2 + 
  theme(
    text = element_text(size = 16, family = "sans"),           
    axis.title = element_text(size = 18, family = "sans"),     
    axis.text = element_text(size = 14, family = "sans"),      
    plot.title = element_text(size = 20, family = "sans"),     
    legend.text = element_text(size = 14, family = "sans"),    
    legend.title = element_text(size = 16, family = "sans")    
  )

ggsave("model_lm_spaic_int_lmer_age.png", 
       plot = p_large, 
       width = 10, 
       height = 7, 
       units = "in",
       dpi = 600)


subject_ids <- unique(df_summary_long$sub)
write.table(subject_ids, 
            file = "/home/data/NDClab/analyses/read-study2-alpha/derivatives/final_subjects_list.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# EM trends ----
## ERN model ----
### Plot 1: First interaction - SPAIC × age----
trends_ern_age <- emtrends(
  model_lm_spaic_int_lmer,
  ~ age_s,
  var = "spaic_s",
  at = list(age_s = c(-1, 0, 1))
)
summary(trends_ern_age, infer = c(TRUE, TRUE))
pairs(trends_ern_age)

### Plot 2: Second interaction - SPAIC × sex----
trends_ern_sex <- emtrends(
  model_lm_spaic_int_lmer,
  ~ sex,
  var = "spaic_s"
)
summary(trends_ern_sex, infer = c(TRUE, TRUE))
pairs(trends_ern_sex)


## Midfrontal theta model----
### First interaction: SPAIC × sex × soc----
trends_theta_sex <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ sex | soc,
  var = "spaic_s"
)
summary(trends_theta_sex, infer = c(TRUE, TRUE))
pairs(trends_theta_sex)

#### Follow up tests at specific contrasts----
trends_theta_by_sex <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ soc | sex,
  var = "spaic_s"
)
summary(trends_theta_by_sex, infer = c(TRUE, TRUE))
pairs(trends_theta_by_sex)

trends_theta_by_soc <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ sex | soc,
  var = "spaic_s"
)
summary(trends_theta_by_soc, infer = c(TRUE, TRUE))
pairs(trends_theta_by_soc)

### Second interaction: SPAIC × age × soc----
trends_theta_age <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ age_s | soc,
  var = "spaic_s",
  at = list(age_s = c(-1, 0, 1))
)
summary(trends_theta_age, infer = c(TRUE, TRUE))
pairs(trends_theta_age)

#### Follow up contrasts for specific contrasts----
trends_theta_older <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ soc,
  var = "spaic_s",
  at = list(age_s = 1)
)
summary(trends_theta_older, infer = c(TRUE, TRUE))
pairs(trends_theta_older)

trends_theta_avg <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ soc,
  var = "spaic_s",
  at = list(age_s = 0)
)
summary(trends_theta_avg, infer = c(TRUE, TRUE))
pairs(trends_theta_avg)

trends_theta_young <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ soc,
  var = "spaic_s",
  at = list(age_s = -1)
)
summary(trends_theta_young, infer = c(TRUE, TRUE))
pairs(trends_theta_young)

trends_theta_social <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ age_s,
  var = "spaic_s",
  at = list(age_s = c(-1, 1), soc = "social")
)
summary(trends_theta_social, infer = c(TRUE, TRUE))
pairs(trends_theta_social)

trends_theta_males <- emtrends(
  model_lm_spaictheta_int_lmer,
  ~ age_s,
  var = "spaic_s",
  at = list(age_s = c(-1, 1), soc = "social")
)
summary(trends_theta_social, infer = c(TRUE, TRUE))
pairs(trends_theta_social)