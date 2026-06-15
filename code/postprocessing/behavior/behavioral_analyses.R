#### Behavioral analysis for READ ####
# Created by Charles Knowlton at the NDC Lab at FIU. Based on mfe_d_face_data_organizer_v1.R by Kianoosh Hosseini                                                                        

# Purpose is to include/exclude participants based on their summary behavioral data.
# New table with the following for each condition:
# 1. accuracy
# 2. RT avg
# 3. Did they have lower than 60% accuracy (excluding multiResp)?
# 4. Did they have an outlier number of fastRTs?
# 5. Is the condition valid (3, 4, 5 are not met)?


######## Start ##########
library(dplyr)
library(stringr)
library(readr)

#Set directory
setwd("/home/data/NDClab/datasets/read-study2-dataset/code/behavior")

#Import data
flanker_summary <- read_csv("/home/data/NDClab/analyses/read-study2-alpha/derivatives/behavior/s1_r1/summary_17_12_2025_14_30_36.csv", na = 'NA')

#Removing participant with partial data that didn't complete the task for preliminary analyses
flanker_summary <- flanker_summary[-16,]

#Initializing variables:
subList <- unique(flanker_summary$sub)

mean_fastRTs_overall <- mean(flanker_summary$fastRTs_overall)
sd_fastRTs_overall <- sd(flanker_summary$fastRTs_overall)

mean_fastRTs_nonsoc <- mean(flanker_summary$invalid_rt_percent_nonsoc)
sd_fastRTs_nonsoc <- sd(flanker_summary$invalid_rt_percent_nonsoc)

mean_fastRTs_soc <- mean(flanker_summary$invalid_rt_percent_soc)
sd_fastRTs_soc <- sd(flanker_summary$invalid_rt_percent_soc)

#Add columns to flanker summary
flanker_summary[ , c("exclude_overall", 
              "exclude_social", 
              "exclude_nonsocial"
              )] <- NA

#Make sure any NaNs are NA
flanker_summary[flanker_summary == "NaN"] <- NA

#Remove participants with under 60% accuracy
#Starting all participant loop
for(i in 1:length(subList)){
  
  # Alone rejections
  
  if (flanker_summary$n_trials_nonsoc[flanker_summary$sub == subList[i]] != 0){
    # if numTrials_condition were greater than zero, then do the below exclusions:
    if (flanker_summary$acc_nonsoc[flanker_summary$sub == subList[i]] <.60 |  # if the mean is less than 60%
        flanker_summary$invalid_rt_percent_nonsoc[flanker_summary$sub == subList[i]]> (3*sd_fastRTs_nonsoc) | # OR an outlier num of fastRTs
        flanker_summary$`6_or_more_err_nonsoc`[flanker_summary$sub == subList[i]] == 0){ # OR less than 6 error trials
      # then the exclude variable for this condition is labeled as 1
      flanker_summary$exclude_nonsocial[flanker_summary$sub == subList[i]] <- 1 
    }else{ # if it doesn't meet exclusions, then exclude = 0
      flanker_summary$exclude_nonsocial[flanker_summary$sub == subList[i]] <- 0
    }
  }else{ # if numTrials_conditon is 0, then exclude it
    flanker_summary$exclude_nonsocial[flanker_summary$sub == subList[i]] <- 1
  }

  # Social rejections
  if (flanker_summary$n_trials_soc[flanker_summary$sub == subList[i]]!=0){
    if (flanker_summary$acc_soc[flanker_summary$sub == subList[i]] < .60 |  # if the mean is less than 60%
        flanker_summary$invalid_rt_percent_soc[flanker_summary$sub == subList[i]] > (3*sd_fastRTs_soc) | # OR an outlier num of fastRTs
        flanker_summary$`6_or_more_err_soc`[flanker_summary$sub == subList[i]] == 0){ #OR less than 6 error trials
      # then the exclude variable for this condition is labeled as 1
      flanker_summary$exclude_social[flanker_summary$sub == subList[i]] <- 1 
    }else{
      flanker_summary$exclude_social[flanker_summary$sub == subList[i]] <- 0
    }
  }else{ # if numTrials_conditon is 0, then exclude it
    flanker_summary$exclude_social[flanker_summary$sub == subList[i]] <- 1
  }
}

write.csv(flanker_summary, "/home/data/NDClab/analyses/read-study2-alpha/derivatives/behavior/s1_r1/prelimBehavioralRejection.csv", row.names = FALSE)