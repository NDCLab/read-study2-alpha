#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
from glob import glob
import numpy as np
from functools import reduce
from datetime import datetime
import re


# In[ ]:


def replace_outliers_with_nan(df, sd_thresh = 3):
    for column in df.columns:
        if pd.api.types.is_numeric_dtype(df[column]):
            mean = df[column].mean()
            std = df[column].std()
            threshold_upper = mean + sd_thresh * std
            threshold_lower = mean - sd_thresh * std

            # Replace outliers with NaN
            df[column] = df[column].apply(lambda x: np.nan if (x > threshold_upper or x < threshold_lower) else x)
    return df

def replace_outliers_with_nan_cols(df, columns_to_check, sd_thresh=3):
    for column in columns_to_check:
        if column in df.columns and pd.api.types.is_numeric_dtype(df[column]):
            mean = df[column].mean()
            std = df[column].std()
            threshold_upper = mean + sd_thresh * std
            threshold_lower = mean - sd_thresh * std

            # Replace outliers with NaN
            df[column] = df[column].apply(lambda x: np.nan if (x > threshold_upper or x < threshold_lower) else x)
    return df


# In[ ]:


analysis_path = "/home/data/NDClab/analyses/read-study2-alpha/"
dataset_path = "/home/data/NDClab/datasets/read-study2-dataset/"

session = "s1_r1"

behavivor_summary_path = f"{analysis_path}/derivatives/behavior/{session}/summary_17_12_2025_14_30_36.csv"


# In[ ]:


behav_data = pd.read_csv(behavivor_summary_path)
behav_id = behav_data["sub"].to_frame()

ern_data = pd.read_csv(f"{analysis_path}/derivatives/csv/{session}/read_erp.csv")
#tf_data = pd.read_csv(f"{analysis_path}/derivatives/csv/{session}/thrive_power_itps.csv")
#icps_data = pd.read_csv(f"{analysis_path}/derivatives/csv/{session}/thrive_icps.csv")

data_frames = [
    behav_id,
    ern_data,
    #tf_data,
    #icps_data,
]

eeg_data = reduce(lambda left, right: pd.merge(left, right, on="sub", how='left'), data_frames)


# In[ ]:


eeg_data_soc = eeg_data[
[i for i in eeg_data.columns if ("_soc" in i or i == "sub")]
]

eeg_data_nonsoc = eeg_data[
[i for i in eeg_data.columns if ("_nonsoc" in i or i == "sub")]
]


# In[ ]:


behavior_df = pd.read_csv(behavivor_summary_path)

behavior_df_nonsoc = behavior_df[[col for col in behavior_df.columns if ("_nonsoc" in col or "sub" in col)]]

behavior_df_nonsoc = behavior_df_nonsoc[behavior_df_nonsoc["acc_nonsoc"] >= 0.6]
behavior_df_nonsoc = replace_outliers_with_nan_cols(behavior_df_nonsoc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="invalid_rt_percent_nonsoc")
behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="skipped_percent_nonsoc")

valid_behavior_nonsoc = behavior_df_nonsoc["sub"].to_frame()

behavior_df_soc = behavior_df[[col for col in behavior_df.columns if ("_soc" in col or "sub" in col)]]

behavior_df_soc = behavior_df_soc[behavior_df_soc["acc_soc"] >= 0.6]
behavior_df_soc = replace_outliers_with_nan_cols(behavior_df_soc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
behavior_df_soc = behavior_df_soc.dropna(subset="invalid_rt_percent_soc")
behavior_df_soc = behavior_df_soc.dropna(subset="skipped_percent_soc")

valid_behavior_soc = behavior_df_soc["sub"].to_frame()


# In[ ]:


valid_eeg_soc = valid_behavior_soc.merge(eeg_data_soc, on="sub", how="left")
valid_eeg_nonsoc = valid_behavior_nonsoc.merge(eeg_data_nonsoc, on="sub", how="left")

merged_valid_eeg_data = valid_eeg_soc.merge(valid_eeg_nonsoc, on="sub", how="outer")


# In[ ]:


behavior_df = pd.read_csv(behavivor_summary_path)

behavior_df_nonsoc = behavior_df[[col for col in behavior_df.columns if ("_nonsoc" in col or "sub" in col)]]

behavior_df_nonsoc = behavior_df_nonsoc[behavior_df_nonsoc["acc_nonsoc"] >= 0.6]
behavior_df_nonsoc = behavior_df_nonsoc[behavior_df_nonsoc["6_or_more_err_nonsoc"] == 1]
behavior_df_nonsoc = replace_outliers_with_nan_cols(behavior_df_nonsoc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="invalid_rt_percent_nonsoc")
behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="skipped_percent_nonsoc")

valid_behavior_nonsoc = behavior_df_nonsoc["sub"].to_frame()

behavior_df_soc = behavior_df[[col for col in behavior_df.columns if ("_soc" in col or "sub" in col)]]

behavior_df_soc = behavior_df_soc[behavior_df_soc["acc_soc"] >= 0.6]
behavior_df_soc = behavior_df_soc[behavior_df_soc["6_or_more_err_soc"] == 1]
behavior_df_soc = replace_outliers_with_nan_cols(behavior_df_soc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
behavior_df_soc = behavior_df_soc.dropna(subset="invalid_rt_percent_soc")
behavior_df_soc = behavior_df_soc.dropna(subset="skipped_percent_soc")

behav_data = replace_outliers_with_nan(
    behavior_df_soc[[i for i in behavior_df_soc if ("sub" in i or\
                                                    ("acc" in i and "con" in i) or\
                                                    ("peri" in i or "pea" in i or "pes" in i))]
    ].merge(
    behavior_df_nonsoc[[i for i in behavior_df_nonsoc if ("sub" in i or\
                                                          ("acc" in i and "con" in i) or\
                                                         ("peri" in i or "pea" in i or "pes" in i))]],
    on="sub", how="outer")
)

# merged_valid_eeg_data = merged_valid_eeg_data.merge(behav_data, on="sub", how="left")


# In[ ]:


age_data = pd.read_csv(iqs_parent_path)
age_data = age_data[[i for i in age_data.columns if "agemos" in i or "record_id" in i]]
age_data = age_data.rename({"record_id" : "sub"}, axis=1)
age_data["sub"] = age_data["sub"] - 80000

age_data = age_data.dropna(how="all", subset = age_data.columns[1:]).reset_index(drop=True)

# to merge Eng and Spanish versions
for i in range(age_data.shape[0]):
    age_data.loc[i, "age_m"] = age_data.iloc[i, 1:].dropna().values[0]

age_data = age_data[["sub", "age_m"]]
if session == "s1_r1":
    age_data["age_m"] = age_data["age_m"] + 9

merged_valid_eeg_data = merged_valid_eeg_data.merge(age_data, on="sub", how="left")


# In[ ]:


full_behavior = pd.read_csv(f"{analysis_path}/derivatives/behavior/{session}/full_df_17_12_2025_14_30_36.csv")
full_behavior.head()
first_soc = pd.DataFrame()

for i, id in enumerate(full_behavior["sub"].unique()):
    id_data = full_behavior[full_behavior["sub"] == id]
    first_soc.loc[i, "sub"] = id
    first_soc.loc[i, "first_soc"] = id_data["condition_soc"].iloc[10]

merged_valid_eeg_data = merged_valid_eeg_data.merge(first_soc, on="sub", how="left")


# In[ ]:


sex_data = pd.read_csv(iqs_parent_path)

sex_data = sex_data.rename({"record_id" : "sub"}, axis=1)
sex_data["sub"] = sex_data["sub"] - 80000

sex_data = sex_data[["sub", f"demo_d_sexbirth_s1_r1_e1", f"demoes_d_sexbirth_s1_r1_e1"]]
sex_data = sex_data.dropna(subset = sex_data.columns[1:], how="all").reset_index(drop=True)

for i in range(sex_data.shape[0]):
    if np.isnan(sex_data.iloc[i, 1]) == 0: # english columns
        sex_data.loc[i, "sex"] = int(sex_data.iloc[i, 1])
    elif np.isnan(sex_data.iloc[i, 2]) == 0: # spanish columns
        sex_data.loc[i, "sex"] = int(sex_data.iloc[i, 2])

sex_data = sex_data[["sub", "sex"]]

merged_valid_eeg_data = merged_valid_eeg_data.merge(sex_data, on="sub", how="left")


# In[ ]:


bbs_p_data = pd.read_csv(f"{dataset_path}/derivatives/preprocessed/redcap/Thrivebbsparents2r1_SCRD_2025-04-09_1555.csv")
bbs_p_data = bbs_p_data.rename({"record_id" : "sub"}, axis=1)
bbs_p_data["sub"] = bbs_p_data["sub"] - 80000
bbs_p_data = bbs_p_data[
    [i for i in bbs_p_data.columns if (i == "sub" or "scrd" in i)]
]

bbs_ch_data = pd.read_csv(f"{dataset_path}/derivatives/preprocessed/redcap/Thrivebbschilds2r1_SCRD_2025-04-09_1555.csv")
bbs_ch_data = bbs_ch_data.rename({"record_id" : "sub"}, axis=1)
bbs_ch_data = bbs_ch_data[
    [i for i in bbs_ch_data.columns if (i == "sub" or "scrd" in i)]
]

merged_valid_eeg_data = merged_valid_eeg_data.merge(bbs_p_data, on="sub", how="left")
merged_valid_eeg_data = merged_valid_eeg_data.merge(bbs_ch_data, on="sub", how="left")


# In[ ]:


spanish_columns = [i for i in merged_valid_eeg_data.columns if ("es_" in i and merged_valid_eeg_data[i].isna().sum() >=100)]

for sp_c in spanish_columns:
    eng_c = "_".join(sp_c.split("es_"))
    assert merged_valid_eeg_data[sp_c].shape[0] == merged_valid_eeg_data[eng_c].shape[0], "ERROR!"
    spanglish_df = merged_valid_eeg_data[["sub", eng_c, sp_c]]
    for row in range(spanglish_df.shape[0]):
        if len(spanglish_df.iloc[row, 1:].dropna().values) > 0:
            merged_valid_eeg_data.loc[row, eng_c] = spanglish_df.iloc[row, 1:].dropna().values[0]
        else:
            merged_valid_eeg_data.loc[row, eng_c] = np.nan

merged_valid_eeg_data = merged_valid_eeg_data.drop(spanish_columns, axis = 1)


# In[ ]:


tf_columns_original = [i for i in merged_valid_eeg_data.columns if not (
    "acc" in i or\
    session in i or\
    "RN_" in i or\
    i=="sub" or\
    i=="first_soc" or\
    i=="sex" or\
    i=="age_m" or\
    "peri" in i or\
    "pes" in i or\
    "pea" in i
)]
tf_columns_renamed = ["_".join(c.split("_err_")) + "_err" if "_err_" in c else "_".join(c.split("_corr_")) + "_corr" if "_corr_" in c else np.nan for c in tf_columns_original]

for i, orig_c in enumerate(tf_columns_original):
    merged_valid_eeg_data.rename({orig_c: tf_columns_renamed[i]}, axis=1, inplace=True)

merged_valid_eeg_data = replace_outliers_with_nan(merged_valid_eeg_data, sd_thresh = 3)

columns = tf_columns_renamed

# Function to isolate pairs
def isolate_pairs(columns):
    pairs = []
    seen = set()

    for col in columns:
        parts = col.split('_')
        measure, condition, window, accuracy = parts[0], parts[1], parts[2], parts[-1]

        # Create a base identifier without the accuracy part
        base_id = '_'.join(parts[:-1])

        if base_id in seen:
            continue

        # Find the corresponding pair
        if accuracy == 'err':
            corr_col = f"{base_id}_corr"
        else:
            corr_col = f"{base_id}_err"

        if corr_col in columns:
            pairs.append((col, corr_col))
            seen.add(base_id)

    return pairs

# Get the pairs
pairs = isolate_pairs(columns)

# Print the pairs
for pair in pairs:
    print(pair)


# In[ ]:


merged_valid_eeg_data = replace_outliers_with_nan(merged_valid_eeg_data, sd_thresh = 3)


# In[ ]:


for s in range(merged_valid_eeg_data.shape[0]):
    for p in pairs:
        if not pd.isnull(merged_valid_eeg_data.loc[s, p[0]]) and not pd.isnull(merged_valid_eeg_data.loc[s, p[1]]):
            merged_valid_eeg_data.loc[s, ("_".join(p[0].split("_")[:-1]) + "_diff")] = merged_valid_eeg_data.loc[s, p[0]] - merged_valid_eeg_data.loc[s, p[1]]
        else:
            merged_valid_eeg_data.loc[s, ("_".join(p[0].split("_")[:-1]) + "_diff")] = np.nan


# In[ ]:


erp_columns = [i for i in merged_valid_eeg_data.columns if ("ERN" in i or "CRN" in i)]
erp_columns


# In[ ]:


for s in range(merged_valid_eeg_data.shape[0]):
    if not pd.isnull(merged_valid_eeg_data.loc[s, "ERN_soc"]) and not pd.isnull(merged_valid_eeg_data.loc[s, "CRN_soc"]):
        merged_valid_eeg_data.loc[s, "ERN_min_CRN_soc"] = merged_valid_eeg_data.loc[s, "ERN_soc"] - merged_valid_eeg_data.loc[s, "CRN_soc"]
    else:
        merged_valid_eeg_data.loc[s, "ERN_min_CRN_soc"] = np.nan
    if not pd.isnull(merged_valid_eeg_data.loc[s, "ERN_nonsoc"]) and not pd.isnull(merged_valid_eeg_data.loc[s, "CRN_nonsoc"]):
        merged_valid_eeg_data.loc[s, "ERN_min_CRN_nonsoc"] = merged_valid_eeg_data.loc[s, "ERN_nonsoc"] - merged_valid_eeg_data.loc[s, "CRN_nonsoc"]
    else:
        merged_valid_eeg_data.loc[s, "ERN_min_CRN_nonsoc"] = np.nan


# In[ ]:


merged_valid_eeg_data = replace_outliers_with_nan_cols(merged_valid_eeg_data, [i for i in merged_valid_eeg_data.columns if ("_diff" in i or "_min_" in i)])


# In[ ]:


state_surveys = pd.read_csv(f"{dataset_path}/derivatives/preprocessed/redcap/Thrivebbschilds2r1_SCRD_2025-04-09_1555.csv")
state_surveys = state_surveys.rename({"record_id": "sub"}, axis=1)

state_surveys = state_surveys[
[i for i in state_surveys.columns if (i == "sub" or "selfnowa" in i or "initstatec" in i or "posttaske" in i or "dyada" in i or "initstated" in i or "posttaskf" in i or "dyadb" in i)\
 and ("timestamp" not in i) and ("_complete" not in i)]
]

state_surveys = state_surveys.merge(first_soc, on="sub", how="left")

new_state_survey_df = pd.DataFrame()
# state_surveys[[i for i in state_surveys.columns if ("initstatec" in i or "first_soc" in i)]]
for c, num_items in zip(["initstatec", "posttaske"], [5, 10]):
    for i in range(state_surveys.shape[0]):
        new_state_survey_df.loc[i, "sub"] = state_surveys.loc[i, "sub"]
        for item in range(1, num_items + 1):
            if state_surveys.loc[i, "first_soc"] == 1:
                new_state_survey_df.loc[i, f"{c}_i{item}_{session}_soc"] = state_surveys.loc[i, f"{c}_i{item}_{session}_e1"]
                new_state_survey_df.loc[i, f"{c}_i{item}_{session}_nonsoc"] = state_surveys.loc[i, f"{c}_i{item}_{session}_e2"]
            elif state_surveys.loc[i, "first_soc"] == 0:
                new_state_survey_df.loc[i, f"{c}_i{item}_{session}_soc"] = state_surveys.loc[i, f"{c}_i{item}_{session}_e2"]
                new_state_survey_df.loc[i, f"{c}_i{item}_{session}_nonsoc"] = state_surveys.loc[i, f"{c}_i{item}_{session}_e1"]
            else:
                new_state_survey_df.loc[i, f"{c}_i{item}_{session}_soc"] = np.nan
                new_state_survey_df.loc[i, f"{c}_i{item}_{session}_nonsoc"] = np.nan

new_state_survey_df = new_state_survey_df.dropna(how="all", subset = new_state_survey_df.columns[1:]).reset_index(drop=True)
state_surveys = new_state_survey_df.merge(state_surveys[[i for i in state_surveys.columns if not ("initstatec" in i or "posttaske" in i or "first_soc" in i)]], on="sub", how="left")


# In[ ]:


#merged_valid_eeg_data = merged_valid_eeg_data.merge(new_state_survey_df, on="sub", how="left")
merged_valid_eeg_data = merged_valid_eeg_data.merge(state_surveys, on="sub", how="left")


# In[ ]:


session = "s1_r1"
merged_valid_eeg_data = pd.read_csv("thrive_mega_df_2025-03-25_21-11-17.csv")
merged_valid_eeg_data.columns = [(i + "_" + session) if (session not in i and i!="sub") else i for i in merged_valid_eeg_data.columns]
merged_valid_eeg_data = merged_valid_eeg_data.merge(pd.read_csv(
    "/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/csv/s2_r1/thrive_mega_df_s2_r1_2025-05-07_19-49-50.csv"), on="sub", how="left")

list(merged_valid_eeg_data.columns)


# In[ ]:


merged_valid_eeg_data.columns = [(i + "_" + session) if (session not in i and i!="sub") else i for i in merged_valid_eeg_data.columns]


# In[ ]:


date_now = datetime.today().strftime('%Y-%m-%d_%H-%M-%S')

# merged_valid_eeg_data.to_csv(f"{analysis_path}/derivatives/csv/{session}/thrive_mega_df_{session}_{date_now}.csv", index=False)
merged_valid_eeg_data.to_csv(f"{analysis_path}/derivatives/csv/s1_r1/thrive_mega_df_s1_s2_{date_now}.csv", index=False)


# In[ ]:


post_error_df = pd.read_csv("/Users/fzaki001/IDENTIFIABLE/ddm_diff_collapsed_both_cond.csv")
post_error_df = post_error_df[[i for i in post_error_df.columns if not\
("ICPS" in i or "power" in i or "ITPS" in i or "dom_hand" in i or "dp_inperson" in i or "task_num" in i or "rd_sda" in i or "sex" in i or "age_m" in i or "first_soc" in i)]]

post_error_df_soc = post_error_df[post_error_df["soc"] == 1]
post_error_df_soc = post_error_df_soc.drop("soc", axis=1).reset_index(drop=True)
post_error_df_soc.columns = [i+"_soc" if not ("sub" in i or "pre_acc" in i) else i for i in post_error_df_soc.columns]
post_error_df_soc = post_error_df_soc[post_error_df_soc["pre_acc"] == 1].drop("pre_acc", axis=1).reset_index(drop=True)
post_error_df_soc = post_error_df_soc[[i for i in post_error_df_soc.columns if ("diff" in i or "sub" in i or "pe" in i or "peri" in i)]]

post_error_df_nonsoc = post_error_df[post_error_df["soc"] == 0]
post_error_df_nonsoc = post_error_df_nonsoc.drop("soc", axis=1).reset_index(drop=True)
post_error_df_nonsoc.columns = [i+"_nonsoc" if not ("sub" in i or "pre_acc" in i) else i for i in post_error_df_nonsoc.columns]
post_error_df_nonsoc = post_error_df_nonsoc[post_error_df_nonsoc["pre_acc"] == 1].drop("pre_acc", axis=1).reset_index(drop=True)
post_error_df_nonsoc = post_error_df_nonsoc[[i for i in post_error_df_nonsoc.columns if ("diff" in i or "sub" in i or "pe" in i or "peri" in i)]]


# In[ ]:


post_error_df = post_error_df_nonsoc.merge(post_error_df_soc, on="sub", how="outer")


# # Merge behav and FNE for LDA

# In[1]:


import os
import pandas as pd
from glob import glob
import numpy as np
from functools import reduce
from datetime import datetime
import re


# In[2]:


def replace_outliers_with_nan(df, sd_thresh = 3):
    for column in df.columns:
        if pd.api.types.is_numeric_dtype(df[column]):
            mean = df[column].mean()
            std = df[column].std()
            threshold_upper = mean + sd_thresh * std
            threshold_lower = mean - sd_thresh * std

            # Replace outliers with NaN
            df[column] = df[column].apply(lambda x: np.nan if (x > threshold_upper or x < threshold_lower) else x)
    return df

def replace_outliers_with_nan_cols(df, columns_to_check, sd_thresh=3):
    for column in columns_to_check:
        if column in df.columns and pd.api.types.is_numeric_dtype(df[column]):
            mean = df[column].mean()
            std = df[column].std()
            threshold_upper = mean + sd_thresh * std
            threshold_lower = mean - sd_thresh * std

            # Replace outliers with NaN
            df[column] = df[column].apply(lambda x: np.nan if (x > threshold_upper or x < threshold_lower) else x)
    return df


# In[ ]:


# merged_valid_eeg_data = merged_valid_eeg_data.merge(behav_data, on="sub", how="left")


# In[48]:


analysis_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"
dataset_path = "/home/data/NDClab/datasets/thrive-dataset/"
dfs = []
for session in ["s1_r1", "s2_r1", "s3_r1"]:
    # list_of_data = []
    print(session)
    pattern = r'(\d{4}-\d{2}-\d{2}_\d{4})'
    most_recent_iqs = max(
        glob(f"{dataset_path}/derivatives/preprocessed/redcap/Thrive*iqschild{"".join(session.split("_"))}*.csv"),
        key=lambda x: re.search(pattern, x).group(1)
    )
    print(most_recent_iqs)
    iqs_ch = pd.read_csv(most_recent_iqs)
    iqs_ch = iqs_ch.rename({"record_id": "sub"}, axis=1)
    iqs_ch = iqs_ch[["sub", f"bfne_b_scrdTotal_{session}_e1"]]
    # print(iqs_ch[f"bfne_b_scrdTotal_{session}_e1"].dropna().shape[0])
    # list_of_data.append(iqs_ch)

    if session == "s1_r1":
        iqs_parent_path = max(
            glob(f"{dataset_path}/sourcedata/checked/redcap/Thrive*iqsparent{"".join(session.split("_"))}*.csv"),
            key=lambda x: re.search(pattern, x).group(1)
        )
        print(iqs_parent_path)
        age_data = pd.read_csv(iqs_parent_path)
        age_data = age_data[[i for i in age_data.columns if "agemos" in i or "record_id" in i]]
        age_data = age_data.rename({"record_id" : "sub"}, axis=1)
        age_data["sub"] = age_data["sub"] - 80000

        age_data = age_data.dropna(how="all", subset = age_data.columns[1:]).reset_index(drop=True)

        # to merge Eng and Spanish versions
        for i in range(age_data.shape[0]):
            age_data.loc[i, "age_m"] = age_data.iloc[i, 1:].dropna().values[0]

        age_data = age_data[["sub", "age_m"]]
        age_data = age_data.rename({"age_m": "age"}, axis=1)

    # if session == "s2_r1":
    #     age_data["age"] = age_data["age"] + 9
    # elif session == "s3_r1":
    #     age_data["age"] = age_data["age"] + 18

        age_data.columns = [i+f"_{session}_e1" if i != "sub" else i for i in age_data.columns]
        # list_of_data.append(age_data)

    # redcap_data = iqs_ch.merge(age_data, on="sub", how="left")

    pattern = r'(\d{2}_\d{2}_\d{4}_\d{2}_\d{2}_\d{2})'
    behavivor_summary_path = max(
        glob(f"{analysis_path}/derivatives/behavior/{session}/summary*{session}*.csv"),
        key=lambda x: re.search(pattern, x).group(1)
    )

    behavior_df = pd.read_csv(behavivor_summary_path)

    behavior_df_soc = behavior_df[[col for col in behavior_df.columns if ("_soc" in col or "sub" in col)]]

    behavior_df_soc = behavior_df_soc[behavior_df_soc["acc_soc"] >= 0.6]
    behavior_df_soc = behavior_df_soc[behavior_df_soc["6_or_more_err_soc"] == 1]
    behavior_df_soc = replace_outliers_with_nan_cols(behavior_df_soc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
    behavior_df_soc = behavior_df_soc.dropna(subset="invalid_rt_percent_soc")
    behavior_df_soc = behavior_df_soc.dropna(subset="skipped_percent_soc")

    behav_data = behavior_df_soc[["sub", "peri_rt_soc"]]
    behav_data.columns = [i+f"_{session}_e1" if i != "sub" else i for i in behav_data.columns]

    full_data = behav_data.merge(iqs_ch, on="sub", how="left")
    if session == "s1_r1":
        full_data = full_data.merge(age_data, on="sub", how="left")
    # list_of_data.append(behav_data)
    # full_data = reduce(lambda left, right: pd.merge(left, right, on='sub', how='left'), list_of_data)
    # full_data = replace_outliers_with_nan(full_data)

    dfs.append(full_data)

merged_df = reduce(lambda left, right: pd.merge(left, right, on='sub', how='left'), dfs)
for session in ["s2_r1", "s3_r1"]:
    if session == "s2_r1":
        merged_df[f"age_{session}_e1"] = merged_df["age_s1_r1_e1"] + 9
    elif session == "s3_r1":
        merged_df[f"age_{session}_e1"] = merged_df["age_s1_r1_e1"] + 18
for age_column in [i for i in merged_df.columns if "age" in i]:
    merged_df[age_column] = merged_df[age_column]/12

merged_df = merged_df.dropna().reset_index(drop=True)
merged_df.columns = [i.replace("_r1_e1", "") if i != "sub" else i for i in merged_df.columns]
merged_df.columns = [i.replace("_b_scrdTotal", "") if "bfne" in i else i for i in merged_df.columns]
merged_df.columns = [i.replace("_rt_soc", "") if "peri" in i else i for i in merged_df.columns]

merged_df


# In[49]:


merged_df.to_csv("LDA_data.csv", index=False)


# In[38]:


dfs[0][dfs[0]["sub"] == 3000123]


# In[24]:


behav_data.merge(iqs_ch[["sub", f"bfne_b_scrdTotal_{session}_e1"]].dropna(), on="sub", how="left").dropna()


# In[ ]:




