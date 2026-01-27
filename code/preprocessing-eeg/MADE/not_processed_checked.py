import os
import re
from glob import glob
from datetime import datetime as dt
import sys

session = "s1_r1"
dataset_path = "/home/data/NDClab/datasets/read-study2-dataset/"
data_path = "derivatives/preprocessed/"
start_date = dt.strptime("18/09/34", "%d/%m/%y")

dir_list = sorted(
    [
        p for p in glob(
            f"{dataset_path}{data_path}sub-*/{session}/eeg/MADE_preprocessing_report*.csv"
        )
    ],
    key=os.path.getmtime,
    reverse=True
)


subs_to_process = []
for path in dir_list:
    # print(path)
    sub_id = re.search(r"sub-(\d+)", path)
    # print(sub_id[1])
    timestamp = os.path.getmtime(path)
    if dt.fromtimestamp(timestamp) < start_date:
        subs_to_process.append(sub_id[1])
subs_to_process = sorted(list(set(subs_to_process)))
# print("/".join(subs_to_process))
print("")
print(f"Subjects already preprocessed in {dataset_path+data_path}: {len(subs_to_process)}")
print("")
print("/".join(subs_to_process))

derivatives = subs_to_process.copy()

dataset_path = "/home/data/NDClab/datasets/read-study2-dataset/"
data_path = "/sourcedata/raw/s1_r1/eeg/"

dir_list = sorted(
    [
        p for p in glob(
    f"{dataset_path}{data_path}/sub-*/sub-*all_eeg_s1_r1_e1*.eeg")
    ], # if "ERROR" not in os.path.basename(p)],
                  key=os.path.getmtime, reverse=True
                 )

subs_to_process = []
for path in dir_list:
    # print(path)
    sub_id = re.search(r"sub-(\d+)", path)
    # print(sub_id[1])
    timestamp = os.path.getmtime(path)
    if dt.fromtimestamp(timestamp) < start_date:
        subs_to_process.append(sub_id[1])
subs_to_process = sorted(list(set(subs_to_process)))
# print("/".join(subs_to_process))
print("")
print(f"Subjects total in {dataset_path+data_path}: {len(subs_to_process)}")
sourcedata = subs_to_process.copy()

difference = sorted(list(set(sourcedata) - set(derivatives)))
print("")
print("/".join(difference))
print("")
print(f"Subjects left to preprocess: {len(difference)}")
