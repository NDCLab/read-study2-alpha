
import pandas as pd
import sys

session = sys.argv[1]

deviation_csv = pd.read_csv(f"subjects_with_deviations_{session}.csv")
stim_counts_csv = pd.read_csv(f"logs/qa_log_eeg_deviation_stim_counts_{session}.csv")
stim_counts_csv = stim_counts_csv.rename({"Subject": "sub"}, axis=1)
merged_csv = deviation_csv.merge(stim_counts_csv, on="sub", how="left")
merged_csv.to_csv(f"merged_deviation_{session}.csv")
