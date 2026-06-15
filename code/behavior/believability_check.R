data <- read_csv("my_lab_shortcut/datasets/read-study2-dataset/sourcedata/checked/redcap/Read2bbsRAs1r1_DATA_2026-02-23_1540.csv")

belief_df <- subset(data)
belief <- na.omit(data$bbsratrk_pdbfq2_s1_r1_e1)
true = 0
false = 0
for (i in 1:length(belief)){
  if (belief[i] == 0){
    true = true + 1
  }
  if (belief[i] == 1){
    false = false + 1
  }
}
belief_df <- data |>
  summarise(
    sub = bbsratrk_acid_s1_r1_e1,
    zoom_q = bbsratrk_pdbfq1_s1_r1_e1,
    unusual = bbsratrk_pdbfq2_s1_r1_e1,
    comments = bbsratrk_pdbfq3_s1_r1_e1,
    doubts_age = bbsratrk_dbfq1_s1_r1_e1,
    doubts_age_comm = bbsratrk_dbfq2_s1_r1_e1,
  )

belief_df <- na.omit(belief_df$unusual)

believability_perc = 100 * (true / sum(true, false))