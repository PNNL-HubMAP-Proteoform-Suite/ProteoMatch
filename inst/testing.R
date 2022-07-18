## Pancreas

run_proteomatch(
  ProteoformFile = "~/Desktop/HuBMAP_Data/HubMAP_Pancreas/ProteoMatch_Inputs/Acinar_islet_collapsed.csv",
  MSFile = "~/Desktop/HuBMAP_Data/HubMAP_Pancreas/ProteoMatch_Inputs/20220228_VU_pancreas_40um_IMAGE-16078AVG.mzML",
  SettingsFile = "~/Desktop/HuBMAP_Data/HubMAP_Pancreas/ProteoMatch_Inputs/Pancreas_Defaults.xlsx",
  Path = "~/Desktop/HuBMAP_Data/Results/Pancreas/"
)

## Rat brain

run_proteomatch(
  ProteoformFile = "~/Desktop/HuBMAP_Data/Rat_Brain/Inputs/Rat_Brain_Proteins.csv",
  MSFile = "~/Desktop/HuBMAP_Data/Rat_Brain/Inputs/PeakData.csv",
  SettingsFile = "~/Desktop/HuBMAP_Data/Rat_Brain/Inputs/RatBrain_Defaults.xlsx",
  Path = "~/Desktop/HuBMAP_Data/Results/RatBrain/"
)

## Soybean

# WT
run_proteomatch(
  ProteoformFile = "~/Desktop/HuBMAP_Data/Soybean_Example/ProteoMatch_Inputs/LFQ_soybean_topdown_new.csv",
  MSFile = "~/Desktop/HuBMAP_Data/Soybean_Example/ProteoMatch_Inputs/MALDI_nodule_WT_DHA_recrystalhistoneprep_15um_IMAGE-avg11567_from100to200min.mzML",
  SettingsFile = "~/Desktop/HuBMAP_Data/Soybean_Example/ProteoMatch_Inputs/Soybean_WT_Defaults.xlsx",
  Path = "~/Desktop/HuBMAP_Data/Results/SoyBean/WT/"
)

# MU
run_proteomatch(
  ProteoformFile = "~/Desktop/HuBMAP_Data/Soybean_Example/ProteoMatch_Inputs/LFQ_soybean_topdown_new.csv",
  MSFile = "~/Desktop/HuBMAP_Data/Soybean_Example/ProteoMatch_Inputs/MALDI_nodule_MUTANT_DHA_recrystalhistoneprep_15um_IMAGE-avg11567_from100to200min.mzML",
  SettingsFile = "~/Desktop/HuBMAP_Data/Soybean_Example/ProteoMatch_Inputs/Soybean_MU_Defaults.xlsx",
  Path = "~/Desktop/HuBMAP_Data/Results/SoyBean/MU/"
)
