## Pancreas

# 1. Get Molecular Formula
Acinar <- data.table::fread("~/Desktop/HuBMAP_Data/HubMAP_Pancreas/Acinar_islet_collapsed.csv")
MolForms <- calculate_molform(
  Proteoform = Acinar$Proteoform,
  Protein = Acinar$Protein,
  Charge = c(1:2)
)

# 2. Get Matched Peaks
ScanMetadata <- pspecterlib::get_scan_metadata("~/Desktop/HuBMAP_Data/HubMAP_Pancreas/20220228_VU_pancreas_40um_IMAGE-16078AVG.mzML")
PeakData <- pspecterlib::get_peak_data(ScanMetadata, 1)


match_proteoform_to_ms1(
  PeakData =
)
