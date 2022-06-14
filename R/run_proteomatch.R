#' Run the ProteoMatch pipeline with a Proteoform, mzML, and Settings file
#'
#' @description A wrapper for quickly running all the steps of the ProteoMatch pipeline
#'
#' @param ProteoformFile Path to a csv with a column of "Proteoforms" and their "Protein"
#'    identifier. Each Proteoform-Protein pairing should be unique. Proteoforms
#'    should be written with proper annotations (i.e. "M.AA`[`Acetyl`]`AA`[`3.2`]`.V"). Required.
#' @param mzMLFile Path to an mzML file containing one spectra. The file must be centroided.
#'    Use the [msconverteR](https://github.com/wilsontom/msconverteR) package and
#'    set msconvert_args = 'peakPicking true 1-'. Required.
#' @param SettingsFile Path to a xlsx file with all user-set parameters. Required.
#' @param Path The base directory of the output trelliscope application. Default is Downloads/Ms1Match.
#'
#' @details The algorithm is run in 4 main steps:
#' 1. Molecular Formulas are calculated with calculate_molform
#' 2. Input spectra is filtered with filter_peaks
#' 3. Isotopic distributions are calculated with the molecular formulas, and matched
#'     to the filtered spectra
#' 4. All high scoring matches are visualized with a trelliscope display
#'
#' @export
run_proteomatch <- function(ProteoformFile,
                            mzMLFile,
                            SettingsFile,
                            Path = file.path(.getDownloadsFolder(), "Ms1Match"),
                            Messages = TRUE) {

  #################
  ## CHECK FILES ##
  #################

  # The proteoform file should be a csv
  if (!grepl(".csv", ProteoformFile)) {
    stop("ProteoformFile should be a csv.")
  }
  Proteoforms <- data.table::fread(ProteoformFile)

  # The proteoform file should contain a proteoform and protein
  if (!all(c("Protein", "Proteoform") %in% colnames(Proteoforms))) {
    stop("ProteoformFile should have a column labeled 'Protein' and another labeled 'Proteoform'.")
  }

  # The MS data should be an mzML file
  if (!grepl(".mzML", mzMLFile)) {
    stop("mzMLFile should be in mzML format.")
  }
  MSData <- pspecterlib::get_scan_metadata(mzMLFile)

  # The mzML file should have one scan
  if (nrow(MSData) > 1) {
    stop("mzMLFile has more than one scan. This is unexpected.")
  }
  PeakData <- pspecterlib::get_peak_data(MSData, 1)

  # The settings file should have be an xlsx
  if (!grepl(".xlsx", SettingsFile)) {
    stop("SettingsFile needs to be an xlsx file.")
  }
  Settings <- suppressWarnings({xlsx::read.xlsx(SettingsFile, 1)})

  # The settings file should have all the required parameters
  RequiredRow <- c("MZRange", "NoiseFilter", "Charges", "MinAbsoluteChange", "CorrelationMinimum",
    "IsotopicPercentage","PPMThreshold", "IsotopeRange", "PlottingWindow", "ProtonMass")
  if (!all(Settings$Parameter %in% RequiredRow)) {
    stop("Settings file is missing: ",
         paste0(RequiredRow[!RequiredRow  %in% Settings$Parameter], ", ", collapse = ""))
  }

  ##################
  ## RUN PIPELINE ##
  ##################

  # 0. Create output directory
  if (!dir.exists(Path)) {dir.create(Path)}

  # 1. Calculate Molecular Formula
  if (Messages) {message("Calculating molecular formulas...")}
  MolForm <- calculate_molform(
    Proteoform = Proteoforms$Proteoform,
    Protein = Proteoforms$Protein,
    Charge = Settings[Settings$Parameter == "Charges", "Default"] %>% strsplit(",") %>% unlist() %>% as.numeric(),
    ProtonMass = Settings[Settings$Parameter == "ProtonMass", "Default"] %>% as.numeric()
  )
  write.csv(MolForm, file.path(Path, "Molecular_Formulas.csv"), row.names = F, quote = F)

  # 2. Filter the data
  if (Messages) {message("Filtering peaks...")}
  FilteredData <- filter_peaks(
    PeakData = PeakData,
    MZRange = Settings[Settings$Parameter == "MZRange", "Default"] %>% strsplit("-") %>% unlist() %>% as.numeric(),
    NoiseFilter = Settings[Settings$Parameter == "NoiseFilter", "Default"] %>% as.numeric()
  )
  write.csv(FilteredData, file.path(Path, "Filtered_Peaks.csv"), row.names = F, quote = F)

  # 3. Match Peaks
  if (Messages) {message("Matching spectra...")}
  MatchedPeaks <- match_proteoform_to_ms1(
    PeakData = FilteredData,
    MolecularFormulas = MolForm,
    IsotopicPercentage = Settings[Settings$Parameter == "IsotopicPercentage", "Default"] %>% as.numeric(),
    PPMThreshold = Settings[Settings$Parameter == "PPMThreshold", "Default"] %>% as.numeric(),
    MinAbsoluteChange = Settings[Settings$Parameter == "MinAbsoluteChange", "Default"] %>% as.numeric(),
    IsotopeRange = Settings[Settings$Parameter == "IsotopeRange", "Default"] %>% strsplit(",") %>% unlist() %>% as.numeric(),
    ProtonMass = Settings[Settings$Parameter == "ProtonMass", "Default"] %>% as.numeric()
  )
  write.csv(MatchedPeaks, file.path(Path, "Matched_Isotope_Distributions.csv"), row.names = F, quote = F)

  # 4. Make the trelliscope display
  if (Messages) {message("Generating trelliscope...")}
  proteomatch_trelliscope(
    PeakData = FilteredData,
    Ms1Match = MatchedPeaks,
    MinCorrelationScore = Settings[Settings$Parameter == "CorrelationMinimum", "Default"] %>% as.numeric(),
    Window = Settings[Settings$Parameter == "PlottingWindow", "Default"] %>% as.numeric()
  )

}
