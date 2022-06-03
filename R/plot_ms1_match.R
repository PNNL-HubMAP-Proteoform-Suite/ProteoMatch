#' Plots the Proteoform Isotope Profile on top of the experimental sequence
#'
#' @description Returns a static plot with identified calculated proteoform peaks
#'     plotted over the experimental spectrum.
#'
#' @param PeakData A pspecterlib peak_data object or data.table with "M/Z" and "Intensity". Required.
#' @param Ms1Match A ProteoMatch_MatchedPeaksclass object from match_full_seq_ms1. Required.
#' @param ID The ID in the ProteoMatch_MatchedPeaks object to plot. Required.
#'
#' @returns A ggplot object
#'
#' @export
plot_ms1_match <- function(PeakData,
                           Ms1Match,
                           ID) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Check that peak data is of the appropriate class
  if ("peak_data" %in% class(PeakData) == FALSE) {
    stop("PeakData must be a pspecterlib peak_data object.")
  }

  # Ms1Match should be a ProteoMatch_MatchedPeaks object
  if ("ProteoMatch_MatchedPeaks" %in% class(Ms1Match) == FALSE) {
    stop("Ms1Match must be a ProteoMatch_MatchedPeaks object from match_proteoform_to_ms1")
  }

  # ID should be length one
  if (length(ID) != 1) {
    stop("ID should be of length 1.")
  }

  # ID must be an acceptable option
  if (ID %in% Ms1Match$ID == FALSE) {
    stop(paste(ID, "is not a recognized ID."))
  }

  ###################################
  ## MAKE A DATAFRAME FOR PLOTTING ##
  ###################################

  browser()

  # Adjust PeakData to be within range
  class(PeakData) <- c("data.table", "data.frame")
  AdjPeakData <- PeakData %>%
    subset(`M/Z` >= min(ms1_match$`M/Z`) - 1 & `M/Z` <= max(ms1_match$`M/Z`) + 1)

  # Fix up ms1 match
  MS1 <- data.table::data.table(
    `M/Z` = c(ms1_match$`M/Z` - 1e-9, ms1_match$`M/Z`, ms1_match$`M/Z` + 1e-9),
    Abundance = c(rep(0, nrow(ms1_match)), ms1_match$`Intensity`, rep(0, nrow(ms1_match)))
  ) %>%
    dplyr::arrange(`M/Z`)

  ###############
  ## MAKE PLOT ##
  ###############
  ggplot(AdjPeakData, aes(x = `M/Z`, y = Intensity)) + geom_line() +
    geom_line(data = MS1, aes(x = `M/Z`, y = Abundance), color = "red") +
    theme_bw() + ylab("Abundance")
}
