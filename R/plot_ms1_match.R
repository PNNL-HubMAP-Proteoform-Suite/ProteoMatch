#' Matches peaks to isotoping profile, given the molecular formula and any mass changes
#'
#' @param PeakData A pspecterlib peak_data object or data.table with "M/Z" and "Intensity". Required.
#' @param ms1_match A ms1_match class object from match_full_seq_ms1. Required.
plot_ms1_match <- function(PeakData, ms1_match) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Check that peak data is of the appropriate class
  if ("peak_data" %in% class(PeakData) == FALSE) {

    # And if not, is a data.table with M/Z and Intensity
    if (colnames(PeakData)[1] != "M/Z" | colnames(PeakData)[2] != "Intensity") {
      stop("PeakData must be of the peak_data class or include two columns: M/Z and Intensity.")
    }

  }

  ###################################
  ## MAKE A DATAFRAME FOR PLOTTING ##
  ###################################

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
