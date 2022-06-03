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
#' @examples
#' \dontrun{
#'
#' # Run two examples with two charge states
#' MolForms_Test <- calculate_molform(
#'    Proteoform = c("M.SS[Methyl]S.V", "M.SS[6]S[7].V"),
#'    Protein = c("Test1", "Test2"),
#'    Charge = 1:2
#' )
#'
#' # Generate some experimental peak data to match
#' PeakData <- pspecterlib::make_peak_data(
#'    MZ = c(294.1296, 295.1325, 296.1343, 297.1369, 298.1390),
#'    Intensity = c(868.3680036, 110.9431876, 18.7179196, 1.7871629, 0.1701294)
#' )
#'
#' # Run algorithm
#' AllMatches <- match_proteoform_to_ms1(
#'    PeakData = PeakData,
#'    MolecularFormula = MolForms_Test,
#' )
#'
#' # Make plot
#' plot_Ms1Match(PeakData = PeakData, Ms1Match = AllMatches, ID = 1)
#'
#'
#' }
#'
#' @export
plot_Ms1Match <- function(PeakData,
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
    subset(`M/Z` >= min(Ms1Match$`M/Z`) - 1 & `M/Z` <= max(Ms1Match$`M/Z`) + 1)

  # Fix up ms1 match
  MS1 <- data.table::data.table(
    `M/Z` = c(Ms1Match$`M/Z` - 1e-9, Ms1Match$`M/Z`, Ms1Match$`M/Z` + 1e-9),
    Abundance = c(rep(0, nrow(Ms1Match)), Ms1Match$`Intensity`, rep(0, nrow(Ms1Match)))
  ) %>%
    dplyr::arrange(`M/Z`)

  ###############
  ## MAKE PLOT ##
  ###############
  ggplot2::ggplot(AdjPeakData, ggplot2::aes(x = `M/Z`, y = Intensity)) + ggplot2::geom_line() +
    ggplot2::geom_line(data = MS1, ggplot2::aes(x = `M/Z`, y = Abundance), color = "red") +
    ggplot2::theme_bw() + ggplot2::ylab("Abundance")

}
