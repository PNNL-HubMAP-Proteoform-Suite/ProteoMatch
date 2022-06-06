# Get downloads folder
#' @export
.getDownloadsFolder <- function(.test_mode = FALSE) {
  if (Sys.info()['sysname'] == "Windows" | .test_mode) {
    folder <- dirname("~")
    return(folder)
  } else {
    folder <- path.expand("~")
    folder <- file.path(folder, "Downloads")
    folder <- paste0(folder, .Platform$file.sep)
    return(folder)
  }
}

#' Generate a trelliscope display for all Proteoform to MS1 Matches
#'
#' @description [Trelliscope](https://github.com/hafen/trelliscopejs) allows
#'    users to filter and sort plots based on cognostic variables. For this dataset,
#'    the plots are generated with the plot_Ms1Match function, and the cognostics
#'    are: Protein, Absolute Relative Error, Correlation, Charge, Proteoform, and ID.
#'
#' @param PeakData A pspecterlib peak_data object or data.table with "M/Z" and "Intensity". Required.
#' @param Ms1Match A ProteoMatch_MatchedPeaks class object from match_full_seq_ms1. Required.
#' @param Path The base directory of the trelliscope application. Default is Downloads/Ms1Match.
#' @param Window The -/+ m/z value on either side of the matched spectra. Default is 2 m/z.
#'
#' @returns An html trelliscope display
#'
#' @export
build_trelliscope_display <- function(PeakData,
                                      Ms1Match,
                                      Path = file.path(.getDownloadsFolder(), "Ms1Match"),
                                      Window = 2,
                                      ...) {

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

  # Check the +/- window range
  if (!is.numeric(Window) | length(Window) != 1) {
    stop("Window must be a single numeric.")
  }
  Window <- abs(Window)





}




