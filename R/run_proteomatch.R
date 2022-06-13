#' Run the ProteoMatch pipeline with a Proteoform, mzML, and Settings file
#'
#' @description A wrapper for quickly running all the steps of the ProteoMatch pipeline
#'
#' @param ProteoformFile Path to a csv with a column of "Proteoforms" and their "Protein"
#'    identifier. Each Proteoform-Protein pairing should be unique. Proteoforms
#'    should be written with proper annotations (i.e. "M.AA`[`Acetyl`]`AA`[`3.2`]`.V"). Required.
#' @param mzMLFile Path to an mzML file containing one spectra. The file must be centroided.
#'    Use the \url{https://github.com/wilsontom/msconverteR} package and
#'    set msconvert_args = 'peakPicking true 1-'. Required.
#' @param SettingsFile Path to a text file with all user-set parameters. Required.
#'
#' @export
run_proteomatch <- function(ProteoformFile,
                            mzMLFile,
                            SettingsFile) {

  ###########
  ##




}
