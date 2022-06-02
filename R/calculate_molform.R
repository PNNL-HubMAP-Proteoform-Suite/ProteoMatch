#' Generate a Molecular Formula Table from proteoforms
#'
#' @description Returns the "ProteoMatch_molform" object with each proteoform's
#'     molecular formula and adjusted mass.
#'
#' @param Proteoform A vector of strings written in proteoform annotation. i.e. "M.AA`[`Acetyl`]`AA`[`3.2`]`.V"
#' @param Protein A vector of identifiers for each protein. Default is NULL.
#' @param Charge The range of charges to test. Default is 1.
#' @param ProtonMass The AMU mass of a proton. Default is 1.00727647.
#'
#' @details
#' The data.table outputted by this function returns 8 columns.
#' \tabular{ll}{
#' Molecular Formula \tab The molecular formula derived from the peptide/protein sequence and its modifications \cr
#' \tab \cr
#' Mass Shift \tab The total mass change indicated in the proteoform \cr
#' \tab \cr
#' Monoisotopic Mass \tab The proteoform's monoisotopic peak based on molecular formula and mass shift \cr
#' \tab \cr
#' Most Abundant Isotope \tab The proteoform's most abundant isotope based on molecular formula and mass shift \cr
#' \tab \cr
#' Protein \tab The provided protein identifier \cr
#' \tab \cr
#' Proteoform \tab The provided proteoform \cr
#' \tab \cr
#' }
#'
#' @returns A data.table containing the molecular formula, mass shift, monoisotopic
#'     mass, most abundant isotope, protein, and proteoform.
#'
#' @examples
#' \dontrun{
#'
#' calculate_molform(
#'    Proteoform = "M.(S)[Acetyl]ATNNIAQARKLVEQLRIEAGIERIKVSKAASDLMSYCEQHARNDPLLVGVPASENPFKDK(KPCIIL)[-52.9879].",
#'    Protein = "O60262",
#'    Charge = 1:3
#' )
#'
#' }
#'
#' @export
calculate_molform <- function(Proteoform,
                              Protein = NULL,
                              Charge = 1,
                              ProtonMass = 1.00727647) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Proteoform should be a string
  if (is.null(Proteoform) || !is.character(Proteoform)) {
    stop("Proteoform must be a vector of characters.")
  }

  # If protein is not NULL...
  if (!is.null(Protein)) {

    # Proteins should be a string
    if (!is.character(Protein)) {
      stop("Protein must be a vector of characters.")
    }

    # Proteoform and Proteins should be the same length
    if (length(Proteoform) != length(Protein)) {
      stop("Proteoform and Protein must be of the same length.")
    }

  }

  # Charge must be numeric and will be rounded to integers
  if (!is.numeric(Charge)) {
    stop("Charge should be numeric.")
  }
  Charge <- unique(round(abs(Charge)))

  # Proton Mass must be a single numeric
  if (length(ProtonMass) != 1 | !is.numeric(ProtonMass)) {
    stop("ProtonMass must be a single numeric.")
  }

  ###################################
  ## MODIFICATION NAMES AND MASSES ##
  ###################################

  # Grab the data within each of the hard brackets []
  Bracketed_Data <- Proteoform %>%
    gsub(pattern = "[", replacement = "^[", fixed = T) %>%
    gsub(pattern= "]", replacement = "]^", fixed = T) %>%
    strsplit("^", fixed = T) %>%
    unlist() %>%
    lapply(function(x) {if (grepl("[", x, fixed = T)) {x}}) %>%
    unlist()

  # Separate out the modifications from the mass changes
  Modifications <- NULL
  MassChanges <- NULL

  for (Mod in Bracketed_Data) {
    Mod <- gsub("\\[|\\]", "", Mod)
    StringTest <- suppressWarnings(is.na(as.numeric(Mod)))
    if (StringTest) {Modifications <- c(Modifications, Mod)} else {MassChanges <- c(MassChanges, as.numeric(Mod))}
  }

  if (!is.null(Modifications)) {
    if (all(Modifications %in% names(Library)) == FALSE) {
      stop(paste("Modification", Modifications, "is not in our library."))
    }
  }

  #######################
  ## CLEAN UP SEQUENCE ##
  #######################

  # First, use the proteoform annotation
  Sequence <- Proteoform

  # Then, remove each of the bracketed pieces of information
  for (Mod in Bracketed_Data) {
    Sequence <- gsub(Mod, "", Sequence, fixed = TRUE)
  }

  # Then, remove any of the remaining parenthesis
  Sequence <- gsub("\\(|\\)", "", Sequence)

  # Get the number of periods in the sequence
  NumPeriods <- str_count(Sequence, "\\.")

  # There shouldn't be more than 2 periods at this point in the pipeline
  if (NumPeriods > 2) {
    stop("There are too many periods in the input sequence.")
  }

  # Split by the periods and take the largest chunk
  SeqSplit <- strsplit(Sequence, "\\.") %>% unlist()
  SeqPosition <- lapply(SeqSplit, nchar) %>% unlist() %>% which.max()
  Sequence <- SeqSplit[SeqPosition]

  ####################################
  ## Generate the Molecular Formula ##
  ####################################

  # Generate a pspecter molecule object
  Formula <- getAtomsFromSeq(Sequence) %>% make_molecule()

  # If there are modifications, add those as well
  if (length(Modifications) > 0) {
    for (Mod in Modifications) {
      Formula <- add_molecules(Formula, Library[[Mod]]$Formula)
    }
  }

  #####################################################
  ## Add monoisotopic mass and most abundant isotope ##
  #####################################################

  # Calculate isotopes
  Isotoping <- Formula$Formula %>% getMolecule()

  # Get monoisotopic mass, which will always be the first one
  MonoMass <- (Isotoping$isotopes[[1]][,1][1] + (Charge * ProtonMass)) / Charge

  # Get the most abundant isotope
  MAI <- (Isotoping$exactmass + (Charge * ProtonMass)) / Charge

  # Add mass changes if they exist
  if (length(MassChanges) > 0) {
    MonoMass <- MonoMass + (sum(MassChanges) / Charge)
    MAI <- MAI + (sum(MassChanges) / Charge)
  }

  return(
    c(
      MolecularFormula = Formula$Formula %>% unlist(),
      MassShift = sum(MassChanges),
      MonoisotopicMass = MonoMass %>% unlist(),
      MostAbundantIsotope = MAI
    )
  )

}
