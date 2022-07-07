#' Matches experimental peaks to calculated isotope profile
#'
#' @description Returns the "ProteoMatch_MatchedPeaks" object with peaks that
#'     have matched to the isotoping profile. This object can be passed to downstream
#'     plotting functions.
#'
#' @param PeakData A pspecterlib peak_data object made with "make_peak_data" or
#'     extracted from a raw or mzML file with "get_peak_data." Use centroided data. Required.
#' @param MolecularFormulas A "ProteoMatch_MolForm" object with Molecular Formulas,
#'     Mass Shifts, and Charges.
#' @param IsotopicPercentage The minimum isotopic percentage (calculated intensity) permitted
#'     to be matched. Default is 1, which is 1%.
#' @param PPMThreshold The maximum m/z error permitted. Default is 10 ppm.
#' @param MinAbundanceChange An abundance (every peak is scaled to the largest peak)
#'     absolute change required to count a subsequent peak as an isotope. Default is 0.1.
#' @param IsotopeRange The minimum and maximum number of isotopes to consider. Default is c(5,20).
#' @param ProtonMass The AMU mass of a proton. Default is 1.00727647.
#'
#' @details
#' The data.table outputted by this function contains 12 columns
#' \tabular{ll}{
#' Protein \tab The provided protein identifier \cr
#' \tab \cr
#' M/Z \tab The calculated M/Z of a particular isotope \cr
#' \tab \cr
#' Intensity \tab The calculated relative intensity of a particular isotope \cr
#' \tab \cr
#' Isotope \tab The isotope number "n" in M+n \cr
#' \tab \cr
#' M/Z Tolerance \tab Given the ppm, the M/Z tolerance is the window in which peaks are searched and matched \cr
#' \tab \cr
#' M/Z Experimental \tab The M/Z of the matched from the experimental peaks \cr
#' \tab \cr
#' Intensity Experimental \tab The intensity of the matched peak from the experimental peaks \cr
#' \tab \cr
#' PPM Error \tab How far off the calculated peak is from the experimental peak. Only peaks within the ppm threshold are kept \cr
#' \tab \cr
#' Absolute Relative Error \tab The sum of the absolute relative difference of calculated and experimental peak intensities. \cr
#' \tab \cr
#' Correlation \tab The cosine correlation of calculated and experimental peak intensities \cr
#' \tab \cr
#' Charge \tab The provided charges \cr
#' \tab \cr
#' Proteoform \tab The provided proteoform \cr
#' \tab \cr
#' ID \tab A unique ID for each Proteoform, Protein, and Charge combination used in plotting functions \cr
#' \tab \cr
#' }
#'
#' @returns A ProteoMatch_MatchedPeaks object, which is a data.table containing the
#'     Protein, M/Z, Intensity, Isotope, M/Z Tolerance, Experimental M/Z,
#'     Experimental Intensity, PPM Error, Absolute Relative Error, Correlation,
#'     Charge, and Proteoform.
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
#' match_proteoform_to_ms1(
#'    PeakData = PeakData,
#'    MolecularFormula = MolForms_Test,
#'    IsotopeRange = c(3,20)
#' )
#'
#' }
#'
#' @export
match_proteoform_to_ms1 <- function(PeakData,
                                    MolecularFormulas,
                                    IsotopicPercentage = 1,
                                    PPMThreshold = 10,
                                    MinAbundanceChange = 0.1,
                                    IsotopeRange = c(5, 20),
                                    ProtonMass = 1.00727647) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Check that peak data is of the appropriate class
  if ("peak_data" %in% class(PeakData) == FALSE) {
    stop("PeakData must be a pspecterlib peak_data object.")
  }

  # Check that molecular formula is a string
  if (inherits(MolecularFormulas, "ProteoMatch_MolForm") == FALSE) {
    stop("MolecularFormula must be a ProteoMatch_MolForm object.")
  }

  # Check that max isotope is a numeric
  if (!is.numeric(IsotopicPercentage)) {
    stop("IsotopicPercentage must be numeric.")
  }
  IsotopicPercentage <- abs(IsotopicPercentage)
  if (IsotopicPercentage > 100) {IsotopicPercentage <- 100}

  # PPM Threshold
  if (!is.numeric(PPMThreshold) || PPMThreshold == 0) {
    stop("PPMThreshold must be a nonzero numeric.")
  }
  PPMThreshold <- abs(PPMThreshold)

  # MinAbundanceChange should be a numeric value
  if (!is.numeric(MinAbundanceChange) || MinAbundanceChange < 0 | MinAbundanceChange > 100) {
    stop("MinAbundanceChange should be a numeric between 0 and 100, inclusive.")
  }

  # Check that max isotope is a numeric
  if (!is.numeric(IsotopeRange) | length(unique(abs(round(IsotopeRange)))) != 2) {
    stop("IsotopeRange must be a numeric with at least two unique values.")
  }
  IsotopeRange <- abs(round(IsotopeRange))

  ##################
  ## RUN ITERATOR ##
  ##################

  .match_proteoform_to_ms1_iterator <- function(PeakData,
                                                MolForm,
                                                MassShift,
                                                MonoisotopicMass,
                                                Charge,
                                                Protein,
                                                Proteoform,
                                                IsotopicPercentage,
                                                PPMThreshold,
                                                MaxIsotopes) {

    ########################
    ## CALCULATE ISOTOPES ##
    ########################

    # Get Isotopes
    Isotopes <- Rdisop::getMolecule(formula = MolForm, maxisotopes = max(IsotopeRange))

    # Pull isotope distribution
    IsoDist <- Isotopes$isotopes[[1]] %>%
      t() %>%
      data.table::data.table() %>%
      dplyr::rename(`M/Z` = V1, Intensity = V2) %>%
      dplyr::mutate(
        `M/Z` = (`M/Z` + (Charge * ProtonMass)) / Charge,
        Isotope = (1:length(Intensity)) - 1
      )

    # Add mass change if it is not NULL
    if (!is.null(MassShift) & MassShift != 0) {
      IsoDist$`M/Z` <- IsoDist$`M/Z` + (MassShift / Charge)
    }

    ####################
    ## MATCH ISOTOPES ##
    ####################

    # Determine the theoertical mz tolerance
    IsoDist$`M/Z Tolerance` <- IsoDist$`M/Z` * (PPMThreshold / 1e6)

    # For each theoretical peak, find the closest index in ms, where ms = theoretical
    LeftIndex <- findInterval(IsoDist$`M/Z`, PeakData$`M/Z`, rightmost.closed = FALSE, all.inside = TRUE)

    # Compute mz differences (absolute) to closest element to each side, smaller to the left and next greater to the right:
    IsoDist$`Left Difference` <- abs(PeakData$`M/Z`[LeftIndex] - IsoDist$`M/Z`)
    IsoDist$`Right Difference` <- abs(PeakData$`M/Z`[LeftIndex + 1] - IsoDist$`M/Z`)
    IsoDist$`Closest Index` <- LeftIndex

    # Set closest index as right side one, if difference is smaller:
    RightIndexBest <- which(IsoDist$`Right Difference` < IsoDist$`Left Difference`)
    IsoDist$`Closest Index`[RightIndexBest] <- IsoDist$`Closest Index`[RightIndexBest] + 1
    IsoDist$`M/Z Difference` <- abs(PeakData$`M/Z`[IsoDist$`Closest Index`] - IsoDist$`M/Z`)

    # Keep only matches within the tolerance
    IsoDist <- IsoDist[which(IsoDist$`M/Z Difference` < IsoDist$`M/Z Tolerance`), ]
    IsoDist$`M/Z Experimental` <- PeakData$`M/Z`[IsoDist$`Closest Index`]
    IsoDist$`Intensity Experimental` <- PeakData$Intensity[IsoDist$`Closest Index`]
    IsoDist$`Abundance Experimental` <- PeakData$Abundance[IsoDist$`Closest Index`]

    # Remove non-necessary rows moving forward
    IsoDist <- IsoDist %>% dplyr::select(-c(`Left Difference`, `Right Difference`, `Closest Index`, `M/Z Difference`))

    # Calculate PPM Error
    IsoDist$`PPM Error` <- ((IsoDist$`M/Z Experimental` - IsoDist$`M/Z`) / IsoDist$`M/Z`) * 1e6

    # Subset down to numbers that are within 1 place of each other
    CloseValues <- IsoDist %>%
      dplyr::select(Isotope) %>%
      dplyr::mutate(
        Order = Isotope - dplyr::lag(Isotope),
        Order = ifelse(is.na(Order), 1, Order),
        Take = Order == dplyr::lag(Order) & Order == 1,
        Take = ifelse(is.na(Take), TRUE, Take),
        Take = ifelse(Take == TRUE, ifelse(dplyr::lag(Take) == FALSE, FALSE, TRUE), Take),
        Take = ifelse(is.na(Take), TRUE, Take)
      )

    IsoDist <- IsoDist[CloseValues$Take,]

    if (nrow(IsoDist) < min(IsotopeRange)) {return(NULL)}

    ################################
    ## FILTER ISOTOPIC PERCENTAGE ##
    ################################

    # Filter by isotopes with an intensity percentage
    IsoDist <- IsoDist[IsoDist$Intensity >= IsotopicPercentage/100,]

    ##########################
    ## CALCULATED ABUNDANCE ##
    ##########################

    # Get max calculated intensity and max measured abundance
    calcInten <- unlist(IsoDist$Intensity)[which.max(unlist(IsoDist$Intensity))]
    maxAbun <- unlist(IsoDist$Abundance)[which.max(unlist(IsoDist$Abundance))]

    # Determine scale and scale intensity
    scalingFactor <- maxAbun / calcInten
    IsoDist$Abundance <- IsoDist$Intensity * scalingFactor

    ######################
    ## ABUNDANCE FILTER ##
    ######################

    # Flag abundance changes that are greater than the threshold
    IsoDist <- IsoDist %>%
      dplyr::mutate(
        `Abundance Diff` = `Abundance Experimental` -
          dplyr::lag(`Abundance Experimental`, default = dplyr::first(`Abundance Experimental`)),
        Flag = abs(`Abundance Diff`) >= MinAbundanceChange | `Abundance Diff` == 0
      )

    # Determine where to subset from
    if (FALSE %in% IsoDist$Flag) {
      IsoDist <- IsoDist[1:(min(which(IsoDist$Flag == FALSE))-1),]
    }

    ######################
    ## CALCULATE SCORES ##
    ######################

    # Add correlation score
    if (nrow(IsoDist) >= min(IsotopeRange)) {

      # Absolute Relative Error and Cosine Correlation
      IsoDist$AbsRelError <- 1/nrow(IsoDist) * sum(abs(IsoDist$Abundance - IsoDist$`Abundance Experimental`) / IsoDist$Abundance)
      IsoDist$Correlation <- lsa::cosine(IsoDist$`Abundance Experimental`, IsoDist$Abundance)[1,1]

      # Figure of merit
      IsoDist$`Figure of Merit` <- nrow(IsoDist) / (sum((IsoDist$Abundance - IsoDist$`Abundance Experimental`)^2 + attributes(PeakData)$pspecter$MinimumAbundance^2))
      IsoDist$`Figure of Merit` <- ifelse(is.infinite(IsoDist$`Figure of Merit`), NA, IsoDist$`Figure of Merit`)

      # Add missing columns and reorder
      IsoDist <- IsoDist %>%
        dplyr::mutate(
          Protein = Protein,
          `Absolute Relative Error` = AbsRelError,
          Charge = Charge,
          Proteoform = Proteoform,
          `Monoisotopic Mass` = MonoisotopicMass,
        ) %>%
        dplyr::select(
          Protein, `M/Z`, `Monoisotopic Mass`, Intensity, Abundance, Isotope, `M/Z Tolerance`, `M/Z Experimental`,
          `Intensity Experimental`, `Abundance Experimental`, `PPM Error`, `Absolute Relative Error`,
          Correlation, `Figure of Merit`, Charge, Proteoform
        ) %>%
        dplyr::mutate(ID = uuid::UUIDgenerate())

      return(IsoDist)

    } else {return(NULL)}

  }

  # Iterate through molecular formulas
  AllMatches <- do.call(rbind, lapply(1:nrow(MolecularFormulas), function(el) {
    .match_proteoform_to_ms1_iterator(
      PeakData = PeakData,
      MolForm = MolecularFormulas$`Molecular Formula`[el] %>% unlist(),
      MassShift = MolecularFormulas$`Mass Shift`[el] %>% unlist(),
      MonoisotopicMass = MolecularFormulas$`Monoisotopic Mass`[el] %>% unlist(),
      Charge = MolecularFormulas$Charge[el] %>% unlist(),
      Protein = MolecularFormulas$Protein[el] %>% unlist(),
      Proteoform = MolecularFormulas$Proteoform[el] %>% unlist(),
      IsotopicPercentage = IsotopicPercentage,
      PPMThreshold = PPMThreshold,
      MaxIsotopes = MaxIsotopes
    )
  }))

  # If no matches, return NULL
  if (is.null(AllMatches)) {
    message("No matches detected.")
    return(NULL)
  }

  # Simplify unique ID
  AllMatches$ID <- as.numeric(as.factor(AllMatches$ID))

  # Add class, and return object
  class(AllMatches) <- c(class(AllMatches), "ProteoMatch_MatchedPeaks")

  return(AllMatches)

}
