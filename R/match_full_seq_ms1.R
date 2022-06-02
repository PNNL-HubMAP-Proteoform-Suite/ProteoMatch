#' Matches peaks to isotoping profile, given the molecular formula and any mass changes
#'
#' @param PeakData A pspecterlib peak_data object or data.table with "M/Z" and "Intensity". Required.
#' @param MolecularFormula The molecular formula written as string with element, numeric.
#'     It should include the modifications' molecular formula calculated in calculate_molform. Required.
#' @param MassShift Any changes to the mass written as a numeric. Default is NULL
#' @param Charge The charge of the mass. Numeric. Default is 1.
#' @param IsotopicPercentage The minimum isotopic percentage (calculated intensity) permitted
#'     to be matched. Default is 0.1, which means 0.1%.
#' @param PPMThreshold The maximum m/z error permitted. Default is 10 ppm.
#' @param MaxIsotopes The maximum number of isotopes to consider. Default is 30.
match_full_seq_ms1 <- function(PeakData,
                               MolecularFormula,
                               MassShift = NULL,
                               Charge = 1,
                               IsotopicPercentage = 1,
                               PPMThreshold = 10,
                               MaxIsotopes = 20) {

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

  # Check that molecular formula is a string
  if (!is.character(MolecularFormula)) {
    stop("MolecularFormula must be a string.")
  }

  # Check that MassShift is numeric, if included
  if (!is.null(MassShift)) {
    if (!is.numeric(MassShift)) {stop("MassShift must be a number.")}
  }

  # Check that charge is a numeric
  if (!is.numeric(Charge)) {
    stop("Charge must be numeric.")
  }
  Charge <- abs(round(Charge))
  if (Charge == 0) {Charge <- 1}

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

  # Check that max isotope is a numeric
  if (!is.numeric(MaxIsotopes)) {
    stop("MaxIsotopes must be numeric.")
  }
  MaxIsotopes <- abs(round(MaxIsotopes))
  if (MaxIsotopes == 0) {MaxIsotopes <- 1}

  ########################
  ## CALCULATE ISOTOPES ##
  ########################

  # Get Isotopes
  Isotopes <- Rdisop::getMolecule(formula = MolecularFormula, maxisotopes = MaxIsotopes)

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
  if (!is.null(MassShift) | MassShift != 0) {
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

  # Remove non-necessary rows moving forward
  IsoDist <- IsoDist %>% dplyr::select(-c(`Left Difference`, `Right Difference`, `Closest Index`, `M/Z Difference`))

  # Calculate PPM Error
  IsoDist$`PPM Error` <- ((IsoDist$`M/Z Experimental` - IsoDist$`M/Z`) / IsoDist$`M/Z`) * 1e6

  # Subset down to numbers that are within 1 place of each other
  CloseValues <- IsoDist %>%
    select(Isotope) %>%
    mutate(
      Order = Isotope - lag(Isotope),
      Order = ifelse(is.na(Order), 1, Order),
      Take = Order == lag(Order) & Order == 1,
      Take = ifelse(is.na(Take), TRUE, Take),
      Take = ifelse(Take == TRUE, ifelse(lag(Take) == FALSE, FALSE, TRUE), Take),
      Take = ifelse(is.na(Take) & Isotope == 0, TRUE, Take)
    )


  IsoDist <- IsoDist[CloseValues$Take,]

  if (nrow(IsoDist) < 3) {return(NULL)}

  ####################
  ## FILTER MATCHES ##
  ####################

  # Filter by intensity
  IsoDist <- IsoDist[IsoDist$Intensity >= IsotopicPercentage/100,]

  # Since there is a high number of matches, stop once the intensity begins to increase at the end
  # of isotope matches
  IsoDist <- IsoDist %>%
    dplyr::mutate(
      `Intensity Diff` = `Intensity Experimental` -
        dplyr::lag(`Intensity Experimental`, default = dplyr::first(`Intensity Experimental`)),
      Flag = `Intensity Diff` >= 0 & `M/Z` > Isotopes$exactmass + (3/Charge)
    )

  # Determine where to subset from
  if (TRUE %in% IsoDist$Flag) {
    IsoDist <- IsoDist[1:(min(which(IsoDist$Flag == TRUE))-1),]
  }

  # Fix intensity to be on the same scale
  IsoDist$Intensity <- IsoDist$Intensity * max(PeakData$Intensity)

  # Add correlation score
  if (nrow(IsoDist) > 3) {

    # Return NULL if there is not at least a 50 between min and max peak
    if (max(IsoDist$`Intensity Experimental`) - min(IsoDist$`Intensity Experimental`) < 50) {
      return(NULL)
    }

    IsoDist$AbsRelError <- 1/nrow(IsoDist) * sum(abs(IsoDist$Intensity - IsoDist$`Intensity Experimental`) / IsoDist$Intensity)
    IsoDist$Correlation <- lsa::cosine(IsoDist$`Intensity Experimental`, IsoDist$Intensity * 1000)[1,1]

  } else {return(NULL)}

  ##########################
  ## RETURN MATCHED PEAKS ##
  ##########################

  return(IsoDist)

}
