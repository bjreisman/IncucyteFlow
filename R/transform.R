#' Transform GatingSet using flowVS-optimized asinh transformation
#'
#' Estimates optimal cofactors using flowVS, then applies asinh transformation.
#' Returns a list containing the transformed GatingSet, estimated cofactors,
#' and a transformerList for proper axis display in plots.
#'
#' @param gs GatingSet to transform
#' @param channels Character vector of channels to transform
#' @param cofactors Named numeric vector of pre-computed cofactors; if NULL, will estimate
#' @param select List of selection criteria for samples used to estimate cofactors
#'
#' @return A list with components:
#'   \item{gs}{Transformed GatingSet}
#'   \item{cofactors}{Named vector of cofactors used}
#'   \item{trans}{transformerList for axis display}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- transform_gs_flowVS(gs,
#'   channels = c("FLR1MeanIntensity", "FLR2MeanIntensity"))
#' gs_transformed <- result$gs
#' }
transform_gs_flowVS <- function(gs,
                                channels = c("FLR1MeanIntensity",
                                             "FLR2MeanIntensity",
                                             "FLR3MeanIntensity"),
                                cofactors = NULL,
                                select = NULL) {

  if (!requireNamespace("flowVS", quietly = TRUE)) {
    stop("Package 'flowVS' is required. Install from Bioconductor.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }

  # Estimate cofactors if not provided
  if (is.null(cofactors)) {
    if (!is.null(select)) {
      # Use CytoExploreR if available for selection
      if (requireNamespace("CytoExploreR", quietly = TRUE)) {
        gs_subset <- do.call(CytoExploreR::cyto_select, c(list(x = gs), select))
        message("Estimating cofactors from ", length(gs_subset), " selected samples")
        cs_subset <- flowWorkspace::gs_pop_get_data(gs_subset, "root")
        fs_subset <- flowWorkspace::cytoset_to_flowSet(cs_subset)
      } else {
        message("CytoExploreR not available, using first sample")
        cs <- flowWorkspace::gs_pop_get_data(gs[1], "root")
        fs_subset <- flowWorkspace::cytoset_to_flowSet(cs)
      }
    } else {
      message("Estimating cofactors from first sample...")
      cs <- flowWorkspace::gs_pop_get_data(gs[1], "root")
      fs_subset <- flowWorkspace::cytoset_to_flowSet(cs)
    }

    cofactors <- flowVS::estParamFlowVS(fs_subset, channels = channels)
    message("Estimated cofactors:")
    print(cofactors)
  }

  # Apply transformation to data
  message("Extracting data...")
  cs <- flowWorkspace::gs_pop_get_data(gs, "root")
  fs <- flowWorkspace::cytoset_to_flowSet(cs)

  message("Applying asinh(x/cofactor) transformation...")
  fs_trans <- flowVS::transFlowVS(fs, channels, cofactors)

  # Save original pData
  original_pdata <- flowWorkspace::pData(gs)

  # Create new GatingSet from transformed data
  message("Creating transformed GatingSet...")
  cs_trans <- flowWorkspace::flowSet_to_cytoset(fs_trans)
  gs_trans <- flowWorkspace::GatingSet(cs_trans)

  # Restore pData
  flowWorkspace::pData(gs_trans) <- original_pdata

  # Create transformerList for axis display
  make_asinh_trans <- function(cofactor) {

    trans_func <- function(x) asinh(x / cofactor)
    inv_func <- function(x) sinh(x) * cofactor

    breaks_func <- function(x) {
      rng <- range(x, na.rm = TRUE, finite = TRUE)
      if (any(!is.finite(rng))) return(c(0, 1, 10, 100))

      max_pow <- ceiling(log10(max(abs(rng[2]), 1)))
      pos_breaks <- c(0, 10^seq(-1, max_pow, by = 1))
      breaks <- pos_breaks[pos_breaks >= rng[1] & pos_breaks <= rng[2]]
      if (length(breaks) < 2) breaks <- c(0, 1, 10, 100)
      return(breaks)
    }

    minor_breaks_func <- function(b, limits, n) {
      rng <- limits
      if (any(!is.finite(rng))) return(numeric(0))

      max_pow <- ceiling(log10(max(abs(rng[2]), 1)))
      minor <- c()
      for (pow in seq(-2, max_pow, by = 1)) {
        decade_minors <- c(2, 3, 4, 5, 6, 7, 8, 9) * 10^pow
        minor <- c(minor, decade_minors)
      }
      minor <- minor[minor >= rng[1] & minor <= rng[2]]
      return(minor)
    }

    scales::trans_new(
      name = "asinh",
      transform = trans_func,
      inverse = inv_func,
      breaks = breaks_func,
      minor_breaks = minor_breaks_func,
      domain = c(-Inf, Inf)
    )
  }

  trans_list <- lapply(channels, function(ch) {
    make_asinh_trans(cofactors[ch])
  })
  names(trans_list) <- channels
  class(trans_list) <- c("transformerList", "list")

  message("Done!")

  list(
    gs = gs_trans,
    cofactors = cofactors,
    trans = trans_list
  )
}
