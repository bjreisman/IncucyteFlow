#' Load Incucyte FCS files into a GatingSet
#'
#' Reads FCS files created by [incucyte_to_fcs()] and creates a GatingSet
#' with parsed metadata from filenames and optional platemap annotation.
#'
#' @param fcs_dir Directory containing FCS files
#' @param platemap Optional platemap data frame from [read_incucyte_platemap()]
#' @param prefix Filename prefix used when creating FCS files
#'
#' @return A GatingSet object with metadata in pData
#' @export
#'
#' @examples
#' \dontrun{
#' gs <- load_incucyte_fcs("fcs_output", platemap, prefix = "experiment1")
#' }
load_incucyte_fcs <- function(fcs_dir,
                              platemap = NULL,
                              prefix = "incucyte") {

  fcs_files <- list.files(fcs_dir, pattern = "\\.fcs$", full.names = TRUE)
  message("Loading ", length(fcs_files), " FCS files...")

  fs <- flowCore::read.flowSet(files = fcs_files)

  filenames <- flowCore::sampleNames(fs)
  basenames <- basename(filenames)

  pattern <- paste0("^", prefix, "_([A-H][0-9]+)_([0-9]+)\\.fcs$")
  parsed <- stringr::str_match(basenames, pattern)

  if (all(is.na(parsed[, 1]))) {
    stop("Filename pattern did not match. Example filename: ", basenames[1])
  }

  metadata <- data.frame(
    name = basenames,
    Wells = parsed[, 2],
    h = as.numeric(parsed[, 3]),
    stringsAsFactors = FALSE
  )

  if (!is.null(platemap)) {
    metadata <- metadata %>%
      dplyr::left_join(platemap, by = c("Wells" = "well"))
  }

  cs <- flowWorkspace::flowSet_to_cytoset(fs)
  gs <- flowWorkspace::GatingSet(cs)

  for (col in colnames(metadata)) {
    if (col != "name") {
      flowWorkspace::pData(gs)[[col]] <- metadata[[col]]
    }
  }

  message("Created GatingSet with ", length(gs), " samples")
  return(gs)
}
