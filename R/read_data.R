#' Read Incucyte CSV data files
#'
#' Reads all CSV files from a folder containing Incucyte object-level data
#' and combines them into a single data frame.
#'
#' @param folder Path to folder containing Incucyte CSV files
#'
#' @return A data frame with all object measurements and a `filename` column
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read_incucyte("path/to/data/folder")
#' }
read_incucyte <- function(folder) {
  files <- list.files(folder, pattern = "*.csv", full.names = TRUE)
  message("Reading ", length(files), " files from ", folder)

  data_list <- lapply(files, function(file) {
    df <- readr::read_csv(file, show_col_types = FALSE)
    df$filename <- basename(file)
    return(df)
  })

  combined_data <- dplyr::bind_rows(data_list)
  message("Read ", nrow(combined_data), " total objects")
  return(combined_data)
}


#' Tidy Incucyte data
#'
#' Parses filenames to extract timepoint information and optionally joins
#' platemap metadata. Expects filenames in format: `dye_dX_Xh00_suffix.csv`
#'
#' @param data Data frame from [read_incucyte()]
#' @param platemap Optional platemap data frame from [read_incucyte_platemap()]
#'
#' @return A tidied data frame with parsed time columns and platemap metadata
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read_incucyte("Data")
#' platemap <- read_incucyte_platemap("experiment.platemap")
#' tidy_data <- tidy_incucyte(data, platemap)
#' }
tidy_incucyte <- function(data, platemap = NULL) {

  tidied <- data %>%
    tidyr::separate_wider_delim(
      "filename",
      delim = "_",
      names = c("dye", "day", "time", "drop")
    ) %>%
    dplyr::select(-"drop", -"dye") %>%
    tidyr::separate_wider_delim(
      "time",
      delim = "h",
      names = c("h", "m")
    ) %>%
    dplyr::mutate(
      day = as.numeric(gsub("d", "", .data$day)),
      h = as.numeric(.data$h)
    ) %>%
    dplyr::select(-"m") %>%
    dplyr::mutate(h = .data$day * 24 + .data$h) %>%
    dplyr::select(-"day") %>%
    dplyr::mutate(
      `Image Column` = stringr::str_pad(.data$`Image Column`, width = 2, side = "left", pad = "0")
    ) %>%
    tidyr::unite("Wells", c("Image Row", "Image Column"), sep = "", remove = FALSE)

  if (!is.null(platemap)) {
    tidied <- tidied %>%
      dplyr::left_join(platemap, by = c("Wells" = "well"))

    platemap_cols <- setdiff(colnames(platemap), "well")

    incomplete_data <- tidied %>%
      dplyr::filter(dplyr::if_any(dplyr::all_of(platemap_cols), is.na)) %>%
      dplyr::distinct(.data$Wells, dplyr::across(dplyr::all_of(platemap_cols)))

    if (nrow(incomplete_data) > 0) {
      incomplete_wells <- unique(incomplete_data$Wells)
      n_rows_affected <- tidied %>%
        dplyr::filter(.data$Wells %in% incomplete_wells) %>%
        nrow()

      message("\n--------------------------------------------------")
      message("Note: ", length(incomplete_wells), " wells have incomplete platemap data")
      message("Affected wells: ", paste(sort(incomplete_wells), collapse = ", "))
      message("Rows affected: ", n_rows_affected, " (",
              round(n_rows_affected / nrow(tidied) * 100, 1), "% of data)")
      message("--------------------------------------------------\n")
    }
  }

  return(tidied)
}
