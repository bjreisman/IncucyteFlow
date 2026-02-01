#' Parse Incucyte .platemap XML to tidy data frame
#'
#' Reads an Incucyte platemap XML file and extracts well-level metadata
#' including cell type, compound, concentration, and growth conditions.
#'
#' @param file_path Path to the .platemap XML file
#'
#' @return A data frame with columns:
#'   \item{well}{Well ID (e.g., "A01", "B12")}
#'   \item{Compound}{Compound name}
#'   \item{Concentration}{Compound concentration}
#'   \item{ConcentrationUnits}{Units for concentration}
#'   \item{GrowthCondition}{Growth condition name}
#'   \item{CellType}{Cell type name}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' platemap <- read_incucyte_platemap("experiment.platemap")
#' }
read_incucyte_platemap <- function(file_path) {
  doc <- XML::xmlParse(file_path)
  well_nodes <- XML::getNodeSet(doc, "//well")

  message("Processing ", length(well_nodes), " wells")

  wells_list <- lapply(seq_along(well_nodes), function(i) {
    node <- well_nodes[[i]]

    # Get row/col from XML attributes
    row <- as.integer(XML::xmlGetAttr(node, "row")) + 1
    col <- as.integer(XML::xmlGetAttr(node, "col")) + 1
    well_id <- paste0(LETTERS[row], sprintf("%02d", col))

    # Initialize
    cell_type <- compound <- concentration <- units <- growth_condition <- NA_character_

    # Extract wellItem elements
    well_items <- XML::getNodeSet(node, ".//wellItem")

    for (item_node in well_items) {
      item_type <- XML::xmlGetAttr(item_node, "type")

      if (item_type == "Compound") {
        concentration <- as.numeric(XML::xmlGetAttr(item_node, "concentration"))
        units <- XML::xmlGetAttr(item_node, "concentrationUnits")

        ref_item <- XML::getNodeSet(item_node, ".//referenceItem")
        if (length(ref_item) > 0) {
          compound <- XML::xmlGetAttr(ref_item[[1]], "displayName")
        }

      } else if (item_type == "GrowthCondition") {
        ref_item <- XML::getNodeSet(item_node, ".//referenceItem")
        if (length(ref_item) > 0) {
          growth_condition <- XML::xmlGetAttr(ref_item[[1]], "displayName")
        }
      } else if (item_type == "CellType") {
        ref_item <- XML::getNodeSet(item_node, ".//referenceItem")
        if (length(ref_item) > 0) {
          cell_type <- XML::xmlGetAttr(ref_item[[1]], "displayName")
        }
      }
    }

    data.frame(
      well = well_id,
      Compound = ifelse(is.null(compound), NA_character_, compound),
      Concentration = concentration,
      ConcentrationUnits = ifelse(is.null(units), "", units),
      GrowthCondition = ifelse(is.null(growth_condition), "", growth_condition),
      CellType = ifelse(is.null(cell_type), "", cell_type),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, wells_list)
}


#' Convert tidy platemap to plater well-list CSV format
#'
#' Exports a platemap data frame to CSV format with one row per well.
#'
#' @param platemap_df Data frame from [read_incucyte_platemap()]
#' @param output_file Path for output CSV file
#'
#' @return The expanded data frame (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' platemap <- read_incucyte_platemap("experiment.platemap")
#' write_plater_csv(platemap, "platemap_welllist.csv")
#' }
write_plater_csv <- function(platemap_df, output_file) {

  # Ensure all 96 wells present (fill blanks)
  all_wells <- expand.grid(row = LETTERS[1:8], col = sprintf("%02d", 1:12))
  all_wells$well <- paste0(all_wells$row, all_wells$col)

  df_plater <- platemap_df %>%
    dplyr::right_join(all_wells, by = "well") %>%
    dplyr::mutate(
      Compound = ifelse(is.na(.data$Compound), "Blank", .data$Compound),
      Concentration = ifelse(is.na(.data$Concentration), NA_real_, .data$Concentration),
      ConcentrationUnits = ifelse(is.na(.data$ConcentrationUnits), "", .data$ConcentrationUnits),
      GrowthCondition = ifelse(is.na(.data$GrowthCondition), "Blank", .data$GrowthCondition),
      CellType = ifelse(is.na(.data$CellType), "Blank", .data$CellType)
    ) %>%
    dplyr::select("well", "Compound", "Concentration", "ConcentrationUnits",
                  "GrowthCondition", "CellType") %>%
    dplyr::arrange(.data$well)

  utils::write.csv(df_plater, output_file, row.names = FALSE)

  message("Wrote ", nrow(df_plater), " wells to: ", output_file)

  invisible(df_plater)
}


#' Convert well-list CSV to plater plate-shaped format
#'
#' Converts a well-list CSV (one row per well) to plater's plate-shaped
#' format where each variable is represented as an 8x12 grid.
#'
#' @param welllist_csv Path to well-list CSV file
#' @param output_file Path for output plater-format CSV
#'
#' @return NULL (called for side effects)
#' @export
#'
#' @examples
#' \dontrun{
#' convert_to_plater_format("platemap_welllist.csv", "platemap_plater.csv")
#' }
convert_to_plater_format <- function(welllist_csv, output_file) {

  df <- utils::read.csv(welllist_csv, stringsAsFactors = FALSE)

  df$row_letter <- substr(df$well, 1, 1)
  df$col_num <- as.integer(substr(df$well, 2, 3))

  column_names <- names(df)[!names(df) %in% c("well", "row_letter", "col_num")]

  con <- file(output_file, "w")
  on.exit(close(con))

  for (i in seq_along(column_names)) {
    col_name <- column_names[i]

    plate_shaped <- df %>%
      dplyr::select("row_letter", "col_num", dplyr::all_of(col_name)) %>%
      tidyr::pivot_wider(
        names_from = "col_num",
        values_from = dplyr::all_of(col_name),
        values_fn = list
      ) %>%
      dplyr::arrange(factor(.data$row_letter, levels = LETTERS[1:8]))

    plate_shaped <- plate_shaped %>%
      dplyr::mutate(dplyr::across(
        -"row_letter",
        ~sapply(., function(x) if(is.list(x)) as.character(x[[1]]) else as.character(x))
      ))

    writeLines(paste(col_name, paste(1:12, collapse = ","), sep = ","), con)

    for (j in 1:nrow(plate_shaped)) {
      row_name <- plate_shaped$row_letter[j]
      row_values <- as.character(plate_shaped[j, -1])
      writeLines(paste(row_name, paste(row_values, collapse = ","), sep = ","), con)
    }

    if (i < length(column_names)) {
      writeLines("", con)
    }
  }

  message("Converted to plater format: ", output_file)
  message("Variables included: ", paste(column_names, collapse = ", "))

  invisible(NULL)
}
