#' Create channel description mapping for Incucyte data
#'
#' Creates a named vector mapping Incucyte channel names to human-readable
#' descriptions for use in FCS file metadata.
#'
#' @param FLR1 Description for FLR1 channel (e.g., "GFP", "Phase")
#' @param FLR2 Description for FLR2 channel (e.g., "TMRM", "mCherry")
#' @param FLR3 Description for FLR3 channel (e.g., "Sytox Deep Red", "Annexin")
#'
#' @return Named vector of channel descriptions
#' @export
#'
#' @examples
#' # Default descriptions
#' desc <- create_channel_desc()
#'
#' # Custom descriptions
#' desc <- create_channel_desc(FLR1 = "GFP", FLR2 = "TMRM", FLR3 = "Sytox Deep Red")
create_channel_desc <- function(FLR1 = "FLR1",
                                FLR2 = "FLR2",
                                FLR3 = "FLR3") {
  c(
    "Area" = "Cell Area",
    "Eccentricity" = "Eccentricity",
    "FLR1IntegratedIntensity" = paste(FLR1, "Integrated"),
    "FLR1MeanIntensity" = paste(FLR1, "Mean"),
    "FLR2IntegratedIntensity" = paste(FLR2, "Integrated"),
    "FLR2MeanIntensity" = paste(FLR2, "Mean"),
    "FLR3IntegratedIntensity" = paste(FLR3, "Integrated"),
    "FLR3MeanIntensity" = paste(FLR3, "Mean"),
    "Texture" = "Texture"
  )
}


#' Convert Incucyte data to FCS files
#'
#' Converts tidied Incucyte data to FCS format files, creating one FCS file
#' per unique combination of grouping variables (default: Wells and timepoint).
#'
#' @param data Tidied Incucyte data frame from [tidy_incucyte()]
#' @param output_dir Directory to save FCS files (created if doesn't exist)
#' @param group_by_cols Columns to group by; one FCS file per unique combination
#' @param channels Vector of channel/column names to include in FCS files
#' @param channel_desc Named vector of channel descriptions from [create_channel_desc()]
#' @param fluor_channels Fluorescence channels for spillover matrix
#' @param filename_prefix Prefix for FCS filenames
#' @param add_spillover Logical; add identity spillover matrix to FCS files
#'
#' @return Character vector of created file paths
#' @export
#'
#' @examples
#' \dontrun{
#' fcs_files <- incucyte_to_fcs(
#'   data = tidy_data,
#'   output_dir = "fcs_output",
#'   group_by_cols = c("Wells", "h"),
#'   channel_desc = create_channel_desc(FLR1 = "GFP", FLR2 = "TMRM"),
#'   filename_prefix = "experiment1"
#' )
#' }
incucyte_to_fcs <- function(data,
                            output_dir = "fcs_output",
                            group_by_cols = c("Wells", "h"),
                            channels = NULL,
                            channel_desc = NULL,
                            fluor_channels = NULL,
                            filename_prefix = "incucyte",
                            add_spillover = TRUE) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Default channels
  if (is.null(channels)) {
    channels <- c(
      "Area",
      "Eccentricity",
      "FLR1IntegratedIntensity",
      "FLR1MeanIntensity",
      "FLR2IntegratedIntensity",
      "FLR2MeanIntensity",
      "FLR3IntegratedIntensity",
      "FLR3MeanIntensity",
      "Texture"
    )
    channels <- channels[channels %in% colnames(data)]
    message("Using channels: ", paste(channels, collapse = ", "))
  }

  # Default channel descriptions
  if (is.null(channel_desc)) {
    channel_desc <- create_channel_desc()
    message("Using default channel descriptions. Use create_channel_desc() for custom labels.")
  }

  # Default fluorescence channels
  if (is.null(fluor_channels)) {
    fluor_channels <- c(
      "FLR1IntegratedIntensity",
      "FLR1MeanIntensity",
      "FLR2IntegratedIntensity",
      "FLR2MeanIntensity",
      "FLR3IntegratedIntensity",
      "FLR3MeanIntensity"
    )
    fluor_channels <- fluor_channels[fluor_channels %in% channels]
    message("Fluorescence channels: ", paste(fluor_channels, collapse = ", "))
  }

  # Create identity spillover matrix for fluorescence channels
  if (add_spillover) {
    n_fluor <- length(fluor_channels)
    spillover_matrix <- diag(n_fluor)
    rownames(spillover_matrix) <- fluor_channels
    colnames(spillover_matrix) <- fluor_channels
    message("Adding identity spillover matrix for ", n_fluor, " fluorescence channels")
  }

  grouped <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by_cols))) %>%
    dplyr::group_split()

  group_keys <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by_cols))) %>%
    dplyr::group_keys()

  message("Creating ", nrow(group_keys), " FCS files...")

  filepaths <- purrr::map2_chr(grouped, seq_len(nrow(group_keys)), function(df, i) {

    key_vals <- group_keys[i, ] %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
      unlist()
    group_label <- paste(key_vals, collapse = "_")
    group_label <- gsub("[^A-Za-z0-9_-]", "_", group_label)
    filename <- paste0(filename_prefix, "_", group_label, ".fcs")

    mat <- as.matrix(df[, channels, drop = FALSE])
    mode(mat) <- "numeric"
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]

    if (nrow(mat) == 0) {
      warning("No valid data for: ", filename)
      return(NA_character_)
    }

    # Get descriptions
    descriptions <- ifelse(channels %in% names(channel_desc),
                           channel_desc[channels],
                           channels)

    params <- Biobase::AnnotatedDataFrame(
      data = data.frame(
        name = channels,
        desc = descriptions,
        range = apply(mat, 2, function(x) diff(range(x, na.rm = TRUE))),
        minRange = apply(mat, 2, min, na.rm = TRUE),
        maxRange = apply(mat, 2, max, na.rm = TRUE),
        row.names = paste0("$P", seq_along(channels))
      )
    )

    ff <- flowCore::flowFrame(exprs = mat, parameters = params)

    # Build keywords
    keywords <- list(
      "$FIL" = filename,
      "$SRC" = "Incucyte",
      "$DATE" = format(Sys.Date(), "%d-%b-%Y"),
      "ORIGINALGUID" = filename,
      "INCUCYTE_FLUOR_CHANNELS" = paste(fluor_channels, collapse = ",")
    )

    for (col in group_by_cols) {
      keywords[[paste0("INCUCYTE_", toupper(col))]] <- as.character(key_vals[col])
    }

    if ("CellType" %in% colnames(df)) {
      keywords[["INCUCYTE_CELLTYPE"]] <- unique(df$CellType)[1]
    }
    if ("Compound" %in% colnames(df)) {
      keywords[["INCUCYTE_COMPOUND"]] <- unique(df$Compound)[1]
    }
    if ("GrowthCondition" %in% colnames(df)) {
      keywords[["INCUCYTE_GROWTHCONDITION"]] <- unique(df$GrowthCondition)[1]
    }

    # Add spillover matrix as keyword
    if (add_spillover) {
      keywords[["SPILL"]] <- spillover_matrix
    }

    flowCore::keyword(ff) <- keywords

    filepath <- file.path(output_dir, filename)
    flowCore::write.FCS(ff, filepath)

    if (i %% 100 == 0) message("  Created ", i, " of ", nrow(group_keys), " files")

    return(filepath)
  })

  message("Done! Created ", sum(!is.na(filepaths)), " FCS files in ", output_dir)
  return(filepaths[!is.na(filepaths)])
}
