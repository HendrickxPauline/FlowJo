# helpers.R
# Utility functions for the Flow Cytometry Shiny app.

#' Load all .fcs files from a directory into a flowSet.
#' Returns NULL (invisibly) if no .fcs files are found or loading fails.
load_fcs_folder <- function(folder_path) {
  fcs_files <- list.files(folder_path, pattern = "\\.fcs$",
                           full.names = TRUE, ignore.case = TRUE)
  if (length(fcs_files) == 0) return(NULL)

  tryCatch(
    flowCore::read.flowSet(fcs_files,
                           transformation      = FALSE,
                           truncate_max_range  = FALSE),
    error = function(e) {
      message("Error loading .fcs files: ", conditionMessage(e))
      NULL
    }
  )
}

#' Extract all events for one channel from a flowSet into a tidy data frame.
#' Returns a data frame with columns: sample (filename) and value (raw intensity).
extract_channel_data <- function(fs, channel) {
  do.call(rbind, lapply(flowCore::sampleNames(fs), function(sname) {
    data.frame(
      sample = sname,
      value  = flowCore::exprs(fs[[sname]])[, channel],
      stringsAsFactors = FALSE
    )
  }))
}
