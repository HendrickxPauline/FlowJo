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

#' Extract two channels from a flowSet for a dot plot.
#' Returns a data frame with columns: sample, x, y.
extract_two_channels <- function(fs, ch_x, ch_y) {
  do.call(rbind, lapply(flowCore::sampleNames(fs), function(sname) {
    m <- flowCore::exprs(fs[[sname]])
    data.frame(
      sample = sname,
      x      = m[, ch_x],
      y      = m[, ch_y],
      stringsAsFactors = FALSE
    )
  }))
}

#' Estimate the 2-D kernel density at each (x, y) point.
#' Returns a numeric vector the same length as x (and y).
#' Uses MASS::kde2d internally; returns zeros on failure.
point_density <- function(x, y, n = 100L) {
  ok  <- is.finite(x) & is.finite(y)
  out <- rep(0, length(x))
  if (sum(ok) < 10L) return(out)
  tryCatch({
    dens    <- MASS::kde2d(x[ok], y[ok], n = n)
    ix      <- pmax(1L, pmin(findInterval(x[ok], dens$x), length(dens$x)))
    iy      <- pmax(1L, pmin(findInterval(y[ok], dens$y), length(dens$y)))
    out[ok] <- dens$z[cbind(ix, iy)]
    out
  }, error = function(e) out)
}
