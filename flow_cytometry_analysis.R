# =============================================================================
# Flow Cytometry Data Analysis
# =============================================================================
# Requirements:
#   BiocManager::install("flowCore")
#   BiocManager::install("ggcyto")
#   install.packages(c("ggplot2", "dplyr", "tidyr"))
# =============================================================================

library(flowCore)
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------

# Path to the folder containing your .fcs files
FCS_DIR <- "data/"

# Output folder for plots and results
OUTPUT_DIR <- "output/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Threshold for "positive" cells (as a percentile of the unstained control,
# or set a fixed value if no control is available)
# Options:
#   "fixed"    -> uses POSITIVE_THRESHOLD_VALUE for all channels
#   "percentile" -> uses the top POSITIVE_PERCENTILE of the first file as cutoff
THRESHOLD_METHOD <- "fixed"
POSITIVE_THRESHOLD_VALUE <- 1000   # adjust per your instrument / panel
POSITIVE_PERCENTILE <- 0.99        # used when THRESHOLD_METHOD = "percentile"

# Channels to EXCLUDE from single-channel plots (e.g. scatter parameters)
SCATTER_CHANNELS <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
                       "Time")

# -----------------------------------------------------------------------------
# 2. READ FCS FILES
# -----------------------------------------------------------------------------

fcs_files <- list.files(FCS_DIR, pattern = "\\.fcs$", full.names = TRUE,
                         ignore.case = TRUE)

if (length(fcs_files) == 0) {
  stop("No .fcs files found in: ", FCS_DIR,
       "\nPlease set FCS_DIR to the folder containing your data.")
}

message("Found ", length(fcs_files), " FCS file(s):")
message(paste(" -", basename(fcs_files), collapse = "\n"))

# Read all files into a flowSet
fs <- read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)
message("\nFlowSet loaded successfully.")
message("Channels available: ", paste(colnames(fs), collapse = ", "))

# -----------------------------------------------------------------------------
# 3. HELPER: extract a tidy data frame from one flowFrame
# -----------------------------------------------------------------------------

fcs_to_df <- function(ff, sample_name) {
  df <- as.data.frame(exprs(ff))
  df$sample <- sample_name
  df
}

# Build a combined data frame (all samples, all channels)
all_data <- lapply(seq_along(fcs_files), function(i) {
  fcs_to_df(fs[[i]], sampleNames(fs)[i])
}) %>% bind_rows()

# -----------------------------------------------------------------------------
# 4. SCATTER PLOTS (FSC vs SSC) — one per sample
# -----------------------------------------------------------------------------

message("\n--- Generating scatter plots (FSC vs SSC) ---")

fsc_col <- grep("^FSC-A", colnames(all_data), value = TRUE)[1]
ssc_col <- grep("^SSC-A", colnames(all_data), value = TRUE)[1]

if (!is.na(fsc_col) && !is.na(ssc_col)) {
  for (sname in unique(all_data$sample)) {
    df_s <- filter(all_data, sample == sname)

    p <- ggplot(df_s, aes(x = .data[[fsc_col]], y = .data[[ssc_col]])) +
      geom_hex(bins = 100) +
      scale_fill_viridis_c(trans = "log1p") +
      labs(title = paste("Scatter plot —", sname),
           x = fsc_col, y = ssc_col) +
      theme_bw()

    out_file <- file.path(OUTPUT_DIR,
                          paste0("scatter_", gsub("[^A-Za-z0-9_]", "_", sname), ".png"))
    ggsave(out_file, p, width = 6, height = 5, dpi = 150)
    message("  Saved: ", out_file)
  }
} else {
  message("  FSC-A / SSC-A columns not found — skipping scatter plots.")
}

# -----------------------------------------------------------------------------
# 5. PER-CHANNEL HISTOGRAMS (fluorescence channels only)
# -----------------------------------------------------------------------------

message("\n--- Generating per-channel histograms ---")

fluor_channels <- setdiff(colnames(fs), SCATTER_CHANNELS)

for (ch in fluor_channels) {
  if (!ch %in% colnames(all_data)) next

  p <- ggplot(all_data, aes(x = .data[[ch]], colour = sample, fill = sample)) +
    geom_density(alpha = 0.3) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
    labs(title = paste("Channel:", ch),
         x = ch, y = "Density") +
    theme_bw() +
    theme(legend.position = "bottom")

  out_file <- file.path(OUTPUT_DIR,
                        paste0("hist_", gsub("[^A-Za-z0-9_]", "_", ch), ".png"))
  ggsave(out_file, p, width = 7, height = 4, dpi = 150)
  message("  Saved: ", out_file)
}

# -----------------------------------------------------------------------------
# 6. MFI (Median Fluorescence Intensity) PER CHANNEL PER SAMPLE
# -----------------------------------------------------------------------------

message("\n--- Calculating MFI ---")

mfi_table <- all_data %>%
  group_by(sample) %>%
  summarise(across(all_of(fluor_channels), median, na.rm = TRUE),
            .groups = "drop")

message("MFI table:")
print(as.data.frame(mfi_table))

write.csv(mfi_table,
          file = file.path(OUTPUT_DIR, "MFI_per_channel.csv"),
          row.names = FALSE)
message("  Saved: ", file.path(OUTPUT_DIR, "MFI_per_channel.csv"))

# -----------------------------------------------------------------------------
# 7. PERCENTAGE POSITIVE CELLS PER CHANNEL PER SAMPLE
# -----------------------------------------------------------------------------

message("\n--- Calculating % positive cells ---")

# Determine thresholds
if (THRESHOLD_METHOD == "percentile") {
  # Use the top N-th percentile of the FIRST sample as cutoff for each channel
  ref_df <- filter(all_data, sample == unique(all_data$sample)[1])
  thresholds <- sapply(fluor_channels, function(ch) {
    quantile(ref_df[[ch]], probs = POSITIVE_PERCENTILE, na.rm = TRUE)
  })
  message("  Thresholds derived from percentile (", POSITIVE_PERCENTILE * 100,
          "%) of sample: ", unique(all_data$sample)[1])
} else {
  thresholds <- setNames(rep(POSITIVE_THRESHOLD_VALUE, length(fluor_channels)),
                         fluor_channels)
  message("  Using fixed threshold: ", POSITIVE_THRESHOLD_VALUE)
}

# Calculate % positive
pct_positive <- all_data %>%
  group_by(sample) %>%
  summarise(
    across(
      all_of(fluor_channels),
      ~ round(mean(. > thresholds[cur_column()], na.rm = TRUE) * 100, 2)
    ),
    .groups = "drop"
  )

message("% Positive cells:")
print(as.data.frame(pct_positive))

write.csv(pct_positive,
          file = file.path(OUTPUT_DIR, "percent_positive.csv"),
          row.names = FALSE)
message("  Saved: ", file.path(OUTPUT_DIR, "percent_positive.csv"))

# -----------------------------------------------------------------------------
# 8. SUMMARY PLOT — MFI heatmap across samples and channels
# -----------------------------------------------------------------------------

message("\n--- Generating MFI summary heatmap ---")

mfi_long <- mfi_table %>%
  pivot_longer(-sample, names_to = "channel", values_to = "MFI")

p_heat <- ggplot(mfi_long, aes(x = channel, y = sample, fill = log1p(MFI))) +
  geom_tile(colour = "white") +
  scale_fill_viridis_c(name = "log1p(MFI)") +
  labs(title = "MFI heatmap", x = "Channel", y = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "MFI_heatmap.png"), p_heat,
       width = max(6, length(fluor_channels) * 0.8), height = max(4, nrow(mfi_table) * 0.6 + 2),
       dpi = 150)
message("  Saved: ", file.path(OUTPUT_DIR, "MFI_heatmap.png"))

# -----------------------------------------------------------------------------
message("\nAnalysis complete. All outputs saved to: ", OUTPUT_DIR)
