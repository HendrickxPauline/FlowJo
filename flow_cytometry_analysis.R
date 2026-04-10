# =============================================================================
# Flow Cytometry Data Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# 1. INSTALL & LOAD PACKAGES  (safe, non-interactive)
# -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c("flowCore", "ggcyto")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "hexbin")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

library(flowCore)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# -----------------------------------------------------------------------------
# 2. CONFIGURATION  — edit these paths and settings
# -----------------------------------------------------------------------------

# Folder containing your .fcs files  (use forward slashes on Windows)
FCS_DIR <- "E:/Exported data/FLOW CYTOMETRY/09042026_MGC_test"

# Where to save plots and CSV results
OUTPUT_DIR <- "E:/Exported data/FLOW CYTOMETRY/09042026_MGC_test/Output_10042026"

# Threshold method for "positive" cells:
#   "fixed"      -> every channel uses POSITIVE_THRESHOLD_VALUE
#   "percentile" -> cutoff = top POSITIVE_PERCENTILE of the first sample
THRESHOLD_METHOD         <- "fixed"
POSITIVE_THRESHOLD_VALUE <- 1000    # adjust to match your instrument / panel
POSITIVE_PERCENTILE      <- 0.99   # only used when THRESHOLD_METHOD = "percentile"

# Scatter / time channels to exclude from fluorescence analysis
SCATTER_CHANNELS <- c("FSC-A", "FSC-H", "FSC-W",
                       "SSC-A", "SSC-H", "SSC-W",
                       "Time")

# -----------------------------------------------------------------------------
# 3. CREATE OUTPUT FOLDER
# -----------------------------------------------------------------------------

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 4. READ FCS FILES
# -----------------------------------------------------------------------------

fcs_files <- list.files(FCS_DIR, pattern = "\\.fcs$",
                         full.names = TRUE, ignore.case = TRUE)

if (length(fcs_files) == 0)
  stop("No .fcs files found in: ", FCS_DIR)

message("Found ", length(fcs_files), " FCS file(s):")
message(paste0("  - ", basename(fcs_files), collapse = "\n"))

fs <- read.flowSet(fcs_files, transformation = FALSE,
                   truncate_max_range = FALSE)

message("\nChannels in data: ", paste(colnames(fs), collapse = ", "))

# Identify fluorescence channels (everything that is not scatter/time)
fluor_channels <- setdiff(colnames(fs), SCATTER_CHANNELS)
message("Fluorescence channels: ", paste(fluor_channels, collapse = ", "))

# -----------------------------------------------------------------------------
# 5. BUILD COMBINED DATA FRAME
# -----------------------------------------------------------------------------

all_data <- do.call(rbind, lapply(seq_along(fs), function(i) {
  df        <- as.data.frame(exprs(fs[[i]]))
  df$sample <- sampleNames(fs)[i]
  df
}))

# -----------------------------------------------------------------------------
# 6. SCATTER PLOTS  (FSC-A vs SSC-A, one per sample)
# -----------------------------------------------------------------------------

message("\n--- Scatter plots (FSC-A vs SSC-A) ---")

fsc_col <- grep("^FSC-A", colnames(all_data), value = TRUE)[1]
ssc_col <- grep("^SSC-A", colnames(all_data), value = TRUE)[1]

if (!is.na(fsc_col) && !is.na(ssc_col)) {

  for (sname in unique(all_data$sample)) {
    df_s <- all_data[all_data$sample == sname, ]

    p <- ggplot(df_s, aes(x = .data[[fsc_col]], y = .data[[ssc_col]])) +
      geom_hex(bins = 80) +
      scale_fill_viridis_c(trans = "log1p", name = "Count") +
      labs(title = paste("Scatter —", sname), x = fsc_col, y = ssc_col) +
      theme_bw(base_size = 13)

    out_file <- file.path(OUTPUT_DIR,
                          paste0("scatter_", gsub("[^A-Za-z0-9_-]", "_", sname), ".png"))
    ggsave(out_file, p, width = 6, height = 5, dpi = 150)
    message("  Saved: ", out_file)
  }

} else {
  message("  FSC-A / SSC-A not found — skipping scatter plots.")
}

# -----------------------------------------------------------------------------
# 7. PER-CHANNEL HISTOGRAMS  (all fluorescence channels, all samples overlaid)
# -----------------------------------------------------------------------------

message("\n--- Channel histograms ---")

for (ch in fluor_channels) {
  if (!ch %in% colnames(all_data)) next

  p <- ggplot(all_data, aes(x = .data[[ch]],
                             colour = sample, fill = sample)) +
    geom_density(alpha = 0.25, linewidth = 0.7) +
    scale_x_continuous(
      trans  = pseudo_log_trans(sigma = 1),
      labels = label_number(scale_cut = cut_short_scale())
    ) +
    labs(title = paste("Channel:", ch), x = ch, y = "Density") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          legend.title    = element_blank())

  out_file <- file.path(OUTPUT_DIR,
                        paste0("hist_", gsub("[^A-Za-z0-9_-]", "_", ch), ".png"))
  ggsave(out_file, p, width = 7, height = 4, dpi = 150)
  message("  Saved: ", out_file)
}

# -----------------------------------------------------------------------------
# 8. MFI  (Median Fluorescence Intensity per channel per sample)
# -----------------------------------------------------------------------------

message("\n--- MFI ---")

mfi_table <- all_data %>%
  group_by(sample) %>%
  summarise(across(all_of(fluor_channels),
                   ~ median(.x, na.rm = TRUE)),
            .groups = "drop")

print(as.data.frame(mfi_table))

write.csv(mfi_table,
          file      = file.path(OUTPUT_DIR, "MFI_per_channel.csv"),
          row.names = FALSE)
message("  Saved: MFI_per_channel.csv")

# -----------------------------------------------------------------------------
# 9. PERCENTAGE POSITIVE CELLS per channel per sample
# -----------------------------------------------------------------------------

message("\n--- % positive cells ---")

# Build threshold vector
if (THRESHOLD_METHOD == "percentile") {
  ref_df     <- all_data[all_data$sample == unique(all_data$sample)[1], ]
  thresholds <- sapply(fluor_channels,
                       function(ch) quantile(ref_df[[ch]],
                                             probs   = POSITIVE_PERCENTILE,
                                             na.rm   = TRUE))
  message("  Using ", POSITIVE_PERCENTILE * 100,
          "th-percentile thresholds from: ", unique(all_data$sample)[1])
} else {
  thresholds <- setNames(rep(POSITIVE_THRESHOLD_VALUE, length(fluor_channels)),
                         fluor_channels)
  message("  Using fixed threshold: ", POSITIVE_THRESHOLD_VALUE)
}

pct_positive <- all_data %>%
  group_by(sample) %>%
  summarise(
    across(
      all_of(fluor_channels),
      ~ round(mean(.x > thresholds[cur_column()], na.rm = TRUE) * 100, 2)
    ),
    .groups = "drop"
  )

print(as.data.frame(pct_positive))

write.csv(pct_positive,
          file      = file.path(OUTPUT_DIR, "percent_positive.csv"),
          row.names = FALSE)
message("  Saved: percent_positive.csv")

# -----------------------------------------------------------------------------
# 10. MFI HEATMAP  (summary overview)
# -----------------------------------------------------------------------------

message("\n--- MFI heatmap ---")

mfi_long <- mfi_table %>%
  pivot_longer(-sample, names_to = "channel", values_to = "MFI")

p_heat <- ggplot(mfi_long,
                 aes(x = channel, y = sample, fill = log1p(MFI))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(MFI, 0)), size = 3, colour = "white") +
  scale_fill_viridis_c(name = "log1p(MFI)") +
  labs(title = "MFI per channel per sample", x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

heatmap_w <- max(6, length(fluor_channels) * 1.1)
heatmap_h <- max(4, nrow(mfi_table) * 0.8 + 2)

ggsave(file.path(OUTPUT_DIR, "MFI_heatmap.png"), p_heat,
       width = heatmap_w, height = heatmap_h, dpi = 150)
message("  Saved: MFI_heatmap.png")

# -----------------------------------------------------------------------------
message("\nDone. All outputs in: ", OUTPUT_DIR)
