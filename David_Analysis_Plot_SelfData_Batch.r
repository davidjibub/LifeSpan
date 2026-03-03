# ============================================================
# README
# Function:
#   Batch scatter plots for all features in master_report_merged.xlsx.
#   Layout follows GAMLSS batch: 5 features per page (5x1).
#   Results are saved into the same 5 group folders created by
#   David_Analysis_fp_robust_SelfData_Batch.r.
#
# Input:
#   D:/Tang_DATA/BrainAge/Results/master_report_merged.xlsx
#   Required base columns (English headers):
#     Subject_ID, Site, Sex, Seq, AgeDays
#
# Output folders (under current working directory):
#   Global
#   Subcortical_All
#   Cortical_Thickness
#   Cortical_Volume
#   Cortical_Area
#   Each folder contains:
#     <FolderName>_Scatter.pdf
# ============================================================

# =========================
# 0) packages
# =========================
pkgs <- c("ggplot2", "readxl", "patchwork")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2)
library(readxl)
library(patchwork)

# =========================
# 1) load data
# =========================
xlsx_path <- "D:/Tang_DATA/BrainAge/Results/master_report_merged.xlsx"
raw_xlsx <- read_excel(xlsx_path)

BASE_COLS <- c("Subject_ID", "Site", "Sex", "Seq", "AgeDays")

# Sex mapping to Grp (M/F)
Grp_raw <- as.character(raw_xlsx[["Sex"]])
Grp_raw[Grp_raw %in% c("\u7537", "Male", "M", "m")] <- "M"
Grp_raw[Grp_raw %in% c("\u5973", "Female", "F", "f")] <- "F"

BASE_DF <- data.frame(
  AgeDays = raw_xlsx[["AgeDays"]],
  Grp = Grp_raw,
  stringsAsFactors = FALSE
)

# =========================
# 2) feature groups (same as GAMLSS batch)
# =========================
GLOBAL_COLS <- c(
  "TIV", "Brain_Parenchyma", "BrainSegVolNotVent", "CortexVol", "SubCortGrayVol",
  "Total_GM", "Total_WM", "CSF",
  "lhCortexVol", "rhCortexVol",
  "lhCerebralWhiteMatterVol", "rhCerebralWhiteMatterVol",
  "SupraTentorialVolNotVent"
)

SUBCORTICAL_ALL_COLS <- c(
  "Left-Hippocampus", "Right-Hippocampus",
  "Left-Amygdala", "Right-Amygdala",
  "Left-Thalamus", "Right-Thalamus",
  "Left-Caudate", "Right-Caudate",
  "Left-Putamen", "Right-Putamen",
  "Left-Pallidum", "Right-Pallidum",
  "Left-Accumbens-area", "Right-Accumbens-area",
  "Left-VentralDC", "Right-VentralDC",
  "Brain-Stem",
  "Left-Lateral-Ventricle", "Right-Lateral-Ventricle",
  "Left-Inf-Lat-Vent", "Right-Inf-Lat-Vent",
  "3rd-Ventricle", "4th-Ventricle", "5th-Ventricle",
  "Left-choroid-plexus", "Right-choroid-plexus",
  "Left-Cerebellum-Cortex", "Right-Cerebellum-Cortex",
  "Left-Cerebellum-White-Matter", "Right-Cerebellum-White-Matter",
  "CC_Anterior", "CC_Mid_Anterior", "CC_Central", "CC_Mid_Posterior", "CC_Posterior",
  "Optic-Chiasm"
)

all_cols <- colnames(raw_xlsx)
GLOBAL_COLS <- intersect(GLOBAL_COLS, all_cols)
SUBCORTICAL_ALL_COLS <- intersect(SUBCORTICAL_ALL_COLS, all_cols)

CORTICAL_THICK_COLS <- all_cols[grepl("_thick$", all_cols)]
CORTICAL_VOL_COLS <- all_cols[grepl("_vol$", all_cols)]
CORTICAL_AREA_COLS <- all_cols[grepl("_area$", all_cols)]

EXCLUDE_COLS <- unique(c(BASE_COLS, GLOBAL_COLS, SUBCORTICAL_ALL_COLS))
CORTICAL_THICK_COLS <- setdiff(CORTICAL_THICK_COLS, EXCLUDE_COLS)
CORTICAL_VOL_COLS <- setdiff(CORTICAL_VOL_COLS, EXCLUDE_COLS)
CORTICAL_AREA_COLS <- setdiff(CORTICAL_AREA_COLS, EXCLUDE_COLS)

GROUPS <- list(
  Global = GLOBAL_COLS,
  Subcortical_All = SUBCORTICAL_ALL_COLS,
  Cortical_Thickness = CORTICAL_THICK_COLS,
  Cortical_Volume = CORTICAL_VOL_COLS,
  Cortical_Area = CORTICAL_AREA_COLS
)

# =========================
# 3) helpers
# =========================
make_scatter <- function(feature_name) {
  df <- data.frame(
    Value = raw_xlsx[[feature_name]],
    AgeDays = BASE_DF$AgeDays,
    Grp = BASE_DF$Grp,
    stringsAsFactors = FALSE
  )
  df <- df[complete.cases(df), ]
  df$Grp <- droplevels(as.factor(df$Grp))

  ggplot(df, aes(x = AgeDays, y = Value, color = Grp, shape = Grp)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_shape_manual(values = c(M = 17, F = 16)) +
    labs(
      title = paste0(feature_name, " - Scatter"),
      x = "Age (days)",
      y = feature_name
    ) +
    theme_minimal()
}

error_plot <- function(feature_name, msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = paste0(feature_name, "\n", msg), hjust = 0) +
    theme_void()
}

# =========================
# 4) batch plotting
# =========================
for (group_name in names(GROUPS)) {
  features <- GROUPS[[group_name]]
  if (length(features) == 0) next

  dir.create(group_name, showWarnings = FALSE)

  chunks <- split(features, ceiling(seq_along(features) / 5))
  pdf_path <- file.path(group_name, sprintf("%s_Scatter.pdf", group_name))
  grDevices::pdf(pdf_path, width = 10, height = 14)

  for (chunk in chunks) {
    plot_list <- list()
    for (feat in chunk) {
      message("Scatter: ", group_name, " / ", feat)
      p <- try(make_scatter(feat), silent = TRUE)
      if (inherits(p, "try-error")) {
        plot_list <- c(plot_list, list(error_plot(feat, "Plot failed")))
      } else {
        plot_list <- c(plot_list, list(p))
      }
    }

    fig <- wrap_plots(plot_list, ncol = 1)
    print(fig)
  }

  grDevices::dev.off()

}

message("=== Done: batch scatter plots saved by group folders ===")
