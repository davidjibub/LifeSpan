# ============================================================
# David_Analysis_Plot_SelfData.r
# Goal:
#   Scatter plot of TIV from master_report_merged.xlsx
#   Male: triangle, Female: circle
#   Colors consistent with David_Analysis_fp_robust_SelfData.r (Grp-based)
# ============================================================

# =========================
# 0) packages
# =========================
pkgs <- c("ggplot2", "readxl")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2)
library(readxl)

# =========================
# 1) load + map self data (same mapping as robust script)
# =========================
xlsx_path <- "D:/Tang_DATA/BrainAge/Results/master_report_merged.xlsx"
raw_xlsx <- read_excel(xlsx_path)

# Expect English headers: Site, Sex, AgeDays, TIV
raw_data <- data.frame(
  TIV = raw_xlsx[["TIV"]],
  Sex = raw_xlsx[["Sex"]],
  AgeDays = raw_xlsx[["AgeDays"]],
  stringsAsFactors = FALSE
)

# Map Sex -> Grp (M/F)
raw_data$Grp <- as.character(raw_data$Sex)
raw_data$Grp[raw_data$Grp %in% c("\u7537", "Male", "M", "m")] <- "M"
raw_data$Grp[raw_data$Grp %in% c("\u5973", "Female", "F", "f")] <- "F"

# Drop missing
raw_data <- raw_data[complete.cases(raw_data[, c("TIV", "AgeDays", "Grp")]), ]
raw_data$Grp <- droplevels(as.factor(raw_data$Grp))

# =========================
# 2) scatter plot
# =========================
# Scale factor to compress y-axis length
# Adjust this value as needed (e.g., 1000, 1000000)
TIV_SCALE <- 6000
raw_data$TIV_scaled <- raw_data$TIV / TIV_SCALE

p <- ggplot(raw_data, aes(x = AgeDays, y = TIV_scaled, color = Grp, shape = Grp)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_shape_manual(values = c(M = 17, F = 16)) +
  labs(
    title = "TIV Scatter (Self Data)",
    x = "Age (days)",
    y = sprintf("TIV / %g", TIV_SCALE)
  ) +
  theme_minimal()

print(p)

# =========================
# 3) save output
# =========================
ggsave("TIV_Scatter_SelfData.png", p, width = 7, height = 5, dpi = 300)

message("=== Done: TIV scatter saved as TIV_Scatter_SelfData.png ===")
