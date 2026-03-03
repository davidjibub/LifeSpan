# ============================================================
# Plot FreeSurfer aparc ROI heatmap using ggsegExtra (RStudio)
# - Reads:   D:/Tang_DATA/BrainAge/ServerData/SHOSU/sub058/stats/lh.aparc.stats
#            D:/Tang_DATA/BrainAge/ServerData/SHOSU/sub058/stats/rh.aparc.stats
# - Creates: random ROI values (for preview)
# - Plots:   ROI-wise heatmap (Desikan-Killiany aparc)
# - Saves:   D:/Tang_DATA/BrainAge/ServerData/SHOSU/sub058/ggseg_aparc_random.png
# ============================================================

# ---------- 0) Install & load packages ----------
pkgs <- c("ggseg", "ggsegExtra", "dplyr", "ggplot2", "viridis")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggseg)
library(ggsegExtra)
library(dplyr)
library(ggplot2)
library(viridis)

# ---------- 1) Your FreeSurfer subject path ----------
subjects_dir <- "D:/Tang_DATA/BrainAge/ServerData/SHOSU"
subject_id   <- "sub058"
stats_dir    <- file.path(subjects_dir, subject_id, "stats")

lh_stats <- file.path(stats_dir, "lh.aparc.stats")
rh_stats <- file.path(stats_dir, "rh.aparc.stats")

if (!file.exists(lh_stats) || !file.exists(rh_stats)) {
  stop(
    "找不到 aparc.stats 文件，请确认路径存在：\n",
    lh_stats, "\n", rh_stats,
    "\n\n（recon-all 的完整输出应包含 stats/ 目录）"
  )
}

# ---------- 2) Read aparc.stats (robust parser) ----------
read_aparc_stats <- function(path, hemi) {
  lines <- readLines(path, warn = FALSE)
  
  # 找到 ColHeaders 行（示例：# ColHeaders StructName NumVert SurfArea ...）
  header_line <- grep("^#\\s*ColHeaders", lines, value = TRUE)
  if (length(header_line) == 0) stop("在文件中找不到 '# ColHeaders'：", path)
  
  colnames_vec <- strsplit(gsub("^#\\s*ColHeaders\\s+", "", header_line[1]), "\\s+")[[1]]
  
  # 用 comment.char="#" 跳过所有以 # 开头的行
  dat <- read.table(path, comment.char = "#", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(dat) != length(colnames_vec)) {
    stop("列数与 ColHeaders 不一致：", path,
         "\n读取到列数: ", ncol(dat),
         " | 头部列数: ", length(colnames_vec))
  }
  colnames(dat) <- colnames_vec
  
  dat$hemi <- hemi
  dat
}

lh_df <- read_aparc_stats(lh_stats, "left")
rh_df <- read_aparc_stats(rh_stats, "right")

aparc_df <- bind_rows(lh_df, rh_df)

# ---------- 3) Make ROI-wise values (random, for preview) ----------
# ggseg 需要列名：region, hemi, value（hemi 用 left/right）
plot_df <- aparc_df %>%
  transmute(
    region = tolower(StructName),
    hemi   = hemi
  ) %>%
  # 去掉常见无效分区（可按需要增减）
  filter(!region %in% c("unknown", "corpuscallosum")) %>%
  distinct(region, hemi) %>%
  mutate(
    value = { set.seed(42); rnorm(n()) }   # ✅ 随机值：先看效果
  )

# ---------- 4) Choose atlas and plot ----------
# FreeSurfer aparc 默认对应 Desikan-Killiany（DK）分区
atlas_to_use <- ggsegExtra::dk

p <- ggseg(
  .data   = plot_df,
  atlas   = atlas_to_use,
  mapping = aes(fill = value),
  colour  = "white",
  size    = 0.2
) +
  scale_fill_viridis_c(option = "C", na.value = "grey90") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = paste0(subject_id, " | FreeSurfer aparc (DK) | random ROI values"),
    fill  = "Random"
  )

print(p)

# ---------- 5) Save figure ----------
out_png <- file.path(subjects_dir, subject_id, "ggseg_aparc_random.png")
ggsave(out_png, plot = p, width = 10, height = 6, dpi = 300)
message("Saved: ", out_png)
