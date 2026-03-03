############################################################
# DK ROI values on FreeSurfer surfaces
#  - Use your own FreeSurfer subject (no fsbrain demo download)
#  - atlas: aparc (Desikan-Killiany / DK)
#  - surface: pial & inflated
#  - output: two PNGs (t4 view mosaics)
############################################################

# 0) 安装包（如已安装可注释掉）
# install.packages("fsbrain")
# install.packages("freesurferformats")

# 1) 导入包
library(fsbrain)
library(freesurferformats)

# 2) Set subjects_dir and subject_id (your FreeSurfer results)
#    subject path: D:/Tang_DATA/BrainAge/Results/brainage0130/brainage0130_1H/sub001
#    SUBJECTS_DIR should be the parent folder containing the subject
# 3) 设置 subjects_dir 与 subject_id（这里用 demo 数据的 subject1）
subjects_dir <- normalizePath("D:/Tang_DATA/BrainAge/Results/brainage0130/brainage0130_1H")
subject_id   <- "sub001"

# 4) 读取 DK(aparc) annotation，拿到 ROI 名称列表
#    注：aparc 对应 DK（Desikan-Killiany）皮层分区
annot_lh <- freesurferformats::read.fs.annot(
  file.path(subjects_dir, subject_id, "label", "lh.aparc.annot")
)
annot_rh <- freesurferformats::read.fs.annot(
  file.path(subjects_dir, subject_id, "label", "rh.aparc.annot")
)

# colortable$struct_names 是 ROI 名称（通常包含 unknown / corpuscallosum 等）
roi_lh <- annot_lh$colortable$struct_names
roi_rh <- annot_rh$colortable$struct_names

# 可选：去掉不想显示/不属于皮层的名字（常见是 unknown、corpuscallosum）
drop_names <- c("unknown", "corpuscallosum")
roi_lh_use <- setdiff(roi_lh, drop_names)
roi_rh_use <- setdiff(roi_rh, drop_names)

# 5) 为每个 ROI 生成一个随机值（可重复：设置随机种子）
set.seed(1234)
vals_lh <- stats::rnorm(length(roi_lh_use))
vals_rh <- stats::rnorm(length(roi_rh_use))

# fsbrain 的 vis.region.values.on.subject() 需要“具名 list”：list(roiName=value, ...)
lh_region_value_list <- as.list(vals_lh); names(lh_region_value_list) <- roi_lh_use
rh_region_value_list <- as.list(vals_rh); names(rh_region_value_list) <- roi_rh_use


# 6) Detect usable pial surface (handle broken symlinks on Windows)
pial_lh <- file.path(subjects_dir, subject_id, "surf", "lh.pial")
pial_rh <- file.path(subjects_dir, subject_id, "surf", "rh.pial")
pial_ok <- all(file.exists(c(pial_lh, pial_rh))) && all(file.info(c(pial_lh, pial_rh))$size > 0)
surface_pial <- if (pial_ok) "pial" else "pial.T1"
# 6) （可选）在无显示器/服务器环境下渲染
#    如果你在本地 RStudio 一般不需要。服务器/无 GUI 时建议打开：
# rgl::rgl.useNULL(TRUE)

# 7) 渲染 + 导出：pial 表面

mk <- list(
  colFn = viridisLite::viridis, # 或 viridisLite::magma/inferno/plasma/cividis
  n = 256
)

out_pial <- file.path(getwd(), "DK_random_values_pial.png")
vis.region.values.on.subject(
  subjects_dir = subjects_dir,
  subject_id   = subject_id,
  atlas        = "aparc",
  lh_region_value_list = lh_region_value_list,
  rh_region_value_list = rh_region_value_list,
  # 关键：表面与视角
  surface = surface_pial,
  views   = c("t4"),
  # 相关设置
  makecmap_options = mk,
  draw_colorbar = TRUE,
  bg = "curv_light",
  border = TRUE,
  # 导出 png（同时也会弹出/创建 rgl 场景）
  rglactions = list(
    "snapshot_png" = out_pial
  )
)

# 8) 渲染 + 导出：inflated 表面
out_inflated <- file.path(getwd(), "DK_random_values_inflated.png")
vis.region.values.on.subject(
  subjects_dir = subjects_dir,
  subject_id   = subject_id,
  atlas        = "aparc",
  lh_region_value_list = lh_region_value_list,
  rh_region_value_list = rh_region_value_list,
  surface = "inflated",
  views   = c("t4"),
  draw_colorbar = TRUE,
  border = TRUE,
  rglactions = list(
    "snapshot_png" = out_inflated
  )
)

cat("Done!\nSaved:\n", out_pial, "\n", out_inflated, "\n")