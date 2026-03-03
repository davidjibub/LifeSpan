# ============================================================
# README
# Function:
#   Batch GAMLSS (GGalt) analysis for all features in
#   master_report_merged.xlsx, excluding base info columns.
#   Results are saved into 5 group folders, with 10 features
#   per figure (10 rows x 2 cols: Trajectories + Variance).
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
#     <FolderName>.pdf  (all pages)
#     <FolderName>_fit_status.csv
# ============================================================

# =========================
# 0) packages
# =========================
pkgs <- c("gamlss", "ggplot2", "readxl", "patchwork")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(gamlss)
library(ggplot2)
library(readxl)
library(patchwork)

# =========================
# 1) source BrainChart repo scripts
# =========================
source("100.common-variables.r")
source("101.common-functions.r")
source("300.variables.r")
source("301.functions.r")

# =========================
# 2) load data
# =========================
xlsx_path <- "D:/Tang_DATA/BrainAge/Results/master_report_merged.xlsx"
raw_xlsx <- read_excel(xlsx_path)

# Base mapping
Study_raw <- as.character(raw_xlsx[["Site"]])
Study_raw[Study_raw == "SHOSU"] <- "S"
Study_raw[grepl("^brainage", Study_raw, ignore.case = TRUE)] <- "E"

Grp_raw <- as.character(raw_xlsx[["Sex"]])
Grp_raw[Grp_raw %in% c("\u7537", "Male", "M", "m")] <- "M"
Grp_raw[Grp_raw %in% c("\u5973", "Female", "F", "f")] <- "F"

TimeTransformed_raw <- raw_xlsx[["AgeDays"]] / 365.25 / 10
INDEX_ID_raw <- paste0(Study_raw, "_", as.character(raw_xlsx[["Subject_ID"]]))

BASE_DF <- data.frame(
  Study = Study_raw,
  Grp = Grp_raw,
  TimeTransformed = TimeTransformed_raw,
  INDEX.ID = INDEX_ID_raw,
  INDEX.OB = 1,
  INDEX.TYPE = "CN",
  stringsAsFactors = FALSE
)

# =========================
# 3) feature groups
# =========================
BASE_COLS <- c("Subject_ID", "Site", "Sex", "Seq", "AgeDays")

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

# Exclude base and previously assigned columns from cortical sets
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
# 4) helpers
# =========================
prep_feature_data <- function(feature_name) {
  df <- data.frame(
    Wand = raw_xlsx[[feature_name]],
    Study = BASE_DF$Study,
    Grp = BASE_DF$Grp,
    TimeTransformed = BASE_DF$TimeTransformed,
    INDEX.ID = BASE_DF$INDEX.ID,
    INDEX.OB = BASE_DF$INDEX.OB,
    INDEX.TYPE = BASE_DF$INDEX.TYPE,
    stringsAsFactors = FALSE
  )

  key_cols <- c("Wand", "Study", "Grp", "TimeTransformed", "INDEX.ID")
  df <- df[complete.cases(df[, key_cols]), ]

  df$Grp <- droplevels(as.factor(df$Grp))
  df$Study <- droplevels(as.factor(df$Study))
  df$INDEX.TYPE <- droplevels(as.factor(df$INDEX.TYPE))

  if (any(is.na(df$Wand))) stop("Wand has NA after filtering.")
  min_y <- min(df$Wand)
  if (min_y <= 0) {
    eps <- abs(min_y) + 1e-6
    df$Wand <- df$Wand + eps
  }

  attr(df, "columns") <- list(Index = c("INDEX.ID", "INDEX.OB", "INDEX.TYPE"))
  attr(df, "tag") <- "SelfData"
  attr(df, "Transformations") <- list(Type = "identity")
  df
}

fit_one <- function(raw_data, mu_formula, sigma_formula, nu_formula, fam_string) {
  HOLDER <- list()
  HOLDER$SUBSET <- raw_data
  HOLDER$MODEL <- list(
    mu = mu_formula,
    sigma = sigma_formula,
    nu = nu_formula,
    family = fam_string,
    covariates = list(
      Y = "Wand", X = "TimeTransformed", BY = "Grp",
      RANEF = "Study", COND = "INDEX.TYPE", ID = "INDEX.ID",
      OTHER = c("INDEX.OB")
    ),
    stratify = NULL,
    inc.fp = FALSE
  )

  fit_debug <- tryCatch(
    {
      gamlssWrapper(
        Model = HOLDER$MODEL,
        Data = HOLDER$SUBSET[, unlist(HOLDER$MODEL$covariates[c("Y","X","BY","OTHER","RANEF")])],
        Ctrl = GAMLSS.CTRL,
        trace = FALSE
      )
    },
    error = function(e) {
      message("gamlssWrapper error: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(fit_debug)) return(NULL)

  EX <- tryCatch(
    {
      Extract.Wrapper(HOLDER, Store.Full = TRUE)
    },
    error = function(e) {
      message("Extract.Wrapper error: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(EX) || is.null(EX$param)) return(NULL)

  if (anyNA(EX$param$mu$fixef)) {
    EX$param$mu$fixef[is.na(EX$param$mu$fixef)] <- 0
  }

  HOLDER$FIT.EXTRACT <- EX
  list(EXTRACT = EX, HOLDER = HOLDER)
}

run_gamlss_for_feature <- function(raw_data, feature_name, fast_preview = FALSE) {
  USE_FP_SEARCH <- TRUE
  FP_NPOLY <- 2
  pw <- c(1, 3)
  if (fast_preview) {
    USE_FP_SEARCH <- FALSE
    pw <- c(1, 3)
  }

  if (USE_FP_SEARCH) {
    fit_fp <- try(
      gamlss(
        Wand ~ fp(TimeTransformed, npoly = FP_NPOLY) + Grp,
        sigma.formula = ~ TimeTransformed,
        nu.formula = ~ 1,
        family = GGalt,
        data = raw_data,
        trace = FALSE
      ),
      silent = TRUE
    )
    if (!inherits(fit_fp, "try-error")) {
      pw_try <- try(getSmo(fit_fp, what = "mu")$power, silent = TRUE)
      if (!inherits(pw_try, "try-error") && !is.null(pw_try)) {
        pw <- pw_try
      }
    }
  }

  neg_tt <- which(raw_data$TimeTransformed < 0)
  if (length(neg_tt) > 0) {
    raw_data <- raw_data[raw_data$TimeTransformed >= 0, ]
    raw_data <- droplevels(raw_data)
  }

  if (any(raw_data$TimeTransformed <= 0) && any(abs(pw - round(pw)) > 1e-8)) {
    pw <- c(1, 3)
  }

  if (fast_preview) {
    models <- list(
      list(
        name = "M2",
        mu = as.formula(sprintf("Wand ~ bfpNA(TimeTransformed, powers=c(%s)) + Grp",
                                paste(pw, collapse=","))),
        sigma = as.formula(~ TimeTransformed),
        nu = as.formula(~ 1),
        fam = "GGalt"
      )
    )
  } else {
    models <- list(
      list(
        name = "M1",
        mu = as.formula(sprintf("Wand ~ bfpNA(TimeTransformed, powers=c(%s)) + Grp + random(Study)",
                                paste(pw, collapse=","))),
        sigma = as.formula(~ TimeTransformed),
        nu = as.formula(~ 1),
        fam = "GGalt"
      ),
      list(
        name = "M2",
        mu = as.formula(sprintf("Wand ~ bfpNA(TimeTransformed, powers=c(%s)) + Grp",
                                paste(pw, collapse=","))),
        sigma = as.formula(~ TimeTransformed),
        nu = as.formula(~ 1),
        fam = "GGalt"
      ),
      list(
        name = "M3",
        mu = as.formula(sprintf("Wand ~ bfpNA(TimeTransformed, powers=c(%s)) + Grp + Study",
                                paste(pw, collapse=","))),
        sigma = as.formula(~ TimeTransformed),
        nu = as.formula(~ 1),
        fam = "GGalt"
      )
    )
  }

  RES <- NULL
  for (m in models) {
    RES <- fit_one(raw_data, m$mu, m$sigma, m$nu, m$fam)
    if (!is.null(RES)) break
  }
  if (is.null(RES)) stop("All models failed for feature: ", feature_name)

  EXTRACT <- RES$EXTRACT
  HOLDER <- RES$HOLDER

  time_seq <- seq(min(raw_data$TimeTransformed), max(raw_data$TimeTransformed), length.out = 200)
  NEWData <- expand.grid(
    TimeTransformed = time_seq,
    Grp = levels(raw_data$Grp)
  )

  CURVE <- Apply.Param(
    NEWData = NEWData,
    FITParam = EXTRACT$param,
    Pred.Set = c("l025" = 0.025, "m500" = 0.5, "u975" = 0.975),
    Add.Moments = TRUE,
    Add.Derivative = FALSE,
    MissingToZero = TRUE,
    NAToZero = TRUE,
    verbose = FALSE
  )

  need <- c("PRED.l025.pop","PRED.m500.pop","PRED.u975.pop","PRED.variance.pop")
  miss <- setdiff(need, names(CURVE))
  if (length(miss) > 0) stop("Apply.Param missing columns for feature: ", feature_name)

  B <- if (fast_preview) 0 else 10
  if (B > 0) {
    Base.Seed <- 12345
    var_mat <- matrix(NA_real_, nrow = nrow(CURVE), ncol = B)

    for (b in 1:B) {
      EX_b <- try(Boot.Function(n = b, Base.Seed = Base.Seed, Holder = HOLDER, Create.INIT = TRUE), silent = TRUE)
      if (inherits(EX_b, "try-error") || is.null(EX_b$param)) next
      if (anyNA(EX_b$param$mu$fixef)) EX_b$param$mu$fixef[is.na(EX_b$param$mu$fixef)] <- 0

      CURVE_b <- try(
        Apply.Param(
          NEWData = NEWData,
          FITParam = EX_b$param,
          Pred.Set = NULL,
          Add.Moments = TRUE,
          Add.Derivative = FALSE,
          MissingToZero = TRUE,
          NAToZero = TRUE,
          verbose = FALSE
        ),
        silent = TRUE
      )
      if (inherits(CURVE_b, "try-error")) next
      if (!("PRED.variance.pop" %in% names(CURVE_b))) next
      var_mat[, b] <- CURVE_b$PRED.variance.pop
    }

    CURVE$VAR.ci.l025 <- apply(var_mat, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
    CURVE$VAR.ci.u975 <- apply(var_mat, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
  } else {
    CURVE$VAR.ci.l025 <- NA_real_
    CURVE$VAR.ci.u975 <- NA_real_
  }

  p_traj <- ggplot(CURVE, aes(x = TimeTransformed, group = Grp, color = Grp, fill = Grp)) +
    geom_ribbon(aes(ymin = PRED.l025.pop, ymax = PRED.u975.pop), alpha = 0.2, colour = NA) +
    geom_line(aes(y = PRED.m500.pop), linewidth = 1) +
    labs(
      title = paste0(feature_name, " - Trajectories"),
      x = "TimeTransformed",
      y = feature_name
    ) +
    theme_minimal()

  p_var <- ggplot(CURVE, aes(x = TimeTransformed, group = Grp, color = Grp, fill = Grp)) +
    {if (B > 0) geom_ribbon(aes(ymin = VAR.ci.l025, ymax = VAR.ci.u975), alpha = 0.2, colour = NA)} +
    geom_line(aes(y = PRED.variance.pop), linewidth = 1) +
    labs(
      title = paste0(feature_name, " - Variance"),
      x = "TimeTransformed",
      y = "Variance"
    ) +
    theme_minimal()

  list(traj = p_traj, var = p_var)
}

error_plot <- function(feature_name, msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = paste0(feature_name, "\\n", msg), hjust = 0) +
    theme_void()
}

# =========================
# 5) batch processing
# =========================
# Preview mode: only generate the first page for one group
PREVIEW_ONLY <- FALSE #确认图大小满意后，把：PREVIEW_ONLY <- FALSE
FAST_PREVIEW <- PREVIEW_ONLY
PREVIEW_GROUP <- "Cortical_Thickness"  # change to any group name
PREVIEW_PAGE <- 1
for (group_name in names(GROUPS)) {
  features <- GROUPS[[group_name]]
  if (length(features) == 0) next

  dir.create(group_name, showWarnings = FALSE)

  chunks <- split(features, ceiling(seq_along(features) / 5))
  page_idx <- 1
  status_rows <- list()
  pdf_path <- file.path(group_name, sprintf("%s.pdf", group_name))
  grDevices::pdf(pdf_path, width = 14, height = 12)

  for (chunk in chunks) {
    plot_list <- list()
    for (feat in chunk) {
      message("Processing: ", group_name, " / ", feat)
      plots <- try({
        raw_data <- prep_feature_data(feat)
        run_gamlss_for_feature(raw_data, feat, fast_preview = FAST_PREVIEW)
      }, silent = TRUE)

      if (inherits(plots, "try-error")) {
        status_rows[[length(status_rows) + 1]] <- data.frame(
          Feature = feat,
          Status = "Failed",
          stringsAsFactors = FALSE
        )
        plot_list <- c(plot_list, list(
          error_plot(feat, "Model failed"),
          error_plot(feat, "Model failed")
        ))
      } else {
        status_rows[[length(status_rows) + 1]] <- data.frame(
          Feature = feat,
          Status = "Success",
          stringsAsFactors = FALSE
        )
        plot_list <- c(plot_list, list(plots$traj, plots$var))
      }
    }

    n_rows <- length(chunk)
    fig <- wrap_plots(plot_list, ncol = 2)
    height_use <- max(8, n_rows * 3.0)
    fig <- fig + plot_layout(heights = rep(1, n_rows))
    print(fig)
    page_idx <- page_idx + 1

    if (PREVIEW_ONLY && group_name == PREVIEW_GROUP && page_idx > PREVIEW_PAGE) break
  }

  grDevices::dev.off()

  if (length(status_rows) > 0) {
    status_df <- do.call(rbind, status_rows)
    write.csv(status_df, file.path(group_name, sprintf("%s_fit_status.csv", group_name)), row.names = FALSE)
  }
  if (PREVIEW_ONLY && group_name == PREVIEW_GROUP) break
}

message("=== Done: batch GAMLSS results saved by group folders ===")
