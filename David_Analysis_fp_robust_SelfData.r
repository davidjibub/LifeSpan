# ============================================================
# README
# Function:
#   Run GAMLSS (GGalt) on self data from master_report_merged.xlsx
#   and generate normative trajectory + variance plots.
#
# Input:
#   D:/Tang_DATA/BrainAge/Results/master_report_merged.xlsx
#   Required columns (English headers):
#     Subject_ID, Site, Sex, AgeDays, Brain_Parenchyma
#
# Output (saved to current working directory):
#   If patchwork installed:
#     BrainChart_Figure1_like_SelfData.png
#   Else:
#     BrainChart_Trajectories_SelfData.png
#     BrainChart_Variance_SelfData.png
#
# Notes:
#   This script only reads the xlsx above.
#   It does NOT read any simulated data files.
# ============================================================

# =========================
# 0) packages
# =========================
pkgs <- c("gamlss", "ggplot2", "readxl")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(gamlss)
library(ggplot2)
library(readxl)

# has_patchwork 是什么？
# has_patchwork 用来检测你是否安装了 patchwork 包。

# 如果有：会把轨迹图和方差图拼成一张图保存 BrainChart_Figure1_like_SelfData.png。
# 如果没有：会分别保存两张图 BrainChart_Trajectories_SelfData.png 和 BrainChart_Variance_SelfData.png。
# 它不影响分析本身，只影响图怎么组合保存。

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

# =========================
# 1) source BrainChart repo scripts (same as robust)
# =========================
source("100.common-variables.r")
source("101.common-functions.r")
source("300.variables.r")
source("301.functions.r")

# =========================
# 2) load + map self data
# =========================
# master_report_merged.xlsx -> columns mapping (English headers)
# Wand              <- Brain_Parenchyma
# Study             <- Site   (SHOSU -> S; brainage* -> E)
# Grp               <- Sex    (Male -> M; Female -> F)
# TimeTransformed   <- AgeDays / 365.25 / 10  (align with simulated Time/10)
# INDEX.ID          <- Subject_ID (prefixed by Study)
# INDEX.OB          <- 1
# INDEX.TYPE        <- CN (constant, if no clinical grouping)

xlsx_path <- "D:/Tang_DATA/BrainAge/Results/master_report_merged.xlsx"
raw_xlsx <- read_excel(xlsx_path)

raw_data <- data.frame(
  Wand = raw_xlsx[["Brain_Parenchyma"]],
  Study = raw_xlsx[["Site"]],
  Grp = raw_xlsx[["Sex"]],
  TimeTransformed = raw_xlsx[["AgeDays"]] / 365.25 / 10,
  INDEX.ID = raw_xlsx[["Subject_ID"]],
  INDEX.OB = 1,
  INDEX.TYPE = "CN",
  stringsAsFactors = FALSE
)

# Study mapping
raw_data$Study <- as.character(raw_data$Study)
raw_data$Study[raw_data$Study == "SHOSU"] <- "S"
raw_data$Study[grepl("^brainage", raw_data$Study, ignore.case = TRUE)] <- "E"

# Grp mapping
raw_data$Grp <- as.character(raw_data$Grp)
raw_data$Grp[raw_data$Grp %in% c("\u7537", "Male", "M", "m")] <- "M"
raw_data$Grp[raw_data$Grp %in% c("\u5973", "Female", "F", "f")] <- "F"

# INDEX.ID: prefix with Study to resemble simulated IDs
raw_data$INDEX.ID <- paste0(raw_data$Study, "_", as.character(raw_data$INDEX.ID))

# drop rows with missing key fields
key_cols <- c("Wand", "Study", "Grp", "TimeTransformed", "INDEX.ID")
raw_data <- raw_data[complete.cases(raw_data[, key_cols]), ]

# factor levels (consistent with robust)
raw_data$Grp <- droplevels(as.factor(raw_data$Grp))
raw_data$Study <- droplevels(as.factor(raw_data$Study))
raw_data$INDEX.TYPE <- droplevels(as.factor(raw_data$INDEX.TYPE))

# GG/GGalt requires Wand > 0
if (any(is.na(raw_data$Wand))) stop("Wand has NA; please handle missing values first.")
min_y <- min(raw_data$Wand)
if (min_y <= 0) {
  eps <- abs(min_y) + 1e-6
  raw_data$Wand <- raw_data$Wand + eps
  message(sprintf("Detected Wand<=0; shifted by +%.6g to satisfy GGalt y>0.", eps))
}

# BrainChart pipeline attributes
attr(raw_data, "columns") <- list(Index = c("INDEX.ID","INDEX.OB","INDEX.TYPE"))
attr(raw_data, "tag") <- "SelfData"
attr(raw_data, "Transformations") <- list(Type="identity")

# =========================
# 3) Step 1: fp search (same as robust)
# =========================
USE_FP_SEARCH <- TRUE
FP_NPOLY <- 2
pw <- c(1, 3)

if (USE_FP_SEARCH) {
  message("=== FP search: using fp() to select powers ===")
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
      message("fp selected powers = ", paste(pw, collapse = ", "))
    } else {
      message("Cannot read fp powers; fallback to powers = 1,3")
    }
  } else {
    message("fp search failed; fallback to powers = 1,3")
  }
}

# TimeTransformed should not be negative
neg_tt <- which(raw_data$TimeTransformed < 0)
if (length(neg_tt) > 0) {
  message(sprintf("Detected TimeTransformed<0 rows: %d, removed", length(neg_tt)))
  raw_data <- raw_data[raw_data$TimeTransformed >= 0, ]
  raw_data <- droplevels(raw_data)
}

# Avoid NaN for non-integer powers if TimeTransformed <= 0
if (any(raw_data$TimeTransformed <= 0) && any(abs(pw - round(pw)) > 1e-8)) {
  message("TimeTransformed<=0 with non-integer fp powers; fallback to powers = 1,3")
  pw <- c(1, 3)
}

# =========================
# 4) helper: fit once (same as robust)
# =========================
fit_one <- function(mu_formula, sigma_formula, nu_formula, fam_string) {
  HOLDER <- list()
  HOLDER$SUBSET <- raw_data
  HOLDER$MODEL <- list(
    mu = mu_formula,
    sigma = sigma_formula,
    nu = nu_formula,
    family = fam_string,
    covariates = list(
      Y="Wand", X="TimeTransformed", BY="Grp",
      RANEF="Study", COND="INDEX.TYPE", ID="INDEX.ID",
      OTHER=c("INDEX.OB")
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
  if (is.null(EX) || is.null(EX$param)) {
    if (!is.null(EX) && !is.null(EX$diagnostics)) {
      message("Extract.Wrapper diagnostics: ", EX$diagnostics)
    }
    return(NULL)
  }

  if (anyNA(EX$param$mu$fixef)) {
    na_names <- names(EX$param$mu$fixef)[is.na(EX$param$mu$fixef)]
    message("mu.fixef has NA; set to 0: ", paste(na_names, collapse = ", "))
    EX$param$mu$fixef[is.na(EX$param$mu$fixef)] <- 0
  }

  HOLDER$FIT.EXTRACT <- EX
  list(EXTRACT = EX, HOLDER = HOLDER)
}

# =========================
# 5) model candidates (same as robust)
# =========================
models <- list(
  list(
    name = "M1: bfpNA(fp powers) + Grp + random(Study)",
    mu   = as.formula(sprintf("Wand ~ bfpNA(TimeTransformed, powers=c(%s)) + Grp + random(Study)",
                              paste(pw, collapse=","))),
    sigma= as.formula(~ TimeTransformed),
    nu   = as.formula(~ 1),
    fam  = "GGalt"
  ),
  list(
    name = "M2: bfpNA(fp powers) + Grp (no random)",
    mu   = as.formula(sprintf("Wand ~ bfpNA(TimeTransformed, powers=c(%s)) + Grp",
                              paste(pw, collapse=","))),
    sigma= as.formula(~ TimeTransformed),
    nu   = as.formula(~ 1),
    fam  = "GGalt"
  ),
  list(
    name = "M3: bfpNA(fp powers) + Grp + Study (fixed Study)",
    mu   = as.formula(sprintf("Wand ~ bfpNA(TimeTransformed, powers=c(%s)) + Grp + Study",
                              paste(pw, collapse=","))),
    sigma= as.formula(~ TimeTransformed),
    nu   = as.formula(~ 1),
    fam  = "GGalt"
  )
)

RES <- NULL
for (m in models) {
  message("Trying model: ", m$name)
  RES <- fit_one(m$mu, m$sigma, m$nu, m$fam)
  if (!is.null(RES)) { message("Model fitted: ", m$name); break }
}
if (is.null(RES)) stop("All models failed. Check column mapping or data issues.")

EXTRACT <- RES$EXTRACT
HOLDER  <- RES$HOLDER

# =========================
# 6) NEWData grid
# =========================
time_seq <- seq(min(raw_data$TimeTransformed), max(raw_data$TimeTransformed), length.out = 200)
NEWData <- expand.grid(
  TimeTransformed = time_seq,
  Grp = levels(raw_data$Grp)
)

# =========================
# 7) Normative trajectories + variance
# =========================
CURVE <- Apply.Param(
  NEWData = NEWData,
  FITParam = EXTRACT$param,
  Pred.Set = c("l025"=0.025, "m500"=0.5, "u975"=0.975),
  Add.Moments = TRUE,
  Add.Derivative = FALSE,
  MissingToZero = TRUE,
  NAToZero = TRUE,
  verbose = FALSE
)

need <- c("PRED.l025.pop","PRED.m500.pop","PRED.u975.pop","PRED.variance.pop")
miss <- setdiff(need, names(CURVE))
if (length(miss) > 0) stop(paste("Apply.Param missing columns:", paste(miss, collapse=", ")))

# =========================
# 8) Normative variance 95% CI (bootstrap)
# =========================
#B <- 200
B <- 10
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

# =========================
# 9) plots
# =========================
p_traj <- ggplot(CURVE, aes(x = TimeTransformed, group = Grp, color = Grp, fill = Grp)) +
  geom_ribbon(aes(ymin = PRED.l025.pop, ymax = PRED.u975.pop), alpha = 0.2, colour = NA) +
  geom_line(aes(y = PRED.m500.pop), linewidth = 1) +
  labs(
    title = "Normative trajectories (PI: 2.5-97.5 centiles)",
    x = "TimeTransformed",
    y = "Wand (Brain_Parenchyma)"
  ) +
  theme_minimal()

p_var <- ggplot(CURVE, aes(x = TimeTransformed, group = Grp, color = Grp, fill = Grp)) +
  geom_ribbon(aes(ymin = VAR.ci.l025, ymax = VAR.ci.u975), alpha = 0.2, colour = NA) +
  geom_line(aes(y = PRED.variance.pop), linewidth = 1) +
  labs(
    title = "Normative variance (CI: bootstrap 95%)",
    x = "TimeTransformed",
    y = "Variance"
  ) +
  theme_minimal()

# =========================
# 10) save outputs
# =========================
if (has_patchwork) {
  fig <- p_traj / p_var
  print(fig)
  ggsave("BrainChart_Figure1_like_SelfData.png", fig, width = 9, height = 7, dpi = 300)
} else {
  print(p_traj); print(p_var)
  ggsave("BrainChart_Trajectories_SelfData.png", p_traj, width = 9, height = 4, dpi = 300)
  ggsave("BrainChart_Variance_SelfData.png", p_var, width = 9, height = 4, dpi = 300)
}

message("=== Done: figures saved in current working directory ===")
