# ============================================================
# David_Analysis.R  (完整可运行版本)
# 目标：
#   1) 用 GAMLSS (GGalt) 做 normative model
#   2) Normative trajectories：PI/centile band (2.5%, 50%, 97.5%)
#   3) Normative variance：方差曲线 + bootstrap 95% CI（这是 CI，不是 PI）
#   4) 支持 “fp 搜索一次幂” -> 固化成 bfpNA(powers=...) 再拟合（更稳）
#   5) 自动兜底：模型不收敛/不可辨识时降级；出现 NA 系数时置 0 防止预测炸
# ============================================================

# =========================
# 0) packages
# =========================
pkgs <- c("gamlss", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(gamlss)
library(ggplot2)

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

# =========================
# 1) source BrainChart repo scripts（你上传的 R 脚本）
# =========================
source("100.common-variables.r")
source("101.common-functions.r")   # 会 source 102.gamlss-recode.r（含 GGalt/bfpNA 等）
source("300.variables.r")
source("301.functions.r")          # Extract.Wrapper / Apply.Param / Boot.Function 等在这里

# =========================
# 2) load + sanitize data
# =========================
# raw_data <- read.csv("Simulated_Brain_Data.csv")  # 原始nature生成的模拟数据
raw_data <- read.csv("Simulated_Brain_Data_chatgpt.csv")  # chatgpt修改后保留必要列的数据

# 关键列类型（建模/预测强依赖 factor levels 一致）
raw_data$Grp       <- droplevels(as.factor(raw_data$Grp))
raw_data$Study     <- droplevels(as.factor(raw_data$Study))
raw_data$INDEX.TYPE<- droplevels(as.factor(raw_data$INDEX.TYPE))

# GG/GGalt 属于正值分布：Wand 必须 > 0
if (any(is.na(raw_data$Wand))) stop("Wand 含 NA：请先处理缺失值。")
min_y <- min(raw_data$Wand)
if (min_y <= 0) {
  eps <- abs(min_y) + 1e-6
  raw_data$Wand <- raw_data$Wand + eps
  message(sprintf("检测到 Wand<=0：已整体平移 +%.6g 以满足 GGalt 的 y>0 要求。", eps))
}

# BrainChart pipeline 需要的 attributes（常用于 ValidateCleanInput / Check.Attributes 等）
attr(raw_data, "columns") <- list(Index = c("INDEX.ID","INDEX.OB","INDEX.TYPE"))
attr(raw_data, "tag") <- "SimulationData"
attr(raw_data, "Transformations") <- list(Type="identity")

# =========================
# 3) Step 1: fp 搜索一次幂（可选但你想要）
#    思路：先用 fp(...) 自动选择幂 -> 读取 power -> 再用 bfpNA 固化幂重拟合
# =========================
USE_FP_SEARCH <- TRUE      # 你想“跑一次 fp”就 TRUE；不想就 FALSE
FP_NPOLY <- 2              # 通常 2 更稳；3 更容易不收敛/过拟合

pw <- c(1, 3)              # 默认幂（如果 fp 搜索失败就回退）

if (USE_FP_SEARCH) {
  message("=== FP 搜索阶段：用 fp() 自动选幂 ===")
  
  # 为了提高 fp 搜索稳定性：先不加 random(Study)（更稳），只搜年龄形状
  # 你如果坚持把 random(Study) 也加进 fp 搜索，可以把公式里的 + random(Study) 放回去
  fit_fp <- try(
    gamlss(
      Wand ~ fp(TimeTransformed, npoly = FP_NPOLY) + Grp,
      sigma.formula = ~ TimeTransformed,
      nu.formula = ~ 1,
      family = GGalt,     # 102.gamlss-recode.r 里定义的 family
      data = raw_data,
      trace = FALSE
    ),
    silent = TRUE
  )
  
  if (!inherits(fit_fp, "try-error")) {
    # getSmo 用于读取 fp 选择的 power
    pw_try <- try(getSmo(fit_fp, what = "mu")$power, silent = TRUE)
    if (!inherits(pw_try, "try-error") && !is.null(pw_try)) {
      pw <- pw_try
      message("fp 选出的 powers = ", paste(pw, collapse = ", "))
    } else {
      message("未能读取 fp 的 powers，回退默认 powers = 1,3")
    }
  } else {
    message("fp 搜索拟合失败，回退默认 powers = 1,3")
  }
}

# TimeTransformed 不应为负值：直接剔除负值行
neg_tt <- which(raw_data$TimeTransformed < 0)
if (length(neg_tt) > 0) {
  message(sprintf("检测到 TimeTransformed<0 的行：%d，已剔除", length(neg_tt)))
  raw_data <- raw_data[raw_data$TimeTransformed >= 0, ]
  raw_data <- droplevels(raw_data)
}

# 如果 TimeTransformed 含非正值，且 fp 给出非整数幂，会导致 bfpNA 产生 NaN
if (any(raw_data$TimeTransformed <= 0) && any(abs(pw - round(pw)) > 1e-8)) {
  message("检测到 TimeTransformed<=0 且 fp 给出非整数幂，回退到 powers = 1,3 以避免数值错误")
  pw <- c(1, 3)
}

# =========================
# 4) helper: 用 BrainChart Extract.Wrapper 拟合一次（带兜底修复 NA 系数）
# =========================
fit_one <- function(mu_formula, sigma_formula, nu_formula, fam_string) {
  HOLDER <- list()
  HOLDER$SUBSET <- raw_data
  HOLDER$MODEL <- list(
    mu = mu_formula,
    sigma = sigma_formula,
    nu = nu_formula,
    family = fam_string,  # 这里用字符串，让 Extract.Wrapper/Apply.Param 走同一套逻辑
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
  
  # 关键兜底：如果 mu 的 fixef 中有 NA（设计矩阵奇异/不可辨识），置 0 防止 Apply.Param 预测炸
  if (anyNA(EX$param$mu$fixef)) {
    na_names <- names(EX$param$mu$fixef)[is.na(EX$param$mu$fixef)]
    message("mu.fixef 出现 NA（不可辨识列），将这些系数置 0：", paste(na_names, collapse = ", "))
    EX$param$mu$fixef[is.na(EX$param$mu$fixef)] <- 0
  }
  
  HOLDER$FIT.EXTRACT <- EX
  list(EXTRACT = EX, HOLDER = HOLDER)
}

# =========================
# 5) 模型候选（从最贴近目标 -> 更稳的降级）
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
  message("尝试拟合：", m$name)
  RES <- fit_one(m$mu, m$sigma, m$nu, m$fam)
  if (!is.null(RES)) { message("拟合成功：", m$name); break }
}
if (is.null(RES)) stop("三套模型都拟合失败。请检查数据列名/分布要求/或把 Extract.Wrapper 的报错贴出。")

EXTRACT <- RES$EXTRACT
HOLDER  <- RES$HOLDER

# =========================
# 6) NEWData 网格（population curves）
# =========================
time_seq <- seq(min(raw_data$TimeTransformed), max(raw_data$TimeTransformed), length.out = 200)
NEWData <- expand.grid(
  TimeTransformed = time_seq,
  Grp = levels(raw_data$Grp)
)

# =========================
# 7) Normative trajectories：PI/centile band (2.5/50/97.5)
#    + Normative variance：点估计（方差）
# =========================
CURVE <- Apply.Param(
  NEWData = NEWData,
  FITParam = EXTRACT$param,
  Pred.Set = c("l025"=0.025, "m500"=0.5, "u975"=0.975),
  Add.Moments = TRUE,       # 需要 variance
  Add.Derivative = FALSE,
  MissingToZero = TRUE,
  NAToZero = TRUE,
  verbose = FALSE
)

need <- c("PRED.l025.pop","PRED.m500.pop","PRED.u975.pop","PRED.variance.pop")
miss <- setdiff(need, names(CURVE))
if (length(miss) > 0) stop(paste("Apply.Param 未生成关键列：", paste(miss, collapse=", "), "\n实际列名：", paste(names(CURVE), collapse=", ")))

# =========================
# 8) Normative variance 的 95% CI（bootstrap over subjects）
#    注意：这是“方差曲线的不确定性 CI”，不是 PI
# =========================
#B <- 200
B <- 10
Base.Seed <- 12345
var_mat <- matrix(NA_real_, nrow = nrow(CURVE), ncol = B)

for (b in 1:B) {
  EX_b <- try(Boot.Function(n = b, Base.Seed = Base.Seed, Holder = HOLDER, Create.INIT = TRUE), silent = TRUE)
  if (inherits(EX_b, "try-error") || is.null(EX_b$param)) next
  
  # 同样兜底：某次 bootstrap 出现 NA fixef 就置 0，防止预测整列 NA
  if (anyNA(EX_b$param$mu$fixef)) EX_b$param$mu$fixef[is.na(EX_b$param$mu$fixef)] <- 0
  
  CURVE_b <- try(
    Apply.Param(
      NEWData = NEWData,
      FITParam = EX_b$param,
      Pred.Set = NULL,          # variance CI 不需要每次算分位数
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
    title = "Normative trajectories (PI: 2.5–97.5 centiles)",
    x = "TimeTransformed",
    y = "Wand"
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
  ggsave("BrainChart_Figure1_like.png", fig, width = 9, height = 7, dpi = 300)
} else {
  print(p_traj); print(p_var)
  ggsave("BrainChart_Trajectories.png", p_traj, width = 9, height = 4, dpi = 300)
  ggsave("BrainChart_Variance.png", p_var, width = 9, height = 4, dpi = 300)
}

message("=== 完成：输出图已保存到当前工作目录 ===")
