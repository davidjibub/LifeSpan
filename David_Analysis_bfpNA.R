# =========================
# 0) packages
# =========================
pkgs <- c("gamlss", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(gamlss)
library(ggplot2)

has_patchwork <- requireNamespace("patchwork", quietly=TRUE)
if (has_patchwork) library(patchwork)

# =========================
# 1) source repo scripts
# =========================
source("100.common-variables.r")
source("101.common-functions.r")   # loads 102.gamlss-recode.r  :contentReference[oaicite:3]{index=3}
source("300.variables.r")
source("301.functions.r")

# =========================
# 2) load + sanitize data
# =========================
raw_data <- read.csv("Simulated_Brain_Data.csv")

raw_data$Grp <- droplevels(as.factor(raw_data$Grp))
raw_data$Study <- droplevels(as.factor(raw_data$Study))
raw_data$INDEX.TYPE <- droplevels(as.factor(raw_data$INDEX.TYPE))

# 确保 Wand > 0 (GGalt/GG 都需要正值)
if (any(is.na(raw_data$Wand))) stop("Wand 含 NA：请先处理缺失值。")
min_y <- min(raw_data$Wand)
if (min_y <= 0) {
  eps <- abs(min_y) + 1e-6
  raw_data$Wand <- raw_data$Wand + eps
  message(sprintf("检测到 Wand<=0：已整体平移 +%.6g 以满足 GGalt 的 y>0 要求。", eps))
}

# attributes needed by pipeline checks  :contentReference[oaicite:4]{index=4}
attr(raw_data, "columns") <- list(Index = c("INDEX.ID","INDEX.OB","INDEX.TYPE"))
attr(raw_data, "tag") <- "SimulationData"
attr(raw_data, "Transformations") <- list(Type="identity")

# =========================
# 3) helper: fit with fallback + fix NA fixef
# =========================
fit_one <- function(mu_formula, sigma_formula, nu_formula, fam) {
  HOLDER <- list()
  HOLDER$SUBSET <- raw_data
  HOLDER$MODEL <- list(
    mu = mu_formula,
    sigma = sigma_formula,
    nu = nu_formula,
    family = fam,
    covariates = list(
      Y="Wand", X="TimeTransformed", BY="Grp",
      RANEF="Study", COND="INDEX.TYPE", ID="INDEX.ID", OTHER=c("INDEX.OB")
    ),
    stratify = NULL,
    inc.fp = FALSE
  )
  
  EX <- try(Extract.Wrapper(HOLDER, Store.Full=TRUE), silent=TRUE)
  if (inherits(EX, "try-error") || is.null(EX$param)) return(NULL)
  
  # 如果 mu fixef 有 NA：先尝试“置零修复”，保证 Apply.Param 可预测
  if (anyNA(EX$param$mu$fixef)) {
    na_names <- names(EX$param$mu$fixef)[is.na(EX$param$mu$fixef)]
    message("mu.fixef 出现 NA（不可辨识列）：将这些系数置 0 以避免预测失败：",
            paste(na_names, collapse=", "))
    EX$param$mu$fixef[is.na(EX$param$mu$fixef)] <- 0
  }
  
  HOLDER$FIT.EXTRACT <- EX
  list(EXTRACT=EX, HOLDER=HOLDER)
}

# =========================
# 4) try models (most similar -> most stable)
# =========================
models <- list(
  list(
    name="M1: bfpNA + Grp + random(Study)",
    mu = as.formula(Wand ~ bfpNA(TimeTransformed, powers=c(1,3)) + Grp + random(Study)),
    sigma = as.formula(~ TimeTransformed),
    nu = as.formula(~ 1),
    fam = "GGalt"
  ),
  list(
    name="M2: bfpNA + Grp (no random)",
    mu = as.formula(Wand ~ bfpNA(TimeTransformed, powers=c(1,3)) + Grp),
    sigma = as.formula(~ TimeTransformed),
    nu = as.formula(~ 1),
    fam = "GGalt"
  ),
  list(
    name="M3: bfpNA + Grp + Study (fixed Study)",
    mu = as.formula(Wand ~ bfpNA(TimeTransformed, powers=c(1,3)) + Grp + Study),
    sigma = as.formula(~ TimeTransformed),
    nu = as.formula(~ 1),
    fam = "GGalt"
  )
)

RES <- NULL
for (m in models) {
  message("尝试拟合：", m$name)
  RES <- fit_one(m$mu, m$sigma, m$nu, m$fam)
  if (!is.null(RES)) { message("拟合成功：", m$name); break }
}
if (is.null(RES)) stop("三套模型都拟合失败：请把 Extract.Wrapper 的完整报错贴出来（通常是数据格式/分布不满足）。")

EXTRACT <- RES$EXTRACT
HOLDER  <- RES$HOLDER

# =========================
# 5) population grid + Apply.Param (PI trajectories + variance point)
# =========================
time_seq <- seq(min(raw_data$TimeTransformed), max(raw_data$TimeTransformed), length.out=200)
NEWData <- expand.grid(TimeTransformed=time_seq, Grp=levels(raw_data$Grp))

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
if (length(miss)>0) stop(paste("预测仍缺列：", paste(miss, collapse=", "), "\n请把 names(CURVE) 发我。"))

# =========================
# 6) Normative variance: 95% CI via bootstrap (CI of variance curve, not PI)
# =========================
B <- 200
Base.Seed <- 12345

var_mat <- matrix(NA_real_, nrow=nrow(CURVE), ncol=B)

for (b in 1:B) {
  EX_b <- try(Boot.Function(n=b, Base.Seed=Base.Seed, Holder=HOLDER, Create.INIT=TRUE), silent=TRUE)
  if (inherits(EX_b, "try-error") || is.null(EX_b$param)) next
  
  # 同样修复：防止 bootstrap 某次出现 NA fixef 直接把整列变 NA
  if (anyNA(EX_b$param$mu$fixef)) EX_b$param$mu$fixef[is.na(EX_b$param$mu$fixef)] <- 0
  
  CURVE_b <- try(Apply.Param(
    NEWData = NEWData,
    FITParam = EX_b$param,
    Pred.Set = NULL,
    Add.Moments = TRUE,
    Add.Derivative = FALSE,
    MissingToZero = TRUE,
    NAToZero = TRUE,
    verbose = FALSE
  ), silent=TRUE)
  
  if (inherits(CURVE_b, "try-error")) next
  if (!("PRED.variance.pop" %in% names(CURVE_b))) next
  var_mat[, b] <- CURVE_b$PRED.variance.pop
}

CURVE$VAR.ci.l025 <- apply(var_mat, 1, function(x) quantile(x, 0.025, na.rm=TRUE))
CURVE$VAR.ci.u975 <- apply(var_mat, 1, function(x) quantile(x, 0.975, na.rm=TRUE))

# =========================
# 7) plots
# =========================
p_traj <- ggplot(CURVE, aes(x=TimeTransformed, group=Grp, color=Grp, fill=Grp)) +
  geom_ribbon(aes(ymin=PRED.l025.pop, ymax=PRED.u975.pop), alpha=0.2, colour=NA) +
  geom_line(aes(y=PRED.m500.pop), linewidth=1) +
  labs(title="Normative trajectories (PI: 2.5–97.5 centiles)",
       x="TimeTransformed", y="Wand") +
  theme_minimal()

p_var <- ggplot(CURVE, aes(x=TimeTransformed, group=Grp, color=Grp, fill=Grp)) +
  geom_ribbon(aes(ymin=VAR.ci.l025, ymax=VAR.ci.u975), alpha=0.2, colour=NA) +
  geom_line(aes(y=PRED.variance.pop), linewidth=1) +
  labs(title="Normative variance (CI: bootstrap 95%)",
       x="TimeTransformed", y="Variance") +
  theme_minimal()

if (has_patchwork) {
  fig <- p_traj / p_var
  print(fig)
  ggsave("BrainChart_Figure1_like.png", fig, width=9, height=7, dpi=300)
} else {
  print(p_traj); print(p_var)
  ggsave("BrainChart_Trajectories.png", p_traj, width=9, height=4, dpi=300)
  ggsave("BrainChart_Variance.png", p_var, width=9, height=4, dpi=300)
}
