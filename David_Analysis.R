# =========================
# 0) packages
# =========================
pkgs <- c("gamlss", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(gamlss)
library(ggplot2)

# patchwork 仅用于拼图；没装也不影响核心结果
has_patchwork <- requireNamespace("patchwork", quietly=TRUE)
if (has_patchwork) library(patchwork)

# =========================
# 1) source repo scripts (你现在已补齐 101/102)
# =========================
source("100.common-variables.r")
source("101.common-functions.r")  # 会 source 102.gamlss-recode.r  :contentReference[oaicite:4]{index=4}
source("300.variables.r")
source("301.functions.r")

# =========================
# 2) load data
# =========================
raw_data <- read.csv("Simulated_Brain_Data.csv")
raw_data$Grp <- as.factor(raw_data$Grp)
raw_data$Study <- as.factor(raw_data$Study)
raw_data$INDEX.TYPE <- as.factor(raw_data$INDEX.TYPE)

# 这些 attributes/INDEX 列是 pipeline 需要的（101 里会检查）:contentReference[oaicite:5]{index=5}
attr(raw_data, "columns") <- list(Index = c("INDEX.ID","INDEX.OB","INDEX.TYPE"))
attr(raw_data, "tag") <- "SimulationData"
attr(raw_data, "Transformations") <- list(Type="identity")

# =========================
# 3) build HOLDER + model
#    用 bfpNA + GGalt（来自 102 的重写）
# =========================
HOLDER <- list()
HOLDER$SUBSET <- raw_data
HOLDER$MODEL <- list(
  mu    = as.formula(Wand ~ bfpNA(TimeTransformed, powers=c(1,3)) + Grp + random(Study)),
  sigma = as.formula(~ TimeTransformed),
  nu    = as.formula(~ 1),
  family = "GGalt",
  covariates = list(
    Y="Wand", X="TimeTransformed", BY="Grp",
    RANEF="Study", COND="INDEX.TYPE", ID="INDEX.ID",
    OTHER=c("INDEX.OB")
  ),
  stratify = NULL,
  inc.fp = FALSE
)

# =========================
# 4) fit
# =========================
EXTRACT <- Extract.Wrapper(HOLDER, Store.Full=TRUE)
# 给 Boot.Function 用（它内部会用 Holder$FIT.EXTRACT$param 做 INIT）
HOLDER$FIT.EXTRACT <- EXTRACT

# （可选）如果 mu 的 fixef 仍有 NA，直接停止：否则 Apply.Param 会生成不了 PRED（你之前就是这里坏的）
if (anyNA(EXTRACT$param$mu$fixef)) stop("mu 的 fixef 出现 NA：模型不可辨识/设计矩阵奇异。请先简化公式排查。")

# =========================
# 5) population grid
# =========================
time_seq <- seq(min(raw_data$TimeTransformed), max(raw_data$TimeTransformed), length.out=200)
NEWData <- expand.grid(TimeTransformed=time_seq, Grp=levels(raw_data$Grp))

# =========================
# 6) population trajectories (PI: 2.5/50/97.5 centiles) + variance curve (point estimate)
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

# sanity
need <- c("PRED.l025.pop","PRED.m500.pop","PRED.u975.pop","PRED.variance.pop")
miss <- setdiff(need, names(CURVE))
if (length(miss)>0) stop(paste("Apply.Param 未生成：", paste(miss, collapse=", ")))

# =========================
# 7) Normative variance 的 95% CI（bootstrap over subjects）
#    这是“variance curve 的 CI”，不是“individual outcome 的 PI”
# =========================
B <- 200                 # bootstrap 次数：想更稳可调大（会更慢）
Base.Seed <- 12345

var_mat <- matrix(NA_real_, nrow=nrow(CURVE), ncol=B)

for (b in 1:B) {
  EX_b <- Boot.Function(n=b, Base.Seed=Base.Seed, Holder=HOLDER, Create.INIT=TRUE)
  
  CURVE_b <- Apply.Param(
    NEWData = NEWData,
    FITParam = EX_b$param,
    Pred.Set = NULL,            # variance 不需要分位数
    Add.Moments = TRUE,
    Add.Derivative = FALSE,
    MissingToZero = TRUE,
    NAToZero = TRUE,
    verbose = FALSE
  )
  
  if (!("PRED.variance.pop" %in% names(CURVE_b))) next
  var_mat[, b] <- CURVE_b$PRED.variance.pop
}

CURVE$VAR.ci.l025 <- apply(var_mat, 1, function(x) quantile(x, 0.025, na.rm=TRUE))
CURVE$VAR.ci.u975 <- apply(var_mat, 1, function(x) quantile(x, 0.975, na.rm=TRUE))

# =========================
# 8) plots
# =========================
p_traj <- ggplot(CURVE, aes(x=TimeTransformed, group=Grp, color=Grp, fill=Grp)) +
  geom_ribbon(aes(ymin=PRED.l025.pop, ymax=PRED.u975.pop), alpha=0.2, colour=NA) +
  geom_line(aes(y=PRED.m500.pop), linewidth=1) +
  labs(title="Normative trajectories (PI: 2.5–97.5 centiles)", x="TimeTransformed", y="Wand") +
  theme_minimal()

p_var <- ggplot(CURVE, aes(x=TimeTransformed, group=Grp, color=Grp, fill=Grp)) +
  geom_ribbon(aes(ymin=VAR.ci.l025, ymax=VAR.ci.u975), alpha=0.2, colour=NA) +
  geom_line(aes(y=PRED.variance.pop), linewidth=1) +
  labs(title="Normative variance (CI: bootstrap 95%)", x="TimeTransformed", y="Variance") +
  theme_minimal()

if (has_patchwork) {
  print(p_traj / p_var)
} else {
  print(p_traj); print(p_var)
}

ggsave("BrainChart_Figure1_like.png",
       if (has_patchwork) (p_traj / p_var) else p_traj,
       width=9, height=6, dpi=300)
