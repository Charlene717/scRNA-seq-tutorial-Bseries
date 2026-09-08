# =====================================================================
# B17_deconvolution.R — B17「反卷積：用單細胞參考拆解 bulk 資料」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 17 集
# 資料：pbmc3k（承接 B10 註釋後物件 output/pbmc_b10_annotated.rds）
#       自建 signature matrix 與 20 個已知比例的假 bulk，
#       NNLS 實跑還原 → 與真值對答案
# 環境：R >= 4.3、Seurat v5.x、dplyr、tidyr、ggplot2、nnls（CRAN）
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# 全程不需下載新資料；沒跑過 B10 的備援方案見第 0 節註解。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(nnls)                        # install.packages("nnls")
set.seed(1234)                       # 全系列固定 seed，結果可重現

if (file.exists("output/pbmc_b10_annotated.rds")) {
  pbmc <- readRDS("output/pbmc_b10_annotated.rds")
} else {
  # 備援：沒跑過 B10 時，退回 B08 的分群物件並手動掛型別標籤。
  # cluster 編號 → 型別的對應來自 B10 的 marker 判讀（0 號起算）。
  # 若你的分群結果編號不同，請先跑 B10 或自行核對 marker 再對照。
  stopifnot("找不到 output/pbmc_b08_clustered.rds，請先跑 R/B08" =
              file.exists("output/pbmc_b08_clustered.rds"))
  pbmc <- readRDS("output/pbmc_b08_clustered.rds")
  Idents(pbmc) <- "RNA_snn_res.0.6"
  new.ids <- c("Naive CD4 T", "CD14 Mono", "Memory CD4 T",
               "B", "CD8 T", "NK", "FCGR3A Mono",
               "Effector CD8 T", "DC", "Platelet")
  names(new.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.ids)
  pbmc$celltype.manual <- Idents(pbmc)
}

## ---- 1. reference -------------------------------------------------
# 反卷積的第一個決策：拆到多細？
# CD4 的兩個亞群、Mono 的兩個亞群 signature 太像（共線性），
# bulk 資料拆不動——先合併成 6 大類。Platelet 只有 14 顆，退場。
merge.map <- c("Naive CD4 T" = "CD4 T", "Memory CD4 T"   = "CD4 T",
               "CD8 T"       = "CD8 T", "Effector CD8 T" = "CD8 T",
               "CD14 Mono"   = "Mono",  "FCGR3A Mono"    = "Mono",
               "B" = "B", "NK" = "NK", "DC" = "DC")
pbmc$ct <- merge.map[as.character(pbmc$celltype.manual)]
pbmc <- subset(pbmc, subset = !is.na(ct))
table(pbmc$ct)
# 預期：B 344 / CD4 T 1178 / CD8 T 311 / DC 32 / Mono 594 / NK 165

## ---- 2. signature -------------------------------------------------
# 2a. 每型別的 marker：只留有區辨力的基因。
#     沒有區辨力的持家基因整列一樣亮，只會把方程變難解。
Idents(pbmc) <- "ct"
mk <- FindAllMarkers(pbmc, only.pos = TRUE,
                     min.pct = 0.25, logfc.threshold = 1)
top <- mk |> group_by(cluster) |>
       slice_max(avg_log2FC, n = 25) |> pull(gene) |> unique()

# 2b. 每型別平均表達（counts 加總後 CPM 化），取 marker 列。
avg <- AggregateExpression(pbmc, group.by = "ct")$RNA
sig <- t(t(avg) / colSums(avg)) * 1e4          # 每欄總和 = 1e4
sig <- as.matrix(sig[intersect(top, rownames(avg)), ])
dim(sig)                                       # ~138 × 6

# 2c. 一行檢查（本集鐵則）：欄與欄的相關 = 共線性體檢。
#     CD4 T vs CD8 T 接近 0.8 —— 這筆「共線性稅」在對答案時會現形。
round(cor(sig), 2)

## ---- 3. pseudobulk ------------------------------------------------
# 自產 ground truth：比例先抽好、存檔，答案在我們手上。
# Dirichlet 抽樣用 rgamma 實作（正規化的 Gamma 即 Dirichlet），
# 濃度參數以 pbmc 的實際比例為中心 × 80，讓 20 個 bulk 各不相同。
types <- colnames(sig)
alpha <- as.numeric(table(pbmc$ct)[types]) / ncol(pbmc) * 80
truth <- t(replicate(20, { g <- rgamma(6, alpha); g / sum(g) }))
colnames(truth) <- types

cnt  <- GetAssayData(pbmc, layer = "counts")
bulk <- sapply(1:20, function(b) {
  cells <- unlist(lapply(types, function(tp)
    sample(colnames(pbmc)[pbmc$ct == tp],
           round(500 * truth[b, tp]), replace = TRUE)))
  Matrix::rowSums(cnt[, cells])      # 500 顆細胞加總 = 一個假 bulk
})
bulk <- t(t(bulk) / colSums(bulk)) * 1e4       # 每個 bulk 也 CPM 化

## ---- 4. nnls ------------------------------------------------------
# 核心解法，透明到可以在白板上講完：
#   min || b - A p ||^2   subject to  p >= 0
# 解完手動正規化成總和 = 1（NNLS 本身不保證總和）。
A <- sig
fit  <- nnls(A, bulk[rownames(sig), 1])
prop <- fit$x / sum(fit$x)
round(setNames(prop, types), 3)      # 第 1 個 bulk 的估計

est <- t(apply(bulk[rownames(sig), ], 2, function(b) {
  x <- nnls(A, b)$x
  x / sum(x)
}))
colnames(est) <- types

# MuSiC 的加權思想（概念）：跨個體越穩定的基因權重越高。
# pbmc3k 只有一位捐贈者，示範不了跨個體變異；多捐贈者參考請直上
# MuSiC 套件（music_prop()），輸入是 SingleCellExperiment + bulk 矩陣。

## ---- 5. evaluate --------------------------------------------------
df <- bind_rows(
  as.data.frame(truth) |> mutate(bulk = 1:20, what = "truth"),
  as.data.frame(est)   |> mutate(bulk = 1:20, what = "est")) |>
  pivot_longer(all_of(types), names_to = "ct") |>
  pivot_wider(names_from = what, values_from = value)

cor(df$truth, df$est)                # 整體相關，預期 ~0.99
sqrt(mean((df$truth - df$est)^2))    # RMSE，預期 ~0.02

# 5a. 估計 vs 真值散點（影片頁 28 的 R 版）
ggplot(df, aes(truth, est, colour = ct)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(size = 2, alpha = .8) +
  coord_equal() +
  labs(x = "true proportion", y = "NNLS estimate", colour = NULL) +
  theme_classic()
ggsave("output/B17_scatter_truth_vs_est.png",
       width = 6, height = 5.2, dpi = 200)

# 5b. 每型別誤差：豐富型別準、稀有型別飄、CD8 吃共線性稅
err <- df |> group_by(ct) |>
  summarise(true.mean = mean(truth),
            mae       = mean(abs(est - truth)),
            rel.err   = mae / true.mean) |>
  arrange(desc(true.mean))
err
ggplot(err, aes(reorder(ct, -true.mean), rel.err)) +
  geom_col(fill = "#2563A8") +
  labs(x = NULL, y = "relative error (MAE / true mean)") +
  theme_classic()
ggsave("output/B17_per_type_error.png", width = 6.5, height = 4, dpi = 200)

## ---- 6. save ------------------------------------------------------
saveRDS(list(signature = sig,
             truth     = truth,
             bulk      = bulk,
             est       = est,
             per.type  = err),
        "output/b17_deconv_results.rds")

sessionInfo()                        # 收進附錄；審稿人的好朋友

# =====================================================================
# [踩雷示範] 以下三段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：reference 缺型別（頁 34–35）
## 把 NK 從 signature 拿掉重解：NK 的份被塞給最像的 CD8 T，
## CD8 平均比例從 ~12% 被灌成 ~22%，全程零報錯。
# A.nonk <- sig[, setdiff(types, "NK")]
# est.nonk <- t(apply(bulk[rownames(sig), ], 2, function(b) {
#   x <- nnls(A.nonk, b)$x
#   x / sum(x)
# }))
# colnames(est.nonk) <- setdiff(types, "NK")
# round(rbind(truth.CD8 = mean(truth[, "CD8 T"]),
#             full.CD8  = mean(est[, "CD8 T"]),
#             nonk.CD8  = mean(est.nonk[, "CD8 T"])), 3)
## 一行檢查：拆真 bulk 前先問「樣本裡的主要型別，參考都有欄嗎？」

## [踩雷示範] 雷二：跨平台直接套（頁 36–37）
## 模擬長度偏誤：把 bulk 依基因長度乘上偏誤因子（長基因訊號偏多），
## 再用未校正的 10x signature 去解——長基因多的型別被系統性高估。
# gene.len <- setNames(runif(nrow(sig), 0.5, 10), rownames(sig))
# bias     <- (gene.len / mean(gene.len))^0.6      # 全長協定的近似
# bulk.fl  <- bulk[rownames(sig), ] * bias
# est.fl <- t(apply(bulk.fl, 2, function(b) {
#   x <- nnls(A, b)$x; x / sum(x)
# }))
# round(colMeans(est.fl) - colMeans(est), 3)       # 系統性偏移方向
## 修法：同平台驗證、或用 CIBERSORTx 的 batch correction。

## [踩雷示範] 雷三：比例差異當細胞數差異（頁 38–39）
## 成分資料：總和鎖死在 1。只讓 CD8「真的」擴增 4 倍，
## 其他型別的比例全部被動下降——但它們的細胞數一顆都沒少。
# n.a <- c("CD4 T" = 450, "CD8 T" = 120, "B" = 130,
#          "NK" = 60, "Mono" = 230, "DC" = 12)
# n.b <- n.a; n.b["CD8 T"] <- 480
# round(rbind(before = n.a / sum(n.a), after = n.b / sum(n.b)), 3)
## 比較請在 CLR 轉換後做（centered log-ratio，一行）：
# clr <- function(p) log(p) - mean(log(p))
# round(rbind(before = clr(n.a / sum(n.a)),
#             after  = clr(n.b / sum(n.b))), 2)
## CLR 後只有 CD8 明顯移動——其他型別的「下降」消失了。
