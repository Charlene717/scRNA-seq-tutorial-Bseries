# =====================================================================
# B15_trajectory.R — B15「軌跡分析與擬時序：細胞狀態的時間軸」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 15 集
# 資料：真演算法模擬的分化資料（一條主幹、兩個分支）
#       ＋ PBMC 3k 反例（承接 output/pbmc_b10_annotated.rds）
# 環境：R >= 4.3、Seurat v5.x、slingshot + SingleCellExperiment、
#       tradeSeq（Bioconductor）、pheatmap
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# 模擬資料不需下載任何檔案；PBMC 反例段需要先跑完 B5–B10。
# 分節標記與影片段落一一對應。
#
# 為什麼用模擬資料？——PBMC 是離散的成熟型別，「不適合軌跡」正是
# 本集要教的重點（雷一）。主線改用自寫的分化模擬（splatter 的
# path 概念：基因程式沿真實進度 t 漸變、分支點後走 A 或 B），
# 好處是每顆細胞的真實進度 t.true 已知，pseudotime 能對答案。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(slingshot)             # BiocManager::install("slingshot")
library(SingleCellExperiment)
set.seed(1234)                 # 全系列固定 seed，結果可重現

dir.create("output", showWarnings = FALSE)

## ---- 1. simulate --------------------------------------------------
# 模擬一個分化過程：1,500 顆細胞、2,000 個基因。
#   t.true ∈ [0,1]：每顆細胞的真實分化進度（模擬已知，等下對答案）
#   branch：t < 0.45 在主幹（0）；之後隨機走分支 A（1）或 B（2）
# 四組基因程式（每組 150 個基因）：
#   p.start  起點程式：隨 t 衰減（幹性、增殖）
#   p.mature 共同成熟：隨 t 上升（兩個分支都會走）
#   p.brA / p.brB 分支特異：分支點之後才啟動
n.cells <- 1500
n.genes <- 2000
tb      <- 0.45                          # 分支點
t.true  <- sort(runif(n.cells))
branch  <- ifelse(t.true < tb, 0L, sample(1:2, n.cells, replace = TRUE))

base <- rgamma(n.genes, shape = 1.1, rate = 1.5) + 0.05  # 基礎表達

make_prog <- function(idx, scale = 8) {
  p <- numeric(n.genes)
  p[idx] <- rgamma(length(idx), shape = 4, rate = 4 / scale)
  p
}
p.start  <- make_prog(1:150)
p.mature <- make_prog(151:300)
p.brA    <- make_prog(301:450)
p.brB    <- make_prog(451:600)

after <- pmax(t.true - tb, 0) / (1 - tb)          # 分支點後的進度
mu <- base +                                       # 細胞 × 基因 的期望值
  outer(1 - t.true, p.start) +
  outer(t.true,     p.mature) +
  outer(after * (branch == 1), p.brA) +
  outer(after * (branch == 2), p.brB)
mu <- mu * rlnorm(n.cells, 0, 0.35)                # 每顆細胞的深度差

counts <- matrix(rnbinom(length(mu), mu = t(mu), size = 2),
                 nrow = n.genes,
                 dimnames = list(sprintf("G%04d", 1:n.genes),
                                 sprintf("cell_%d", 1:n.cells)))
# 負二項取樣（size = 2 是 scRNA 常見的離散度量級）；
# t(mu)：rnbinom 按列填，轉置後才是 基因 × 細胞

## ---- 2. preprocess ------------------------------------------------
# B5–B9 的標準流程原封不動；軌跡分析接在「分群之後」。
sim <- CreateSeuratObject(counts = counts, project = "sim_traj")
sim$t.true <- t.true                      # 標準答案收進 metadata
sim$branch <- factor(branch, labels = c("trunk", "A", "B"))

sim <- NormalizeData(sim) |> FindVariableFeatures() |>
       ScaleData() |> RunPCA(npcs = 20, verbose = FALSE) |>
       RunUMAP(dims = 1:10, verbose = FALSE)
sim <- FindNeighbors(sim, dims = 1:10) |>
       FindClusters(resolution = 0.4)
table(sim$seurat_clusters)

# 確認 root 的生物學依據：起點程式（G0001–G0150）分數最高的群。
# 真實資料裡這一步是「查幹性 marker」；模擬資料裡我們自己算。
sim <- AddModuleScore(sim, features = list(rownames(sim)[1:150]),
                      name = "start.score")
root.clus <- names(which.max(
  tapply(sim$start.score1, sim$seurat_clusters, mean)))
root.clus                                  # 這一集後面都用它當 root

## ---- 3. slingshot -------------------------------------------------
# 3a. getLineages：cluster 質心 + MST，決定路徑骨架與分岔。
#     重點參數：
#       reducedDim / rd —— 在「PCA 空間」算（雷三的反面）
#       start.clus      —— root，由上面的生物學依據指定
sce <- as.SingleCellExperiment(sim)
rd  <- reducedDim(sce, "PCA")[, 1:2]

lin <- getLineages(rd, clusterLabels = as.character(sim$seurat_clusters),
                   start.clus = root.clus)
slingLineages(lin)          # 預期兩條 lineage：主幹 → A、主幹 → B

# 3b. getCurves：每條 lineage 迭代擬合 principal curve
crv <- getCurves(lin)

# 3c. pseudotime：細胞 × lineage 的矩陣；NA = 不在該 lineage 上
pt <- slingPseudotime(crv)
head(pt, 3)
sim$pt1 <- pt[colnames(sim), 1]
sim$pt2 <- pt[colnames(sim), 2]

FeaturePlot(sim, "pt1", reduction = "umap") +
  scale_colour_viridis_c(na.value = "grey85") +
  ggtitle("Lineage1 pseudotime")

# 3d. 對答案：pseudotime vs 真實進度（模擬資料才有的奢侈）
ok <- !is.na(sim$pt1)
cor(sim$pt1[ok], sim$t.true[ok], method = "spearman")
# 預期 ~0.99：順序抓對了。注意 pseudotime 的「數值」仍然
# 不是物理時間——它是曲線上的弧長，是排序不是時鐘。

# 一行檢查（雷一的解藥）：密度要連續、沒有斷崖
hist(sim$pt1, breaks = 40, main = "hist(pseudotime)")

## ---- 4. root-experiment -------------------------------------------
# 換 root 重跑：唯一的改動是 start.clus 改成 Lineage1 的末端群。
# MST 不變、曲線不變，只有「方向」整條反轉。
end.clus <- tail(slingLineages(lin)$Lineage1, 1)

lin2 <- getLineages(rd, clusterLabels = as.character(sim$seurat_clusters),
                    start.clus = end.clus)
crv2 <- getCurves(lin2)
pt2  <- slingPseudotime(crv2)

cor(pt[, 1], pt2[, 1], method = "spearman", use = "complete.obs")
# 預期 ~ -1：同一批細胞，順序完全倒過來。
# 後果示範：成熟程式的基因（G0151–G0300）沿正確方向「上升」；
# root 錯了會被讀成「下降」——生物學結論整段反向。
mature.gene <- "G0180"
df <- data.frame(pt.good = pt[, 1], pt.bad = pt2[, 1],
                 expr = LayerData(sim, layer = "data")[mature.gene, ])
op <- par(mfrow = c(1, 2))
plot(df$pt.good, df$expr, pch = 16, cex = .3,
     xlab = "pseudotime (root 正確)", ylab = mature.gene)
plot(df$pt.bad, df$expr, pch = 16, cex = .3,
     xlab = "pseudotime (root 亂選)", ylab = mature.gene)
par(op)

## ---- 5. tradeseq --------------------------------------------------
# 沿軌跡的基因動態。tradeSeq 的完整版（fitGAM 全基因）在筆電上
# 要跑數分鐘到數十分鐘；影片示範用 HVG 前 500 個基因。
# 沒裝 tradeSeq 的話，可以跳到 5b 的簡化版（loess smoother 熱圖）。
if (requireNamespace("tradeSeq", quietly = TRUE)) {
  library(tradeSeq)
  hvg <- head(VariableFeatures(sim), 500)
  gam <- fitGAM(counts = counts[hvg, ], sds = crv, nknots = 6)
  assoc <- associationTest(gam)          # 檢定：曲線是不是水平線
  head(assoc[order(assoc$pvalue), ], 5)
  top <- rownames(assoc)[order(assoc$pvalue)][1:50]
} else {
  # 備援：直接拿三組已知程式的基因畫熱圖（教學用）
  top <- rownames(sim)[c(1:17, 151:167, 451:467)]
}

# 5b. 沿軌跡表達熱圖：Lineage1 細胞依 pseudotime 排序、
#     每個基因先做 loess smoother 再畫（讀 smoother，不讀點）
ord  <- order(sim$pt1, na.last = NA)
expr <- as.matrix(LayerData(sim, layer = "data")[top, ord])
ptx  <- sort(sim$pt1[!is.na(sim$pt1)])
sm <- t(apply(expr, 1, function(y)
  predict(loess(y ~ ptx, span = 0.3), ptx)))
sm <- t(scale(t(sm)))                     # 每列 z 分數
if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(sm, cluster_cols = FALSE, cluster_rows = TRUE,
                     show_colnames = FALSE,
                     main = "genes x pseudotime (Lineage1)")
}

## ---- 6. pbmc-counterexample ---------------------------------------
## [踩雷示範] 雷一：離散型別硬拉軌跡。
## PBMC 的 T / B / Mono / NK 是成熟型別，彼此之間沒有分化連續體；
## 但 slingshot 照樣回傳一條「路徑」。整段預設註解掉——
## 打開跑一次、看完 hist() 的斷崖之後，請把它關回去。
# pbmc <- readRDS("output/pbmc_b10_annotated.rds")   # B10 的產物
# sce.p <- as.SingleCellExperiment(pbmc)
# sce.p <- slingshot(sce.p,
#                    clusterLabels = as.character(pbmc$celltype.manual),
#                    reducedDim = "PCA")
# pt.p <- slingPseudotime(sce.p)
# hist(pt.p[, 1], breaks = 40,
#      main = "PBMC forced trajectory: canyons")     # 密度斷崖
# # 對照第 3 節模擬資料的 hist(sim$pt1)：連續 vs 斷崖，一眼判生死。

## [踩雷示範] 雷三：在 UMAP 座標上算 pseudotime。
## 只要把 reducedDim 換成 "UMAP"，pseudotime 就建立在
## 被壓縮與拉伸過的距離上（B9）。跑完與 PCA 版比排序：
# sce.u <- slingshot(sce, clusterLabels = as.character(sim$seurat_clusters),
#                    reducedDim = "UMAP", start.clus = root.clus)
# pt.u <- slingPseudotime(sce.u)
# cor(pt[, 1], pt.u[, 1], method = "spearman", use = "complete.obs")
# # 大方向像、局部被攪亂；扭曲多嚴重取決於那一次 UMAP 的佈局。

## ---- 7. save ------------------------------------------------------
# 把 pseudotime 與 lineage 資訊都收在 Seurat 物件裡存檔
sim$root.clus <- root.clus
saveRDS(sim, "output/sim_b15_trajectory.rds")

sessionInfo()
