# =====================================================================
# B06_normalization.R — B6「Normalization、HVG、Scaling：兩條路線實測」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 6 集
# 資料：PBMC 3k（承接 B5 的 QC 完成物件 output/pbmc_b05_qc.rds）
# 環境：R >= 4.3、Seurat v5.x、sctransform、patchwork
#
# 使用方式：從專案根目錄 source() 或逐段執行（SCT 段約 1–3 分鐘）。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
set.seed(1234)                       # 全系列固定 seed，結果可重現

qc_rds <- "output/pbmc_b05_qc.rds"
stopifnot("找不到 B5 的輸出，請先跑完 R/B05_qc.R" = file.exists(qc_rds))

pbmc <- readRDS(qc_rds)
pbmc
# 預期：13714 features across 2638 samples，1 layer present: counts

## ---- 1. normalize -------------------------------------------------
# 目的：把「深度」除掉。CP10K（除以總量 ×10,000）+ log1p，一行完成。
# 兩個參數都是預設值，寫出來是讓你知道它們存在（影片頁 11）。
pbmc <- NormalizeData(pbmc,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)

# 前後對比：counts 原封不動，data 是「新增」的 layer（影片頁 13）
LayerData(pbmc, layer = "counts")["LYZ", 1:3]
LayerData(pbmc, layer = "data")["LYZ", 1:3]
Layers(pbmc)                          # "counts" "data"

## ---- 2. hvg -------------------------------------------------------
# 目的：從 1.3 萬個基因挑出 2000 個「變異超出技術預期」的基因。
# vst：對 mean-variance 期望曲線的殘差排名（吃 counts，不是 data）。
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",
                             nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
top10   # PPBP / LYZ / GNLY…——HVG 名單本身就有生物學（影片頁 16）

# mean-variance 圖 + 前 10 名標籤
p <- VariableFeaturePlot(pbmc)
LabelPoints(plot = p, points = top10, repel = TRUE)
ggsave("output/B06_hvg_plot.png", width = 8, height = 5.5, dpi = 200)

## ---- 3. scale -----------------------------------------------------
# 目的：每個基因 z-score（量級拉平）+ clip ±10（防單一離群細胞）。
# 預設只 scale HVG 那 2000 個——scale.data 是稠密矩陣，全做很貴。
pbmc <- ScaleData(pbmc)
LayerData(pbmc, layer = "scale.data")["LYZ", 1:3]   # 有正有負、個位數

# 需要全基因（例如 DoHeatmap 的基因不在 HVG 裡）時再補跑：
#   pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# vars.to.regress：本系列預設不用（影片頁 23 的建議）。
# 確認干擾存在、且它不是你的生物學之後，才回頭加，並寫進方法段：
#   pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

## ---- 4. sct -------------------------------------------------------
# 路線二：SCTransform（v2）。一行替代三步，輸出住在新的 SCT assay。
# install.packages("sctransform")   # 第一次需要
pbmc.sct <- SCTransform(pbmc, vst.flavor = "v2", verbose = FALSE)
pbmc.sct
# 預期：Active assay 變成 SCT（~1.2 萬 features、3000 HVG），RNA 原封不動

DefaultAssay(pbmc.sct)               # "SCT"
head(VariableFeatures(pbmc.sct), 10)

# 鐵則：SCT 的 data = log1p(校正後 counts)，出廠即 log，絕不再 log！
LayerData(pbmc.sct, assay = "SCT", layer = "data")["LYZ", 1:3]

## ---- 5. compare ---------------------------------------------------
# 檢查點一：HVG 名單重疊（SCT 預設 3000，公平起見取前 2000）
hvg.log <- VariableFeatures(pbmc)
hvg.sct <- head(VariableFeatures(pbmc.sct), 2000)
length(intersect(hvg.log, hvg.sct))  # PBMC 3k 上約六成（影片頁 31）

# 檢查點二：兩條路線各自跑到 UMAP（PCA/UMAP 參數 B7、B9 詳講）
pbmc     <- RunPCA(pbmc, verbose = FALSE)
pbmc     <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
pbmc.sct <- RunPCA(pbmc.sct, verbose = FALSE)
pbmc.sct <- RunUMAP(pbmc.sct, dims = 1:10, verbose = FALSE)

p1 <- DimPlot(pbmc)     + ggtitle("LogNormalize")
p2 <- DimPlot(pbmc.sct) + ggtitle("SCTransform")
p1 + p2
ggsave("output/B06_umap_two_routes.png", width = 11, height = 5, dpi = 200)

# 本集保命一行：PC1 不該是深度軸（接近 1 = 中招；影片頁 39）
cor(Embeddings(pbmc, "pca")[, 1], pbmc$nCount_RNA)

## ---- 6. save ------------------------------------------------------
saveRDS(pbmc,     "output/pbmc_b06_lognorm.rds")  # B7 開場直接讀這份
saveRDS(pbmc.sct, "output/pbmc_b06_sct.rds")      # 對照用，B11 前還會提到
sessionInfo()                         # 收進附錄；審稿人的好朋友

# =====================================================================
# [踩雷示範] 以下三段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：跳過 normalization 直接 PCA（影片頁 38–39）
## z-score 救不了細胞間的深度差異，PC1 會變成定序深度軸：
# pbmc.bad1 <- readRDS(qc_rds)
# pbmc.bad1 <- FindVariableFeatures(pbmc.bad1)      # vst 吃 counts，跑得動
# pbmc.bad1 <- ScaleData(pbmc.bad1)                 # 注意：沒跑 NormalizeData
# pbmc.bad1 <- RunPCA(pbmc.bad1, verbose = FALSE)
# cor(Embeddings(pbmc.bad1, "pca")[, 1], pbmc.bad1$nCount_RNA)
# # 接近 1 → PC1 = 深度。對照正常流程的同一行檢查（5. compare 的結尾）

## [踩雷示範] 雷二：HVG 只抓 200 個（影片頁 40–41）
## 稀有亞型的 marker 擠不進前 200，該族群在 PC 空間失去座標：
# pbmc.bad2 <- NormalizeData(readRDS(qc_rds))
# pbmc.bad2 <- FindVariableFeatures(pbmc.bad2, nfeatures = 200)
# pbmc.bad2 <- ScaleData(pbmc.bad2) |> RunPCA(verbose = FALSE) |>
#              RunUMAP(dims = 1:10, verbose = FALSE)
# DimPlot(pbmc.bad2)   # 與 B06_umap_two_routes.png 對照：小群糊掉/消失

## [踩雷示範] 雷三：SCT 之後又手動 log（影片頁 42–43）
## SCT 的 data 已在 log 空間，再 log 一次 → 對比塌掉、FeaturePlot 灰掉：
# sct.data <- LayerData(pbmc.sct, assay = "SCT", layer = "data")["LYZ", ]
# summary(sct.data)              # 正常：陽性群平均 ~3
# summary(log1p(sct.data))       # 雙重轉換：塌到 ~1.4，差距縮一半
# hist(sct.data); hist(log1p(sct.data))   # 兩張分布圖並排看
