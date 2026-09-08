# =====================================================================
# B07_pca.R — B7「PCA 與維度選擇：到底取幾個 PC」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 7 集
# 資料：PBMC 3k（承接 B6 的 normalization 後物件 output/pbmc_b06_lognorm.rds；
#       若你沒跑 B6，本腳本會自動退回 B5 的 QC 後物件
#       output/pbmc_b05_qc.rds，重跑 LogNormalize 路線三行）
# 環境：R >= 4.3、Seurat v5.x、ggplot2
#
# 使用方式：從專案根目錄 source() 或逐段執行（全程約 2–3 分鐘）。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(ggplot2)
set.seed(1234)                       # 全系列固定 seed，結果可重現

# 承接點：優先讀 B6 的產物（已含 NormalizeData / HVG / ScaleData）。
# B6 若用了不同檔名，改下面第一個路徑即可，不影響本集其他內容。
if (file.exists("output/pbmc_b06_lognorm.rds")) {
  pbmc <- readRDS("output/pbmc_b06_lognorm.rds")
} else {
  # 備援：從 B5 的 QC 後物件起手，重跑 LogNormalize 路線三行
  stopifnot("找不到 output/pbmc_b05_qc.rds，請先跑完 R/B05_qc.R" =
              file.exists("output/pbmc_b05_qc.rds"))
  pbmc <- readRDS("output/pbmc_b05_qc.rds")
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
                               nfeatures = 2000)
  pbmc <- ScaleData(pbmc)
}
pbmc
# 預期：13714 features across ~2638 samples、2000 variable features

## ---- 1. runpca ----------------------------------------------------
# PCA 跑在 2000 個 HVG 的 scale.data 上，不是全部基因。
# npcs 預設 50：先多算存著，取幾個是下游 dims 的事，不用回頭重算。
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 三種看輸出的方式（影片頁 11–14）：
# (a) 文字版：每個 PC 正負兩端的 top 基因——PC 的「配方表」
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# 預期：PC_1 正端 CST3/TYROBP/LST1/AIF1/FTL（骨髓系），
#       負端 MALAT1/LTB/IL32/IL7R/CD2（淋巴系）

# (b) 圖版 loadings：同一份資訊畫成條狀圖
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave("output/B07_viz_loadings.png", width = 10, height = 5, dpi = 200)

# (c) 散點：細胞落在 PC1 x PC2 平面上的位置
DimPlot(pbmc, reduction = "pca") + NoLegend()
ggsave("output/B07_pca_dimplot.png", width = 6.5, height = 5.5,
       dpi = 200)

## ---- 2. depth-check -----------------------------------------------
# 一行檢查：PC1 是生物學，還是深度殘留？（雷二的預防針）
# |r| 小（<0.3 上下）→ 深度沒有主宰 PC1；r 很高 → 回 B6 檢查
cor(pbmc$nCount_RNA, Embeddings(pbmc, "pca")[, 1])

# 上色看一眼：深度應該「散在各群裡」，而不是沿著某個 PC 排成梯度
FeaturePlot(pbmc, features = "nCount_RNA", reduction = "pca")
ggsave("output/B07_pca_depth.png", width = 6.5, height = 5.5, dpi = 200)

# 順手檢查 PC1 兩端的 top loadings 有沒有被核糖體基因洗版
top.pc1 <- names(sort(abs(Loadings(pbmc, "pca")[, 1]),
                      decreasing = TRUE))[1:20]
mean(grepl("^RP[LS]", top.pc1))   # 比例接近 1 就要警惕（雷二）

## ---- 3. dimheatmap ------------------------------------------------
# 每個 PC 一張小熱圖：分數最極端的 500 顆細胞 x top loading 基因。
# balanced = TRUE 正負兩端各取一半，兩個方向的程式都看得到。
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# 一次看 15 個：訊號從區塊分明 → 逐漸變成鹽粒雜訊
png("output/B07_dimheatmap_15.png", width = 1800, height = 2200,
    res = 150)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

## ---- 4. elbow -----------------------------------------------------
# 方法一：ElbowPlot。y 軸是每個 PC 的標準差，找曲線變平的位置
ElbowPlot(pbmc, ndims = 50)
ggsave("output/B07_elbowplot.png", width = 7, height = 4.5, dpi = 200)

# 方法二：把「拐點」變成數字（量化 cutoff，寫得進方法段）
pct  <- pbmc[["pca"]]@stdev / sum(pbmc[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1  <- which(cumu > 90 & pct < 5)[1]     # 累積 90% 且單 PC < 5%
co2  <- sort(which(diff(pct) < -0.1), decreasing = TRUE)[1] + 1
c(cumulative90 = co1, drop0.1 = co2)
# 判讀：co2（相鄰差 < 0.1% 的位置）通常小很多；
# 兩個數字框出一個區間，拐點目測值落在裡面就放心。

## ---- 5. stability -------------------------------------------------
# 方法三：穩定度法。nPC 掃一輪，各自跑下游分群，看結果動不動。
# （FindNeighbors/FindClusters 的細節 B8 講；這裡先當黑盒子用）
for (d in c(5, 10, 15, 20, 30)) {
  pbmc <- FindNeighbors(pbmc, dims = 1:d, verbose = FALSE)
  pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)
  pbmc@meta.data[[paste0("npc", d)]] <- Idents(pbmc)
}
sapply(pbmc@meta.data[paste0("npc", c(5, 10, 15, 20, 30))],
       \(x) length(unique(x)))
# 預期：nPC=5 群數偏少；10 以上群數與成員幾乎不再變動（平原區）

# 跨 run 對照永遠用交叉表（B8 雷二的預習）：
# 對角線乾淨 = 兩個 nPC 給出幾乎相同的分群
table(npc10 = pbmc$npc10, npc30 = pbmc$npc30)

# JackStraw 的歷史定位（v5 淡出，教學保留概念即可）：
# 把資料隨機打亂再跑 PCA，看真 PC 是否顯著強於隨機——
# 統計上漂亮，但慢、對大資料過於敏感，且 SCT 流程不支援。
# Seurat v5 官方教學已改用 ElbowPlot + 經驗法則。想跑可以：
#   pbmc <- JackStraw(pbmc, num.replicate = 100)   # 只支援 LogNormalize
#   pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#   JackStrawPlot(pbmc, dims = 1:20)

## ---- 6. finalize --------------------------------------------------
# 定案：PBMC 3k 乾淨、族群分明，三法同指 8–12 區間 → dims = 1:10
#（與 Seurat 官方教學一致；B8 的 FindNeighbors 承接這個決定）。
# 方法段三要素：RunPCA(npcs = 50)、ElbowPlot + 量化 cutoff
#（co1/co2）、nPC 5–30 穩定度掃描；複雜組織建議往 20–30 取。
npc.final <- 10

# 用定案的 nPC 重建正式的分群欄位，清掉掃描用的暫存欄
pbmc <- FindNeighbors(pbmc, dims = 1:npc.final)
pbmc <- FindClusters(pbmc, resolution = 0.5)
for (d in c(5, 10, 15, 20, 30)) pbmc@meta.data[[paste0("npc", d)]] <- NULL

saveRDS(pbmc, "output/pbmc_b07_pca.rds")   # 下游從這份接手
sessionInfo()                              # 收進附錄；審稿人的好朋友

# =====================================================================
# [踩雷示範] 以下兩段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：nPC 抓太少，亞群直接消失（影片頁 32–33）
## dims = 1:5 只帶進前五個表現程式；住在 PC 6 以後的亞群差異
## 從此進不了下游——分群、UMAP、註釋全都看不到它：
# pbmc.bad1 <- FindNeighbors(pbmc, dims = 1:5)
# pbmc.bad1 <- FindClusters(pbmc.bad1, resolution = 0.5)
# table(npc5 = Idents(pbmc.bad1), npc10 = pbmc$seurat_clusters)
## 交叉表會看到 nPC=10 的兩個群在 nPC=5 併成同一格。
## 一行檢查：ElbowPlot 上 5 明顯還在陡坡上，就不該停在 5。

## [踩雷示範] 雷二：把 PC1 當生物學，不檢查就開始詮釋（頁 34–36）
## 跳過 normalization 直接 scale + PCA，PC1 會變成深度軸、
## top loadings 被核糖體/管家基因洗版：
# pbmc.bad2 <- readRDS("output/pbmc_b05_qc.rds")
# pbmc.bad2 <- FindVariableFeatures(pbmc.bad2)   # 故意跳過 Normalize
# pbmc.bad2 <- ScaleData(pbmc.bad2)
# pbmc.bad2 <- RunPCA(pbmc.bad2)
# cor(pbmc.bad2$nCount_RNA, Embeddings(pbmc.bad2, "pca")[, 1])
## r 會衝到 0.9 上下——那條軸量的是「誰測得深」，不是細胞身分。
## 一行檢查（正式流程每次都做）：cor(nCount, PC1) + 看 top loadings。
