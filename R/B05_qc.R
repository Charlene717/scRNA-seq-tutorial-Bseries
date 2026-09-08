# =====================================================================
# B05_qc.R — B5「讀入資料與 QC：閾值到底怎麼定」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 5 集
# 資料：PBMC 3k（10x Genomics 官方公開資料，filtered, hg19）
# 環境：R >= 4.3、Seurat v5.x、DoubletFinder（GitHub 版）
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup ----------------------------------------------------
library(Seurat)
library(ggplot2)
set.seed(1234)                       # 全系列固定 seed，結果可重現

# 資料還沒下載的話，先跑 data/download_data.R（見 repo README）
data_dir <- "data/pbmc3k/filtered_gene_bc_matrices/hg19/"
stopifnot("找不到資料資料夾，請先執行 data/download_data.R" =
            dir.exists(data_dir))
dir.create("output", showWarnings = FALSE)

## ---- 1. load ------------------------------------------------------
# Read10X 吃「資料夾」路徑（內含 matrix.mtx / genes.tsv / barcodes.tsv）
pbmc.data <- Read10X(data.dir = data_dir)

# min.cells / min.features 只是建物件時的粗篩，不是 QC 本體（影片頁 11）
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc
# 預期：13714 features across 2700 samples

## ---- 2. metrics ---------------------------------------------------
# 人類粒線體基因是 MT- 開頭；小鼠要改 pattern = "^mt-"（影片頁 12 的坑）
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)

## ---- 3. inspect ---------------------------------------------------
# 鐵則：先看分布，再定閾值。順序不能反（影片雷二）。
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)
ggsave("output/B05_qc_violins.png", width = 12, height = 5, dpi = 200)

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave("output/B05_scatter_mt.png", width = 6, height = 5, dpi = 200)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("output/B05_scatter_nfeature.png", width = 6, height = 5, dpi = 200)

# 過濾前先留一份原始物件——你的後悔藥（影片頁 23、37）
saveRDS(pbmc, "output/pbmc_raw.rds")

## ---- 4. thresholds ------------------------------------------------
# 路線一：MAD 法（percent.mt 直接算；R 的 mad() 已含 1.4826 校正）
mt.upper <- median(pbmc$percent.mt) + 3 * mad(pbmc$percent.mt)
mt.upper                              # PBMC 3k 約 4.9

# 右偏嚴重的指標（nCount / nFeature）建議 log 後再算 MAD（影片頁 21）
log.nf   <- log10(pbmc$nFeature_RNA)
nf.upper <- 10 ^ (median(log.nf) + 3 * mad(log.nf))
nf.lower <- 10 ^ (median(log.nf) - 3 * mad(log.nf))
c(lower = nf.lower, upper = nf.upper)

# 路線二：看圖的固定閾值（本集採用，與 Seurat 官方教學一致）
#   nFeature 200–2500、percent.mt < mt.upper
# 兩派可混用；重點是理由寫得出來（影片頁 19、22）。

n.before <- ncol(pbmc)
pbmc <- subset(pbmc,
               subset = nFeature_RNA > 200 &
                        nFeature_RNA < 2500 &
                        percent.mt   < mt.upper)
n.after <- ncol(pbmc)

# 過濾紀錄：這幾個數字要進分析筆記與論文方法段（影片頁 24）
cat(sprintf("QC filter: %d -> %d cells (removed %d, %.1f%%)\n",
            n.before, n.after, n.before - n.after,
            100 * (n.before - n.after) / n.before))

## ---- 5. doublets --------------------------------------------------
# DoubletFinder 需要標準前處理完成後的 PC 空間。
# 這五步是 B6–B9 的主角，此處先照抄（影片頁 29）。
# 安裝：remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# pK 掃描：用掃的，不要用猜的（影片頁 29）
sweep.res  <- paramSweep(pbmc, PCs = 1:10, sct = FALSE)
sweep.stat <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn      <- find.pK(sweep.stat)
pk.best <- as.numeric(as.character(
  bcmvn$pK[which.max(bcmvn$BCmetric)]))
pk.best

# 預期 doublet 數：回收 ~2,600 顆 → 約 2.3%（每千顆 +0.8% 經驗法則）
# 進階：同型 doublet 校正（McGinnis et al. 2019）——
#   homotypic.prop <- modelHomotypic(<分群標籤>)   # B8 之後才有標籤
#   nExp.adj <- round(nExp * (1 - homotypic.prop))
nExp <- round(0.023 * ncol(pbmc))
pbmc <- doubletFinder(pbmc, PCs = 1:10, pN = 0.25,
                      pK = pk.best, nExp = nExp)

df.col <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)
table(pbmc@meta.data[[df.col]])

# 本集採「標記觀察」路線：doublet 標籤留在 metadata，B8 分群後回頭檢查。
# 若要直接移除（出最終結果前建議）：
#   pbmc <- subset(pbmc, cells = colnames(pbmc)[pbmc@meta.data[[df.col]] == "Singlet"])
DimPlot(pbmc, group.by = df.col) + ggtitle("DoubletFinder classifications")
ggsave("output/B05_doublets_umap.png", width = 7, height = 5.5, dpi = 200)

## ---- 6. save ------------------------------------------------------
saveRDS(pbmc, "output/pbmc_b05_qc.rds")   # B6 開場直接讀這份
sessionInfo()                              # 收進附錄；審稿人的好朋友

# =====================================================================
# [踩雷示範] 以下三段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：percent.mt 一律砍 5%（影片頁 34–35）
## 在高代謝組織（心肌、腫瘤、肝）會整型細胞消失。PBMC 剛好倖免——
## 這正是「運氣好」的例子：
# pbmc.bad1 <- subset(readRDS("output/pbmc_raw.rds"), percent.mt < 5)
# ncol(pbmc.bad1)   # 對 PBMC 差異不大；對心肌資料是團滅

## [踩雷示範] 雷二：先過濾再看分布（影片頁 36–37）
## 截斷後的分布無法回推原始的線該畫哪：
# pbmc.bad2 <- subset(readRDS("output/pbmc_raw.rds"), percent.mt < 8)
# VlnPlot(pbmc.bad2, features = "percent.mt")   # 尾巴消失，分布被切平

## [踩雷示範] 雷三：QC 一輪定案，不回頭（影片頁 38–39）
## 分群後的殘渣群檢查（B8 會正式做）：
# qc.by.cluster <- aggregate(pbmc@meta.data[, c("percent.mt", "nFeature_RNA")],
#                            by = list(cluster = pbmc$seurat_clusters), FUN = median)
# qc.by.cluster   # 某群 mt 高 + nFeature 低 → 殘渣群，回 QC 修線重跑
