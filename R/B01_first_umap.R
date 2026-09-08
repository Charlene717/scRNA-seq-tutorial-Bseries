# =====================================================================
# B01_first_umap.R — B1「環境建置與你的第一張 UMAP」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 1 集
# 資料：PBMC 3k（10x 官方 filtered_gene_bc_matrices, hg19；
#       先執行 R/B00_setup.R 下載並解壓）
# 環境：R >= 4.3、Seurat v5.x、tidyverse（只用到 ggplot2）、renv
#
# 使用方式：從專案根目錄 source() 或逐段執行（全程約 2–3 分鐘）。
# 本集是 quickstart：每一步只做「跑通」，參數全部照預設或教學慣例；
# 每段註解都標了「哪一集講透」；分節與影片段落一一對應。
# =====================================================================

## ---- 0. setup ------------------------------------------------------
# 安裝檢查：缺什麼裝什麼。第一次執行可能要幾分鐘。
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")            # CRAN 上即為 v5.x
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
# renv 建議在 Console 手動跑一次，不放在腳本裡自動執行：
#   renv::init()      # 建立專案獨立套件庫
#   renv::snapshot()  # 套件裝齊後，把版本寫進 renv.lock

library(Seurat)
set.seed(1234)   # 全系列固定的隨機種子；同參數重跑，結果才會一致

# 版本驗證：Seurat 要 5.x、R 要 4.3 以上。求救之前先跑這兩行。
packageVersion("Seurat")
R.version.string
stopifnot(packageVersion("Seurat") >= "5.0.0")

## ---- 1. project ----------------------------------------------------
# 專案資料夾（已存在就跳過）。路徑一律相對於專案根目錄——
# 本系列所有腳本都不用 setwd()，開 .Rproj 工作目錄自動就是專案根。
for (d in c("data", "R", "output")) {
  dir.create(d, showWarnings = FALSE)
}

## ---- 2. data-check -------------------------------------------------
# 資料還沒下載的話，先跑 R/B00_setup.R（約 7.6 MB）。
mtx.dir <- "data/pbmc3k/filtered_gene_bc_matrices/hg19"
stopifnot("找不到 pbmc3k，請先執行 R/B00_setup.R" =
            dir.exists(mtx.dir))
list.files(mtx.dir)   # 應有 barcodes.tsv / genes.tsv / matrix.mtx

## ---- 3. load -------------------------------------------------------
# Read10X 吃的是「資料夾」路徑（不是單一檔案），它把三個檔案
# 組回一個稀疏矩陣。細節 B3、B5 會講。
pbmc.data <- Read10X(data.dir = mtx.dir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc   # 13714 features across 2700 samples

## ---- 4. qc-rough ---------------------------------------------------
# 粗略 QC：閾值今天照教學慣例抄（200–2500、mt < 5%）。
# 這些線該怎麼由「你自己資料的分布」決定，是 B5 整集的主題。
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 &
                 nFeature_RNA < 2500 & percent.mt < 5)
pbmc   # 2700 → 2638 顆

## ---- 5. preprocess -------------------------------------------------
# 四行都用預設值；每個預設值背後的假設，B6（前三行）與 B7（RunPCA）
# 會逐一拆解。今天先不求甚解。
pbmc <- NormalizeData(pbmc)            # 校正定序深度      → B6
pbmc <- FindVariableFeatures(pbmc)     # 挑 2000 個 HVG    → B6
pbmc <- ScaleData(pbmc)                # 壓到同一尺度      → B6
pbmc <- RunPCA(pbmc)                   # 降到主成分空間    → B7
ElbowPlot(pbmc)                        # 「取幾個 PC」→ B7

## ---- 6. cluster-umap -----------------------------------------------
# dims = 1:10 與 resolution = 0.5 都是教學慣例值；
# 怎麼選 dims 看 B7，resolution 掃描與 clustree 看 B8，UMAP 參數看 B9。
pbmc <- FindNeighbors(pbmc, dims = 1:10)        # 建鄰居圖  → B8
pbmc <- FindClusters(pbmc, resolution = 0.5)    # 找社群    → B8
pbmc <- RunUMAP(pbmc, dims = 1:10)              # 壓成二維  → B9
table(Idents(pbmc))   # 9 群（0–8）

## ---- 7. first-umap -------------------------------------------------
# 你的第一張 UMAP。這些群是誰（T 細胞？B 細胞？），B10 教你回答。
DimPlot(pbmc, reduction = "umap", label = TRUE)
ggplot2::ggsave("output/B01_first_umap.png",
                width = 7, height = 5.5, dpi = 200)

## ---- 8. save -------------------------------------------------------
# 收工三件套：物件存 rds、圖已存 output、sessionInfo 記版本。
saveRDS(pbmc, "output/pbmc_b01_first_run.rds")
sessionInfo()

## =====================================================================
## [踩雷示範] 以下皆為「會出錯／不該用」的寫法，預設整段註解。
## 想看後果可逐行解開執行，看完記得改回來。
## =====================================================================

## [踩雷示範] 雷一：v4 時代的教學碼，在 Seurat v5 直接報錯 -----------
## v4 的矩陣放在 @counts slot；v5 改放 layers，舊寫法會炸：
# counts.v4 <- pbmc@assays$RNA@counts
#   # Error: no slot of name "counts" for this object of class "Assay5"
## v5 的正確寫法：
# counts.v5 <- LayerData(pbmc, assay = "RNA", layer = "counts")
## 一行檢查：你的物件是 v5 結構嗎？layers 有哪些？
# Layers(pbmc)      # 應回 "counts" "data" "scale.data"

## [踩雷示範] 雷二：setwd() + 絕對路徑，換台電腦第一行就爆 -----------
# setwd("C:/Users/ming/Desktop/analysis_final_v2")
#   # Error in setwd(...) : cannot change working directory
## 正解：RStudio Project + 相對路徑（本腳本全程如此），
## 檢查目前工作目錄是不是專案根：
# getwd()           # 應以 /scrna-course 之類的專案名結尾
