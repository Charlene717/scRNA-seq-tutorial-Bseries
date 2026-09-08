# =====================================================================
# B18_ecosystem.R — B18「生態系與收官：多組學、空間、Python 世界與下一步」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 18 集（收官）
# 資料：既有物件——WNN 段用 SeuratData 的 bmcite（骨髓 CITE-seq），
#       h5ad 段用 B10 的產物 output/pbmc_b10_annotated.rds
# 環境：R >= 4.3、Seurat v5.x、SeuratData（bmcite）、SeuratDisk
#       Python 端（選配）：scanpy、anndata —— 驗收碼以註解收錄
#
# 使用方式：從專案根目錄 source() 或逐段執行。
#   - WNN 段需要 bmcite（約幾百 MB），預設 run.wnn <- FALSE，
#     InstallData("bmcite") 之後改 TRUE 再跑。
#   - h5ad 段需要 SeuratDisk 與 B10 產物；兩者缺一該段自動跳過。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
set.seed(1234)                 # 全系列固定 seed，結果可重現

dir.create("output", showWarnings = FALSE)

run.wnn <- FALSE               # 需要 SeuratData::bmcite（CITE-seq 物件）；
                               # SeuratData::InstallData("bmcite") 後改 TRUE

## ---- 1. wnn -------------------------------------------------------
# CITE-seq 下游與 WNN（影片頁 10–12）。
# bmcite：人類骨髓 CITE-seq，約 3 萬顆細胞、25 個抗體（ADT assay 已附）。
# 若你手上是自己的 CITE-seq 資料，從 Read10X 起手的裝配寫法見影片頁 10：
#   inputs <- Read10X("data/citeseq/filtered_feature_bc_matrix/")
#   obj <- CreateSeuratObject(counts = inputs[["Gene Expression"]])
#   obj[["ADT"]] <- CreateAssayObject(counts = inputs[["Antibody Capture"]])
if (run.wnn) {
  library(SeuratData)
  bm <- LoadData("bmcite")
  bm <- UpdateSeuratObject(bm)             # 舊版物件升級到 v5

  # 1a. 兩個 assay 各自 normalize：RNA 照 B6；ADT 用 CLR（margin = 2 對細胞）
  DefaultAssay(bm) <- "RNA"
  bm <- NormalizeData(bm)
  bm <- NormalizeData(bm, assay = "ADT",
                      normalization.method = "CLR", margin = 2)

  # 1b. 兩個模態各自降維
  DefaultAssay(bm) <- "RNA"
  bm <- bm |> FindVariableFeatures() |> ScaleData() |>
        RunPCA(verbose = FALSE)                       # -> "pca"
  DefaultAssay(bm) <- "ADT"
  # ADT 只有 25 個 feature：不做 HVG，全部拿去 scale + PCA
  VariableFeatures(bm) <- rownames(bm[["ADT"]])
  bm <- bm |> ScaleData() |>
        RunPCA(reduction.name = "apca", verbose = FALSE)

  # 1c. WNN：兩組 KNN -> 逐細胞模態權重 -> 加權鄰居圖
  #   reduction.list / dims.list —— 兩個模態各用哪個空間、各取幾維
  #   （RNA 30 維承 B7 的判斷；ADT 25 個抗體取 18 維）
  bm <- FindMultiModalNeighbors(bm,
    reduction.list = list("pca", "apca"),
    dims.list      = list(1:30, 1:18))

  # 1d. 下游照舊，只是輸入換成加權圖
  bm <- RunUMAP(bm, nn.name = "weighted.nn",
                reduction.name = "wnn.umap")
  bm <- FindClusters(bm, graph.name = "wsnn",
                     resolution = 0.8, verbose = FALSE)

  # 1e. 每顆細胞的 RNA 權重（1 - ADT 權重）——畫回 UMAP 看哪些族群靠誰
  summary(bm$RNA.weight)
  p.w <- FeaturePlot(bm, features = "RNA.weight",
                     reduction = "wnn.umap")
  ggplot2::ggsave("output/b18_wnn_rna_weight.png", p.w,
                  width = 7, height = 6, dpi = 150)

  saveRDS(bm, "output/bm_b18_wnn.rds")
}

## ---- 2. h5ad-export -----------------------------------------------
# Seurat -> h5ad（影片頁 24）。兩步：SaveH5Seurat -> Convert。
# 預設行為：X = DefaultAssay 的 data 層（normalized）、raw = counts、
#           scale.data 不帶、KNN 圖不帶、多的 assay 要另外處理（坑三）。
has.disk <- requireNamespace("SeuratDisk", quietly = TRUE)
has.pbmc <- file.exists("output/pbmc_b10_annotated.rds")

if (has.disk && has.pbmc) {
  library(SeuratDisk)
  pbmc <- readRDS("output/pbmc_b10_annotated.rds")   # B10 的成品

  SaveH5Seurat(pbmc, filename = "output/pbmc_b10.h5Seurat",
               overwrite = TRUE)
  Convert("output/pbmc_b10.h5Seurat", dest = "h5ad",
          overwrite = TRUE)
  # 產出 output/pbmc_b10.h5ad —— 接下來第 3、4 節做三關驗收
} else {
  message("[skip] h5ad 段：需要 SeuratDisk 與 output/pbmc_b10_annotated.rds")
}

## ---- 3. checks ----------------------------------------------------
# 三個坑的 R 端檢查碼（影片頁 25、32、34、36）。
if (has.pbmc) {
  pbmc <- readRDS("output/pbmc_b10_annotated.rds")

  # 3a. 坑一（方向）：轉出前記下「基因 × 細胞」的正確 shape
  dim(pbmc)                    # 基因數在前、細胞數在後（R 端慣例）
  #  -> Python 端 adata.shape 應該是「細胞數在前」且兩數對調
  head(colnames(pbmc), 3)      # barcode —— 應對應 Python 端 obs_names
  head(rownames(pbmc), 3)      # 基因名 —— 應對應 Python 端 var_names

  # 3b. 坑二（型別）：逐欄記下 meta.data 的型別，轉入後與 dtypes 並排對
  sapply(pbmc@meta.data, class)
  #  factor 欄位（如 seurat_clusters）到 Python 應成為 category；
  #  層級順序要緊的欄位，先把順序另存一份：
  lvl.backup <- lapply(Filter(is.factor, pbmc@meta.data), levels)
  str(lvl.backup)

  # 3c. 坑三（圖層）：盤點哪些東西「不會」跟著 h5ad 過去
  Layers(pbmc)                 # scale.data 在這裡，但預設不轉
  names(pbmc@graphs)           # KNN 圖：一定不轉，Python 端重建
  Assays(pbmc)                 # 第二個以上的 assay（如 ADT）要另外轉
}

## ---- 4. python-side -----------------------------------------------
# Python 端驗收碼（影片頁 25、34）。以註解收錄，開個終端機照打：
#
#   import scanpy as sc
#   adata = sc.read_h5ad("output/pbmc_b10.h5ad")
#   adata                      # 三關：shape / obs 欄位 / obsm
#   adata.shape                # (2638, 13714) —— 細胞數在前才對
#   adata.obs_names[:3]        # barcode（對照 3a 的 colnames）
#   adata.var_names[:3]        # 基因名（對照 3a 的 rownames）
#   adata.obs.dtypes           # 對照 3b 的 sapply(..., class)
#   adata.obs["seurat_clusters"] = (
#       adata.obs["seurat_clusters"].astype("category"))  # 掉了就補
#   adata.obsm.keys()          # X_pca / X_umap 應該都在
#   # 要重建的：sc.pp.scale()（scale.data）、sc.pp.neighbors()（KNN 圖）

## ---- 5. save ------------------------------------------------------
# 轉檔紀錄：把「轉了什麼、掉了什麼」寫成一張表，跟檔案一起交接
if (has.disk && has.pbmc) {
  handoff <- data.frame(
    item   = c("X (= data layer)", "raw (= counts)", "obs (= meta.data)",
               "obsm (pca/umap)", "scale.data", "graphs", "extra assays"),
    status = c("carried", "carried", "carried",
               "carried", "DROPPED - rescale in Python",
               "DROPPED - rebuild with sc.pp.neighbors",
               "DROPPED - export separately"))
  write.csv(handoff, "output/b18_h5ad_handoff.csv", row.names = FALSE)
  print(handoff)
}

## [踩雷示範] 坑一：手動搬矩陣把方向弄反。
## 自己用矩陣拼 AnnData（例如經 reticulate / 手寫 h5）時，
## 忘了轉置就是經典災難——R 的矩陣直接丟過去是「基因 × 細胞」。
## 整段預設註解掉；打開看 dim() 的對比就好，不要真的拿去分析。
# mat.wrong <- LayerData(pbmc, layer = "counts")     # 基因 × 細胞
# dim(mat.wrong)                                     # 13714 × 2638
# mat.right <- Matrix::t(mat.wrong)                  # 細胞 × 基因
# dim(mat.right)                                     # 2638 × 13714
# # AnnData 要的是 mat.right 這個方向（cells × genes）

## [踩雷示範] 坑二：factor 降級成整數碼。
## 有些繞路的轉法會把 factor 存成底層整數（1,2,3...），
## 型別檢查抓得到：
# f <- pbmc$seurat_clusters
# head(as.integer(f))          # 底層整數碼——存成這個就回不去了
# head(as.character(f))        # 轉檔前先轉字串，才保得住標籤本身

sessionInfo()
