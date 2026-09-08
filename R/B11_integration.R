# =====================================================================
# B11_integration.R — B11「多樣本整合：batch effect 與 Harmony」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 11 集
# 資料：ifnb（Kang et al. 2018, Nat Biotechnol；GEO GSE96583）
#       對照（CTRL, 6,548 顆）vs IFN-β 刺激（STIM, 7,451 顆）的人類 PBMC
#       由 SeuratData 套件安裝，本集不承接 B10 的 pbmc 物件
# 環境：R >= 4.3、Seurat v5.x、SeuratData、harmony（CRAN）、
#       dplyr、ggplot2、patchwork
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# CCA 整合段最耗時（約 3–8 分鐘，視機器而定），全程約 10–15 分鐘。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
library(patchwork)
set.seed(1234)                       # 全系列固定 seed，結果可重現

dir.create("output", showWarnings = FALSE)

## ---- 1. load ------------------------------------------------------
# 第一次使用要先安裝 ifnb 資料集（約 100 MB，之後有快取）：
# InstallData("ifnb")
ifnb <- LoadData("ifnb")
ifnb <- UpdateSeuratObject(ifnb)     # SeuratData 物件是舊格式，升到 v5
ifnb
# 預期：14053 features across 13999 samples；1 layer: counts
# metadata 自帶兩個關鍵欄位：
#   stim               條件標籤（CTRL / STIM）
#   seurat_annotations 官方發表的細胞型別註釋（本集當「已知答案」用）
head(ifnb@meta.data, 3)

## ---- 2. split -----------------------------------------------------
table(ifnb$stim)                     # CTRL 6548 / STIM 7451，量級相當
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
# 預期：4 layers: counts.CTRL, counts.STIM, data.CTRL, data.STIM
# 沒有任何細胞被移動或刪除——同一個 assay 內分抽屜而已

## ---- 3. preprocess ------------------------------------------------
# split 狀態下，這四步對每個 layer「各自」執行（看 console 訊息）。
# 各自算的理由：合併算 HVG / scaling 會被批次差異污染，
# 干擾素反應基因會混進 HVG，反過來放大批次方向的訊號。
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)   # 每 layer 各 2000，物件層級彙整
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)

## ---- 4. baseline --------------------------------------------------
# 鐵則（雷一的解藥）：整合之前，先留一張未整合的圖當基線。
# reduction 與 cluster 都取獨立名字，之後才能與整合版本並排對照。
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 0.5,
                     cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca",
                reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated",
        group.by = c("stim", "unintegrated_clusters"))
ggsave("output/B11_umap_unintegrated.png", width = 11, height = 5, dpi = 200)
# 預期畫面：兩個條件各自成島；同一型別（用 seurat_annotations 對照）
# 在兩座島各出現一次
DimPlot(ifnb, reduction = "umap.unintegrated",
        group.by = "seurat_annotations", split.by = "stim")

## ---- 5. integrate -------------------------------------------------
# v5 統一介面：換方法只換 method 參數。三種都跑，reduction 便宜。
#   orig.reduction  吃哪個未整合的降維
#   new.reduction   校正結果（一個新 reduction；表達值不會被動）
# CCA：最強、最慢、最易過度校正。進階參數 k.anchor / k.weight 用預設；
#      若某樣本細胞數 < 100，k.weight 要調小，否則報錯。
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration,
                        orig.reduction = "pca",
                        new.reduction  = "integrated.cca",
                        verbose = FALSE)

# Harmony：快、省記憶體、可擴展。theta（預設 2）控制混合力道：
#          金絲雀檢查顯示校過頭時，第一個往下調的旋鈕就是它。
# install.packages("harmony")        # 第一次使用先安裝
ifnb <- IntegrateLayers(object = ifnb, method = HarmonyIntegration,
                        orig.reduction = "pca",
                        new.reduction  = "harmony",
                        verbose = FALSE)

# RPCA：保守、快，適合同平台多樣本、批次差異溫和的場合
ifnb <- IntegrateLayers(object = ifnb, method = RPCAIntegration,
                        orig.reduction = "pca",
                        new.reduction  = "integrated.rpca",
                        verbose = FALSE)
Reductions(ifnb)
# 預期："pca" "umap.unintegrated" "integrated.cca" "harmony"
#       "integrated.rpca"

## ---- 6. downstream ------------------------------------------------
# 與 B8/B9 完全同一套動作，唯一差別：reduction 指向校正後的空間。
# resolution 沿用 0.5 起步；正式分析請回 B8 流程掃 resolution + clustree。
ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.5,
                     cluster.name = "harmony_clusters")
ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:30,
                reduction.name = "umap.harmony")
DimPlot(ifnb, reduction = "umap.harmony",
        group.by = c("stim", "harmony_clusters"))
ggsave("output/B11_umap_harmony.png", width = 11, height = 5, dpi = 200)

# 作業：CCA 版本的下游，跑完與 harmony 版本並排比較
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.5,
                     cluster.name = "cca_clusters")
ifnb <- RunUMAP(ifnb, reduction = "integrated.cca", dims = 1:30,
                reduction.name = "umap.cca")
table(ifnb$harmony_clusters, ifnb$cca_clusters)   # 兩法分群的交叉表

# 之後的分析以 harmony 版本為正本
Idents(ifnb) <- "harmony_clusters"

## ---- 7. join ------------------------------------------------------
# 時機：整合與分群做完、DE / marker 檢定之前。
# 校正住在 reduction 裡，join 動的是表達值抽屜——兩者無關，不會互毀。
ifnb <- JoinLayers(ifnb)
ifnb                                  # 預期：2 layers: data, counts

## ---- 8. evaluate --------------------------------------------------
# 檢查一：split 圖。同一群應出現在左右兩張圖的同一位置
DimPlot(ifnb, reduction = "umap.harmony", split.by = "stim")
ggsave("output/B11_umap_split.png", width = 11, height = 5, dpi = 200)

# 檢查二：每群的條件組成。margin = 1 → 每列（每群）各自轉比例。
# 全體基準約 0.47 / 0.53；大群應接近對半，重度偏斜的群去查 marker
prop.table(table(ifnb$harmony_clusters, ifnb$stim), margin = 1) |>
  round(2)

# 檢查三：已知 marker 應跨條件亮同一個位置（幾何用整合、數值用 RNA）
FeaturePlot(ifnb, features = c("CD14", "CD3D", "MS4A1"),
            reduction = "umap.harmony", split.by = "stim")
ggsave("output/B11_marker_alignment.png", width = 10, height = 9, dpi = 200)

# 金絲雀檢查（雷二的解藥）：干擾素反應是刺激「應該」誘導的生物訊號。
# ISG 分數在 STIM 應明顯高於 CTRL；被整合抹平 = 過度整合的警訊
isg <- list(c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1"))
ifnb <- AddModuleScore(ifnb, features = isg, name = "ISG.score")
VlnPlot(ifnb, features = "ISG.score1", group.by = "stim", pt.size = 0)
# 用官方註釋對照：活化 T（"T activated"）等刺激相關狀態應仍可辨識
table(ifnb$seurat_annotations, ifnb$stim)

## ---- 9. save ------------------------------------------------------
saveRDS(ifnb, "output/ifnb_b11_integrated.rds")   # B12 開場直接讀這份

sessionInfo()                        # 收進附錄；審稿人的好朋友

# =====================================================================
# [踩雷示範] 以下三段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：跳過未整合的基線（頁 38–39）
## 「整合永遠會給你一張混得很好的圖」——沒有基線，你無從判斷
## 它修了什麼、修得對不對。示範：從頭來過、直接整合：
# ifnb2 <- LoadData("ifnb") |> UpdateSeuratObject()
# ifnb2[["RNA"]] <- split(ifnb2[["RNA"]], f = ifnb2$stim)
# ifnb2 <- NormalizeData(ifnb2) |> FindVariableFeatures() |>
#          ScaleData() |> RunPCA()
# ifnb2 <- IntegrateLayers(ifnb2, method = HarmonyIntegration,
#                          orig.reduction = "pca",
#                          new.reduction = "harmony")
# ifnb2 <- RunUMAP(ifnb2, reduction = "harmony", dims = 1:30)
# DimPlot(ifnb2, group.by = "stim")   # 混得很勻、很漂亮——然後呢？
## 「整合改變了什麼？」——這個物件裡沒有任何東西能回答。
## 解藥就是主流程第 4 節那幾行：reduction 取名、基線常駐物件。

## [踩雷示範] 雷二：把本來就同批的資料硬整合（頁 40–41）
## 把 CTRL 隨機切成兩個假「批次」（其實毫無批次差異），
## 硬跑 CCA 整合——校正力道全部花在雜訊上，稀有族群被拉散：
# ctrl <- subset(ifnb, subset = stim == "CTRL")
# ctrl$fake.batch <- sample(c("A", "B"), ncol(ctrl), replace = TRUE)
# ctrl[["RNA"]] <- split(ctrl[["RNA"]], f = ctrl$fake.batch)
# ctrl <- NormalizeData(ctrl) |> FindVariableFeatures() |>
#         ScaleData() |> RunPCA()
# ctrl <- IntegrateLayers(ctrl, method = CCAIntegration,
#                         orig.reduction = "pca",
#                         new.reduction = "integrated.fake",
#                         k.weight = 50)
# ctrl <- RunUMAP(ctrl, reduction = "integrated.fake", dims = 1:30)
# DimPlot(ctrl, group.by = "seurat_annotations", label = TRUE)
## 對照未整合版本：pDC、Mk 這類小族群的邊界被攪糊。
## 教訓：沒有證據顯示存在批次差異，就不要整合——整合不是免費的。

## [踩雷示範] 雷三：在校正後的數值上做 DE（頁 42–43）
## v5 的 IntegrateLayers 只給 reduction，這個雷主要埋在 v4 舊物件
## （IntegrateData 產生的 integrated assay）與「忘記 JoinLayers」裡。
## (a) 忘記 join：v5 的防呆會直接報錯，請你先 JoinLayers——
# ifnb.split <- ifnb; ifnb.split[["RNA"]] <-
#   split(ifnb.split[["RNA"]], f = ifnb.split$stim)
# FindMarkers(ifnb.split, ident.1 = "STIM", group.by = "stim")
## (b) v4 風格的錯誤（拿到別人的舊物件時最常見）：
# DefaultAssay(old.obj) <- "integrated"      # ← 校正值，不是測量值
# FindMarkers(old.obj, ident.1 = "STIM", group.by = "stim")
## 校正壓縮了變異、又在細胞間互相借值——檢定假設全毀，p 值不可信。
## 正解：DE 一律回 RNA 的表達值（JoinLayers 後的 counts / data）。
## 至於「回 RNA 之後怎麼做才對」（pseudobulk + DESeq2）→ B12。
