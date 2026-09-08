# =====================================================================
# B09_umap.R — B9「UMAP / t-SNE：參數、隨機性與可重現性」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 9 集
# 資料：PBMC 3k（承接 B8 分群後物件 output/pbmc_b08_clustered.rds，
#       其中已含 pca reduction 與 resolution 0.6 的 10 群——見 R/B08_clustering.R）
# 環境：R >= 4.3、Seurat v5.x（UMAP 引擎為 R 原生 uwot，無需 Python）、
#       patchwork、ggplot2
#
# 使用方式：從專案根目錄 source() 或逐段執行（全程約 2–3 分鐘，
# 掃描迴圈最花時間）。分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
set.seed(1234)   # 第一層（全域）：管一般亂數；RunUMAP 內部另有 seed.use

stopifnot("找不到 output/pbmc_b08_clustered.rds，請先跑完 R/B08_clustering.R" =
            file.exists("output/pbmc_b08_clustered.rds"))
pbmc <- readRDS("output/pbmc_b08_clustered.rds")   # 2,638 顆、10 群
table(Idents(pbmc))                                # RNA_snn_res.0.6

## ---- 1. umap -------------------------------------------------------
# dims = 1:10 不是預設值——是 B7 決定的 nPC。原則：分群（FindNeighbors）
# 用幾維，視覺化就用幾維，兩邊看的才是同一個空間。
# 其餘全用預設：umap.method = "uwot"、n.neighbors = 30、
# min.dist = 0.3、metric = "cosine"、seed.use = 42。
# 第一次跑會出現「default method ... changed to UWOT」Warning：
# 是資訊不是錯誤，v5 用 R 原生引擎，每個 session 只提醒一次。
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("output/B09_umap_default.png", width = 7, height = 5.5, dpi = 200)

## ---- 2. reproducibility --------------------------------------------
# 可重現性驗證：同參數重跑，座標應逐點一致（功臣：seed.use = 42）。
u1 <- Embeddings(pbmc, "umap")
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
u2 <- Embeddings(pbmc, "umap")
all.equal(u1, u2)
# 預期：TRUE（同機器、同版本、同參數）

# 想看「另一種長相」：只換 seed，存到另一個 reduction，不蓋正式版。
# 佈局會重排，但每一群的成員一顆都不變。
pbmc <- RunUMAP(pbmc, dims = 1:10, seed.use = 123,
                reduction.name = "umap_alt",
                reduction.key  = "UMAPalt_", verbose = FALSE)
DimPlot(pbmc, reduction = "umap") | DimPlot(pbmc, reduction = "umap_alt")

# seed 的三個層級（影片頁 14）：
#   1. 全域：set.seed(1234)——管不進 RunUMAP（它內部自設 seed）
#   2. 函式層：seed.use = 42（預設）——同機同版本重跑必相同
#   3. 實作層：uwot 預設單執行緒 SGD → 同 seed 必同圖；
#      跨版本／跨平台不保證位元一致 → renv.lock + sessionInfo() 記版本。
#      若用 umap.method = "umap-learn" 或開多執行緒，同 seed 也可能漂。

## ---- 3. sweep-neighbors --------------------------------------------
# n.neighbors：拿幾個鄰居定義「局部」。小 → 局部細節；大 → 全域關係。
# 結構性旋鈕（真的改 KNN 圖）。每個結果各存一個 reduction，互不覆蓋。
for (nn in c(5, 30, 100)) {
  pbmc <- RunUMAP(pbmc, dims = 1:10, n.neighbors = nn,
                  reduction.name = paste0("umap_nn", nn),
                  reduction.key  = paste0("UMAPnn", nn, "_"),
                  verbose = FALSE)
}
p.nn <- lapply(c(5, 30, 100), \(nn)
  DimPlot(pbmc, reduction = paste0("umap_nn", nn)) + NoLegend() +
    ggtitle(paste0("n.neighbors = ", nn)))
wrap_plots(p.nn, ncol = 3)
ggsave("output/B09_sweep_nneighbors.png", width = 13, height = 4.5,
       dpi = 200)

## ---- 4. sweep-mindist ----------------------------------------------
# min.dist：二維圖上點跟點最近可以貼多近。純外觀旋鈕（不碰 KNN 圖）。
# 壓小 → 硬塊＋大片留白「分超開」；調大 → 鬆散雲朵。
for (md in c(0.01, 0.3, 0.8)) {
  tag <- gsub("\\.", "", as.character(md))    # "001" "03" "08"
  pbmc <- RunUMAP(pbmc, dims = 1:10, min.dist = md,
                  reduction.name = paste0("umap_md", tag),
                  reduction.key  = paste0("UMAPmd", tag, "_"),
                  verbose = FALSE)
}
p.md <- lapply(c("001", "03", "08"), \(tag)
  DimPlot(pbmc, reduction = paste0("umap_md", tag)) + NoLegend())
wrap_plots(p.md, ncol = 3)
ggsave("output/B09_sweep_mindist.png", width = 13, height = 4.5,
       dpi = 200)

## ---- 5. tsne --------------------------------------------------------
# t-SNE 對照。perplexity（預設 30）之於 t-SNE ＝ n.neighbors 之於 UMAP。
# RunTSNE 也有 seed.use（預設 1）：t-SNE 一樣隨機，一樣要釘 seed。
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne", label = TRUE) |
  DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("output/B09_tsne_vs_umap.png", width = 12, height = 5, dpi = 200)

## ---- 6. pubfig ------------------------------------------------------
# 出版級 DimPlot：色盲友善配色（Okabe-Ito 延伸 10 色）、label + repel、
# NoLegend；匯出時「先定紙面尺寸、再回推字級」，縮排後文字 >= 8 pt。
cb10 <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
          "#56B4E9", "#D55E00", "#F0E442", "#999999",
          "#882255", "#44AA99")

p.pub <- DimPlot(pbmc, reduction = "umap",
                 label = TRUE, repel = TRUE, label.size = 4,
                 cols = cb10, pt.size = 0.3) +
  NoLegend() + ggtitle(NULL) +
  theme(axis.title = element_text(size = 9),
        axis.text  = element_text(size = 8))

ggsave("output/B09_umap_final.pdf", p.pub,
       width = 12, height = 10, units = "cm")            # 投稿用向量檔
ggsave("output/B09_umap_final.png", p.pub,
       width = 12, height = 10, units = "cm", dpi = 600) # 簡報／預覽用
# 圖說模板：「UMAP: Seurat v5.x, RunUMAP(dims = 1:10, seed.use = 42),
# 其餘參數為預設；套件版本見 renv.lock。」——圖跟參數是一組的。

## ---- 7. save --------------------------------------------------------
# 物件裡此時有 umap（正式）、umap_alt、六個掃描 reduction、tsne。
# 正式分析只用 umap；其他留著教學對照用。B10 讀這一份。
saveRDS(pbmc, "output/pbmc_b09_umap.rds")

## ---- 8. session -----------------------------------------------------
sessionInfo()   # 第三層可重現性：把版本釘進紀錄；收進附錄

# =====================================================================
# [踩雷示範] 以下兩段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：拿 UMAP 座標做下游統計（頁 31–32）
## 群間距離的排序會隨 seed 翻轉——建立在 UMAP 距離上的統計是流沙：
# cent <- function(emb) {
#   t(sapply(split(as.data.frame(emb), Idents(pbmc)), colMeans))
# }
# d1 <- dist(cent(Embeddings(pbmc, "umap")))       # seed 42
# d2 <- dist(cent(Embeddings(pbmc, "umap_alt")))   # seed 123
# cor(as.vector(d1), as.vector(d2), method = "spearman")
## 排序相關性通常低得嚇人（本片示意資料實測 rho = -0.04）。
## 一行自保：換個 seed 重跑，你想引用的「距離結論」還在嗎？
## 真要量族群相似度 → 回 PC 空間：dist(cent(Embeddings(pbmc, "pca")))

## [踩雷示範] 雷二：調參數調到群分開為止（頁 33–34）
## min.dist 壓小 + n.neighbors 壓小 + 排斥力加大 →
## 任何資料都能「分得很開」，連續過渡也能被壓出假留白：
# pbmc <- RunUMAP(pbmc, dims = 1:10, n.neighbors = 5,
#                 min.dist = 0.001, repulsion.strength = 2,
#                 reduction.name = "umap_fake",
#                 reduction.key  = "UMAPfake_")
# DimPlot(pbmc, reduction = "umap") |
#   DimPlot(pbmc, reduction = "umap_fake")
## 右圖「涇渭分明」是參數捏的，不是生物學。
## 分不分得開的證據永遠是 marker（FindMarkers）與 PC 空間結構。
