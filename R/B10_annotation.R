# =====================================================================
# B10_annotation.R — B10「細胞註釋：marker、FindAllMarkers、
#                    SingleR / Azimuth」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 10 集
# 資料：PBMC 3k（承接 B8 分群後物件 output/pbmc_b08_clustered.rds，
#       其中已含 res 0.2–1.4 掃描結果、定案 res 0.6 十群、UMAP）
# 環境：R >= 4.3、Seurat v5.x、dplyr、ggplot2、
#       SingleR + celldex（Bioconductor）、Azimuth（GitHub）
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# SingleR / Azimuth 首次執行會下載參考資料（共約 1 GB），
# 全程約 10–15 分鐘；之後有快取會快很多。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
set.seed(1234)                       # 全系列固定 seed，結果可重現

stopifnot("找不到 output/pbmc_b08_clustered.rds，請先跑完 R/B08_clustering.R" =
            file.exists("output/pbmc_b08_clustered.rds"))
pbmc <- readRDS("output/pbmc_b08_clustered.rds")
# 亦可改讀 B9 的 output/pbmc_b09_umap.rds：分群欄位完全相同，
# 只多了 B9 參數實驗的幾個 UMAP reduction。本集以 B8 定案物件為正本。

# 鐵則：不靠默契，明確指定用哪一組分群標籤（B8 定案 res 0.6，十群）
Idents(pbmc) <- "RNA_snn_res.0.6"
table(Idents(pbmc))

## ---- 1. markers ---------------------------------------------------
# 每群輪流「一群 vs 其餘全部」找差異基因（= 傳喚證人）。
#   only.pos = TRUE       只要正向 marker（註釋要的是「有什麼」）
#   min.pct  = 0.25       至少一邊 25% 細胞表現才檢定（濾雜訊）
#   logfc.threshold = 0.25 倍數差太小的直接跳過（提速）
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
dim(pbmc.markers)

# 逐欄判讀（影片頁 11）：看表順序 = avg_log2FC → pct.1/pct.2 → p_val_adj
# p_val_adj 只當把關，不當排序——先分群再檢定是 double dipping，
# p 值天生偏小。
pbmc.markers |> filter(cluster == 3) |> slice_head(n = 4)
# 預期：CD79A / MS4A1 / TCL1A / CD79B——B 細胞受體零件，3 號 = B

# 每群 top5，排序鍵是 log2FC 不是 p 值（雷一）
top5 <- pbmc.markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 5) |>
  ungroup()
top5 |> filter(cluster == 9) |>
  select(gene, avg_log2FC, pct.1, pct.2)
# 預期：PPBP / PF4 / GNG11 / ...，pct.2 ~ 0 → 9 號 = 血小板

write.csv(pbmc.markers, "output/B10_findallmarkers.csv",
          row.names = FALSE)         # 證人名單存檔，審稿補充資料用

## ---- 2. viz -------------------------------------------------------
# 三件套：FeaturePlot 看「亮在哪」、VlnPlot 看「整群高還是少數高」、
# DotPlot 看「多基因 × 多群總覽」（呼應 A5 的判讀課）。
FeaturePlot(pbmc, features = "MS4A1")
ggsave("output/B10_featureplot_ms4a1.png", width = 6, height = 5,
       dpi = 200)

VlnPlot(pbmc, features = "MS4A1", pt.size = 0)
ggsave("output/B10_vlnplot_ms4a1.png", width = 8, height = 4, dpi = 200)

DotPlot(pbmc, features = c("MS4A1", "CD79A", "CD3E", "NKG7")) +
  RotatedAxis()
ggsave("output/B10_dotplot_mini.png", width = 7, height = 4.5, dpi = 200)

## ---- 3. manual ----------------------------------------------------
# 3a. PBMC 經典 marker 面板（影片頁 19 的表）。
#     注意：CD4 mRNA dropout 嚴重，CD4 T 靠 IL7R 佐證。
panel <- c("CD3E", "IL7R", "CCR7", "S100A4",   # T 與 CD4 亞群
           "CD8A",                             # CD8 T
           "MS4A1",                            # B
           "CD14", "LYZ",                      # 古典單核球
           "FCGR3A", "MS4A7",                  # 非古典單核球
           "NKG7", "GNLY",                     # NK 與胞殺程式
           "FCER1A", "CST3",                   # DC
           "PPBP")                             # 血小板

# 3b. 四大譜系門牌先劃大陸塊，再處理細節
FeaturePlot(pbmc, features = c("CD3E", "MS4A1", "CD14", "NKG7"))
ggsave("output/B10_featureplot_lineage.png", width = 10, height = 8,
       dpi = 200)

# 3c. 全景 DotPlot：找「沿對角線一塊塊亮」；
#     7 號同時亮 T 底（CD3E/CD8A）與胞殺（NKG7/GNLY）→ 記在案上
DotPlot(pbmc, features = panel) + RotatedAxis()
ggsave("output/B10_dotplot_panel.png", width = 11, height = 5, dpi = 200)

# 3d. 掛名字。RenameIdents 只改 active identity，
#     所以立刻存進 metadata（celltype.manual）；cluster 編號永遠保留。
new.ids <- c("Naive CD4 T", "CD14 Mono", "Memory CD4 T",
             "B", "CD8 T", "NK", "FCGR3A Mono",
             "Effector CD8 T", "DC", "Platelet")
names(new.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.ids)
pbmc$celltype.manual <- Idents(pbmc)
table(pbmc$celltype.manual)

DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()
ggsave("output/B10_umap_manual.png", width = 7, height = 5.5, dpi = 200)

## ---- 4. singler ---------------------------------------------------
# BiocManager::install(c("SingleR", "celldex"))
library(SingleR)
library(celldex)

# 參考集選擇是一個「決定」：涵蓋要廣、組織要對口。
# 第一輪粗篩用 HPCA（37 種主標籤）。
ref.hpca <- HumanPrimaryCellAtlasData()

# clusters = 逐「群」註釋：快、穩、輸出好判讀。
# Seurat v5 取 log-normalized 值用 layer = "data"。
pred <- SingleR(test = GetAssayData(pbmc, layer = "data"),
                ref = ref.hpca, labels = ref.hpca$label.main,
                clusters = pbmc$RNA_snn_res.0.6)
pred$labels
# 預期（群 0–9）：T_cells Monocyte T_cells B_cell T_cells
#                 NK_cell Monocyte NK_cell DC Platelets
# → 7 號被標 NK_cell：cold open 的來源

# score 熱圖：看的不是「誰最高」，是「甩開第二名多少」
plotScoreHeatmap(pred)
ggsave("output/B10_singler_heatmap.png", width = 8, height = 5,
       dpi = 200)

# delta = 最高分 − 次高分：低 = 工具在猜（影片頁 27）
data.frame(cluster = rownames(pred),
           label   = pred$labels,
           delta   = round(pred$delta.next, 3))
# 預期：7 號 delta ~0.02，全場最低，比其他群低一個數量級

# 一行檢查（本集的「一行抓雷」）：把 delta 攤開，低分群送回 marker 複審
plot(pred$delta.next, type = "h",
     xlab = "cluster", ylab = "delta.next")

# 群層級標籤攤回每顆細胞，之後交叉驗證表用
pbmc$celltype.singler <-
  pred$labels[as.integer(as.character(pbmc$RNA_snn_res.0.6)) + 1]

## ---- 5. azimuth ---------------------------------------------------
# remotes::install_github("satijalab/azimuth")
# 首次執行會下載 pbmcref 參考地圖（數百 MB）。
library(Azimuth)
pbmc <- RunAzimuth(pbmc, reference = "pbmcref")

# l1 / l2 / l3 由粗到細；prediction score 是「參考細胞投票一致度」，
# 不是「正確率」——參考集沒有的型別照樣有分數（雷二）。
head(pbmc@meta.data[, c("predicted.celltype.l2",
                        "predicted.celltype.l2.score")], 3)

# 每群的 l2 多數標籤（進交叉驗證表的那一欄）
azi.tab <- pbmc@meta.data |>
  group_by(cluster = RNA_snn_res.0.6) |>
  count(predicted.celltype.l2) |>
  slice_max(n, n = 1)
azi.tab
# 預期：0 CD4 Naive / 1 CD14 Mono / 2 CD4 TCM / 3 B naive /
#       4 CD8 Naive / 5 NK / 6 CD16 Mono / 7 CD8 TEM / 8 cDC2 /
#       9 Platelet

## ---- 6. crossval --------------------------------------------------
# 6a. 三方交叉驗證表（影片頁 32）。判定原則：
#     粗細之差（T_cells vs CD4 Naive）不算衝突；方向相反才開庭。
crossval <- pbmc@meta.data |>
  group_by(cluster = RNA_snn_res.0.6) |>
  summarise(n       = n(),
            manual  = first(as.character(celltype.manual)),
            singler = first(celltype.singler),
            azimuth = names(sort(table(predicted.celltype.l2),
                                 decreasing = TRUE))[1])
crossval
write.csv(crossval, "output/B10_crossval_table.csv", row.names = FALSE)

# 6b. 唯一吵架的 7 號：手動=效應 CD8 T、SingleR=NK、Azimuth=CD8 TEM。
#     開庭：跟公認的 NK（5 號）當面對質。
Idents(pbmc) <- "RNA_snn_res.0.6"
vs.nk <- FindMarkers(pbmc, ident.1 = 7, ident.2 = 5, min.pct = 0.25)
vs.nk[c("CD3E", "CD8A", "GNLY", "NKG7"),
      c("avg_log2FC", "pct.1", "pct.2")]
# 判讀：CD3E/CD8A 在 7 號一側衝高（TCR 零件，NK 沒有）→ T 細胞；
# GNLY/NKG7 兩邊都高（共享胞殺程式）→ SingleR 被整體相似度騙了。
# 裁決：canonical marker（一級證據）壓過相似度（二級）→ 效應 CD8 T。

## ---- 7. save ------------------------------------------------------
# 最終標籤 = 手動 + 裁決結果；存明確欄位，cluster 編號原樣保留
pbmc$celltype.final <- pbmc$celltype.manual
Idents(pbmc) <- "celltype.final"

DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()
ggsave("output/B10_umap_final.png", width = 7, height = 5.5, dpi = 200)

saveRDS(pbmc, "output/pbmc_b10_annotated.rds")   # B11 之後接續用

sessionInfo()                        # 收進附錄；審稿人的好朋友

# =====================================================================
# [踩雷示範] 以下三段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：按 p 值挑 marker（頁 37–38）
## 兩千多顆細胞讓持家基因的 p 值也小到爆表；
## 按 p 排序你會選到 B2M，按 FC 與 pct 差你會選到 NKG7：
# nk.mk <- FindMarkers(pbmc, ident.1 = 5,
#                      min.pct = 0.1, logfc.threshold = 0)
# nk.mk[c("B2M", "NKG7"),
#       c("p_val_adj", "avg_log2FC", "pct.1", "pct.2")]
## B2M：pct 1.00 vs 0.99、FC ~0.6 ——人人都有，不是 marker。

## [踩雷示範] 雷二：參考集沒有的型別被硬塞最近標籤（頁 39–40）
## DICE 只有 T/B/NK/單核球——沒有 DC、沒有血小板：
# ref.dice <- DatabaseImmuneCellExpressionData()
# pred.d <- SingleR(test = GetAssayData(pbmc, layer = "data"),
#                   ref = ref.dice, labels = ref.dice$label.main,
#                   clusters = pbmc$RNA_snn_res.0.6)
# pred.d$labels[9:10]              # DC 與血小板 → 都變 "Monocytes"
# round(pred.d$delta.next[9:10], 3) # delta 趴在地上：唯一的警報
## 鐵則：標籤照收、delta 必查；低分群回 canonical marker 複審。

## [踩雷示範] 雷三：subcluster 之後還用大類 marker（頁 41–42）
## 把 T 細胞單獨抓出來再分群，然後拿 CD3E 去「區分」亞群：
# tcell <- subset(pbmc, subset = RNA_snn_res.0.6 %in% c(0, 2, 4, 7))
# tcell <- FindNeighbors(tcell, dims = 1:10) |>
#          FindClusters(resolution = 0.4)
# VlnPlot(tcell, features = "CD3E", pt.size = 0)   # 全亮，資訊量 = 0
# VlnPlot(tcell, features = c("CCR7", "S100A4", "CD8A", "GZMB"),
#         pt.size = 0, ncol = 4)                   # 亞群面板才分得開
## 原則：marker 解析度要配得上群的解析度，切一層、換一層面板。
