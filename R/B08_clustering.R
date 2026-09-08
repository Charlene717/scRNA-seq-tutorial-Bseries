# =====================================================================
# B08_clustering.R — B8「分群：resolution 掃描與 clustree」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 8 集
# 資料：PBMC 3k（承接 B5 QC 後物件 output/pbmc_b05_qc.rds，
#       其中已含 NormalizeData / HVG / ScaleData / PCA / UMAP 與
#       DoubletFinder 標籤——見 R/B05_qc.R）
# 環境：R >= 4.3、Seurat v5.x、clustree、dplyr、ggplot2
#
# 使用方式：從專案根目錄 source() 或逐段執行（全程約 2–3 分鐘）。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
set.seed(1234)                       # 全系列固定 seed，結果可重現

stopifnot("找不到 output/pbmc_b05_qc.rds，請先跑完 R/B05_qc.R" =
            file.exists("output/pbmc_b05_qc.rds"))
pbmc <- readRDS("output/pbmc_b05_qc.rds")   # 承接 B5–B7 的物件
pbmc
# 預期：13714 features across ~2638 samples（含 pca 與 umap reductions）

## ---- 1. neighbors -------------------------------------------------
# 在 PC 空間（不是基因空間）建 KNN 圖，再加權成 SNN 圖。
# dims = 1:10 承接 B7 的 nPC 決定；k.param 預設 20，沒特殊理由不動。
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# 物件裡多了兩張圖：RNA_nn（KNN）與 RNA_snn（SNN，分群用的那張）
names(pbmc@graphs)

## ---- 2. clusters --------------------------------------------------
# 在 SNN 圖上跑社群偵測。algorithm = 1（Louvain）是預設；
# Leiden 用 algorithm = 4（需 Python 套件 leidenalg）。
pbmc <- FindClusters(pbmc, resolution = 0.5)
# 預期 console：Number of communities: 9（Louvain、10 random starts）

# 編號 0 起跳、按群大小排序——編號是任意的（雷二會收回來用）
table(Idents(pbmc))

# UMAP 是 B9 的主角，今天先借 B5 算好的座標看結果
DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("output/B08_umap_res0.5.png", width = 7, height = 5.5, dpi = 200)

## ---- 3. sweep -----------------------------------------------------
# 鐵則：不挑一個值，先看整條光譜。每個 resolution 一個 metadata
# 欄位（RNA_snn_res.X），互不覆蓋，全部留著給 clustree。
for (r in seq(0.2, 1.4, by = 0.2)) {
  pbmc <- FindClusters(pbmc, resolution = r, verbose = FALSE)
}

# 掃出來的欄位（注意 1.0 的欄名是 "RNA_snn_res.1"）
res.cols <- grep("^RNA_snn_res", colnames(pbmc@meta.data),
                 value = TRUE)
res.cols

# 每個 resolution 各分幾群：同一份資料，6 到 16 群都做得出來
sapply(pbmc@meta.data[res.cols], \(x) length(unique(x)))

# 注意：掃完之後 active identity 停在「最後跑的」1.4，
# 後面每一步都要先確認 Idents 指向哪個欄位。

## ---- 4. clustree --------------------------------------------------
# install.packages("clustree")   # CRAN 版即可
library(clustree)

# prefix 對準欄位名的共同開頭：每個 resolution 疊成一層，
# 節點 = 群（大小 = 細胞數），箭頭 = 細胞在層間的流向。
p.tree <- clustree(pbmc, prefix = "RNA_snn_res.")
p.tree
ggsave("output/B08_clustree.png", plot = p.tree,
       width = 9, height = 7, dpi = 200)

# 判讀三問（影片頁 21–22）：
#   1) 這一層穩定嗎（連續幾個 resolution 結構不變 → 定案候選）
#   2) 分裂是單向的嗎（乾淨一裂為二、無回流 → 拿 marker 去驗）
#   3) 有沒有交換（箭頭交叉、多個上游 → 亂流層，退回去）

## ---- 5. markers ---------------------------------------------------
# 5a. 先用「交叉表」找出相鄰 resolution 之間被切開的是誰。
#     跨標籤對照永遠用 table()，不要用編號硬對（雷二）。
table(pbmc$RNA_snn_res.0.4, pbmc$RNA_snn_res.0.6)
# 預期：0.4 的 1 號群在 0.6 被切成 1 號（~436 顆）與 6 號（~158 顆）

# 5b. 那一刀的兩邊差在哪些基因？
Idents(pbmc) <- "RNA_snn_res.0.6"
split.mk <- FindMarkers(pbmc, ident.1 = 1, ident.2 = 6,
                        min.pct = 0.25)
head(split.mk, 10)
# 預期：S100A8 / S100A9 / LGALS2 / CD14 在 1 號側衝高；
# 反方向（tail 或 ident 對調）是 FCGR3A / MS4A7 ——
# 教科書級的古典 vs 非古典單核球。
# 注意：同一份資料先分群再檢定是 double dipping，
# p 值只當排序用；判斷靠 log2FC、pct 差與已知生物學。

# 5c. 拿到圖上對質
FeaturePlot(pbmc, features = c("CD14", "FCGR3A"))
ggsave("output/B08_featureplot_cut.png", width = 10, height = 4.5,
       dpi = 200)
VlnPlot(pbmc, features = c("CD14", "LYZ", "FCGR3A", "MS4A7"),
        idents = c(1, 6), ncol = 4, pt.size = 0)
ggsave("output/B08_vlnplot_cut.png", width = 12, height = 4, dpi = 200)

# 5d. 對照組：0.6 → 0.8 多出來的那一刀（實跑會看到 top 基因
#     全是核糖體蛋白、log2FC < 0.5、pct 兩邊都 ~99% → 不切）。
table(pbmc$RNA_snn_res.0.6, pbmc$RNA_snn_res.0.8)

## ---- 6. percluster-qc ---------------------------------------------
# 回收 B5 的兩個鉤子：分群後每群 QC 檢查 + 每群 doublet 比例。

# 6a. 每群 QC 中位數：要抓「mt 一枝獨秀 + nFeature 墊底」的殘渣群
qc.tab <- pbmc@meta.data |>
  group_by(cluster = RNA_snn_res.0.6) |>
  summarise(n  = n(),
            mt = median(percent.mt),
            nF = median(nFeature_RNA))
arrange(qc.tab, desc(mt))
# 這份資料：沒有殘渣群（B5 的過濾線畫對了）。
# 若有：不是刪掉該群，是回 B5 修過濾線、從 QC 重跑（雷三）。

# 視覺化版本：一眼看所有群的分布
VlnPlot(pbmc, features = c("percent.mt", "nFeature_RNA"),
        group.by = "RNA_snn_res.0.6", pt.size = 0, ncol = 2)
ggsave("output/B08_percluster_qc.png", width = 12, height = 4.5,
       dpi = 200)

# 6b. 每群 doublet 比例：B5 留在 metadata 的 DoubletFinder 標籤
df.col <- grep("^DF.classifications",
               colnames(pbmc@meta.data), value = TRUE)
round(prop.table(table(pbmc$RNA_snn_res.0.6,
                       pbmc@meta.data[[df.col]]), 1), 3)
# 判讀：各群都貼著整體預期率（~2.3%）→ 過關。
# 要抓的是「某群過半被標 doublet、又卡在兩大群之間」——
# 那不是群裡有 doublet，是整群就是 doublet，整群移除並記錄。

## ---- 7. finalize --------------------------------------------------
# 定案：clustree 穩定帶 + marker 證據 → resolution 0.6。
# 方法段三要素：參數、掃描範圍、選擇理由（影片頁 34）。
Idents(pbmc) <- "RNA_snn_res.0.6"
pbmc$seurat_clusters <- pbmc$RNA_snn_res.0.6   # 對齊慣用欄位

# 理由記錄（進分析筆記與論文方法段）：
#   FindNeighbors(dims = 1:10, k.param = 20) + Louvain(algorithm = 1)
#   resolution 0.2–1.4 掃描；0.6 起結構穩定，0.4→0.6 那一刀
#   有 CD14/LYZ vs FCGR3A/MS4A7 支持；0.6→0.8 之後無 marker 證據。
#   每群 QC 與 doublet 比例檢查通過（無殘渣群、無 doublet 群）。

saveRDS(pbmc, "output/pbmc_b08_clustered.rds")   # B9/B10 直接讀這份

## ---- 8. session ---------------------------------------------------
sessionInfo()                        # 收進附錄；審稿人的好朋友

# =====================================================================
# [踩雷示範] 以下三段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：調 resolution 直到出現「想要的」群數（頁 37–38）
## 同一份資料，4 / 8 / 16 群任你點菜——群數本身不是證據：
# for (r in c(0.05, 0.4, 1.6)) {
#   pbmc <- FindClusters(pbmc, resolution = r, verbose = FALSE)
#   cat(sprintf("res %.2f -> %d clusters\n", r,
#               length(unique(Idents(pbmc)))))
# }
## 用結論挑參數、再用參數產生結論 = 循環論證。

## [踩雷示範] 雷二：跨 run 用編號硬比（頁 39–40）
## 換 dims 重跑，「cluster 3」很可能已是另一群細胞：
# pbmc.alt <- FindNeighbors(pbmc, dims = 1:20)
# pbmc.alt <- FindClusters(pbmc.alt, resolution = 0.6)
# table(old = pbmc$RNA_snn_res.0.6, new = Idents(pbmc.alt))
## 對角線亂掉 = 編號洗牌。跨 run 對照只能用交叉表，不能看編號。

## [踩雷示範] 雷三：把殘渣群當新細胞型別（頁 41）
## 殘渣群的三個特徵：高 mt、低 nFeature、marker 全是
## 粒線體/熱休克基因（MT-*, HSPA1A, ...）。檢查法：
# Idents(pbmc) <- "RNA_snn_res.0.6"
# junk.mk <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25)
# junk.mk |> group_by(cluster) |> slice_head(n = 5) |>
#   select(cluster, gene, avg_log2FC)
## 某群 top 基因全是 MT-/HSP 開頭 → 回 B5 修過濾線重跑，
## 不是把它寫成新發現，也不是默默刪掉。
