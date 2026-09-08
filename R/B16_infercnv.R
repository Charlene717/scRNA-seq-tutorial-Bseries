# =====================================================================
# B16_infercnv.R — B16「腫瘤場景：inferCNV 與惡性細胞判定」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 16 集
# 資料：infercnv 套件內建的 oligodendroglioma 範例（Smart-seq，
#       下採樣版；經典 chr1p/19q 共缺失事件），不需下載任何外部檔案。
# 環境：R >= 4.3、infercnv（Bioconductor）；HMM 段需要系統先裝 JAGS
#       （Ubuntu: sudo apt install jags；macOS: brew install jags）。
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# 主跑約 3–8 分鐘；HMM 段另需數倍時間，趕時間可跳過 5. hmm。
# 分節標記與影片段落一一對應。
#
# 注意：本集主角不是 Seurat 物件——inferCNV 直接吃三個檔案。
# 你自己的資料要從 Seurat 物件匯出 counts 與 annotations，
# 寫法收在 4.5 節註解。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
# BiocManager::install("infercnv")
library(infercnv)
set.seed(1234)                 # 全系列固定 seed，結果可重現

dir.create("output", showWarnings = FALSE)

## ---- 1. inputs ----------------------------------------------------
# inferCNV 的三件輸入：counts 矩陣、annotations、gene order file。
# 內建範例用 system.file() 取得路徑——這些檔案隨套件安裝，不用下載。
counts_f <- system.file("extdata",
  "oligodendroglioma_expression_downsampled.counts.matrix.gz",
  package = "infercnv")
anno_f <- system.file("extdata",
  "oligodendroglioma_annotations_downsampled.txt",
  package = "infercnv")
gene_f <- system.file("extdata",
  "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
  package = "infercnv")
# gene order 檔名就在提醒你：這是教學用的下採樣版，
# 正式分析要用完整的 GENCODE 座標檔（infercnv wiki 有下載連結）。

# 先看 annotations 長什麼樣：兩欄（cell_id、group），tab 分隔、無表頭
anno <- read.table(anno_f, sep = "\t",
                   col.names = c("cell", "group"))
table(anno$group)
# reference 兩組：Microglia/Macrophage、Oligodendrocytes (non-malignant)
# 其餘 malignant_* 各組來自不同病人的腫瘤細胞（觀察組）

## ---- 2. object ----------------------------------------------------
# ref_group_names 必須「一字不差」對上 annotation 的組名。
# 沒列進 reference 的組，全部自動當觀察組。
cnv0 <- CreateInfercnvObject(
  raw_counts_matrix = counts_f,
  annotations_file  = anno_f,
  gene_order_file   = gene_f,
  delim             = "\t",
  ref_group_names   = c("Microglia/Macrophage",
                        "Oligodendrocytes (non-malignant)"))

## ---- 3. run -------------------------------------------------------
# 第一輪：連續模式（HMM = FALSE），先看訊號長相。
#   cutoff：平均計數低於此值的基因剔除。Smart-seq 用 1；10x 用 0.1。
#   window_length 預設 101 個基因——本集內文解釋為什麼是這個量級。
#   cluster_by_groups：觀察組各自聚類（保留樣本分組看 subclone）。
#   denoise：把 reference 波動範圍內的值壓回白色。
cnv <- infercnv::run(cnv0,
  cutoff            = 1,
  out_dir           = "output/b16_infercnv",
  cluster_by_groups = TRUE,
  denoise           = TRUE,
  HMM               = FALSE)

## ---- 4. outputs ---------------------------------------------------
# out_dir 裡的重點檔案：
#   infercnv.png                     最終 heatmap（先看這張）
#   infercnv.preliminary.png         denoise 前的版本（對照用）
#   infercnv.references.txt          reference 的平滑矩陣
#   infercnv.observations.txt        觀察組的平滑矩陣
#   run.final.infercnv_obj           完整結果物件（下一節讀回）
list.files("output/b16_infercnv")

## ---- 4.5 scores ---------------------------------------------------
# 從結果物件計算 per-cell CNV 分數，並做三角驗證式的粗判定。
obj  <- readRDS("output/b16_infercnv/run.final.infercnv_obj")
expr <- obj@expr.data                    # 基因 × 細胞，中性值在 1 附近
dim(expr)

cnv_score <- colMeans((expr - 1)^2)      # 每顆細胞偏離基線的程度
ref_cells <- unlist(obj@reference_grouped_cell_indices)
summary(cnv_score[ref_cells])            # reference 的分數該貼近 0
summary(cnv_score[-ref_cells])           # 觀察組整體應該高一截

# 與「惡性平均 profile」的相關：抓的是「事件的位置對不對得上」
mal_mean <- rowMeans(expr[, -ref_cells])
corr_mal <- apply(expr, 2, cor, y = mal_mean)

# 粗判定：CNV 分數 × profile 相關。閾值以 reference 的分布為準
# （95 分位當上界），兩路證據一致才下標籤，不一致標 ambiguous。
# 第三路證據（分群、marker）在你自己的 Seurat 流程裡補上。
s_hi <- quantile(cnv_score[ref_cells], 0.95)
lab  <- ifelse(cnv_score > s_hi & corr_mal > 0.4, "malignant",
        ifelse(cnv_score <= s_hi & corr_mal < 0.2, "normal",
               "ambiguous"))
table(lab, ifelse(seq_along(lab) %in% ref_cells, "ref", "obs"))

# —— 用在你自己的資料時（Seurat 物件 pbmc）——
# counts 與 annotations 從物件匯出，gene order 用完整 GENCODE 檔：
#   mat  <- GetAssayData(pbmc, layer = "counts")
#   anno <- data.frame(cell  = colnames(pbmc),
#                      group = as.character(Idents(pbmc)))
#   write.table(anno, "output/b16_anno.txt", sep = "\t",
#               quote = FALSE, row.names = FALSE, col.names = FALSE)
# 跑完 inferCNV 後把標籤掛回 metadata：
#   pbmc <- AddMetaData(pbmc, metadata = lab[colnames(pbmc)],
#                       col.name = "malignancy")

## ---- 5. hmm -------------------------------------------------------
# 第二輪：HMM 模式，把連續訊號翻譯成離散狀態。
#   HMM_type = "i3"：三狀態（缺失／中性／增加）
#   HMM_type = "i6"：六狀態（含程度，輸出整數拷貝數等級）
# 注意：要用「還沒跑過的原始物件」cnv0 重跑；需要 JAGS。
# 跑很久（數十分鐘量級），預設註解掉，需要時再打開。
# cnv_hmm <- infercnv::run(cnv0,
#   cutoff            = 1,
#   out_dir           = "output/b16_infercnv_hmm",
#   cluster_by_groups = TRUE,
#   denoise           = TRUE,
#   HMM               = TRUE,
#   HMM_type          = "i6")
# HMM 輸出（out_dir 裡）：
#   17_HMM_pred*.png / .txt          每區段的離散狀態
#   map_metadata_from_infercnv.txt   可直接併回 metadata 的表

## ---- 6. traps -----------------------------------------------------
## [踩雷示範] 雷一：reference 選錯——把某個惡性組當基線。
## 後果：該組被「歸零」，正常組反而長出假「增加」，方向整個反轉。
## 跑完把 output/b16_trap_wrongref/infercnv.png 和正確版並排看。
# cnv_wrong <- CreateInfercnvObject(
#   raw_counts_matrix = counts_f,
#   annotations_file  = anno_f,
#   gene_order_file   = gene_f,
#   delim             = "\t",
#   ref_group_names   = c("malignant_MGH36"))   # 錯：拿腫瘤組當正常
# cnv_wrong <- infercnv::run(cnv_wrong,
#   cutoff = 1, out_dir = "output/b16_trap_wrongref",
#   cluster_by_groups = TRUE, denoise = TRUE, HMM = FALSE)

## [踩雷示範] 雷二的「一行檢查」：疑似事件落在高危險基因簇時，
## 先看那一段是不是 HLA（chr6p21）或 IG 基因座（chr2/14/22）。
## 把該區段基因抓出來看名字，一眼就能識破：
# gene_pos <- read.table(gene_f,
#   col.names = c("gene", "chr", "start", "stop"))
# head(subset(gene_pos, chr == "chr6" &
#             start > 29e6 & stop < 33e6))     # HLA 區段的基因名單
## 名單裡一排 HLA-*：那段「擴增」多半是表達，不是拷貝數。

## 雷三沒有程式碼可跑——它是推論方向的錯誤：
## CNV 乾淨 ≠ 不是癌（近二倍體血癌、部分兒科腫瘤）。
## 缺 CNV 訊號時，回到突變、融合基因、marker 與臨床資訊。

## ---- 7. save ------------------------------------------------------
saveRDS(obj, "output/b16_infercnv_obj.rds")
sessionInfo()
