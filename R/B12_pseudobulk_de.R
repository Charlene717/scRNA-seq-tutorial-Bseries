# =====================================================================
# B12_pseudobulk_de.R — B12「差異表達的正確做法：pseudobulk + DESeq2」
#                       隨集腳本（B 系列收官集）
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 12 集
# 資料：ifnb（Kang et al. 2018，GEO GSE96583）
#       承接 B11 整合後物件 output/ifnb_b11_integrated.rds
# 環境：R >= 4.3、Seurat v5.x、dplyr、ggplot2、
#       DESeq2 + apeglm（Bioconductor）
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# 第 2 節會從 GitHub 下載 demuxlet 的 donor 指派檔（兩個小文字檔），
# 需要網路連線；下載後會存到 data/，之後離線可跑。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
set.seed(1234)                       # 全系列固定 seed，結果可重現

stopifnot("找不到 output/ifnb_b11_integrated.rds，請先跑完 R/B11_integration.R" =
            file.exists("output/ifnb_b11_integrated.rds"))
ifnb <- readRDS("output/ifnb_b11_integrated.rds")

# B11 把 RNA assay 依 stim 拆成多個 layer 做整合；
# DE 要用完整的 counts，先把 layer 合回一份。
ifnb <- JoinLayers(ifnb, assay = "RNA")
table(ifnb$stim)                     # CTRL 6548 / STIM 7451
table(ifnb$seurat_annotations)

## ---- 1. naive -----------------------------------------------------
# 天真做法（基線，留著等下對帳）：FindMarkers 直接比條件，
# 預設 wilcoxon，把每顆細胞當一個獨立觀察值。
# 注意：這不是本集推薦的做法——它就是雷本體，跑它是為了看見後果。
Idents(ifnb) <- "seurat_annotations"
mono <- subset(ifnb, idents = "CD14 Mono")   # 鎖定同一型細胞
Idents(mono) <- "stim"

naive <- FindMarkers(mono, ident.1 = "STIM", ident.2 = "CTRL")
head(naive, 3)
# 預期：ISG15 / IFIT3 / IFI6，p_val = 0（小到浮點數存不下）
sum(naive$p_val_adj < 0.05)          # ~3,000+ 個「顯著」基因

## ---- 2. donor -----------------------------------------------------
# ifnb 的 metadata 只有 stim / ctrl 兩個條件，表面上 n = 1 vs 1。
# 但這批資料是 8 位捐贈者混樣上機（multiplexing），
# Kang et al. 用 demuxlet 依遺傳變異把每顆細胞指回捐贈者。
# 我們下載該論文公開的指派結果，把 donor_id 加回 metadata。
dir.create("data/demuxlet", showWarnings = FALSE, recursive = TRUE)
base <- paste0("https://raw.githubusercontent.com/yelabucsf/",
               "demuxlet_paper_code/master/fig3/")
for (f in c("ye1.ctrl.8.10.sm.best", "ye2.stim.8.10.sm.best")) {
  dest <- file.path("data/demuxlet", f)
  if (!file.exists(dest)) download.file(paste0(base, f), dest)
}
ctrl.info <- read.table("data/demuxlet/ye1.ctrl.8.10.sm.best",
                        header = TRUE, stringsAsFactors = FALSE)
stim.info <- read.table("data/demuxlet/ye2.stim.8.10.sm.best",
                        header = TRUE, stringsAsFactors = FALSE)
info <- rbind(ctrl.info, stim.info)

# BEST 欄：SNG-1016 = 單顆細胞、判給捐贈者 1016；DBL- 開頭 = doublet。
# 只留單顆（SNG），並把 barcode 對回 ifnb 的細胞名。
info <- info[grepl("^SNG", info$BEST), c("BARCODE", "BEST")]
info$donor_id <- sub("^SNG-", "", info$BEST)
# 兩個 lane 的 barcode 可能重複；重複的一律丟掉，避免張冠李戴
dup <- info$BARCODE[duplicated(info$BARCODE)]
info <- info[!info$BARCODE %in% dup, ]

# ifnb 細胞名形如 AAACATACATTTCC.1；demuxlet 檔用 AAACATACATTTCC-1
cells <- sub("-1$", "", sub("\\.1$", "", colnames(ifnb)))
idx <- match(cells, sub("-1$", "", info$BARCODE))
ifnb$donor_id <- info$donor_id[idx]

table(ifnb$stim, ifnb$donor_id, useNA = "ifany")
# 預期：8 位 donor（101, 1015, 1016, 1039, 107, 1244, 1256, 1488），
# 每格數百到數千顆；NA = demuxlet 沒把握或 doublet
ifnb <- subset(ifnb, subset = donor_id %in% unique(na.omit(ifnb$donor_id)))

## ---- 3. pseudobulk ------------------------------------------------
# AggregateExpression：依 型別 × 條件 × donor 把 RNA counts 逐基因加總。
# 鐵則：加總的是 RNA assay 的原始 counts——
#       不是 data、不是 scale.data、更不是 integrated 值。
pseudo <- AggregateExpression(ifnb, assays = "RNA",
            group.by = c("seurat_annotations", "stim", "donor_id"))
counts <- pseudo$RNA
dim(counts)                          # 基因 × (型別×條件×donor) 欄
colnames(counts)[1:4]                # 例："CD14 Mono_CTRL_101" ...

# 取出 CD14 Mono 的 16 欄（8 CTRL + 8 STIM），整理樣本表
cd14 <- as.matrix(counts[, grepl("^CD14 Mono_", colnames(counts))])
parts <- strsplit(colnames(cd14), "_")
coldata <- data.frame(
  condition = factor(sapply(parts, `[`, 2), levels = c("CTRL", "STIM")),
  donor     = sapply(parts, `[`, 3),
  row.names = colnames(cd14))
table(coldata$condition)             # CTRL 8 / STIM 8

# 一行檢查（本集的「一行抓雷」）：pseudobulk 矩陣必須是整數 counts。
# 若這裡不是 0，代表你加總到了錯的資料層（雷一）。
sum(cd14 != round(cd14))             # 應為 0

## ---- 4. deseq2 ----------------------------------------------------
# 低表現基因先剔除：全部樣本加起來不到 10 個 counts 的基因
# 沒有檢定價值，留著只會拖累多重檢定校正。
keep <- rowSums(cd14) >= 10
# design 提醒：ifnb 是配對設計（同一批 donor 同時出現在兩個條件），
# 寫 ~ donor + condition 會先吸掉個體差異、功效更好；
# 本集為教學簡潔用 ~ condition。你的實驗若一個 donor 只屬一種條件
# （常見病例對照），donor 與 condition 完全混淆，不能放。
dds <- DESeqDataSetFromMatrix(countData = cd14[keep, ],
                              colData   = coldata,
                              design    = ~ condition)
dds <- DESeq(dds)
# 訊息依序：size factor → 離散度（gene-wise → 趨勢 → 收縮）→ 檢定

res <- results(dds, contrast = c("condition", "STIM", "CTRL"),
               alpha = 0.05)
summary(res)                         # 上調 + 下調 ~2,000 內外

# lfcShrink：把低表現、高變異基因的誇張 log2FC 往 0 收，
# 報告與排序一律用 shrink 後的值。
res.shr <- lfcShrink(dds, coef = "condition_STIM_vs_CTRL",
                     type = "apeglm")
head(res.shr[order(res.shr$padj), ], 6)
# 預期榜首：ISG15 / IFIT1 / IFIT3 / CXCL10 等干擾素刺激基因

## ---- 5. compare ---------------------------------------------------
# 對帳：naive vs pseudobulk 的顯著基因數
n.naive <- sum(naive$p_val_adj < 0.05)
n.pb    <- sum(res.shr$padj < 0.05, na.rm = TRUE)
c(naive = n.naive, pseudobulk = n.pb)   # 差一個數量級

# 逐基因對照 p 值：只在 naive 顯著、pseudobulk 不顯著的那群，
# 多半是 donor 間變異被誤讀成條件效應
common <- intersect(rownames(naive), rownames(res.shr))
audit <- data.frame(
  gene    = common,
  p.naive = naive[common, "p_val_adj"],
  p.pb    = res.shr[common, "padj"])
with(audit, table(naive.sig = p.naive < 0.05,
                  pb.sig    = !is.na(p.pb) & p.pb < 0.05))

## ---- 6. viz -------------------------------------------------------
# 6a. 火山圖（shrink 後的 log2FC；雙閾值：padj 與 |log2FC|）
df <- as.data.frame(res.shr)
df$sig <- !is.na(df$padj) & df$padj < 0.05 &
          abs(df$log2FoldChange) > 1
ggplot(df, aes(log2FoldChange, -log10(padj), colour = sig)) +
  geom_point(size = 0.6) +
  scale_colour_manual(values = c("grey70", "#C0392B")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_classic() + theme(legend.position = "none")
ggsave("output/B12_volcano_cd14.png", width = 7, height = 5, dpi = 200)

# 6b. 發表前最後一關：top 基因的「每樣本一點」圖。
#     真差異要在樣本層級站得住；假訊號在這張圖上原形畢露。
plot_gene <- function(gene) {
  d <- data.frame(expr = log2(t(cd14)[, gene] /
                              colSums(cd14) * 1e6 + 1),
                  condition = coldata$condition)
  ggplot(d, aes(condition, expr, colour = condition)) +
    geom_jitter(width = 0.12, size = 2.4) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4) +
    scale_colour_manual(values = c("#2563A8", "#C0392B")) +
    labs(title = gene, y = "log2 CPM (per sample)") +
    theme_classic() + theme(legend.position = "none")
}
plot_gene("ISG15")
ggsave("output/B12_dots_isg15.png", width = 4.5, height = 4, dpi = 200)

## ---- 7. save ------------------------------------------------------
saveRDS(list(naive      = naive,
             pseudobulk = as.data.frame(res.shr),
             coldata    = coldata,
             audit      = audit),
        "output/ifnb_b12_de_results.rds")

## [踩雷示範] 雷一：在 integrated / 校正值上做 DE（頁 37–38）
## corrected 值連批次帶條件差一起抹、可為負、非整數——
## DESeq2 會直接報錯 "some values in assay are negative"：
# int.mat <- GetAssayData(ifnb, assay = "integrated.cca",
#                         layer = "data")   # 依 B11 的整合名稱調整
# bad <- AggregateExpression(ifnb, assays = "integrated.cca",
#          group.by = c("seurat_annotations", "stim", "donor_id"))
# summary(as.numeric(bad[[1]][, 1]))        # 有負值、非整數
# DESeqDataSetFromMatrix(countData = round(bad[[1]]), ...)  # 錯誤示範
## 鐵則：DE 一律回 RNA assay 的 counts（第 3 節的一行檢查抓得到）。

## [踩雷示範] 雷二：沒有生物重複也硬做（頁 39–40）
## 只依條件聚合 → CTRL 一欄、STIM 一欄，n = 1 vs 1：
# pb1 <- AggregateExpression(ifnb, assays = "RNA",
#          group.by = c("seurat_annotations", "stim"))$RNA
# cd14.1 <- as.matrix(pb1[, grepl("^CD14 Mono_", colnames(pb1))])
# cold1 <- data.frame(condition = factor(c("CTRL", "STIM")),
#                     row.names = colnames(cd14.1))
# dds1 <- DESeqDataSetFromMatrix(cd14.1, cold1, ~ condition)
# dds1 <- DESeq(dds1)
## 新版 DESeq2 直接報錯：組內零自由度，變異無從估計。
## 統計救不了設計：n 要在實驗設計時就決定（每條件 >= 3）。

## [踩雷示範] 雷三：火山圖只看 p 不看效應與表現量（頁 41–42）
## 未 shrink 的結果按 p 排序，低表現高變異基因會混進榜單：
# res.raw <- results(dds, contrast = c("condition", "STIM", "CTRL"))
# cand <- rownames(subset(as.data.frame(res.raw),
#                         baseMean < 5 & abs(log2FoldChange) > 3 &
#                         padj < 0.05))
# head(cand)                        # 低表現卻「效應巨大」的可疑名單
# if (length(cand)) print(plot_gene(cand[1]))  # 每樣本一點：原形畢露
## 對策：報告用 lfcShrink 後的 log2FC，榜單基因逐一畫每樣本點圖。

## ---- 8. sessioninfo ------------------------------------------------
sessionInfo()
