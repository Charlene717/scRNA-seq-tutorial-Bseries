# =====================================================================
# B13_enrichment.R — B13「DE 之後：功能富集與基因集分析」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 13 集
# 資料：ifnb（Kang et al. 2018，GEO GSE96583）
#       承接 B12 的 DE 結果 output/ifnb_b12_de_results.rds
#       與 B11 的整合物件 output/ifnb_b11_integrated.rds
# 環境：R >= 4.3、Seurat v5.x、dplyr、ggplot2、
#       msigdbr（CRAN）、fgsea + clusterProfiler + org.Hs.eg.db
#       （Bioconductor）；AUCell 為選配（第 6 節）
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# 全程離線可跑（msigdbr 的基因集內建在套件裡，不需下載）。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(1234)                       # 全系列固定 seed，結果可重現

stopifnot("找不到 output/ifnb_b12_de_results.rds，請先跑完 R/B12_pseudobulk_de.R" =
            file.exists("output/ifnb_b12_de_results.rds"))
de <- readRDS("output/ifnb_b12_de_results.rds")
names(de)                            # naive / pseudobulk / coldata / audit
pb <- de$pseudobulk                  # DESeq2 lfcShrink 後的結果表
head(pb[order(pb$padj), ], 3)        # 榜首應為 ISG15 / IFIT 家族

## ---- 1. ranks -----------------------------------------------------
# 從 DE 表建 ranked list：全部檢定過的基因、由「方向 × 顯著程度」排序。
# 兩派排序指標：
#   (a) DESeq2 的 Wald stat（results() 的 stat 欄）——B12 存的是
#       lfcShrink 後的表，沒有 stat 欄，所以本集用 (b)；
#   (b) sign(log2FC) × -log10(pvalue)，與 (a) 高度一致。
# 注意：排序用「未 shrink 的統計量或 p 值」；shrink 後的 log2FC
#       是拿來報告效應量的，不是拿來排序的。
pb <- pb[!is.na(pb$pvalue) & !is.na(pb$log2FoldChange), ]
pb$metric <- sign(pb$log2FoldChange) * -log10(pmax(pb$pvalue, 1e-300))
ranks <- sort(setNames(pb$metric, rownames(pb)), decreasing = TRUE)

# 交給 fgsea 前的三個檢查：長度、兩端、乾淨度
length(ranks)                        # ~11,000：檢定過的全部基因
head(ranks, 3); tail(ranks, 3)       # 頭該是 ISG15 等 IFN 基因
stopifnot(!anyNA(ranks), !anyDuplicated(names(ranks)))

## ---- 2. hallmark --------------------------------------------------
# MSigDB Hallmark：50 個精選、去冗餘的基因集，富集分析的起手式。
# msigdbr 把基因集打包成 data frame；split 成 fgsea 要的 list。
hm <- msigdbr(species = "Homo sapiens", collection = "H")
# 舊版 msigdbr（< 10.0）參數名是 category = "H"，內容相同
hallmark <- split(hm$gene_symbol, hm$gs_name)
length(hallmark)                     # 50
sapply(hallmark[1:3], length)

## ---- 3. fgsea -----------------------------------------------------
# fgsea：fast GSEA（multilevel 演算法）。
# minSize / maxSize：太小的集估不穩、太大的集什麼都沾一點；
# eps = 0：p 值算到底，不設 1e-10 的下限（要準確 padj 就設 0）。
set.seed(1234)
gsea <- fgsea(pathways = hallmark,
              stats    = ranks,
              minSize  = 15,
              maxSize  = 500,
              eps      = 0)

# 已知答案驗證法：IFN-β 刺激的單核球，
# INTERFERON_ALPHA / GAMMA_RESPONSE 必須是榜首，否則先懷疑自己。
gsea %>% arrange(padj) %>%
  select(pathway, NES, padj, size) %>% head(8)
gsea %>% arrange(NES) %>%
  select(pathway, NES, padj) %>% head(3)   # 下調端（負 NES）

## ---- 4. plots -----------------------------------------------------
# 4a. enrichment plot：報告 GSEA 結果的標準圖
plotEnrichment(hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE, ranks) +
  labs(title = "HALLMARK_INTERFERON_ALPHA_RESPONSE")
ggsave("output/B13_gsea_ifna.png", width = 6, height = 4, dpi = 200)

# 4b. leading edge：峰值之前的成員基因——真正扛起訊號的那批
ia <- which(gsea$pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE")
le.ifna <- gsea$leadingEdge[[ia]]
head(le.ifna, 10)                    # ISG15 / IFIT1 / MX1 / OAS1 ...
length(le.ifna)

# 4c. NES 總覽長條圖（顯著通路，按 NES 排）
gsea %>% filter(padj < 0.05) %>% arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway)) %>%
  ggplot(aes(NES, pathway, fill = NES > 0)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("#2563A8", "#C0392B")) +
  labs(y = NULL) + theme_classic()
ggsave("output/B13_gsea_nes.png", width = 7, height = 5, dpi = 200)

## ---- 5. ora -------------------------------------------------------
# ORA 路線：清單 + 背景。
# 清單：上調且效應夠大（padj < 0.05 且 log2FC > 1）。
# 背景（universe）：這次「真的檢定過」的基因——不是全基因組！
sig.up   <- rownames(pb)[!is.na(pb$padj) &
                         pb$padj < 0.05 & pb$log2FoldChange > 1]
universe <- rownames(pb)
length(sig.up); length(universe)

ego <- enrichGO(gene     = sig.up,
                universe = universe,      # 關鍵參數，沒有之一
                OrgDb    = org.Hs.eg.db,
                keyType  = "SYMBOL",
                ont      = "BP",
                pAdjustMethod = "BH")
nrow(as.data.frame(ego))
# simplify()：依 GO 語義相似度合併高度重疊的父子／兄弟條目
ego <- simplify(ego)
head(as.data.frame(ego)[, c("Description", "GeneRatio",
                            "p.adjust")], 6)
# 預期主題：defense response to virus / response to type I IFN 等

## ---- 6. modulescore -----------------------------------------------
# 細胞層級：把 Hallmark IFN-α response 的分數畫回 UMAP。
stopifnot("找不到 output/ifnb_b11_integrated.rds，請先跑完 R/B11_integration.R" =
            file.exists("output/ifnb_b11_integrated.rds"))
ifnb <- readRDS("output/ifnb_b11_integrated.rds")
ifnb <- JoinLayers(ifnb, assay = "RNA")

ifn.genes <- hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE
sum(ifn.genes %in% rownames(ifnb))   # 97 個成員，資料裡找得到幾個
ifnb <- AddModuleScore(ifnb, features = list(ifn.genes),
                       name = "IFNa_score")
# 預設 nbin = 24（表現量分箱）、ctrl = 100（每 bin 抽的對照基因數）
# 分數欄自動加編號：IFNa_score1

FeaturePlot(ifnb, "IFNa_score1", split.by = "stim")
ggsave("output/B13_score_umap.png", width = 9, height = 4, dpi = 200)

# 每型別 × 條件的分數摘要：CD14 Mono 的 STIM-CTRL 落差應最醒目
ifnb@meta.data %>%
  group_by(seurat_annotations, stim) %>%
  summarise(score = median(IFNa_score1), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = stim, values_from = score) %>%
  mutate(diff = STIM - CTRL) %>% arrange(desc(diff))

# [選配] AUCell：以「細胞內排名」計分，對深度與 normalization 穩健。
# 計算量較大，先用小抽樣試跑即可：
# BiocManager::install("AUCell")
# library(AUCell)
# rk  <- AUCell_buildRankings(GetAssayData(ifnb, layer = "counts"))
# auc <- AUCell_calcAUC(list(IFNa = ifn.genes), rk)
# ifnb$IFNa_auc <- as.numeric(getAUC(auc)["IFNa", colnames(ifnb)])

## [踩雷示範] 雷一：背景用全基因組（頁 39–40）
## enrichGO 不給 universe → 背景變成整個 org.Hs.eg.db 註釋（~2 萬基因）。
## scRNA-seq 偵測得到的基因（偏高表現的管家類）整批變「富集」：
# ego.bad <- enrichGO(gene = sig.up, OrgDb = org.Hs.eg.db,
#                     keyType = "SYMBOL", ont = "BP")
# c(correct = nrow(as.data.frame(ego)),
#   genome  = nrow(as.data.frame(ego.bad)))   # 條目數暴增
# head(as.data.frame(ego.bad)$Description, 15)
## 榜單混入 translation / ribosome / metabolic process 類條目。
## 一行檢查：結果表若滿是管家類主題，先查 universe 給了沒。

## [踩雷示範] 雷二：拿自家 marker 清單當基因集（頁 41–42）
## cluster 的 top marker 拿回同一份資料打分——結論早就藏在輸入裡：
# Idents(ifnb) <- "seurat_annotations"
# mono.markers <- FindMarkers(ifnb, ident.1 = "CD14 Mono",
#                             only.pos = TRUE) %>%
#   slice_min(p_val_adj, n = 50) %>% rownames()
# ifnb <- AddModuleScore(ifnb, features = list(mono.markers),
#                        name = "circular")
# VlnPlot(ifnb, "circular1", group.by = "seurat_annotations")
## CD14 Mono「顯著最高分」——這不是發現，是同義反覆。
## 基因集必須來自獨立於這份資料的先驗知識（MSigDB／文獻）。

## [踩雷示範] 雷三：只報 p 最小的 20 條 GO（頁 43–44）
## 未 simplify 的結果按 p 取前 20，父子項重複計數：
# ego.raw <- enrichGO(gene = sig.up, universe = universe,
#                     OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
#                     ont = "BP")
# top20 <- head(as.data.frame(ego.raw), 20)
# top20$Description        # 十幾條是「干擾素／病毒」的不同說法
## 對策 = 第 5 節的 simplify() + 報告時按主題收斂、附代表項與
## leading edge，而不是挑二十條裡最合心意的那條。

## ---- 7. save ------------------------------------------------------
saveRDS(list(ranks        = ranks,
             gsea         = gsea,
             leading_ifna = le.ifna,
             ora          = as.data.frame(ego),
             sig_up       = sig.up,
             universe     = universe),
        "output/ifnb_b13_gsea.rds")

## ---- 8. sessioninfo ------------------------------------------------
sessionInfo()
