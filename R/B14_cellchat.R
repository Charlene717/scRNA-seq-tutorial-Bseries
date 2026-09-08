# =====================================================================
# B14_cellchat.R — B14「細胞通訊：配體–受體配對與 CellChat」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 14 集
# 資料：ifnb（Kang et al. 2018, Nat Biotechnol；GEO GSE96583）
#       承接 B11 的整合物件 output/ifnb_b11_integrated.rds
#       （含 stim 條件標籤與 seurat_annotations 型別標籤）
# 環境：R >= 4.3、Seurat v5.x、CellChat（GitHub: jinworks/CellChat）、
#       patchwork、future
#       安裝：remotes::install_github("jinworks/CellChat")
#
# 使用方式：從專案根目錄 source() 或逐段執行。
# computeCommunProb 兩個條件合計約 10–20 分鐘（視機器與 workers）。
# 分節標記與影片段落一一對應。
# =====================================================================

## ---- 0. setup -----------------------------------------------------
library(Seurat)
library(CellChat)
library(patchwork)
set.seed(1234)                       # 全系列固定 seed，結果可重現

# 平行化（可選）：加速 computeCommunProb 與 permutation
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 4 * 1024^3)

dir.create("output", showWarnings = FALSE)

## ---- 1. load ------------------------------------------------------
# 承接 B11：整合物件（含未整合 / CCA / Harmony 三套 reduction）。
# 通訊分析完全不用 reduction，只用 RNA 表達值 + 型別標籤——
# 整合的貢獻是「讓跨批次的型別標籤可靠」（B10–B11 的工作）。
ifnb <- readRDS("output/ifnb_b11_integrated.rds")
ifnb <- JoinLayers(ifnb)             # 把 split 的 RNA layers 合回來
Idents(ifnb) <- "seurat_annotations" # 型別標籤 = 全集的承重牆
sort(table(Idents(ifnb)), decreasing = TRUE)
# 預期（官方註釋，13 型別、13,999 顆）：
#   CD14 Mono 4362 / CD4 Naive T 2504 / CD4 Memory T 1762 /
#   CD16 Mono 1044 / B 978 / CD8 T 814 / T activated 633 / NK 619 /
#   DC 472 / B Activated 388 / Mk 236 / pDC 132 / Eryth 55
# 注意小群：Eryth、pDC——踩雷區主角

## ---- 2. create ----------------------------------------------------
# 條件比較的鐵則：兩個條件「各自」建物件，最後 mergeCellChat。
ctrl <- subset(ifnb, subset = stim == "CTRL")   # 6,548 顆
cc.ctrl <- createCellChat(object = ctrl,
                          group.by = "ident", assay = "RNA")
# 吃的是 RNA assay 的 data 層（log-normalized）；
# console 會列出納入分析的細胞群，掃一眼確認 13 型別都在
cc.ctrl@DB <- CellChatDB.human       # 小鼠請換 CellChatDB.mouse
showDatabaseCategory(CellChatDB.human)
# 三大類：Secreted Signaling（大宗）/ ECM-Receptor / Cell-Cell Contact
# 只想分析分泌訊號時：
# cc.ctrl@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")

## ---- 3. infer -----------------------------------------------------
cc.ctrl <- subsetData(cc.ctrl)       # 裁到資料庫基因；必跑
cc.ctrl <- identifyOverExpressedGenes(cc.ctrl)
cc.ctrl <- identifyOverExpressedInteractions(cc.ctrl)
nrow(cc.ctrl@LR$LRsig)               # 進入計分的互動數（依資料而異）

# 本集最重要的一行。兩個參數要能對人解釋：
# * type = "triMean"（預設）：Tukey 三均值，約需 25% 細胞表達才非零，
#   保守、抓到的互動「少而強」。研究小比例細胞的訊號時可改
#   type = "truncatedMean", trim = 0.1（約需 10% 細胞表達），並寫進論文。
# * population.size = FALSE（預設）：不把群大小乘進分數。
#   問「組織總流量」時設 TRUE；解離資料的組成通常已失真，預設較穩。
cc.ctrl <- computeCommunProb(cc.ctrl, type = "triMean",
                             population.size = FALSE)
# permutation 檢定內建於上一步（預設 nboot = 100 → p 值下限 0.01）

cc.ctrl <- filterCommunication(cc.ctrl, min.cells = 10)
# 少於 10 顆細胞的群不參與——為什麼非濾不可，見 [踩雷示範] 雷一

## ---- 4. aggregate -------------------------------------------------
cc.ctrl <- computeCommunProbPathway(cc.ctrl)  # L-R → pathway（netP）
cc.ctrl <- aggregateNet(cc.ctrl)              # 型別×型別 count / weight
groupSize <- as.numeric(table(cc.ctrl@idents))

## ---- 5. visualize -------------------------------------------------
# 5.1 circle plot：count 與 weight 是兩張不同的圖，報告要講清楚
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cc.ctrl@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cc.ctrl@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction strength")
par(mfrow = c(1, 1))

# 5.2 bubble plot：circle 給方向，bubble 給名字
netVisual_bubble(cc.ctrl,
                 sources.use = c("CD14 Mono", "CD16 Mono"),
                 targets.use = c("NK", "CD8 T", "CD4 Naive T", "B"),
                 remove.isolate = TRUE)
# 判讀：點大小 = permutation p；顏色 = 通訊機率；空格 = 沒過門檻，不是零

# 5.3 訊號角色：誰在送、誰在收
cc.ctrl <- netAnalysis_computeCentrality(cc.ctrl, slot.name = "netP")
netAnalysis_signalingRole_scatter(cc.ctrl)
# 預期：單核球在右下（主要發送者），B / pDC 偏左上（主要接收者）

## ---- 6. stim ------------------------------------------------------
# 包成函式的理由：保證兩個條件參數「一模一樣」，比較才有意義
run_cellchat <- function(obj) {
  cc <- createCellChat(obj, group.by = "ident", assay = "RNA")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean",
                          population.size = FALSE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  netAnalysis_computeCentrality(cc, slot.name = "netP")
}
cc.stim <- run_cellchat(subset(ifnb, subset = stim == "STIM"))

cc <- mergeCellChat(list(CTRL = cc.ctrl, STIM = cc.stim),
                    add.names = c("CTRL", "STIM"))

## ---- 7. compare ---------------------------------------------------
# 7.1 最粗：總互動數與總強度
compareInteractions(cc, show.legend = FALSE, group = c(1, 2)) +
  compareInteractions(cc, show.legend = FALSE, group = c(1, 2),
                      measure = "weight")

# 7.2 中層：每條 pathway 的資訊流（= 通訊機率總和）
rankNet(cc, mode = "comparison", stacked = TRUE, do.stat = TRUE)
# 預期：CXCL 等干擾素下游軸在 STIM 大幅上升——與 B12 DE（CXCL10 上調）、
# B13 GSEA（IFN response）交叉印證；三個角度指向同一件事

# 7.3 細層：差異網路（紅 = STIM 強，藍 = STIM 弱）
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cc, weight.scale = TRUE)
netVisual_diffInteraction(cc, weight.scale = TRUE, measure = "weight")
par(mfrow = c(1, 1))

# 7.4 鎖定單一 pathway 看細節（例：CXCL）
# pathways.show <- "CXCL"
# netVisual_aggregate(cc.ctrl, signaling = pathways.show,
#                     layout = "circle")
# netVisual_aggregate(cc.stim, signaling = pathways.show,
#                     layout = "circle")

# 報告措辭（雷二的解藥）：通訊機率（強度）與 p 值（穩定性）成對報告；
# 動詞用「提示 / 預測」，不用「證明」。

## ---- 8. save ------------------------------------------------------
saveRDS(cc, "output/ifnb_b14_cellchat.rds")   # 合併物件（含兩條件）

sessionInfo()                        # 收進附錄；審稿人的好朋友

# =====================================================================
# [踩雷示範] 以下三段預設註解掉。想體驗災難，取消註解單獨執行。
# 千萬不要把這些結果用在正式分析。
# =====================================================================

## [踩雷示範] 雷一：min.cells 太小的型別照算（頁 39–40）
## Eryth 在 CTRL 只有 20 多顆。把 min.cells 調成 3 讓它入榜，
## 再用重抽樣看它的通訊分數穩不穩：
# cc.loose <- filterCommunication(cc.ctrl, min.cells = 3)
# ## 重抽樣：對接收方細胞 bootstrap 30 次，重算 CTRL 網路，
# ## 記錄「CD14 Mono → Eryth」的聚合強度
# probs <- replicate(30, {
#   cells <- unlist(lapply(split(colnames(ctrl), Idents(ctrl)),
#                          function(x) sample(x, length(x),
#                                             replace = TRUE)))
#   cc.b <- createCellChat(ctrl[, cells], group.by = "ident",
#                          assay = "RNA")
#   cc.b@DB <- CellChatDB.human
#   cc.b <- subsetData(cc.b)
#   cc.b <- identifyOverExpressedGenes(cc.b)
#   cc.b <- identifyOverExpressedInteractions(cc.b)
#   cc.b <- computeCommunProb(cc.b)
#   cc.b <- aggregateNet(cc.b)
#   cc.b@net$weight["CD14 Mono", "Eryth"]
# })
# quantile(probs, c(0.05, 0.5, 0.95)) # 90% 區間 vs 中位數：小群大幅抖動
## 對照大群（如 CD4 Naive T）做同樣的事，區間窄得多。
## 教訓：min.cells = 10 是底線；幾十顆的群，結論照樣要打折。

## [踩雷示範] 雷二：把通訊分數當效應量報（頁 41–42）
## 把每條互動的 (通訊機率, p 值) 攤開來看——p 觸底的互動一大排，
## 機率卻從 0.05 到 0.9 都有。只報 p 值 = 把它們寫成一樣強：
# prob <- cc.ctrl@net$prob; pval <- cc.ctrl@net$pval
# df <- data.frame(prob = as.vector(prob), pval = as.vector(pval))
# df <- df[df$prob > 0, ]
# plot(df$prob, -log10(df$pval + 1e-3), pch = 16,
#      col = rgb(0.15, 0.39, 0.66, 0.4),
#      xlab = "communication probability (effect size)",
#      ylab = "-log10 permutation p")
## 一行檢查：p 最小的互動裡，機率排名墊底的有幾條？
# small <- df[df$pval <= 0.01, ]
# mean(small$prob < 0.1)              # 「顯著但微弱」互動的比例
## 教訓：p 值管真假、機率管強弱；報告永遠成對念。

## [踩雷示範] 雷三：cluster 粒度錯配——8 群 vs 14 群翻盤（頁 43–44）
## 用 Seurat clusters 取代型別標籤：resolution 0.1 → 約 8 群，
## 0.8 → 約 14 群（實際群數依資料微動）。同一條流程各跑一次：
# ctrl2 <- FindNeighbors(ctrl, reduction = "integrated.cca",
#                        dims = 1:30)
# ctrl2 <- FindClusters(ctrl2, resolution = 0.1,
#                       cluster.name = "coarse")   # 併群：T 全在一起
# ctrl2 <- FindClusters(ctrl2, resolution = 0.8,
#                       cluster.name = "fine")     # 活化 T 單獨成群
# table(ctrl2$coarse); table(ctrl2$fine)
# cc.coarse <- run_cellchat(SetIdent(ctrl2, value = "coarse"))
# cc.fine   <- run_cellchat(SetIdent(ctrl2, value = "fine"))
# par(mfrow = c(1, 2))
# netVisual_circle(cc.coarse@net$weight, weight.scale = TRUE,
#                  title.name = "coarse (~8 clusters)")
# netVisual_circle(cc.fine@net$weight, weight.scale = TRUE,
#                  title.name = "fine (~14 clusters)")
## 一行檢查（跑通訊之前就該做）：你關心的受體專一在哪個子群？
# VlnPlot(ctrl, features = "CXCR3", group.by = "seurat_annotations")
## CXCR3 集中在 T activated / NK——粗粒度把它平均進大 T 群，軸線消失。
## 教訓：粒度是結論的一部分；重要結論換一個粒度重跑，並寫進方法段。
