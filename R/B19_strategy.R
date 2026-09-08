# =====================================================================
# B19 · 複雜樣本的分析策略設計：不同腫瘤、不同戰法
# ---------------------------------------------------------------------
# 本腳本不是一條資料管線，是一套「策略模板」：
#   1. strategy_checklist()：五個決策點的問答式檢查表
#   2. strategy_presets()  ：六種場景的參數起點（範例值，依資料調整）
#   3. 模擬示範一：病人效應大小 → cell-level DE 的假陽性
#   4. 模擬示範二：CNV 訊號強弱 → 惡性判定的可行性
# 兩個模擬都可從頭 source() 重現投影片上的數字。
# 依賴：base R + stats（模擬不需要 Seurat；本系列環境為 Seurat v5 + renv）
# =====================================================================

set.seed(1234)                       # 全系列固定，結果才能重現

## ---- 1. strategy_checklist：動手前的五個問題 ----
# 用法：每開一個新專案，先跑一次，把答案寫進分析筆記。
# 答案的格式跟 QC 閾值一樣：有數字、有理由、有紀錄。

strategy_checklist <- function() {
  qs <- c(
    "統計單位：這個比較裡，獨立的觀察值是誰？（預設：病人；細胞當 n＝假重複）",
    "整合強度：哪些細胞該對齊、哪些不能硬拉？（預設：TME 對齊、惡性 per-patient）",
    "惡性判定：這種腫瘤，哪一種證據拿得到？（CNV 大→inferCNV；乾淨→突變+表型+參照系）",
    "比較設計：批次跟生物效應拆得開嗎？配對了嗎？（能配對就配對；混批在收樣本前決定）",
    "下游路線：資料撐得起這條下游嗎？（先驗收最小格的 MIN_CELLS / MIN_PAIRS，再選路）"
  )
  cat("— B19 策略檢查表（動手前逐題作答並記錄）—\n")
  for (i in seq_along(qs)) cat(sprintf(" [%d] %s\n", i, qs[i]))
  invisible(qs)
}

strategy_checklist()

## ---- 2. strategy_presets：六場景參數起點 ----
# 重要：所有數字都是「範例值，依資料調整」。
# preset 是策略的起點，不是可以照抄的答案——
# 每一個值都要回到你自己的資料分布重新驗證（B5 的鐵則）。

strategy_presets <- function() {
  list(
    ## -- GBM（實體瘤原型）------------------------------------------
    GBM = list(
      unit        = "病人；狀態比例與分數以病人為 n",
      integration = "TME 用整合空間註釋；惡性回未整合空間 per-patient",
      malignancy  = "inferCNV 為主（chr7+/10- 為錨），CNV 分數×相關×marker 三角驗證",
      design      = "core vs margin 同病人配對；margin 先驗 MIN_CELLS",
      downstream  = "Neftel 四狀態分數；缺氧/解離簽名先排除；TME 通訊",
      params = list(                       # 範例值，依資料調整
        qc_mt_rule    = "median + 3*MAD（雙峰時看谷；GBM 常落在 10-20%）",
        infercnv_cutoff = 0.1,             # 10x 用 0.1；Smart-seq2 用 1
        window_length = 101,               # 整臂事件用預設大視窗
        min_cells = 20, min_pairs = 2      # pseudobulk 門檻（B12）
      )
    ),
    ## -- 黑色素瘤（免疫腫瘤學）--------------------------------------
    melanoma = list(
      unit        = "病人；responder 定義收案前寫死",
      integration = "TME 對齊是重點：T 細胞取出後亞群重跑 HVG/PCA/分群",
      malignancy  = "譜系 marker（MLANA/PMEL）+ CNV 佐證即可，難度較低",
      design      = "responder vs non-responder：組成與表現分開問",
      downstream  = "exhaustion 打分數；TCR clonotype 追蹤；PD-1 軸通訊",
      params = list(                       # 範例值，依資料調整
        tcell_resolution = "亞群重跑時 resolution 0.6-1.2 掃描（B8）",
        min_cells = 20, min_pairs = 2
      )
    ),
    ## -- 乳癌（亞型驅動）--------------------------------------------
    breast = list(
      unit        = "病人，且亞型內比較（ER/HER2/TNBC 不混）",
      integration = "同亞型內整合 TME；跨亞型只做描述性對比",
      malignancy  = "上皮腫瘤 EPCAM+ 幾乎即腫瘤；HER2 focal amp 當錨",
      design      = "配對正常乳腺組織（切除檢體），設計期就要留",
      downstream  = "CAF 亞群、空間驗證；亞型各自 DE",
      params = list(                       # 範例值，依資料調整
        window_length = 31,                # focal amp 要小視窗才看得到
        min_cells = 20, min_pairs = 2
      )
    ),
    ## -- 肝癌（高代謝＋病毒）----------------------------------------
    liver = list(
      unit        = "病人",
      integration = "TME 對齊；肝細胞認 zonation 梯度後再談惡性分群",
      malignancy  = "CNV（1q+/8q+ 常見）+ HBV 整合證據",
      design      = "腫瘤 vs 癌旁配對；纖維化程度記進 metadata",
      downstream  = "zonation 打分數；CAF/免疫抑制軸",
      params = list(                       # 範例值，依資料調整
        qc_mt_rule = "肝細胞 mt 中位數可達 15-20%；MAD 重算，勿抄 5%",
        reference  = "GRCh38 + HBV 序列自建 mkref（B2）",
        min_cells = 20, min_pairs = 2
      )
    ),
    ## -- AML（血液腫瘤；策略幾乎全反）-------------------------------
    AML = list(
      unit        = "病人（唯一不反的一條）",
      integration = "與正常骨髓參照圖譜共嵌入：偏離參照系處才是病",
      malignancy  = "分化阻滯+突變/融合+免疫表型；CNV 退居輔助",
      design      = "診斷/緩解/復發縱向配對是主線",
      downstream  = "沿分化階層的組成分析、幹性分數；軌跡有意義",
      params = list(                       # 範例值，依資料調整
        reference  = "正常骨髓圖譜（必需品，不是選配）",
        cnv_note   = "近二倍體時 CNV 分數 AUC 可低到 0.5-0.6（見模擬二）",
        min_cells = 20, min_pairs = 2
      )
    ),
    ## -- 多區域／縱向 -----------------------------------------------
    multiregion = list(
      unit        = "病人（區域/時間點是病人內的重複測量）",
      integration = "病人內各區域可整合（共享主幹 CNV）；病人間分開",
      malignancy  = "主幹事件=早期、分支事件=晚期；subclone 家系樹",
      design      = "同病人配對；凍存湊批、混批處理，設計期防 confound",
      downstream  = "subclone 演化、區域組成、治療前後配對比較",
      params = list(                       # 範例值，依資料調整
        batch_rule = "每個定序批次都要同時含前與後（或多區域混排）",
        min_cells = 20, min_pairs = 2
      )
    )
  )
}

presets <- strategy_presets()
names(presets)
presets$AML$malignancy

## ---- 3. 模擬示範一：病人效應大小 → cell-level DE 假陽性 ----
# 問題：為什麼「病人是統計單位」在腫瘤加倍重要？
# 做法：300 個「完全沒有真效應」的基因，兩組各 3 位病人、
#       每位 120 顆細胞。病人隨機效應（log2 SD）從 0 掃到 0.6，
#       比較 cell-level wilcoxon 與 pseudobulk t 檢定的假陽性數。
# 對應投影片：頁 43（曲線圖 b19_15）。約跑 1-2 分鐘。

sim_de_fp <- function(sd_patient, n_gene = 300, n_donor = 3,
                      cells_per = 120) {
  base <- rgamma(n_gene, 1.6, scale = 1.1) + 0.05
  X <- NULL; cond <- c(); donor <- c()
  for (cc in 0:1) {
    for (dd in seq_len(n_donor)) {
      eff   <- rnorm(n_gene, 0, sd_patient)      # 病人隨機效應（log2）
      mu    <- base * 2^eff
      depth <- rlnorm(cells_per, 0, 0.25)
      m     <- outer(depth, mu) * 0.6            # cells x genes
      size  <- 1.8                               # NB 離散參數
      Xi <- matrix(rnbinom(length(m), mu = m, size = size),
                   nrow = cells_per)
      X     <- rbind(X, Xi)
      cond  <- c(cond,  rep(cc, cells_per))
      donor <- c(donor, rep(cc * n_donor + dd, cells_per))
    }
  }
  ## cell-level wilcoxon（雷三的做法）
  p_cell <- apply(X, 2, function(g)
    wilcox.test(g[cond == 1], g[cond == 0])$p.value)
  ## pseudobulk + Welch t（正確做法的最小版；正式分析用 DESeq2，B12）
  ids <- sort(unique(donor))
  pb  <- t(sapply(ids, function(d) colSums(X[donor == d, , drop = FALSE])))
  cpm <- log2(pb / rowSums(pb) * 1e6 + 1)
  cc2 <- sapply(ids, function(d) cond[donor == d][1])
  p_pb <- sapply(seq_len(n_gene), function(j)
    t.test(cpm[cc2 == 1, j], cpm[cc2 == 0, j])$p.value)
  c(fp_cell = sum(p_cell * n_gene < 0.05),       # Bonferroni 簡化版
    fp_pb   = sum(p_pb   * n_gene < 0.05))
}

sds <- c(0, 0.15, 0.30, 0.45, 0.60)
fp  <- t(sapply(sds, sim_de_fp))
rownames(fp) <- paste0("SD=", sds)
print(fp)
# 預期形狀（種子不同數字略異）：fp_cell 隨 SD 單調上升到數十個，
# fp_pb 全程貼 0。腫瘤的病人效應天生大（共享 CNV），你永遠在右端。

## ---- 4. 模擬示範二：CNV 訊號強弱 → 惡性判定可行性 ----
# 問題：為什麼「惡性判定看癌別」？近二倍體腫瘤為何換證據？
# 做法：3,000 個基因沿 22 條染色體排列，惡性細胞帶不同大小的
#       CNV 事件。跑 inferCNV 式原理鏈（log 化→減參考→滑動平均），
#       算每顆細胞的 CNV 分數，看惡性/正常兩組分數的 AUC。
# 對應投影片：頁 27（小提琴＋AUC 圖 b19_10）。

chr_ng <- round(c(2058, 1309, 1078, 752, 886, 1048, 979, 677, 786, 733,
                  1298, 1034, 327, 830, 613, 873, 1197, 270, 1472, 544,
                  234, 488) / 18262 * 3000)
chr_b  <- c(0, cumsum(chr_ng)); G <- sum(chr_ng)

cnv_profile <- function(events) {           # events: list(c(chr, f0, f1, ratio))
  prof <- rep(1, G)
  for (e in events) {
    s <- chr_b[e[1]] + 1; n <- chr_ng[e[1]]
    idx <- (s + floor(e[2] * n)):(s + ceiling(e[3] * n) - 1)
    prof[idx] <- e[4]
  }
  prof
}

movavg <- function(x, w) as.numeric(stats::filter(x, rep(1 / w, w),
                                                  sides = 2)) |>
  (\(v) ifelse(is.na(v), x, v))()           # 邊緣視窗不足時退回原值

sim_cnv_auc <- function(profile, n_mal = 60, n_ref = 60, w = 51) {
  base  <- rgamma(G, 2.0, scale = 1.6) + 0.05
  make  <- function(ratio, n) {
    depth <- rlnorm(n, 0, 0.3)
    m <- outer(depth, base * ratio)
    matrix(rnbinom(length(m), mu = m, size = 2), nrow = n)
  }
  X    <- rbind(make(rep(1, G), n_ref), make(profile, n_mal))
  expr <- log2(X / rowSums(X) * 1e4 + 1)
  rel  <- sweep(expr, 2, colMeans(expr[1:n_ref, ]))   # 減掉參考平均
  rel  <- pmin(pmax(rel, -3), 3)
  sm   <- rel
  for (ci in 1:22) {                                  # 逐染色體平滑
    idx <- (chr_b[ci] + 1):chr_b[ci + 1]
    sm[, idx] <- t(apply(rel[, idx, drop = FALSE], 1, movavg, w = w))
  }
  score <- rowMeans(sm^2)                             # 每細胞 CNV 分數
  ref_s <- score[1:n_ref]; mal_s <- score[-(1:n_ref)]
  ## AUC = Mann-Whitney U / (n1*n2)
  u <- wilcox.test(mal_s, ref_s, alternative = "greater")$statistic
  as.numeric(u) / (n_mal * n_ref)
}

auc_gbm   <- sim_cnv_auc(cnv_profile(list(c(7, 0, 1, 1.5),   # chr7 整條增
                                          c(10, 0, 1, 0.5))))# chr10 整條缺
auc_small <- sim_cnv_auc(cnv_profile(list(c(13, 0, 0.6, 1.5))))
auc_dipl  <- sim_cnv_auc(cnv_profile(list()))          # 近二倍體：無事件
cat(sprintf("AUC — GBM 式: %.2f｜單一小事件: %.2f｜近二倍體: %.2f\n",
            auc_gbm, auc_small, auc_dipl))
# 預期：GBM 式 ~0.95+；小事件 ~0.6-0.8；近二倍體 ~0.5（＝擲硬幣）。
# 工具沒有變，是訊號不在——這時惡性判定要換證據
# （突變/融合、免疫表型、對照正常參照系的分化阻滯；見 B16 與本集 4-5）。

## [踩雷示範] 雷一：把 PBMC 教學管線的參數原封不動套進腫瘤
## 後果見投影片頁 2 與頁 38-39；勿在真實分析中執行這種「照抄」。
# pbmc_params <- list(mt_cut = 5, resolution = 0.5, dims = 1:10)
# tumor <- run_pipeline(tumor_data, params = pbmc_params)   # 錯誤示範

## [踩雷示範] 雷二：跨腫瘤抄參數（GBM 的 mt<20% 套到 PBMC）
## 一行檢查：先看分布再定線——分布圖會直接拆穿抄來的數字。
# hist(pbmc$percent.mt, breaks = 60)   # 若主體在 2-4%，20% 線形同虛設

## ---- 5. session info ----
sessionInfo()
