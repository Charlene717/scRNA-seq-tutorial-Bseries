# =====================================================================
# B00_setup.R — B 系列練習腳本 0：套件安裝、專案結構、資料下載
#
# ⚠ 尚未實跑驗證：本檔的安裝清單依各集腳本的 library() 反推、並比照
#   Q 系列 00_setup.R 的寫法整理。請在乾淨環境完整跑過一次、確認無誤
#   之後，把這段警語刪掉再發布。
#
# 對應影片：B1「環境建置與你的第一張 UMAP」
# 執行方式：在 RStudio 開啟專案（.Rproj），從專案根目錄逐段執行。
#           第一次執行約 20–40 分鐘（安裝套件為主；資料只有約 8 MB）。
# 之後每一支 B0x / B1x 腳本都假設本檔已經跑過。
# =====================================================================

## ---- 1. 套件 -------------------------------------------------------
# CRAN
cran <- c("Seurat",        # 全系列主力，需 v5
          "ggplot2", "dplyr", "tidyr", "patchwork",   # 繪圖與資料整理
          "clustree",      # B8 解析度掃描
          "msigdbr",       # B13 基因集
          "nnls",          # B17 反卷積
          "remotes")       # 裝 GitHub 套件用
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c("SingleCellExperiment",           # 物件轉換
          "SingleR", "celldex",             # B10 自動註釋
          "DESeq2",                         # B12 pseudobulk
          "fgsea", "clusterProfiler",       # B13 GSEA / ORA
          "org.Hs.eg.db",                   # B13 基因 ID 轉換
          "slingshot",                      # B15 軌跡
          "infercnv")                       # B16 惡性判定（另需系統層級 JAGS，見 §2）
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)

# GitHub
# 常見錯誤「HTTP error 401 Bad credentials」= 環境變數 GITHUB_PAT 裡有過期／無效的 token。
# 下面的 helper 先用匿名下載（不需 token）；匿名每小時 60 次，超過再設定有效的 PAT。
install_gh <- function(repo, pkg = basename(repo)) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  old <- Sys.getenv("GITHUB_PAT"); Sys.unsetenv("GITHUB_PAT")
  on.exit(if (nzchar(old)) Sys.setenv(GITHUB_PAT = old))
  ok <- tryCatch({ remotes::install_github(repo, upgrade = "never"); TRUE },
                 error = function(e) { message("  ! ", repo, " 安裝失敗：", conditionMessage(e)); FALSE })
  if (!ok) message("    → 稍後手動執行：remotes::install_github(\"", repo, "\")")
  invisible(ok)
}
install_gh("chris-mcginnis-ucsf/DoubletFinder")            # B5 doublet 偵測（GitHub only）
install_gh("satijalab/seurat-data", "SeuratData")          # B11–B14 的 ifnb 資料
install_gh("satijalab/azimuth", "Azimuth")                 # B10 映射式註釋（相依較多）
install_gh("jinworks/CellChat")                            # B14 細胞通訊（約 5–10 分鐘）

## ---- 2. 系統層級相依：JAGS（只有 B16 需要）------------------------
# rjags 只是介面，R 裝不了 JAGS 本體；裝完要重開 R。
#   macOS  : brew install jags
#   Ubuntu : sudo apt-get install jags
#   Windows: https://sourceforge.net/projects/mcmc-jags/
cat("infercnv 可載入：", requireNamespace("infercnv", quietly = TRUE), "\n")
# FALSE = JAGS 還沒裝好，B16_infercnv.R 會跑不動；其餘各集不受影響。

## ---- 3. 版本驗證 ---------------------------------------------------
library(Seurat)
cat("R      ", R.version.string, "\n")
cat("Seurat ", as.character(packageVersion("Seurat")), "\n")
stopifnot(packageVersion("Seurat") >= "5.0.0")   # 全系列 Seurat v5
set.seed(1234)                                   # 全系列固定 seed，結果可重現

## ---- 4. 專案結構 ---------------------------------------------------
# 路徑一律相對於專案根目錄；全系列不用 setwd()。
for (d in c("data", "R", "output", "output/figs", "output/rds", "output/tables")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ---- 5. 資料 (1)：PBMC 3k（B1、B5–B10 主線、B17 參考）--------------
# 10x Genomics 官方公開資料：2,700 顆 PBMC，filtered_gene_bc_matrices（hg19），約 7.6 MB。
# 若直連失效：到 10x Genomics Datasets 搜尋 "PBMC 3k"，下載同名檔案放進 data/。
pbmc3k.url <- paste0("https://cf.10xgenomics.com/samples/cell/",
                     "pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz")
pbmc3k.tar <- "data/pbmc3k.tar.gz"
if (!dir.exists("data/pbmc3k/filtered_gene_bc_matrices/hg19")) {
  options(timeout = 600)                  # 網速慢時不要中途放棄
  download.file(pbmc3k.url, pbmc3k.tar, mode = "wb")
  untar(pbmc3k.tar, exdir = "data/pbmc3k")
  file.remove(pbmc3k.tar)                 # 解壓完即可移除壓縮檔
}
list.files("data/pbmc3k", recursive = TRUE)
# 預期三個檔案：
#   filtered_gene_bc_matrices/hg19/barcodes.tsv
#   filtered_gene_bc_matrices/hg19/genes.tsv
#   filtered_gene_bc_matrices/hg19/matrix.mtx

## ---- 6. 資料 (2)：ifnb（B11–B14 整合、DE、富集、通訊）-------------
# 由 SeuratData 提供，不是直接下載。原始出處：Kang et al. 2018（GEO: GSE96583）。
if (requireNamespace("SeuratData", quietly = TRUE)) {
  if (!"ifnb" %in% SeuratData::AvailableData()$Dataset[SeuratData::AvailableData()$Installed]) {
    SeuratData::InstallData("ifnb")       # 約 200 MB
  }
}

## ---- 7. 資料 (3)：PBMC 1k FASTQ（B2–B3 上游，選配）----------------
# 只有要自己跑 Cell Ranger 才需要，檔案較大（數 GB）且需要 Linux 環境。
# 下載位置與指令見 R/B02_mkref_mkfastq.sh 與 R/B03_count.sh；
# 跑不動的話那兩集提供現成輸出，可以跳過運算直接看判讀。

## ---- 8. 完成 -------------------------------------------------------
cat("\n環境建置完成。接著執行 R/B01_first_umap.R。\n")
sessionInfo()
