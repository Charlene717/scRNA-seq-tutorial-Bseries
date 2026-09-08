# =====================================================================
# data/download_data.R — 下載本系列使用的公開資料集
#
# 從專案根目錄執行：source("data/download_data.R")
# 原始資料不進版本控制；刪掉重跑本腳本即可還原。
# =====================================================================

dir.create("data", showWarnings = FALSE)

## ---- PBMC 3k（B1、B5–B10 主線資料）---------------------------------
# 10x Genomics 官方公開資料：2,700 顆 PBMC，
# filtered_gene_bc_matrices（hg19），約 7.6 MB。
pbmc3k.url <- paste0("https://cf.10xgenomics.com/samples/cell/",
                     "pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz")
pbmc3k.tar <- "data/pbmc3k.tar.gz"

if (!dir.exists("data/pbmc3k/filtered_gene_bc_matrices/hg19")) {
  options(timeout = 600)                  # 網速慢時不要中途放棄
  download.file(pbmc3k.url, pbmc3k.tar, mode = "wb")
  untar(pbmc3k.tar, exdir = "data/pbmc3k")
  file.remove(pbmc3k.tar)                 # 解壓完即可移除壓縮檔
}

# 驗證：三個檔案都在才算成功
list.files("data/pbmc3k", recursive = TRUE)
# [1] "filtered_gene_bc_matrices/hg19/barcodes.tsv"
# [2] "filtered_gene_bc_matrices/hg19/genes.tsv"
# [3] "filtered_gene_bc_matrices/hg19/matrix.mtx"

## ---- ifnb（B11–B12 整合與 DE）--------------------------------------
# 由 SeuratData 提供，於 B11 首次使用時安裝：
#   install.packages("SeuratData",
#     repos = c("https://seurat.nygenome.org", getOption("repos")))
#   SeuratData::InstallData("ifnb")
# 原始出處：Kang et al. 2018（GEO: GSE96583）
