#!/usr/bin/env bash
## =====================================================================
## B03_count.sh — B3「cellranger count 實戰與 web_summary 判讀」隨集腳本
## ---------------------------------------------------------------------
## 用法：不建議整檔一次執行；跟著影片「逐段」複製到終端機執行。
##       分節標記與影片段落一一對應。
## 前置：B2 已下載參考基因組；Cell Ranger 8.x 已安裝並在 PATH 上。
## 注意：[踩雷示範] 區塊預設全部註解掉，看懂後果即可，不要照跑。
## =====================================================================
set -euo pipefail

## ---- 0. setup ----------------------------------------------------------
## 專案根目錄下執行（cellranger 會「在哪裡跑就在哪裡生資料夾」）。
## 路徑一律相對於專案根目錄，不用絕對路徑。
REF=refdata-gex-GRCh38-2020-A          # B2 下載的預建參考（資料夾）
FASTQ_DIR=fastq/pbmc_1k_v3_fastqs      # FASTQ 所在資料夾
SAMPLE=pbmc_1k_v3                      # FASTQ 檔名的樣本前綴（一字不差）
RUN_ID=pbmc1k_count                    # 輸出資料夾名（自訂）

cellranger --version                   # 驗證安裝：印出 8.x 版號才繼續

## ---- 1. download --------------------------------------------------------
## PBMC 1k v3 官方公開資料（約 5.2 GB）。原始資料不進 repo。
## 來源頁：10x Genomics 網站 Datasets → "1k PBMCs from a Healthy Donor (v3)"
mkdir -p fastq
# curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
# tar -xf pbmc_1k_v3_fastqs.tar -C fastq/

## 【跑不動的備案】機器不夠力時，直接下載官方跑好的輸出，
## 從第 4 節接著做，判讀內容完全一樣能跟上：
# curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_web_summary.html
# curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz
# curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_raw_feature_bc_matrix.tar.gz

## ---- 2. preflight -------------------------------------------------------
## 開跑前的兩個檢查（影片頁 12–13：--sample 必須對上檔名前綴）。
ls -lh "$FASTQ_DIR" | head            # 親眼看檔名：樣本名_S1_L001_R1_001.fastq.gz
cellranger sitecheck | head -30       # 機器體檢：核心數、記憶體、ulimit

## ---- 3. count -----------------------------------------------------------
## 本集主角。8 核 64 GB 下 PBMC 1k 約 1–2 小時。
## 中斷（斷電、被砍）不用從頭來：同一目錄重跑同一行指令即可續跑。
cellranger count \
  --id="$RUN_ID" \
  --transcriptome="$REF" \
  --fastqs="$FASTQ_DIR" \
  --sample="$SAMPLE" \
  --create-bam=true \
  --localcores=8 \
  --localmem=64
## 成功的判準只有一行：最後印出 "Pipestance completed successfully!"

## ---- 4. inspect ---------------------------------------------------------
## 輸出導覽（影片頁 29–32）。
ls "$RUN_ID"/outs/                    # 報告、矩陣、molecule_info、BAM、cloupe
# open  "$RUN_ID"/outs/web_summary.html     # macOS
# xdg-open "$RUN_ID"/outs/web_summary.html  # Linux

## 矩陣三件套長什麼樣：
ls "$RUN_ID"/outs/filtered_feature_bc_matrix/
##   barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz

## filtered 的細胞數（barcodes 行數 = 欄數），要跟 web_summary 一致：
zcat "$RUN_ID"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | wc -l

## ---- 5. sanity check ----------------------------------------------------
## 「一行檢查」：metrics_summary.csv 自動驗收（影片頁 37 的雷一防線）。
## 回收率 = Estimated cells / 上機數，應約 0.4–0.7；此例上機 1,600 顆。
LOADED=1600
python3 - "$RUN_ID/outs/metrics_summary.csv" "$LOADED" <<'PY'
import csv, sys
row = next(csv.DictReader(open(sys.argv[1])))
cells = int(row["Estimated Number of Cells"].replace(",", ""))
ratio = cells / int(sys.argv[2])
print(f"cells={cells:,}  loaded={sys.argv[2]}  recovery={ratio:.2f}")
print("OK" if 0.4 <= ratio <= 0.7 else "!! recovery 超出 0.4–0.7，先看 knee plot")
for k in ["Valid Barcodes", "Reads Mapped Confidently to Transcriptome",
          "Sequencing Saturation", "Median Genes per Cell"]:
    print(f"{k}: {row[k]}")
PY

## =====================================================================
## [踩雷示範] 以下三段預設註解。每一段都對應影片踩雷區的一個實跑後果。
## =====================================================================

## [踩雷示範] trap 1 — --sample 對不上 FASTQ 前綴（影片頁 13）
## 後果：preflight 30 秒內報錯
##   "error: No input FASTQs were found for the requested parameters."
# cellranger count --id=trap_sample \
#   --transcriptome="$REF" --fastqs="$FASTQ_DIR" \
#   --sample=pbmc1k --create-bam=false \
#   --localcores=4 --localmem=32

## [踩雷示範] trap 2 — 把 raw 當 filtered 餵給 R（影片頁 32）
## 後果：dim() 顯示約 680 萬「顆細胞」；建 Seurat 物件會吃爆記憶體。
## （在 R 內執行；此處以 heredoc 形式收錄，方便對照）
# Rscript - <<'RS'
# library(Seurat)
# raw  <- Read10X("pbmc1k_count/outs/raw_feature_bc_matrix/")
# filt <- Read10X("pbmc1k_count/outs/filtered_feature_bc_matrix/")
# print(dim(raw))   # 36601 x 6794880  ← 幾乎全是空滴
# print(dim(filt))  # 36601 x    1187  ← cell calling 的倖存者
# RS

## [踩雷示範] trap 3 — 用 --force-cells 硬灌細胞數（影片頁 36–37）
## 後果：把 ambient 平台整段收進矩陣，web_summary 的
##   Mean Reads per Cell 與 Median Genes per Cell 被同步稀釋。
# cellranger count --id=trap_force \
#   --transcriptome="$REF" --fastqs="$FASTQ_DIR" \
#   --sample="$SAMPLE" --create-bam=false \
#   --force-cells=8000 --localcores=8 --localmem=64

## ---- 6. session info ----------------------------------------------------
## 可重現性記錄：把版本與機器狀態留在分析筆記裡。
cellranger --version
uname -a
date
