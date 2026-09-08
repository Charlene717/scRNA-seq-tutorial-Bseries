#!/usr/bin/env bash
# ============================================================
# B04_aggr_multi.sh — B4「多樣本與多模態：aggr、multi、Feature Barcode」隨集腳本
#
# 用法：
#   本腳本設計成「逐段執行」：每一節可以單獨複製到終端機跑，
#   也可以整支 `bash R/B04_aggr_multi.sh`（需要 Cell Ranger 與範例資料在位）。
#   沒有環境時，把它當「可複製貼上的範本庫」：三張 CSV 的 heredoc
#   範本（aggregation / multi config / feature reference）都在這裡。
#
# 對應影片：分節標記與影片段落一一對應。
# 慣例：
#   - 路徑一律相對於專案根目錄（在專案根目錄執行本腳本）
#   - 踩雷示範以 ## [踩雷示範] 標明並預設註解掉，避免照跑
# ============================================================
set -euo pipefail

## ---- 0. setup ----------------------------------------------------------
# 路徑變數集中在這裡改。RUNS_DIR 底下應有每個樣本的 count 輸出：
#   runs/<sample>/outs/molecule_info.h5
#   runs/<sample>/outs/filtered_feature_bc_matrix/
REF="refs/refdata-gex-GRCh38-2024-A"   # B2 準備的參考基因組
RUNS_DIR="runs"                        # 各樣本 count 輸出的上層目錄
OUT_DIR="output/b04"
mkdir -p "$OUT_DIR"

command -v cellranger >/dev/null 2>&1 \
  || echo "[提醒] 找不到 cellranger——本腳本可以只當範本讀"

## ---- 1. aggr-csv --------------------------------------------------------
# aggregation CSV：兩欄。sample_id 決定 barcode 字尾（-1、-2…按列序），
# molecule_h5 指向 count 輸出的 molecule_info.h5（不是矩陣資料夾！）
cat > "$OUT_DIR/aggr.csv" <<'EOF'
sample_id,molecule_h5
ctrl_1,runs/ctrl_1/outs/molecule_info.h5
ctrl_2,runs/ctrl_2/outs/molecule_info.h5
treat_1,runs/treat_1/outs/molecule_info.h5
treat_2,runs/treat_2/outs/molecule_info.h5
EOF
echo "[1] aggregation CSV → $OUT_DIR/aggr.csv"

## ---- 2. aggr-run --------------------------------------------------------
# --normalize 的兩種寫法（影片頁 11–12）：
#   mapped（預設）：把每個樣本的 reads 抽到最淺樣本的深度——會丟 reads
#   none          ：不抽樣，全部保留；深度差留給下游 normalization（B6）
#                   與整合（B11）處理。拿去正式分析建議用這個。
# cellranger aggr --id=all_samples \
#   --csv="$OUT_DIR/aggr.csv" \
#   --normalize=none
#
# 跑完第一件事：開 all_samples/outs/web_summary.html 的 aggregation 表，
# 看每個樣本被抽掉的比例（normalize=none 時應為 0%）。

## ---- 3. multi-config ----------------------------------------------------
# multi config CSV：一個檔案、數個 [區塊]（影片頁 18–19）。
# 範例：GEX + CellPlex（CMO）多重化，四個 donor 混一個 channel。
#   [gene-expression]  全域資源與開關
#   [libraries]        每組 FASTQ 是哪種 library（第三欄 feature_types）
#   [samples]          標籤 ↔ 樣本 對照表——漏了它就拆不了樣本
cat > "$OUT_DIR/multi_config.csv" <<'EOF'
[gene-expression]
reference,refs/refdata-gex-GRCh38-2024-A
create-bam,false

[libraries]
fastq_id,fastqs,feature_types
pbmc_gex,fastq/gex,Gene Expression
pbmc_cmo,fastq/cmo,Multiplexing Capture

[samples]
sample_id,cmo_ids
donor_A,CMO301
donor_B,CMO302
donor_C,CMO303
donor_D,CMO304
EOF
echo "[3] multi config CSV → $OUT_DIR/multi_config.csv"

## ---- 4. multi-run -------------------------------------------------------
# cellranger multi --id=pbmc_multi --csv="$OUT_DIR/multi_config.csv"
#
# 跑完看輸出結構（影片頁 20）。下游要讀的是各樣本的
# sample_filtered_feature_bc_matrix：
# ls pbmc_multi/outs/per_sample_outs/
# ls pbmc_multi/outs/per_sample_outs/donor_A/count/

## ---- 5. fb-csv ----------------------------------------------------------
# CITE-seq / hashing 走 count 的兩張表（影片頁 24–25）。
# (a) libraries CSV：兩組 FASTQ 各是什麼 library
cat > "$OUT_DIR/libraries.csv" <<'EOF'
fastqs,sample,library_type
fastq/gex,pbmc_gex,Gene Expression
fastq/adt,pbmc_adt,Antibody Capture
EOF

# (b) feature reference CSV：抗體的「字典」。
#     sequence 欄只能複製貼上（來源＝試劑廠商官方對照表），
#     read / pattern 照廠商文件填。ADT 與 HTO 都寫 Antibody Capture。
cat > "$OUT_DIR/feature_ref.csv" <<'EOF'
id,name,read,pattern,sequence,feature_type
CD3,CD3_TotalB,R2,5PNNNNNNNNNN(BC),AACAAGACCCTTGAG,Antibody Capture
CD19,CD19_TotalB,R2,5PNNNNNNNNNN(BC),CTGGGCAATTACTCG,Antibody Capture
HTO1,Hashtag1,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
HTO2,Hashtag2,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
EOF
echo "[5] libraries.csv 與 feature_ref.csv → $OUT_DIR/"

## ---- 6. fb-count --------------------------------------------------------
# cellranger count --id=pbmc_citeseq \
#   --transcriptome="$REF" \
#   --libraries="$OUT_DIR/libraries.csv" \
#   --feature-ref="$OUT_DIR/feature_ref.csv" \
#   --create-bam=false
#
# 跑完第一件事：web_summary 的 Antibody 段看 Antibody Reads Usable。
# 正常的 FB library 通常過半；個位數百分比＝紅燈（先查 feature_ref）。

## ---- 7. r-handoff -------------------------------------------------------
# R 端承接（影片頁 29–30）：讀入雙 assay ＋ HTODemux。
# 需要 Seurat v5；把矩陣路徑換成你的輸出。逐行內容與投影片一致。
# Rscript - <<'RSCRIPT'
# library(Seurat)
# set.seed(1234)
#
# d <- Read10X("pbmc_citeseq/outs/filtered_feature_bc_matrix/")
# names(d)                        # "Gene Expression" "Antibody Capture"
# pbmc <- CreateSeuratObject(counts = d[["Gene Expression"]])
# pbmc[["HTO"]] <- CreateAssay5Object(counts = d[["Antibody Capture"]])
#
# pbmc <- NormalizeData(pbmc, assay = "HTO",
#                       normalization.method = "CLR")
# pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
# print(table(pbmc$HTO_classification.global))
#
# ## 常規下一步：只留 singlet（但 negative 丟之前先做雷三的檢查！）
# # pbmc <- subset(pbmc, HTO_classification.global == "Singlet")
# saveRDS(pbmc, "output/b04/pbmc_demux.rds")
# RSCRIPT

## ==========================================================================
## [踩雷示範] 以下三段全部預設註解掉。要看後果，逐段解開。
## ==========================================================================

## [踩雷示範] 雷一：aggr 預設 normalize=mapped 把深樣本抽掉近半 reads
##（影片頁 33–34）。跑完比較兩個輸出的 web_summary：
##   mapped 版的 aggregation 表會列出每個樣本被丟棄的比例，
##   最深樣本的 median genes/cell 明顯低於它自己 count 時的值。
# cellranger aggr --id=all_mapped --csv="$OUT_DIR/aggr.csv" \
#   --normalize=mapped
# cellranger aggr --id=all_none   --csv="$OUT_DIR/aggr.csv" \
#   --normalize=none
## 一行檢查（兩個輸出的細胞數與深度對照）：
# grep -o '"median_genes_per_cell":[0-9]*' \
#   all_mapped/outs/metrics_summary_json.json 2>/dev/null || true

## [踩雷示範] 雷二：feature reference 序列打錯一個鹼基 → 整條 ADT 全掛
##（影片頁 35–36）。下面這張表把 CD3 的序列第 14 碼 A 改成 C：
# cat > "$OUT_DIR/feature_ref_typo.csv" <<'EOF'
# id,name,read,pattern,sequence,feature_type
# CD3,CD3_TotalB,R2,5PNNNNNNNNNN(BC),AACAAGACCCTTGCG,Antibody Capture
# EOF
# cellranger count --id=pbmc_typo --transcriptome="$REF" \
#   --libraries="$OUT_DIR/libraries.csv" \
#   --feature-ref="$OUT_DIR/feature_ref_typo.csv" --create-bam=false
## 後果：流程正常結束、零錯誤；web_summary 的 Antibody Reads Usable
## 暴跌、CD3 counts 幾乎歸零。預防：sequence 欄永遠複製貼上＋交叉核對。

## [踩雷示範] 雷三：hashing 的 negative 一律丟掉（影片頁 37–38）。
## negative 是否集中在特定族群，要在 R 端用交叉表檢查——
## 這行檢查是本集「一行檢查」的主角：
# Rscript - <<'RSCRIPT'
# pbmc <- readRDS("output/b04/pbmc_demux.rds")
# ## 假設已有粗略註釋欄 celltype（B10 教怎麼來）：
# ## negative 均勻散在各群＝技術性，可丟；
# ## 集中在單一族群＝生物性訊號，先查染色再決定去留。
# print(table(pbmc$celltype, pbmc$HTO_classification.global))
# RSCRIPT

echo "B04 腳本走完。三張 CSV 範本在 $OUT_DIR/，Cell Ranger 指令見註解。"

## ---- session info（bash 版）---------------------------------------------
# 可重現性紀錄：把工具版本留在輸出目錄。
{ echo "date: $(date -Iseconds)"
  command -v cellranger >/dev/null 2>&1 && cellranger --version || true
  bash --version | head -1
} > "$OUT_DIR/session_info.txt" 2>&1 || true
