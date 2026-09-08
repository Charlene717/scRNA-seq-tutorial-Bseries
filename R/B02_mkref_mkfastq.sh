#!/usr/bin/env bash
# =====================================================================
# B02_mkref_mkfastq.sh — B2「參考基因組與 FASTQ：mkref 與 mkfastq」隨集腳本
#
# 對應影片：scRNA-seq 教學影片系列 · B 系列 · 第 2 集
# 資料：PBMC 1k FASTQ（10x Genomics 官方公開資料，v3 chemistry）
#
# 系統需求（跑不動就只看不跑，B3 起有現成輸出可接）：
#   - Linux x86-64（macOS / Windows 不支援 Cell Ranger）
#   - CPU >= 8 核（官方建議 16）
#   - RAM >= 64 GB（mkref / count 對人類基因組的硬需求）
#   - 硬碟可用空間 >= 數百 GB（BCL、FASTQ、reference、輸出）
#   - Cell Ranger 下載需先在 10x 官網註冊帳號
#
# 使用方式：逐段複製執行（建議），或 bash B02_mkref_mkfastq.sh 全跑。
# 路徑一律相對於專案根目錄；分節標記與影片段落一一對應。
# 版本註記：以 cellranger-9.0.1 為例；請以你下載當下的版本為準。
# =====================================================================
set -euo pipefail

## ---- 0. install ---------------------------------------------------
# （影片頁 11–12）下載、解壓、加 PATH、驗收。
# 10x 官網註冊後會給你一段帶簽名的 curl 指令（連結有時效），形如：
# curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/..."

# tar -xzf cellranger-9.0.1.tar.gz
# export PATH=$PWD/cellranger-9.0.1:$PATH   # 常駐請寫進 ~/.bashrc

cellranger --version          # 預期輸出：cellranger cellranger-9.0.1

# 驗收安裝：用內建迷你資料集跑完整流程（約 5–10 分鐘）
# 你要看的是最後一行 "Pipestance completed successfully!"
# cellranger testrun --id=check_install

## ---- 1. fastq-download --------------------------------------------
# （影片頁 16）抓 PBMC 1k 的官方 FASTQ（約 5 GB，網路差就先跳過本節）
mkdir -p data output
# curl -o data/pbmc_1k_v3_fastqs.tar \
#   https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
# tar -xf data/pbmc_1k_v3_fastqs.tar -C data/

## ---- 2. fastq-anatomy ---------------------------------------------
# （影片頁 16–17）R1 = 16 bp barcode + 12 bp UMI；R2 = cDNA。
FQDIR=data/pbmc_1k_v3_fastqs

# 看 R1 的第一筆紀錄：四行一筆，第二行應為 28 個鹼基
zcat "$FQDIR"/pbmc_1k_v3_S1_L001_R1_001.fastq.gz | head -4

# 讀長檢查（前 1000 筆）：全部 28 → v3/v4；全部 26 → v2
zcat "$FQDIR"/pbmc_1k_v3_S1_L001_R1_001.fastq.gz | head -4000 |
  awk 'NR % 4 == 2 { print length($0) }' | sort | uniq -c

# R2 是 cDNA（約 90 bp），這才是拿去比對基因組的讀段
zcat "$FQDIR"/pbmc_1k_v3_S1_L001_R2_001.fastq.gz | head -4000 |
  awk 'NR % 4 == 2 { print length($0) }' | sort | uniq -c

## ---- 3. fastq-qc --------------------------------------------------
# （影片頁 17）快速品質檢查：筆數是否一致、檔案是否截斷
# R1 與 R2 的 read 數必須完全相同（gzip 壞檔會在這裡露餡）
for f in "$FQDIR"/pbmc_1k_v3_S1_L001_R[12]_001.fastq.gz; do
  echo -n "$f : "
  zcat "$f" | wc -l | awk '{ print $1/4 " reads" }'
done
# gzip 完整性檢查（截斷的下載檔會報錯）
gzip -t "$FQDIR"/*.fastq.gz && echo "gzip OK"
# 想看逐鹼基品質曲線可另跑 FastQC（非必需）：fastqc *.fastq.gz

## ---- 4. mkfastq ---------------------------------------------------
# （影片頁 22–24）多數人拿到的已是 FASTQ，本節僅示範、預設不執行。
# samplesheet 三欄：Lane,Sample,Index（Index 填 10x index「組名」）
cat > output/samplesheet_demo.csv <<'CSV'
Lane,Sample,Index
1,pbmc_1k,SI-TT-A1
CSV

# mkfastq 是 Illumina bcl2fastq 的包裝；10x 已建議改用 BCL Convert。
# 需要整個 BCL run 資料夾（數百 GB），因此註解掉：
# cellranger mkfastq \
#   --run=/data/runs/240801_A00228_0279_BHFWFVDMXX \
#   --csv=output/samplesheet_demo.csv \
#   --output-dir=fastq_out

## ---- 5. ref-prebuilt ----------------------------------------------
# （影片頁 26）預建 reference：人／小鼠、無轉基因時的正解
# curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
# tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
# ls refdata-gex-GRCh38-2024-A
#   fasta/ genes/ star/ reference.json
# reference.json 記著版本來源——方法段直接引用：
# cat refdata-gex-GRCh38-2024-A/reference.json

## ---- 6. mkgtf-mkref -----------------------------------------------
# （影片頁 28–30）自建 reference：以斑馬魚 GRCz11 為例。
# 鐵則：FASTA 與 GTF 必須同一個 Ensembl release（此處都是 110）。
# wget https://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
# wget https://ftp.ensembl.org/pub/release-110/gtf/danio_rerio/Danio_rerio.GRCz11.110.gtf.gz
# gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz \
#        Danio_rerio.GRCz11.110.gtf.gz

# 第一步 mkgtf：biotype 白名單過濾（照 10x 預建的配方）
# cellranger mkgtf \
#   Danio_rerio.GRCz11.110.gtf \
#   Danio_rerio.GRCz11.110.filtered.gtf \
#   --attribute=gene_biotype:protein_coding \
#   --attribute=gene_biotype:lncRNA \
#   --attribute=gene_biotype:IG_C_gene \
#   --attribute=gene_biotype:TR_C_gene

# 第二步 mkref：建 STAR 索引（人類約 1–2 小時、64 GB RAM）
# cellranger mkref \
#   --genome=GRCz11_cr \
#   --fasta=Danio_rerio.GRCz11.dna.primary_assembly.fa \
#   --genes=Danio_rerio.GRCz11.110.filtered.gtf \
#   --nthreads=8 --memgb=64

## ---- 7. transgene -------------------------------------------------
# （影片頁 31–32）加轉基因（GFP）：FASTA 與 GTF「兩邊都要加」。
# 假設 GFP.fa 是一條 720 bp 的序列，FASTA header 為 >GFP

# ① 序列接到基因組 FASTA 之後（多一條迷你染色體）
# cat refdata-gex-GRCh38-2024-A/fasta/genome.fa GFP.fa \
#   > genome_plus_gfp.fa

# ② GTF 加上 GFP 的 exon 條目（gene_id / transcript_id 一致）
# zcat refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz > genes_plus_gfp.gtf
# echo -e 'GFP\tcustom\texon\t1\t720\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";' \
#   >> genes_plus_gfp.gtf

# ③ 重建 reference
# cellranger mkref --genome=GRCh38_gfp \
#   --fasta=genome_plus_gfp.fa --genes=genes_plus_gfp.gtf

# ④ 一行驗證（建完必做、開 count 之前做）：GTF 裡查得到 GFP 嗎？
# grep -c 'gene_id "GFP"' genes_plus_gfp.gtf    # 預期 >= 1

## ---- 8. traps ------------------------------------------------------
# （影片頁 35–40）三個雷的示範，全部預設註解掉，不要照跑。

# [踩雷示範] 雷一：把 R2 餵給 R1 的位置（R1/R2 對調）
# 做法：把 R2 改名成 R1、R1 改名成 R2 再跑 count。
# 後果：valid barcodes 掉到 ~0%，cellranger 報錯中止，訊息點名
#       barcode 對不上 whitelist。
# cp fastq_swap_demo 的建立與執行略——不要在正式資料夾玩這個。

# [踩雷示範] 雷二：檔名不合 bcl2fastq 規範
# mkdir -p fastq_renamed
# cp "$FQDIR"/pbmc_1k_v3_S1_L001_R1_001.fastq.gz fastq_renamed/pbmc1k.R1.fq.gz
# cp "$FQDIR"/pbmc_1k_v3_S1_L001_R2_001.fastq.gz fastq_renamed/pbmc1k.R2.fq.gz
# cellranger count --id=demo_badname --sample=pbmc1k \
#   --fastqs=fastq_renamed --transcriptome=refdata-gex-GRCh38-2024-A \
#   --create-bam=false
# 後果：error: No input FASTQs were found for the requested parameters.
# 解法：改回標準命名（--sample 要與樣本名段一致）
# mv fastq_renamed/pbmc1k.R1.fq.gz fastq_renamed/pbmc1k_S1_L001_R1_001.fastq.gz
# mv fastq_renamed/pbmc1k.R2.fq.gz fastq_renamed/pbmc1k_S1_L001_R2_001.fastq.gz

# [踩雷示範] 雷三：GFP 只加了 FASTA、忘了 GTF
# 做法：mkref 時 --fasta 用 genome_plus_gfp.fa、--genes 用原版 GTF。
# 後果：mkref 與 count 全程無錯誤，GFP 在矩陣裡永遠是 0（安靜的零）。
# 一行檢查（上面第 7 節的 ④）就能在開跑前抓到。

## ---- 9. session ----------------------------------------------------
# 記錄工具版本（可重現性；相當於 R 腳本結尾的 sessionInfo()）
cellranger --version
bash --version | head -1
uname -srm
date
