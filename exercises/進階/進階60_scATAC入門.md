# 進階 60 · 沒有表達值的單細胞：scATAC 入門

**難度**：★★★ ｜ **預估時間**：2 個工作天 ｜ **對應集數**：B7、B8、B18、B19 ｜ **資料集**：10x「10k Human PBMCs, ATAC v1.1, Chromium X」（本題獨用）

## 背景與研究主題

這一題把「表達值」整個拿走。scATAC-seq 量的是染色質開放性——每顆細胞核裡，基因組的哪些區域是打開的。你面對的不再是基因 × 細胞的表達矩陣，而是 peak（基因組區間）× 細胞的事件計數：**比 RNA 再稀疏一個量級**（一顆二倍體細胞每個位置最多兩份拷貝，開或不開近乎二元），而且 feature 不是基因名，是 `chr1-1234-5678` 這種座標字串。整條 B 系列流程要逐站重審：QC 換一套語言（TSS enrichment、nucleosome signal）、正規化與降維換 LSI（TF-IDF + SVD，B7 的 PCA 直覺只能部分平移——第一成分丟不丟是本題的招牌決策）、註解沒有 marker 基因可直接看，要靠 gene activity 近似再借你在基礎01 練出的 PBMC RNA 知識當先驗。做完這題再去進階29（multiome），你會知道同細胞 RNA 是多大的奢侈。你的研究主題：

1. ATAC 的稀疏性與近二元性實際長什麼樣？它如何改寫 QC、正規化、降維的每個決策？
2. 只靠染色質開放性，PBMC 的型別結構能分到多細？哪些型別在 ATAC 空間反而比 RNA 空間更清楚（提示：分化譜系的調控狀態）？
3. gene activity 近似的極限在哪——哪些型別它註解得動、哪些它含糊，為什麼？

## 資料集

- **10x Datasets「10k Human PBMCs, ATAC v1.1, Chromium X」**（scATAC-seq v1.1；9,030 nuclei）。
- 下載（10x datasets 頁搜尋名稱，免註冊）三個檔：`fragments.tsv.gz`（**2.6 GB**）與其 **.tbi index**（頁面有就一起抓，沒有就自己 tabix 建——見真實數據關卡）、`filtered_peak_bc_matrix.h5`（130 MB）、`singlecell.csv`（per-barcode QC 指標）。磁碟預留 15–20 GB（資料集總覽）。
- 需另裝 **Signac**（資料集總覽有註明）＋ EnsDb 人類基因註解（見提示）。

## 真實數據關卡

真實資料在你跑第一行 Signac 之前就開始出題。這一節照做並記錄：

1. **fragments 不是矩陣，是原始事件流**：`zcat fragments.tsv.gz | head` 偷看——每列一個片段：染色體、起、訖、barcode、讀數。peak 矩陣只是把這些事件對 peak 區間做的一種摘要；TSS enrichment、nucleosome signal、gene activity 全都要回 fragments 現算，所以 2.6 GB 逃不掉。用 `wget -c` 斷點續傳；下載完驗檔案大小。**index 是生死線**：`.tbi` 要跟 fragments 同資料夾同名，沒有 index 的話 Signac 大半功能直接罷工（自建：`tabix -p bed fragments.tsv.gz`，需要 htslib；把你走的路記下來）。
2. **h5 解剖與稀疏性測量**：讀入 peak 矩陣後親手量化：非零元素比例是多少？跟你基礎01 的 RNA 矩陣比差幾倍？非零值的分布長怎樣（大多是 1 和 2——近二元性的直接證據）？把這兩個數字寫進筆記——它們是你之後每個「為什麼 ATAC 不能照抄 RNA 流程」論述的證據。peak 有幾萬個、rownames 是座標字串——feature 數比 RNA 多得多，這對 HVF 選擇（`FindTopFeatures` 的邏輯 vs RNA 的 HVG）意味著什麼？
3. **ATAC QC 是全新一套**：`singlecell.csv` 是 Cell Ranger ATAC 的 per-barcode 指標表——先逐欄搞懂再用。核心 QC：**TSS enrichment**（訊噪比——片段在轉錄起始位點附近該富集）、**nucleosome signal**（片段長度該有核小體週期性）、fragments in peaks 比例、blacklist 比例。每個指標查 Signac 文件弄懂它量什麼、為什麼壞值代表什麼技術失敗，閾值三件事照寫——nCount 與 mito% 的直覺在這裡只夠用一半。
4. **決策日誌**：開 `decisions.md`，凡是「資料逼你做的決定」（index 路線、稀疏性測量結果、QC 閾值與依據、LSI 成分取捨）都記一筆：發現了什麼 → 選項 → 你選了什麼 → 理由。這是繳交物，也是論文 Methods 的草稿。

## 任務

### 階段 A：建物件、QC、LSI（B7）

1. 建 ChromatinAssay 物件（掛 fragments 路徑與基因註解），算齊 QC 指標（TSS enrichment、nucleosome signal、fragments in peaks），畫片段長度分佈圖與 TSS enrichment 圖——這兩張是 ATAC 的「QC 名片」，要能對別人解釋圖上每個特徵（核小體週期的鋸齒從哪來）。
2. 設閾值過濾，前後細胞數進 decisions.md。
3. `RunTFIDF` → `FindTopFeatures` → `RunSVD`（合稱 LSI）。**招牌決策**：`DepthCor` 檢查各成分與測序深度的相關性——第 1 成分丟不丟？用你自己的圖回答，不是抄慣例；並跟 B7 的 elbow 思維對照：LSI 選成分數的邏輯哪裡像 PCA、哪裡不像？

### 階段 B：分群與 gene activity 註解（B8、B18）

4. LSI 空間分群 + UMAP（resolution 三件事照 B8 的紀律）。先誠實記錄：這時候你手上只有「群 0、群 1…」，連一個基因名都還沒有。
5. `GeneActivity` 算基因活性近似（把基因體＋上游區間的片段數當表達的代理），建成第二個 assay。用你在基礎01 練熟的 PBMC marker（CD3E、MS4A1、CD14、NKG7…那一套先驗）在 gene activity 上註解各群——至少分出 T（CD4/CD8 分不分得動？誠實寫）、NK、B、單核球、DC。
6. 評估 gene activity 的極限：哪些 marker 在 activity 層級訊號清楚、哪些糊掉（想想這個近似假設了什麼——開放 ≠ 表達）？挑一個含糊案例寫 150 字分析。

### 階段 C：ATAC 才看得到的東西（發表導向；B19）

7. 對你有把握的兩個型別做 differentially accessible regions（`FindMarkers` 於 peak assay，注意 ATAC 建議的檢定設定，查 Signac 文件），挑代表 peak 畫 `CoveragePlot`；找 2–3 個「開放性差異明顯、但基礎01 的 RNA 資料裡對應基因表達差異平平」的區域——調控層先於（或獨立於）表達層的候選。
8. 用 label transfer 借 RNA 資料當外援：拿一份帶註解的 PBMC RNA 參考（你基礎01 的物件即可）對 gene activity 做 Seurat label transfer，跟你的手動註解對照——不同細胞、不同 modality 的橋接，先天限制是什麼？寫進筆記（這正是進階40 的主題預告）。
9. 穩健性（換 LSI 成分數、換 `FindTopFeatures` 門檻、換 resolution）＋寫「發表路徑」評估（見下節，300–500 字）。

## 繳交物

1. 可重跑的 R 專案（renv + set.seed(1234)）＋ `decisions.md` 決策日誌（含稀疏性測量數字）。
2. 圖：片段長度分佈與 TSS enrichment 圖、DepthCor 圖、UMAP（註解前後）、gene activity dotplot、2–3 張 CoveragePlot、label transfer 對照表（英文標籤）。
3. 分析筆記：各階段回答＋每個參數三件事＋gene activity 極限的 150 字。
4. 發表路徑評估（300–500 字）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- 你的「開放性差異但表達平平」候選屬於指南第二節的哪一類？它在換 LSI 成分數與 peak 檢定設定後還在嗎？
- **驗證設計**：天然下一步就在題庫——**進階29**（10x Multiome，RNA+ATAC 同細胞）：你的候選區域在同細胞資料裡可以直接驗「開放但不表達」是真解耦還是橋接誤差；**進階40**（婦科雙癌 scRNA + scATAC 非同細胞）則是本題橋接手法的疾病版壓力測試。
- **novelty 定位**：健康 PBMC scATAC 是工具論文的標準展示品，生物 novelty 天花板低——合理的發表形態是方法比較或教學資源型（指南第四節第三層）；查 Signac/ArchR 等工具文獻確認你的比較角度沒被官方 vignette 做完。
- 若要成文：缺什麼（同細胞 multiome 驗證？疾病對照？多樣本）？誠實寫出把這套流程搬到哪個疾病場景才有生物故事。

## 自我檢核點

- [ ] fragments + index 齊備，稀疏性與近二元性親手量化過且數字進了筆記
- [ ] ATAC QC 指標能各自說出「量什麼、壞值代表什麼」，不是抄 RNA 的
- [ ] LSI 第 1 成分的去留用自己的 DepthCor 圖回答
- [ ] 註解順序誠實：先分群、後 gene activity、再對照 label transfer
- [ ] 開放性 vs 表達的候選有 CoveragePlot 證據，且知道要去進階29 驗什麼
- [ ] 發表路徑評估誠實面對健康 PBMC 的 novelty 天花板

## 提示（卡關再看）

<details><summary>提示 1：建物件與註解資源</summary>
`counts <- Read10X_h5("filtered_peak_bc_matrix.h5")`；`CreateChromatinAssay(counts, sep=c(":","-") 或 c("-","-") 依 rownames 實際格式, fragments="fragments.tsv.gz", annotation=...)`。annotation 用 `GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)`，染色體命名風格不合時 `seqlevelsStyle(annotation) <- "UCSC"`——這是經典翻車點。`singlecell.csv` 用 read.csv 讀進來後依 barcode 對齊塞進 meta.data。
</details>

<details><summary>提示 2：LSI 與第一成分</summary>
`RunTFIDF()` → `FindTopFeatures(min.cutoff="q0")` → `RunSVD()` → `DepthCor(obj)`。第 1 成分通常與 nCount 相關性極高（TF-IDF 沒能完全消掉深度效應），慣例用 2:30——但你的任務是畫圖、看數字、寫三件事，把「慣例」變成「證據支持的決策」。順手想：RNA 的 PC1 有沒有類似問題？為什麼 B7 沒教你丟 PC1？
</details>

<details><summary>提示 3：label transfer 的正確期待</summary>
`FindTransferAnchors(reference=RNA物件, query=ATAC物件, reduction="cca")` 時 query 端用 gene activity assay——兩邊「假裝」在同一個特徵空間。它給的是機率性標籤，粗型別通常可靠、細亞群（naive vs memory）常糊。跟你的手動註解列聯表對照：兩邊都心虛的群，就標「未定」，這比硬給名字誠實——參考 mapping 的謙卑是 B18 圖譜時代的核心素養。
</details>

## 進階挑戰

- 用 chromVAR 加 motif 分析（進階29 的提示 3 可借），看各型別的 TF motif activity——不靠 gene activity 的第二條註解路線，兩條路線互相印證。
- 從 fragments 自己重新 call peak（MACS2，按型別分開 call），跟 Cell Ranger 的統一 peak set 比：型別特異的 peak 被統一 call 吃掉多少？這是 ATAC 分析的已知系統性議題。
- 做完進階29 後回來重看本題第 7 步的候選區域清單，寫 200 字：同細胞證據改判了幾個？「非同細胞橋接」的誤差率給你什麼教訓？

## 參考文獻

- 10x Genomics. 10k Human PBMCs, ATAC v1.1, Chromium X dataset.（10x datasets 頁）
- Stuart T, et al. Single-cell chromatin state analysis with Signac. *Nature Methods* (2021).（工具）
- Cusanovich DA, et al. Multiplex single-cell profiling of chromatin accessibility by combinatorial cellular indexing. *Science* (2015).（scATAC 方法源流）
- 10x Genomics. PBMC granulocytes removed (10k) Multiome dataset.（同細胞驗證場，見進階29）
