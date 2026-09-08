# 進階 29 · 同一顆細胞的兩本帳：Multiome RNA + ATAC

**難度**：★★★ ｜ **預估時間**：2 個工作天 ｜ **對應集數**：B4、B18、B19 ｜ **資料集**：10x Multiome「PBMC granulocytes removed (10k)」（本題獨用）

## 背景與研究主題

10x Multiome 在**同一顆細胞核**裡同時量基因表達（RNA）與染色質開放性（ATAC）。這不是把兩份資料對齊——是天生成對，因此可以問單模態問不了的問題：一個基因要表達，它的調控區通常要先開放；那麼「開放但還沒表達」的位置，就是**調控先行**的候選現場。本題走完 Signac 的完整流程：ATAC 有一套跟 RNA 完全不同的 QC 語言、LSI 取代 PCA、然後用 WNN 把兩個 modality 接起來——你在**基礎05** 用 CITE-seq 練過的 WNN，這次把蛋白換成染色質，體會同一框架吃不同 modality 的威力與陷阱。最後用 linkage 分析把 peak 連到基因、用 motif activity 對照 TF 表達，抓「開放性先行於表達」的候選調控。你的研究主題：

1. RNA-only、ATAC-only、WNN 三種分群，各自看到（與看不到）PBMC 的哪些結構？
2. 哪些 peak 與哪些基因的表達顯著連動（linkage）？這些連動在細胞型別間怎麼分佈？
3. 哪些 TF 呈現「motif activity 與自身 RNA 表達脫鉤」的模式？其中有沒有「開放先行、表達未至」的候選調控事件？

## 資料集

- **10x Datasets「PBMC from a healthy donor - granulocytes removed through cell sorting (10k)」**（Multiome ATAC + Gene Expression；11,898 cells 估計）。
- 下載（10x datasets 頁搜尋名稱，填 email 免費）至少三個檔：`filtered_feature_bc_matrix.h5`（RNA + ATAC peak counts 同一個 h5）、`atac_fragments.tsv.gz`、以及它的 **index（.tbi）**——fragments 沒有 index，Signac 大半功能直接罷工。
- 需另裝 **Signac**（資料集總覽有註明）＋註解資源（EnsDb 人類註解、JASPAR motif 資料庫，見提示）。

## 真實數據關卡

真實資料在你跑第一行 Seurat 之前就開始出題。這一節照做並記錄：

1. **一個 h5、兩個 modality**：`Read10X_h5` 讀進來的是一個 list（`Gene Expression` 與 `Peaks` 兩個矩陣）——先各看 `dim()` 與 rownames 長相：RNA 是 gene symbol，ATAC 是 `chr-start-end` 的 peak 座標字串。兩個矩陣的欄（barcode）是同一批細胞——驗證這件事（交集比例），這是「同細胞雙模態」的立足點。fragments 檔用 zcat 偷看前幾列：每列是一個片段（染色體、起訖、barcode、讀數），這才是 ATAC 的原始層，peak 矩陣只是它的一種摘要。
2. **ATAC 的 QC 是另一種語言**：nCount 與 mito% 的直覺只夠用一半。ATAC 核心 QC 是 **TSS enrichment**（片段在轉錄起始位點附近的富集度——訊噪比指標）與 **nucleosome signal**（單核小體 vs 次核小體片段長度比——片段長度分佈該有核小體週期性）；再加 fragments in peaks 比例。每個指標查 Signac 文件搞懂它量什麼、為什麼低（或高）是壞事，閾值三件事照寫——**照抄 RNA QC 的人這題直接死在起點**。
3. **雙過濾的交集**：RNA 過濾與 ATAC 過濾各自會殺掉一批細胞——最終物件只能留**兩邊都過**的交集。記錄三個數字：只死於 RNA 標準的、只死於 ATAC 標準的、雙雙陣亡的。一個 RNA 完美但 ATAC 稀爛的 barcode 說明了什麼（核質量？技術層面哪一步失敗）？想一想寫進筆記。
4. **決策日誌**：開 `decisions.md`，凡是「資料逼你做的決定」（h5 解剖發現、ATAC QC 閾值與依據、交集過濾數字、註解版本選擇）都記一筆：發現了什麼 → 選項 → 你選了什麼 → 理由。這是繳交物，也是論文 Methods 的草稿。

## 任務

### 階段 A：建物件與雙 QC（B4）

1. 建含 RNA assay + ChromatinAssay 的 Seurat 物件（ChromatinAssay 掛 fragments 路徑與基因註解），算齊 ATAC QC 指標（TSS enrichment、nucleosome signal、fragments in peaks）與 RNA QC 指標。
2. 執行雙過濾與交集，三個死亡數字進 decisions.md。
3. 各 modality 各自處理：RNA 走標準流程到 PCA；ATAC 走 TF-IDF → `FindTopFeatures` → SVD（合稱 LSI），**檢查第 1 個 LSI 成分與測序深度的相關性，高相關就丟棄**（`DepthCor` 圖為證）——這是 B4「不同 modality 有不同統計性質」在染色質上的具體化。

### 階段 B：WNN 聯合分群與註解（B18；呼應基礎05）

4. `FindMultiModalNeighbors`（PCA + LSI，各自 dims 有理由）→ wsnn 分群 → WNN UMAP。三套分群（RNA-only、ATAC-only、WNN）列聯表互相對照。
5. 註解 WNN 分群：RNA marker 為主、ATAC 的 gene activity（`GeneActivity`）為輔，至少分出 CD4/CD8 T、NK、B（naive/memory）、單核球（classical/non-classical）、DC。與 Azimuth PBMC reference 的自動註解對答案（自行查 Azimuth 用法——PBMC 是型別標準答案最齊的組織，這是本題選它的原因）。
6. 畫各群的 modality weight（呼應基礎05 的 RNA.weight 分析）：哪些型別更靠染色質分、哪些更靠表達分？跟基礎05 的 CITE-seq 版結論並排，寫 200 字：「第二 modality 換成 ATAC 之後，買到的解析度跟蛋白買到的有什麼質的不同？」

### 階段 C：調控先行的偵探工作（發表導向；B19）

7. **Linkage**：`RegionStats` → `LinkPeaks`，對你感興趣的型別 marker 基因集找顯著連動的 peak；挑 2–3 個案例畫 `CoveragePlot`（軌跡圖上同時看開放性、連動弧線與表達）。
8. **Motif activity vs TF 表達**：跑 chromVAR 得每細胞 motif activity，挑 10–20 個免疫關鍵 TF，對每個 TF 畫「motif activity vs 自身 RNA 表達」的型別層級對照。分類：兩者同步的、activity 高但表達低的（開放先行？母系 TF 家族 motif 混淆？）、表達高但 activity 平的。挑一個「開放先行」候選講出一個可檢驗的調控故事。
9. 穩健性檢查（換 LSI dims、換 linkage 距離窗、motif 家族混淆的討論）＋寫「發表路徑」評估（見下節，300–500 字）。

## 繳交物

1. 可重跑的 R 專案（renv + set.seed(1234)、相對路徑）＋ `decisions.md` 決策日誌（含雙過濾死亡數字）。
2. 圖：ATAC QC 圖組（TSS enrichment、片段長度分佈）、DepthCor 圖、三套分群 UMAP + 列聯表、modality weight 圖、2–3 張 CoveragePlot、TF activity vs 表達對照圖（英文標籤）。
3. 分析筆記：各階段回答＋每個參數三件事＋基礎05 對照短文。
4. 發表路徑評估（300–500 字）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- 你的「開放先行」候選屬於指南第二節的哪一類？它在換參數（LSI dims、linkage 窗）後還在嗎？motif 家族混淆排除到什麼程度？
- **驗證設計**：型別註解的驗證用 **Azimuth PBMC reference**（自行查證用法與版本）；WNN 框架本身的行為驗證回到**基礎05**（bmcite CITE-seq）——同一套 WNN、不同第二 modality，權重結構的異同就是現成的方法學對照。調控候選的外部驗證可查公開 PBMC multiome/ATAC 資料（自行 GEO 查證），健康 PBMC 的資源不缺。
- **novelty 定位**：健康 PBMC 的 multiome 是方法學示範常客，生物 novelty 天花板低——本題更可能的發表形態是方法比較或教學資源型（指南第四節第三層）；若你想往生物發現走，誠實評估「健康單一供體」的限制。
- 若要成文：缺哪一塊（疾病對照組？多供體？擾動實驗佐證調控假說？）——寫出把這條 pipeline 搬到哪個疾病資料集才有故事。

## 自我檢核點

- [ ] fragments + index 齊備，h5 雙矩陣與 barcode 交集驗證過
- [ ] ATAC QC 三指標能各自說出「量什麼、為什麼是這個閾值」，不是抄 RNA 的
- [ ] LSI 第 1 成分與深度的相關性檢查過且有處置
- [ ] 雙過濾交集的三個死亡數字進了 decisions.md
- [ ] linkage 與 motif 結論有 CoveragePlot／表達層證據，母系 motif 混淆有討論
- [ ] 發表路徑評估誠實面對健康單供體的 novelty 天花板

## 提示（卡關再看）

<details><summary>提示 1：建 ChromatinAssay</summary>
`counts <- Read10X_h5(...)`; RNA 用 `counts$'Gene Expression'` 建 assay；ATAC 用 `CreateChromatinAssay(counts = counts$Peaks, sep = c("-","-"), fragments = "atac_fragments.tsv.gz", annotation = ...)`。annotation 從 `EnsDb.Hsapiens.v86` 取（`GetGRangesFromEnsDb`），染色體命名風格（chr1 vs 1）不合時用 `seqlevelsStyle` 統一——這是本題經典翻車點之一。.tbi 要跟 fragments 同資料夾同名。
</details>

<details><summary>提示 2：ATAC 流程骨架</summary>
`RunTFIDF()` → `FindTopFeatures(min.cutoff = "q0")` → `RunSVD()`；`DepthCor(obj)` 看各成分與 nCount 的相關性，慣例上丟棄成分 1、用 2:30——但把「你自己的 DepthCor 圖」當證據，不是抄慣例。WNN 呼叫形如 `FindMultiModalNeighbors(obj, reduction.list = list("pca","lsi"), dims.list = list(1:30, 2:30))`。
</details>

<details><summary>提示 3：chromVAR 與 motif 混淆</summary>
`AddMotifs`（JASPAR2020，人類 collection）→ `RunChromVAR`。解讀最大陷阱：同家族 TF（如 FOS/JUN 家族、GATA 家族）motif 幾乎相同，activity 分不出是哪個成員在做事——所以才要對照「該 TF 自身的 RNA 表達」：家族裡只有一個成員在該型別表達，指認才有底氣。把這個推理過程寫進筆記。
</details>

## 進階挑戰

- 對一條你有把握的分化軸（如 naive → memory B）檢查 peak 開放與基因表達的先後：沿偽時間（或型別序）畫兩者的曲線，開放是否真的先行？
- 用 `FindMarkers` 對 ATAC assay 找型別特異 peak（differentially accessible regions），跟 RNA marker 的資訊量比一比。
- 想把這套流程用在疾病上？進階05（AML）的惡性分化阻滯是天然候選——寫 200 字設計一個「AML multiome 能回答、純 RNA 回答不了」的問題。

## 參考文獻

- 10x Genomics. PBMC from a healthy donor - granulocytes removed through cell sorting (10k), Single Cell Multiome ATAC + Gene Expression dataset.（10x datasets 頁）
- Stuart T, et al. Single-cell chromatin state analysis with Signac. *Nature Methods* (2021).（工具）
- Schep AN, et al. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nature Methods* (2017).（工具）
- Hao Y, et al. Integrated analysis of multimodal single-cell data. *Cell* (2021).（WNN 方法；Azimuth reference）
