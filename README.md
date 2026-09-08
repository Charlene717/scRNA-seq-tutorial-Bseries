# scRNA-seq B 系列 · 實作篇

> **🚧 這個倉庫正在建置中。** 投影片、隨集腳本、練習題庫與參考庫都已經齊備、可以直接使用；影片正在製作，上架後會把連結陸續補進下面的表格。
> 內容也會隨著課程改版持續更新。有任何疑問、勘誤或建議，歡迎開 [Issue](../../issues) 告訴我。

> 十九集課程、983 頁投影片、十九支隨集腳本、80 題實作練習、四座參考庫——從 FASTQ 與 Cell Ranger 一路做到富集、通訊、軌跡、惡性判定與反卷積，把一份 scRNA-seq 分析做到能寫進論文。

![code: MIT](https://img.shields.io/badge/code-MIT-blue.svg)
![content: CC BY--NC 4.0](https://img.shields.io/badge/content-CC%20BY--NC%204.0-lightgrey.svg)
![R ≥ 4.3](https://img.shields.io/badge/R-%E2%89%A54.3-276DC3.svg)
![Seurat v5](https://img.shields.io/badge/Seurat-v5-1a6b5a.svg)
![語言: 中英雙語](https://img.shields.io/badge/%E8%AA%9E%E8%A8%80-%E4%B8%AD%E8%8B%B1%E9%9B%99%E8%AA%9E-1a6b5a.svg)

這是「怎麼做」的完整縱貫線。主線用 10x 官方的 PBMC 公開資料貫穿，每一步都可以在自己的電腦上重跑；下游應用段各自帶最合適的公開資料集。課程的核心主張是：**每個參數都是一個假設，要有數字、有理由、有紀錄**——`set.seed(1234)`，所有腳本從專案根目錄相對路徑執行，不用 `setwd()`。

**兩件別的課通常不做的事**

1. **每一集都刻意做錯一次。** 「⑤ 踩雷區」用實跑的輸出示範「這個參數設錯會看到什麼」——不是講原理，是把災難跑給你看。
2. **每一集都有「⑤+ 複雜樣本策略室」。** 同一個步驟，資料換成腫瘤、發炎組織、多區域取樣或血液腫瘤時，決策會怎麼變。標準流程的教學多半預設一份乾淨的樣本；真實資料不是那樣，而破掉的往往正是這一步依賴的假設。

## 課程一覽

### B-I 上游：Cell Ranger（B1–B4）

| # | 標題 | 頁數 | 時長 | 影片 | 投影片 | 腳本 |
|---|---|---|---|---|---|---|
| B1 | 環境建置與你的第一張 UMAP | 47 | ~28 分 | 製作中 | [中文](slides/B1_環境建置與第一張UMAP_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B1_Setup_And_First_UMAP_Slides_EN.pdf) | [`B01`](R/B01_first_umap.R) |
| B2 | 參考基因組與 FASTQ：`mkref` 與 `mkfastq` | 51 | ~36 分 | 製作中 | [中文](slides/B2_參考基因組與FASTQ_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B2_Reference_And_FASTQ_Slides_EN.pdf) | [`B02`](R/B02_mkref_mkfastq.sh) |
| B3 | `cellranger count` 實戰與 `web_summary` 判讀 | 51 | ~40 分 | 製作中 | [中文](slides/B3_cellranger_count實戰_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B3_CellRanger_Count_Slides_EN.pdf) | [`B03`](R/B03_count.sh) |
| B4 | 多樣本與多模態：`aggr`、`multi`、Feature Barcode | 50 | ~36 分 | 製作中 | [中文](slides/B4_多樣本與多模態_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B4_Aggr_Multi_FeatureBarcode_Slides_EN.pdf) | [`B04`](R/B04_aggr_multi.sh) |

### B-II 下游：Seurat 核心流程（B5–B10）

| # | 標題 | 頁數 | 時長 | 影片 | 投影片 | 腳本 |
|---|---|---|---|---|---|---|
| B5 | 讀入資料與 QC：閾值到底怎麼定 | 52 | ~39 分 | 製作中 | [中文](slides/B5_讀入資料與QC_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B5_Load_And_QC_Slides_EN.pdf) | [`B05`](R/B05_qc.R) |
| B6 | Normalization、HVG、Scaling：LogNormalize vs SCTransform 實測 | 55 | ~38 分 | 製作中 | [中文](slides/B6_Normalization與HVG_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B6_Normalization_HVG_Scaling_Slides_EN.pdf) | [`B06`](R/B06_normalization.R) |
| B7 | PCA 與維度選擇：到底取幾個 PC | 48 | ~32 分 | 製作中 | [中文](slides/B7_PCA與維度選擇_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B7_PCA_Dimensionality_Slides_EN.pdf) | [`B07`](R/B07_pca.R) |
| B8 | 分群：resolution 掃描與 clustree | 53 | ~38 分 | 製作中 | [中文](slides/B8_分群與resolution_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B8_Clustering_Resolution_Slides_EN.pdf) | [`B08`](R/B08_clustering.R) |
| B9 | UMAP / t-SNE：參數、隨機性與可重現性 | 46 | ~31 分 | 製作中 | [中文](slides/B9_UMAP與可重現性_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B9_UMAP_Reproducibility_Slides_EN.pdf) | [`B09`](R/B09_umap.R) |
| B10 | 細胞註釋：marker、`FindAllMarkers`、SingleR / Azimuth | 55 | ~42 分 | 製作中 | [中文](slides/B10_細胞註釋_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B10_Cell_Annotation_Slides_EN.pdf) | [`B10`](R/B10_annotation.R) |

### B-III 整合與差異表達（B11–B12）

| # | 標題 | 頁數 | 時長 | 影片 | 投影片 | 腳本 |
|---|---|---|---|---|---|---|
| B11 | 多樣本整合與批次效應診斷 | 56 | ~40 分 | 製作中 | [中文](slides/B11_資料整合_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B11_Integration_Slides_EN.pdf) | [`B11`](R/B11_integration.R) |
| B12 | 差異表達的正確做法：`FindMarkers` 的陷阱與 pseudobulk + DESeq2 | 55 | ~43 分 | 製作中 | [中文](slides/B12_差異表達_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B12_Differential_Expression_Slides_EN.pdf) | [`B12`](R/B12_pseudobulk_de.R) |

### B-IV 下游應用與策略（B13–B19）

| # | 標題 | 頁數 | 時長 | 影片 | 投影片 | 腳本 |
|---|---|---|---|---|---|---|
| B13 | DE 之後：功能富集與基因集分析 | 54 | ~41 分 | 製作中 | [中文](slides/B13_功能富集與基因集分析_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B13_Enrichment_GSEA_Slides_EN.pdf) | [`B13`](R/B13_enrichment.R) |
| B14 | 細胞通訊：配體–受體配對與 CellChat | 55 | ~41 分 | 製作中 | [中文](slides/B14_細胞通訊_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B14_Cell_Communication_Slides_EN.pdf) | [`B14`](R/B14_cellchat.R) |
| B15 | 軌跡分析與擬時序：細胞狀態的時間軸 | 52 | ~40 分 | 製作中 | [中文](slides/B15_軌跡分析與擬時序_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B15_Trajectory_Pseudotime_Slides_EN.pdf) | [`B15`](R/B15_trajectory.R) |
| B16 | 腫瘤場景：inferCNV 與惡性細胞判定 | 51 | ~38 分 | 製作中 | [中文](slides/B16_inferCNV與惡性細胞_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B16_InferCNV_Malignant_Slides_EN.pdf) | [`B16`](R/B16_infercnv.R) |
| B17 | 反卷積：用單細胞參考拆解 bulk 資料 | 50 | ~39 分 | 製作中 | [中文](slides/B17_反卷積_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B17_Deconvolution_Slides_EN.pdf) | [`B17`](R/B17_deconvolution.R) |
| B18 | 生態系與收官：多組學、空間、Python 世界與下一步 | 51 | ~33 分 | 製作中 | [中文](slides/B18_生態系與收官_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B18_Ecosystem_Finale_Slides_EN.pdf) | [`B18`](R/B18_ecosystem.R) |
| B19 | 複雜樣本的分析策略設計：不同腫瘤、不同戰法 | 51 | ~36 分 | 製作中 | [中文](slides/B19_複雜樣本策略設計_投影片_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/B19_Complex_Sample_Strategy_Slides_EN.pdf) | [`B19`](R/B19_strategy.R) |
| — | 總結卡與重點卡（十九集的一頁式整理） | 36 | — | — | [中文](slides/B系列_總結卡與重點卡_ZH.pdf)&nbsp;｜&nbsp;[EN](slides/SeriesB_Summary_Cards_EN.pdf) | — |

投影片**中英雙語各一套**（PDF），十九集共 983 頁，另有 36 頁總結卡。全系列約 11.8 小時。

> **建議觀看順序**：B17 → B19 → B18。B19 是策略課，放在反卷積之後、收官之前；B18 保持系列收尾。

## 適合誰

| 你想要 | 看 | 做 | 動手時間預估 |
|---|---|---|---|
| 從原始資料自己跑一份完整分析 | B1、B5–B10 | 腳本 01、05–10 | 2–4 天 |
| 兩組比較、做到可發表的差異表達 | 加 B11–B13 | 加腳本 11–13 | 再 2–3 天 |
| 論文裡那些下游圖（通訊、軌跡、CNV、反卷積） | 加 B14–B18 | 加腳本 14–18 | 再 4–7 天 |
| 自己跑上游、處理多樣本與多模態 | 加 B2–B4 | 加 `B02`–`B04` 指令稿 | 視資料與機器而定 |
| 面對腫瘤等複雜樣本設計自己的策略 | 全系列 + B19 | 全部腳本 + [`exercises/`](exercises/) | 挑一條主線做完，數週 |

**時間預估怎麼看**：估的是動手的時間，不含看影片，也假設你會停下來想每個參數的理由並寫下「數字、理由、紀錄」。只求「跑完不報錯」會快很多，但那不是這門課要教的。套件安裝與資料下載的等待另計。

前置需求：跑過幾行 R、會開 RStudio。不需要單細胞經驗——沒有的話，先看 [A 系列 · 入門篇](https://github.com/Charlene717/scRNA-seq-tutorial-Aseries) 會順很多。

## 隨集腳本（R/）

十九支腳本，與投影片逐段對應：開頭註明對應集數、資料與環境需求，分節標記對得上影片段落。

| 腳本 | 做什麼 | 集 |
|---|---|---|
| [`B01_first_umap.R`](R/B01_first_umap.R) | 15 行從 `Read10X` 到 `DimPlot`——先給一次完整的成功經驗 | B1 |
| [`B02_mkref_mkfastq.sh`](R/B02_mkref_mkfastq.sh) | 參考基因組自建與 FASTQ 產生（Cell Ranger 指令稿） | B2 |
| [`B03_count.sh`](R/B03_count.sh) | `cellranger count` 與 `web_summary` 判讀 checklist | B3 |
| [`B04_aggr_multi.sh`](R/B04_aggr_multi.sh) | `aggr` / `multi` / Feature Barcode 的設定檔與指令 | B4 |
| [`B05_qc.R`](R/B05_qc.R) | 建物件 → 三大 QC 指標 → 分布定閾值（含 MAD）→ DoubletFinder | B5 |
| [`B06_normalization.R`](R/B06_normalization.R) | LogNormalize 三步 vs SCTransform 一行，兩條路線對決 | B6 |
| [`B07_pca.R`](R/B07_pca.R) | PCA、ElbowPlot / JackStraw / 下游穩定度三種選 nPC 的方式 | B7 |
| [`B08_clustering.R`](R/B08_clustering.R) | `FindNeighbors` / `FindClusters`、resolution 掃描與 clustree | B8 |
| [`B09_umap.R`](R/B09_umap.R) | UMAP / t-SNE 參數、seed 與可重現性檢查 | B9 |
| [`B10_annotation.R`](R/B10_annotation.R) | marker 面板 → `FindAllMarkers` → SingleR / Azimuth → 交叉驗證 | B10 |
| [`B11_integration.R`](R/B11_integration.R) | `IntegrateLayers` 與 Harmony、正負對照與 LISI 診斷 | B11 |
| [`B12_pseudobulk_de.R`](R/B12_pseudobulk_de.R) | cell-level DE 的假陽性示範 → pseudobulk + DESeq2 | B12 |
| [`B13_enrichment.R`](R/B13_enrichment.R) | ORA vs GSEA、Hallmark / GO / KEGG、NES 熱圖與 dotplot | B13 |
| [`B14_cellchat.R`](R/B14_cellchat.R) | CellChat 六種圖、兩條件比較與假陽性控制 | B14 |
| [`B15_trajectory.R`](R/B15_trajectory.R) | Slingshot + tradeSeq、起點選擇與軌跡上的差異檢定 | B15 |
| [`B16_infercnv.R`](R/B16_infercnv.R) | inferCNV 參考細胞選擇、CNV 熱圖判讀與惡性判定 | B16 |
| [`B17_deconvolution.R`](R/B17_deconvolution.R) | 建參考矩陣 → NNLS / MuSiC 反卷積 → 與已知比例對答案 | B17 |
| [`B18_ecosystem.R`](R/B18_ecosystem.R) | 多模態、空間、Python 生態的接口示範 | B18 |
| [`B19_strategy.R`](R/B19_strategy.R) | 兩個小型模擬：同一套流程套在不同資料結構上，哪一步先失效 | B19 |

資料下載腳本在 [`data/download_data.R`](data/download_data.R)。**本倉庫不含任何資料檔與分析輸出**（見 `.gitignore`）。

建議做法：新建一個 RStudio Project，把 `R/` 與 `data/` 複製進去，先跑 `data/download_data.R`，再從 `B01_first_umap.R` 開始。全程不要用 `setwd()`。

## 資料

| 資料 | 內容 | 用在 |
|---|---|---|
| PBMC 3k（10x 官方） | 一份健康捐贈者的 PBMC，2,700 顆 | B1、B5–B10、B17 |
| PBMC 1k FASTQ（10x 官方） | 小檔 FASTQ，筆電跑得動 | B2–B3 |
| ifnb（Kang et al. 2018） | IFN-β 刺激 vs 對照，雙條件 | B11–B14 |
| 模擬分化資料 | 用真實演算法產生、有已知答案的軌跡 | B15 |
| inferCNV 內建範例 | 含正常參考細胞的腫瘤資料 | B16 |
| 兩個小型模擬 | 不同資料結構的策略對照 | B19 |

## 自測題庫

114 題（每集 6 題），中英對照，含單選、多選與是非。互動題庫在 [`quiz/`](quiz/)：把整個資料夾下載到電腦後，點擊 `index.html` 即可開啟作答——即點即答、附解析、可切換中英文、可隨機抽題、進度自動保存。不需要架站也不需要連網。

## 實作練習題（exercises/）

課程之外的動手關卡，共 **80 題**，全部放在 [`exercises/`](exercises/)：

- **基礎 20 題**——把課程裡的某一段流程換一份資料獨立走完，半天到一天一題。分兩圈：第一圈練標準流程（換組織、換物種、換規模、有標準答案的分群），第二圈練核心技能（雙條件比較、軌跡、CITE-seq、doublet 與 ambient、多技術整合、hashing、snRNA、h5ad 轉檔）。
- **進階 60 題**——**60 個各不相同的公開資料集**，從讀檔到臨床驗證走完整條線，一題約 2–3 週。涵蓋二十餘種癌別、非腫瘤疾病、發育與圖譜、空間與多模態。**挑一題你想做的做完就好。**

進階題比基礎題多兩節：**真實數據關卡**（讀檔前的格式解剖、尺度確認、`decisions.md` 決策日誌）與**發表路徑**（訊號分類、驗證設計、novelty 查證、目標期刊層級）。每張題卡都標明對應課程的哪一段與哪一支腳本；資料集的 accession、規模、格式與下載連結整理在 [`exercises/資料集總覽.md`](exercises/資料集總覽.md)（另有 `.xlsx` 表格版）。

題卡之間織了一張**驗證網**：每張進階卡的「發表路徑」都指名可以互相驗證的其他題目（例如 GBM 01↔02↔03、肝癌 08↔09↔22、鱗癌 p-EMT 12↔18↔23），做完一條主線就有交叉驗證的材料。建議路線見 [`exercises/README.md`](exercises/README.md)。

## 參考庫（references/）

分析卡住的時候，缺的往往不是程式碼，而是「該去哪裡查」。[`references/`](references/) 收了四座參考庫，每座都分**精選**（日常分析從這裡選就夠）與**延伸**（什麼情況值得跳出精選），每個資源都寫清楚定位、規模、存取方式、適用與不適用的情境，以及怎麼引用；另附決策樹與 `.xlsx` 表格版。所有資源都逐一打開官方頁查證過，並標了查證日期。

| 參考庫 | 收什麼 | 精選／延伸 | 對應 |
|---|---|---|---|
| [細胞標注](references/細胞標注參考庫/) | marker 字典、自動標注工具與參考集、參考圖譜、惡性判定與細胞狀態 | 16 / 20 | B10、B16 |
| [富集分析](references/富集分析參考庫/) | 基因集資料庫、GSEA / ORA 工具鏈、單細胞活性推斷、網頁工具 | 16 / 13 | B12、B13 |
| [單細胞與空間資料庫](references/單細胞與空間資料庫參考庫/) | 要找新資料時該去哪個儲存庫、圖譜站或空間資料庫 | 16 / 11 | B1、B18 |
| [工具目錄與學習資源](references/工具目錄與學習資源參考庫/) | 做某種分析該用什麼工具、去哪裡繼續學、benchmark 證據在哪 | 12 / 7 | B18 |

> [Q 系列 · 快速上手篇](https://github.com/Charlene717/scRNA-seq-tutorial-Qseries) 的 `references/` 是這四座庫的**精簡版**（精選全收、延伸只挑 Q 系列用得到的）。要完整版就是這裡。

## 環境需求

- R ≥ 4.3 + RStudio，Seurat v5。其餘套件在各腳本開頭列出。
- `B16_infercnv.R` 需要**系統層級**的 [JAGS 4.x](https://sourceforge.net/projects/mcmc-jags/)（rjags 只是介面，R 裝不了 JAGS 本體；裝完要重開 R）。
- B2–B4 需要 [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger)（Linux；記憶體與磁碟需求見官方文件）。跑不動的話，這三集提供現成輸出的下載位置，可以跳過運算直接看判讀。
- 下游工具改版頻繁：B13–B19 的片尾附**工具版本卡**，寫明錄製當下的套件版本與檢查日期。腳本裡的參考輸出都是實跑驗證過的；換版本後結果數字可能不同，判讀原則不變。

## 互動式教學網站

課程之外，同一套教材有三個可以直接點著玩的網站，**中英雙語、免安裝、手機也開得起來**。跑腳本之前先玩過對應的那一章，會比較知道自己在調什麼：

| 網站 | 內容 |
|---|---|
| [單細胞 RNA-seq 互動教學](https://charlene717.github.io/scrna-interactive-tutorial/) | 分析流程十二章，每章都有互動模擬、R 與 Python 兩版程式碼、文獻頁與自測題——與本系列 B5–B15 幾乎逐章對應 |
| [scRNA-seq 進階數學](https://charlene717.github.io/scrna-advanced-math/) | M1–M6 六個模組、六十多個互動示範（含 3D PCA），把降維、分群與統計推論的數學畫出來 |
| [生資互動式教學入口](https://charlene717.github.io/bioinfo-interactive-tutorial-portal/) | 整個生資教學模組的地圖：Git、Linux、R、生物統計、生資概論、單細胞、空間轉錄體、AI 代理人…… |

**互動網站 ↔ 本系列對照**

| 互動頁 | 對應集數 | 建議用法 |
|---|---|---|
| [`qc.html`](https://charlene717.github.io/scrna-interactive-tutorial/qc.html) | B5 | 先在網頁上拉一次閾值看細胞掉多少，再回來跑 `B05_qc.R` |
| [`normalization.html`](https://charlene717.github.io/scrna-interactive-tutorial/normalization.html)、[`variable-features.html`](https://charlene717.github.io/scrna-interactive-tutorial/variable-features.html)、[`scaling.html`](https://charlene717.github.io/scrna-interactive-tutorial/scaling.html) | B6 | 三頁對應 B6 講的三步；SCTransform 的取捨影片講得比較細 |
| [`pca.html`](https://charlene717.github.io/scrna-interactive-tutorial/pca.html) | B7 | 配合 [進階數學 M3](https://charlene717.github.io/scrna-advanced-math/#M3) 的 3D PCA 一起看 |
| [`clustering.html`](https://charlene717.github.io/scrna-interactive-tutorial/clustering.html) | B8 | resolution 拉一遍，感受群數怎麼跟著跑 |
| [`umap.html`](https://charlene717.github.io/scrna-interactive-tutorial/umap.html) | B9 | `n.neighbors` / `min.dist` 的直覺，網頁比講的快 |
| [`annotation.html`](https://charlene717.github.io/scrna-interactive-tutorial/annotation.html) | B10 | marker 面板與自動註釋的對照 |
| [`differential-expression.html`](https://charlene717.github.io/scrna-interactive-tutorial/differential-expression.html) | B12 | 網頁有一張「該用哪種檢定」的決策樹；影片把 pseudobulk + DESeq2 從頭跑一次 |
| [`integration.html`](https://charlene717.github.io/scrna-interactive-tutorial/integration.html) | B11 | 整合前後的對照圖 |
| [`cellchat.html`](https://charlene717.github.io/scrna-interactive-tutorial/cellchat.html) | B14 | 六種圖的讀法先在網頁上認一遍 |
| [`trajectory.html`](https://charlene717.github.io/scrna-interactive-tutorial/trajectory.html) | B15 | 擬時序的直覺 |
| [空間轉錄體互動教學](https://charlene717.github.io/spatial-transcriptomics-interactive-tutorial/) | B18 | B18 只用一段帶過空間轉錄體，想往下走從這裡開始 |
| [生資 R 入門](https://charlene717.github.io/r-bioinformatics-interactive-tutorial/)、[生資 Linux 入門](https://charlene717.github.io/linux-bioinfo-interactive-tutorial/) | B1、B2–B3 | R 還不熟、或第一次在伺服器上跑 Cell Ranger 的前置補課 |

## 系列導覽

這是四個並行的 scRNA-seq 課程系列，可以各自獨立看，也可以互相補位：

| 系列 | 定位 | 集數 | 適合 | 倉庫 |
|---|---|---|---|---|
| **A · 入門篇** | 觀念與判讀，不寫程式 | 6 | 濕實驗背景、要跟生資合作者對話 | [Aseries](https://github.com/Charlene717/scRNA-seq-tutorial-Aseries) |
| **B · 實作篇** | 從 FASTQ 到下游應用的完整縱貫線 | 19 | 要自己做完整分析、要發表 | 本倉庫 |
| **C · 數學篇** | 每一步在算什麼、哪個假設會先破 | 19 | 要審稿、開發方法、解釋參數選擇 | [Cseries](https://github.com/Charlene717/scRNA-seq-tutorial-Cseries) |
| **Q · 快速上手篇** | 三集速成，以 GBM 腫瘤資料貫穿 | 3 | 趕時間、要在幾天內跑完一份分析 | [Qseries](https://github.com/Charlene717/scRNA-seq-tutorial-Qseries) |

**怎麼選路**

- 沒有單細胞背景 → 先看 **A 系列**，再回來。
- 只有幾天、而且做的是腫瘤 → 先走 **Q 系列**，之後要補完整流程再回到 B。
- 跑得出來但不知道為什麼 → 對照著看 **C 系列**：C 的每一集都對應 B 的一個環節（C1–C3 對應 B5–B6、C4–C6 對應 B7–B9、C7–C8 對應 B11–B12、C9–C15 對應 B13–B15、C16–C18 對應 B16–B19）。

## 目錄結構

```
scRNA-seq-tutorial-Bseries/
├── README.md
├── LICENSE                  # 程式碼：MIT
├── slides/                  # 投影片 PDF，中英各二十份
├── R/                       # 十九支隨集腳本
├── data/                    # 資料下載腳本（不放原始資料）
├── exercises/               # 實作練習題 80 題
│   ├── 基礎/                # 20 題：換一份資料把課程流程走完
│   ├── 進階/                # 60 題：發表導向的完整再分析專案
│   └── 資料集總覽.md        # 全部資料集的 accession 與下載方式
├── references/              # 四座參考庫：標注、富集、資料入口、工具與學習資源
└── quiz/                    # 互動自測題庫（index.html + 題目資料）
```

## 授權

- **程式碼**（`R/`、`data/`、`quiz/`）：[MIT License](LICENSE)——可自由使用、修改、再散布。
- **教材**（`slides/` 的投影片、`exercises/` 的題卡與 `references/` 的參考庫，及其中的圖表文字）：[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/deed.zh-hant)——註明出處、非商業使用；商業授權請聯絡作者。

課堂使用（含大學課程、實驗室內部訓練）屬於非商業使用，歡迎直接拿去用，請保留出處。

## 引用

使用本課程材料發表或授課時，請註明：

> Charlene717. *scRNA-seq 教學影片系列 · B 系列 · 實作篇.* https://github.com/Charlene717/scRNA-seq-tutorial-Bseries

用課程材料做出來的分析發表時，也請引用對應的原始資料與工具：

- **資料**：10x Genomics 公開資料集；Kang et al. (2018) *Nature Biotechnology*（ifnb / GSE96583）；各練習題的資料集見 `exercises/資料集總覽.md`。
- **主要工具**：Cell Ranger（10x Genomics）、Seurat v5（Hao et al. 2024）、Harmony（Korsunsky et al. 2019）、DoubletFinder（McGinnis et al. 2019）、SoupX、SingleR、Azimuth、clustree、DESeq2（Love et al. 2014）、edgeR / limma、fgsea、msigdbr、clusterProfiler（Wu et al. 2021）、CellChat v2（Jin et al. 2021）、Slingshot（Street et al. 2018）、tradeSeq、Monocle 3、inferCNV（Broad Institute）、MuSiC（Wang et al. 2019）。

---

問題、勘誤與建議都歡迎開 [Issue](../../issues)。祝分析順利——記住：每個參數是一個假設，有數字、有理由、有紀錄。
