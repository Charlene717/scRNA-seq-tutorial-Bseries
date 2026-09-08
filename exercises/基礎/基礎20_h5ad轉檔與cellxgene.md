# 基礎 20 · 圖譜資料拿得下來：h5ad 轉檔與 cellxgene

**難度**：★★ ｜ **預估時間**：1 個工作天 ｜ **對應集數**：B1、B10 ｜ **前置**：建議先完成基礎02（胰臟 marker 是本題的 sanity check；做過基礎09/17 更佳）

## 背景與研究主題

單細胞領域最大的公共資源已經不在 GEO，而在圖譜平台：CZ CELLxGENE Discover 上躺著上千個 curated 資料集，統一格式（h5ad）、統一命名（Cell Ontology）、附作者註解與 embedding，點一下就能下載。問題是它們是 Python 世界的 anndata 格式，進 R 要過轉檔這一關——而轉檔是會掉東西的：counts 還是 normalized？基因名是 symbol 還是 Ensembl ID？metadata 與 UMAP 座標有沒有跟過來？這一題用 Tabula Sapiens 的胰臟資料練完整的一套取用技能：線上導航與探索、下載 h5ad、**兩條轉檔路線都走一遍**（zellkonverter 與 SeuratDisk）、轉檔後逐項驗證，最後用你在基礎02/09/17 攢下的胰臟知識檢查生物學有沒有在搬運中折損。選胰臟不是巧合——對照組織你懂，轉檔壞掉你才看得出來。你的研究主題：

1. cellxgene 上怎麼找資料？線上 explorer 能回答什麼問題、什麼時候才需要下載？
2. h5ad 進 R 的兩條路各自把什麼帶過來、什麼掉了、哪裡會報錯？你怎麼系統性驗證？
3. Tabula Sapiens 的胰臟註解與你熟悉的 Baron/Muraro 圖景一致嗎？差異出在哪一層？

## 資料集

- **CZ CELLxGENE「Tabula Sapiens」collection 的 Pancreas dataset**：14,140 cells，多平台（TS 標準流程），h5ad 格式。
- 下載：cellxgene.cziscience.com 搜「Tabula Sapiens」，collection 頁找到 Pancreas dataset 點 Download 取得 h5ad 永久連結，免登入。
- 格式注意：cellxgene 的 h5ad 裡 `X` 通常是 normalized 值、raw counts 另放在 `raw` 或 layers；`var` 的索引是 Ensembl ID、基因符號在 `feature_name` 欄——轉檔後 rownames 是誰、值是不是整數，一律親眼確認，不要猜。

## 任務

### 階段 A：導航、線上探索與下載（對應 B1）

1. 先別下載。到 cellxgene.cziscience.com 找到 Tabula Sapiens collection，用線上 explorer 玩至少 10 分鐘：按 cell type／組織上色、查 INS 與 GCG 的表現、篩選 Pancreas。列一張兩欄清單：「線上就能回答的問題」vs「非下載不可的問題」——這張清單決定你以後什麼時候值得花這個下載與轉檔成本。
2. 下載 Pancreas dataset 的 h5ad，記錄下載連結、檔案大小與日期（Data availability 的習慣，總覽引用節的格式照用）。

### 階段 B：兩條轉檔路線與逐項驗證（對應 B1）

3. 路線一 **zellkonverter**：`readH5AD()` 得到 SingleCellExperiment，檢查 `assayNames()`、`colData` 與 `reducedDims()`，再轉成 Seurat 物件（`as.Seurat()` 或手動抽矩陣與 metadata 自建）。
4. 路線二 **SeuratDisk**：`Convert("xx.h5ad", dest = "h5seurat")` → `LoadH5Seurat()`。兩條路都要走；其中一條對新版 anndata 檔報錯是常態——**把錯誤訊息與你的處置記下來**，這比一次跑通更接近真實。
5. 轉檔驗證清單，兩條路線各驗一次做成對照表：（a）細胞數＝14,140？（b）metadata 欄位數與關鍵欄（型別註解、donor、assay/平台）在不在？（c）UMAP 等 embedding 有沒有跟過來？（d）表達值是 counts 還是 normalized（看是不是整數、看 max）？raw counts 找得到嗎？（e）基因名是 Ensembl 還是 symbol，要不要用 `feature_name` 換回來？逐項打勾或打叉。

### 階段 C：生物驗證與圖景對照（對應 B10）

6. 選一條路線的產物，用作者附的型別註解畫 UMAP（優先用帶過來的 embedding；帶不過來就自己重跑並註記），再畫 GCG/INS/SST/PPY 加 acinar/ductal marker 的 DotPlot——用你的胰臟知識當裁判：標籤與 marker 對得上，轉檔才算真的成功；對不上，先懷疑轉檔（值的尺度、基因名對錯位）再懷疑註解。
7. 與基礎02 的 Baron 14 型別圖景對照：TS 用 Cell Ontology 命名，粒度與名字差在哪（哪些型別被合併、哪些拆得更細、哪些名字換了說法）？寫 150–300 字筆記：公共圖譜資料「拿來就用」之前，你至少要驗證哪五件事？

## 繳交物

1. 可重跑的 R 專案（兩條轉檔路線各一支 script，含報錯與處置紀錄）。
2. 轉檔驗證對照表（兩路線 × 五驗證項）。
3. 圖：作者註解 UMAP、胰臟 marker DotPlot（英文標籤）。
4. 「線上 vs 下載」問題清單＋圖譜資料使用前檢查筆記（150–300 字）。

## 自我檢核點

- [ ] 有先線上探索再下載，兩欄清單能具體說出 explorer 的能與不能
- [ ] 兩條轉檔路線都實際跑過，報錯有紀錄與處置，不是只寫成功的那條
- [ ] 驗證對照表五項齊全，「X 是 counts 還是 normalized」有證據（不是用猜的）
- [ ] 基因名（Ensembl vs symbol）的處理有明確決定，marker 查得到不是靠運氣
- [ ] marker 與型別標籤的一致性有圖為證，對不上時先排查了轉檔層
- [ ] 與 Baron 圖景的對照分得清「粒度不同」「命名不同」與「真的不一致」

## 提示（卡關再看）

<details><summary>提示 1：zellkonverter 第一次跑很慢</summary>
zellkonverter 走 Bioconductor，底層用 basilisk 自建一個小 Python 環境，第一次執行會下載安裝、等幾分鐘是正常的。大檔可用 `readH5AD(..., use_hdf5 = TRUE)` 讓矩陣留在硬碟上省記憶體。轉 Seurat 前先在 SCE 層把該看的看完——`assayNames(sce)` 告訴你有哪幾層值。
</details>

<details><summary>提示 2：SeuratDisk 報錯</summary>
常見死法：anndata 版本太新、`uns`/`raw` 結構不合它預期。可試只載部分內容（`LoadH5Seurat(..., assays = ...)`）；還是不行就記錄錯誤、靠路線一完成任務——「兩條路線準備一條會斷」正是要你同時學兩條的原因。
</details>

<details><summary>提示 3：找 raw counts 與換基因名</summary>
SCE 版看 `assayNames()`（常見 `X` 與 `counts`，或藏在 `altExp`／h5ad 的 `raw`）；值是不是 counts 用「非整數比例」與 `max()` 判斷。換 symbol：拿 `rowData(sce)$feature_name`（或轉檔後 meta.features 的對應欄）重建 rownames，注意重複的 symbol 要 `make.unique()`。若真找不到 raw counts，就用 normalized 值做展示層面的圖並把限制寫進筆記——這是一個要記錄的決定，不是一個可以沉默帶過的細節。
</details>

## 進階挑戰

- 有 Python 環境的話，用 scanpy 打開同一個 h5ad 看 `adata.obs`／`adata.var`／`adata.raw`，跟你在 R 端看到的逐項對——理解 anndata 原生結構後，轉檔掉東西時你能一眼指出掉在哪。
- 拿你基礎02 的 Baron 物件與轉好的 TS 胰臟做 label transfer 互驗，看兩份獨立圖譜的註解互相認得多少。
- 這套 h5ad 功夫的正式續集：進階30（HLCA，5.6 GB h5ad 的圖譜級 mapping）與進階41（SCLC atlas，同樣從 cellxgene 取 h5ad）；進階13 的 Zenodo h5ad 也走同一條轉檔路。

## 參考文獻

- The Tabula Sapiens Consortium. The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans. *Science* (2022).
- CZ CELLxGENE Discover（cellxgene.cziscience.com；資料平台）。
- Zappia L, Lun A. zellkonverter: Conversion between scRNA-seq objects（Bioconductor 套件文件）。
