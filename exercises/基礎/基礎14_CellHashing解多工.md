# 基礎 14 · 條碼上的條碼：Cell Hashing 解多工

**難度**：★★ ｜ **預估時間**：1 個工作天 ｜ **對應集數**：B4、B5

## 背景與研究主題

一次上機只跑一個樣本很貴，把八個人的細胞混在一起跑便宜得多——但混完怎麼知道哪顆細胞是誰的？Cell hashing 的答案：上機前給每位 donor 的細胞掛一個抗體偶聯的 **HTO**（hashtag oligo），像行李吊牌；定序後每顆細胞讀它的 HTO 訊號就能認領回家。這是 B4「feature barcode——同一顆細胞可以同時讀多種分子條碼」的經典應用。更妙的是它送你一份禮物：**跨 donor 的 doublet 會同時掛兩種吊牌**，直接現形——基礎10 靠演算法猜的事，這裡有實驗證據。本題照 Seurat hashing vignette 的流程，用 Stoeckius 等人的 8-donor PBMC 資料把解多工完整走一遍。你的研究主題：

1. HTO 訊號怎麼從「連續的計數」變成「這顆細胞屬於 donor X」的離散判定？
2. singlet／doublet／negative 三類各占多少？doublet 率與理論預期對得上嗎？
3. 解多工的品質怎麼驗證？你憑什麼相信 HTODemux 的分類？

## 資料集

- **GSE108313**（Stoeckius；Drop-seq 型平台 + 8 個 HTO，8 位 donor 的 PBMC 混樣）。論文報 14,002 singlets + 2,974 multiplets；vignette 版矩陣約 ~17,916 cells（載入後以 `dim()` 為準）。
- 下載：**用 Seurat hashing vignette 提供的 RNA + HTO 兩個矩陣**（vignette 頁面附 Dropbox 連結）；GEO 上是 raw 資料，本題不用。
- 格式注意：兩個矩陣是分開的檔，cell barcode 集合不完全相同——建物件前要先取交集對齊，這是本題第一個坑。

## 任務

### 階段 A：兩種條碼裝進同一個物件（對應 B4）

1. 分別讀入 RNA 與 HTO 矩陣，`intersect()` 兩邊的 barcode 取交集，各自 subset 到同一組細胞、同一個順序。記錄：兩邊各有多少 barcode、交集剩多少。
2. 用 RNA 建 Seurat 物件，HTO 用 `CreateAssayObject()` 掛成第二個 assay——一個物件、兩種 modality，正是 B4 feature barcode 的資料結構長相。
3. HTO 正規化：`NormalizeData(assay = "HTO", normalization.method = "CLR")`。想一下為什麼 HTO 跟 ADT 一樣用 CLR 而不是 LogNormalize（基礎05 遇過同樣的問題——這類 panel 型計數的統計性質跟全轉錄組不同）。

### 階段 B：HTODemux 與三分類（對應 B4、B5）

4. 跑 `HTODemux()`（預設 positive.quantile = 0.99 起步），看 `table()` 的三分類結果：Singlet／Doublet／Negative 各多少顆、各占幾 %。
5. 用圖驗證分類合理：(a) `RidgePlot()` 畫幾個 HTO 在三分類與各 donor 間的訊號分布——singlet 應該「一個 HTO 高、其餘低」；(b) `HTOHeatmap()` 總覽；(c) 挑一對 HTO 畫 `FeatureScatter()`——doublet 應該落在「兩軸都高」的角落。
6. 看 nCount：doublet 的 RNA nCount 是否系統性偏高（兩顆細胞的內容物）？negative 呢——是空液滴、爛細胞，還是吊牌沒掛好的正常細胞？用 QC 指標佐證你的猜測（B5 的分布思維在這裡完全適用）。

### 階段 C：doublet 率對帳與清理後分析（對應 B4、B5）

7. 對帳：算你抓到的 doublet 率，跟兩個參照比——(a) 總覽給的論文數字（14,002 singlets + 2,974 multiplets，自己算出比例）；(b) 理論盲區：8 個 donor 等比例混樣時，**同 donor 互撞的 doublet 掛同一種吊牌、抓不到**，可辨識的只有跨 donor 的 7/8——所以 HTO 給的 doublet 率是低估，低估多少自己推。寫 150–250 字對帳筆記。
8. `subset()` 出 singlets，快速跑一輪標準流程到 UMAP（參數三件事照基礎01 的標準），以 donor 上色：八位 donor 在各免疫族群裡混得如何？最後回答：如果沒有 HTO，這 2,974 顆 multiplets 混在資料裡，會對下游分析造成什麼樣的污染？

## 繳交物

1. 可重跑的 R 專案（barcode 對齊到解多工全程）。
2. 圖：HTO RidgePlot、HTOHeatmap、HTO pair scatter、三分類 nCount 比較、singlets 的 donor 上色 UMAP（英文標籤）。
3. 分析筆記：barcode 對齊紀錄＋參數三件事＋分類驗證觀察。
4. doublet 率對帳筆記（150–250 字）。

## 自我檢核點

- [ ] barcode 交集對齊有做且有紀錄，兩個 assay 的細胞順序一致
- [ ] 說得出 HTO 為什麼用 CLR 正規化（而不是照抄 RNA 的做法）
- [ ] 三分類結果至少用兩種圖獨立驗證過，不是 `HTODemux` 說了算
- [ ] doublet 率對帳同時比了論文數字與 7/8 理論盲區，知道 HTO doublet 率是低估
- [ ] 能具體說出 negative 那群的可能組成與證據

## 提示（卡關再看）

<details><summary>提示 1：barcode 對不齊</summary>
先各自 `head(colnames())` 看格式——有沒有 `-1` 後綴、大小寫、前綴差異，必要時先清理再 `intersect()`。交集後記得兩個矩陣都要用**同一個向量**做欄 subset，順序才會一致；`identical(colnames(rna), colnames(hto))` 驗證過再往下走。
</details>

<details><summary>提示 2：HTODemux 在做什麼</summary>
它對每個 HTO 擬合背景分布、訂一條「這個 HTO 算陽性」的線（positive.quantile 控制），然後數每顆細胞陽性幾種：0 種 = Negative、1 種 = Singlet、≥2 種 = Doublet。理解這個邏輯你就知道兩件事：quantile 調鬆 doublet 變多、調嚴 negative 變多——可以掃兩三個值看分類穩不穩。
</details>

<details><summary>提示 3：某個 HTO 的 ridge 分不出雙峰</summary>
八種吊牌的染色效率不會一樣好，總有一兩個 HTO 訊噪比較差——這是實驗現實。看它的 negative 率是否特別高、singlet 是否特別少；把「哪個 HTO 品質差」寫進筆記，這正是拿到自己實驗的 hashing 資料時第一個要檢查的事。
</details>

## 進階挑戰

- 用 `MULTIseqDemux()` 重新分類一次，與 HTODemux 做交叉表：兩個演算法對哪些細胞有分歧？分歧細胞的 HTO 訊號長什麼樣？
- 對 singlets 補跑 `scDblFinder`（基礎10 的工具）：它還能抓到 HTO 看不見的**同 donor doublet** 嗎？被它標記的細胞 HTO 分類是什麼？兩種證據（實驗吊牌 vs 演算法）的互補關係想清楚。
- 基礎03 的 ifnb 資料用的是 demuxlet——靠天然遺傳變異解多工、不用掛吊牌。比較兩種策略的適用場景（同一人的多條件混樣，hashing 行、demuxlet 不行——為什麼？）。
- hashing + CITE-seq 全套玩到底的版本在進階58（Hao Azimuth PBMC 參考集，228 個 ADT + hashing）等你。

## 參考文獻

- Stoeckius M, et al. Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. *Genome Biology* (2018). GEO: GSE108313。
- Satija Lab. Demultiplexing with hashtag oligos (HTOs)（Seurat hashing vignette；資料下載來源）。
