# 基礎 19 · 量的是核不是細胞：snRNA-seq 初體驗

**難度**：★★ ｜ **預估時間**：1 個工作天 ｜ **對應集數**：B5、B10 ｜ **前置**：建議先完成基礎07（小鼠腦的層次註解——本題要跨平台對照）

## 背景與研究主題

有些組織拿不到完整的單細胞：冷凍檢體、脂肪、心肌，還有腦——成年神經元在解離時大量死亡，活下來的組成早已失真。single-nucleus RNA-seq（snRNA-seq）改抽細胞核，繞過解離的屠殺，代價是你量到的東西變了：核裡沒有粒線體，**percent.mt 的意義整個反轉**——它不再是「細胞快死了」的指標，而是「胞質污染混進來了」的指標，該接近 0；核內轉錄本大量是未剪接的 pre-mRNA，定序讀到的 intron 也被算進 counts；每顆核的 UMI 與基因數比全細胞低一截；解離壓力基因（Fos/Jun 那一票）反而乾乾淨淨。這一題用 10x 官方的成年小鼠腦核資料走一遍完整流程，重點不是新工具——流程一模一樣——而是**每一步的解讀都要換腦袋**。做完這題，你在基礎07 練的小鼠腦註解功力可以直接沿用，還多了一雙看得懂 sn 資料的眼睛。你的研究主題：

1. snRNA 的 QC 分布長什麼樣？percent.mt 為什麼反轉、閾值該怎麼訂？
2. 腦核資料能註解出哪些型別？與基礎07 的全細胞圖譜比，誰變多了、誰縮水了？
3. 哪些訊號是「核」的招牌（Malat1、低 counts），哪些是全細胞才有的（解離壓力、粒線體）？

## 資料集

- **10x Datasets「5k Adult Mouse Brain Nuclei (Nuclei Isolation Kit)」**：7,377 nuclei detected，10x 3' v3.1（snRNA），Cell Ranger 輸出（filtered h5 / MEX）。
- 下載：10x datasets 頁搜尋名稱，免註冊；優先抓 filtered feature-barcode matrix（HDF5 或 MEX）。
- 格式注意：h5 用 `Read10X_h5()` 讀（需要 `hdf5r` 套件）；小鼠粒線體前綴是小寫 `mt-`（基礎07 踩過的坑，這裡再踩一次就不應該了）。

## 任務

### 階段 A：QC 換腦袋（對應 B5）

1. 讀入、`dim()` 對規模。算 `percent.mt`（`pattern = "^mt-"`）並畫分布：預期大多數核貼著 0。回答兩個問題寫進筆記——核裡為什麼不該有粒線體讀數？percent.mt 偏高的「核」物理上是什麼（破核帶著胞質殘渣、或 ambient RNA 汙染）？據此訂閾值（三件事照規矩；會比你做過的任何 scRNA 題嚴格得多）。
2. nCount/nFeature 分布與你在基礎01/基礎07 看過的全細胞資料並排比較（拿舊圖來對即可），描述量級差異。畫 Malat1 的表現分布——核滯留 lncRNA 是 sn 資料的招牌，順手記下它佔每顆核 counts 的比例大概多少。
3. 標準流程到 UMAP＋分群（參數三件事照基礎01 的標準）。

### 階段 B：腦型別註解（對應 B10）

4. 用基礎07 那套小鼠腦 marker 起手（Snap25、Slc17a7、Gad1/Gad2、Aqp4、Plp1/Mbp、Pdgfra、Cx3cr1/C1qb、Cldn5；不足自查，來源記進筆記）。注意本題是**全腦**核萃取、不是皮質＋海馬，可能出現基礎07 沒見過的區域型別（視丘、紋狀體、顆粒細胞……），註解不出來的群誠實標 unknown 並寫下排除過程。
5. 產出註解後 UMAP、marker DotPlot 與「cluster × 型別 × 證據 marker」對照表。

### 階段 C：sn vs sc 跨平台對照（對應 B5、B10）

6. 型別組成對照：把本題的型別比例與你基礎07（Zeisel，全細胞）的結果做成並排表。預期方向：神經元在 sn 版佔比高得多（它們在解離中活不下來的問題消失了）。但先別急著下結論——兩份資料除了 sn/sc 之外還差了區域（全腦 vs 皮質＋海馬）、平台、年代，哪些組成差異能歸給 sn/sc、哪些不能？把「可歸因」與「混雜」分開寫，這是實驗設計對照組思維的練習。
7. 解離壓力訊號：查 immediate-early genes（Fos、Jun、Hspa1a 等），在你的 sn 資料裡看它們的表現——理論上應該乾淨。寫 150–300 字筆記總結 sn 與 sc 的系統性差異：QC 邏輯、counts 性質、組成偏差、各自適合的場景。這份筆記做進階25/27 時要拿出來用。

## 繳交物

1. 可重跑的 R 專案。
2. 圖：percent.mt 與 nCount/nFeature 分布（標注你的閾值）、Malat1 分布或 FeaturePlot、註解後 UMAP、marker DotPlot、sn vs sc 型別比例並排圖（英文標籤）。
3. 分析筆記：參數三件事＋marker 查證來源＋組成對照的「可歸因／混雜」清單。
4. sn vs sc 差異筆記（150–300 字）。

## 自我檢核點

- [ ] percent.mt 的閾值來自這份資料自己的分布，且說得出它在 sn 情境下「量的是什麼」
- [ ] nCount/nFeature 與全細胞資料的量級差異有圖有數字，不是一句「比較低」
- [ ] 全腦冒出來的陌生群有誠實的 unknown 與排除紀錄，沒有硬套皮質型別
- [ ] 組成比較分得清「sn/sc 效應」與「區域／平台混雜」，沒有把一切都歸給 sn
- [ ] 解離壓力基因有實際查證，筆記能列出至少三條 sn 與 sc 的系統性差異

## 提示（卡關再看）

<details><summary>提示 1：h5 讀不進來</summary>
`Read10X_h5()` 需要 `hdf5r`（`install.packages("hdf5r")`）；裝不了就改抓 MEX 三件套用 `Read10X()`。讀入後 `grep("^mt-", rownames(obj), value = TRUE)` 確認抓到粒線體基因——抓到 0 個先檢查大小寫再檢查物種。
</details>

<details><summary>提示 2：percent.mt 閾值訂不下手</summary>
分布大概率擠在 0–2%，別套 scRNA 的 10–20%。畫 log 尺度或放大 0–10% 區間找「主峰之外的小尾巴」，尾巴就是胞質污染候選。訂 5% 還是 2% 沒有標準答案——數字、理由、紀錄，照三件事走，並記下你因此丟了幾顆核。
</details>

<details><summary>提示 3：組成對照怎麼比才誠實</summary>
只比「大類方向」（神經元 vs 膠細胞的相對佔比、微膠佔比），不要比精確百分比——區域不同，精確比例本來就該不同。若想控制區域變因，可只取你能辨識的皮質型別子集來比。「這個對照不完美」寫出來就是分數，假裝完美才是扣分。
</details>

## 進階挑戰

- 用 label transfer 把你的基礎07 Zeisel 註解轉移到這份核資料上，看跨 sn/sc 的 anchor 品質——哪些型別轉得動、哪些垮掉？
- 查 10x 文件了解 Cell Ranger 的 `include-introns` 設定對 snRNA 定量的影響，寫 100 字說明「intron reads 算不算」為什麼對核資料特別關鍵。
- 對這份核資料跑一次 doublet 偵測（如 DoubletFinder），想想核的 doublet 與細胞的 doublet 形成機制有何不同。
- 這題是 snRNA 疾病資料的前置：人類腦的阿茲海默 snRNA 見進階25，同研究 sn vs sc 正面對決見進階27，腎臟 snRNA 見進階49。

## 參考文獻

- 10x Genomics Datasets: 5k Adult Mouse Brain Nuclei Isolated with Chromium Nuclei Isolation Kit（資料來源頁）。
- Bakken TE, et al. Single-nucleus and single-cell transcriptomes compared in matched cortical cell types. *PLoS ONE* (2018).（sn vs sc 系統性比較）
- Zeisel A, et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. *Science* (2015). GEO: GSE60361.（基礎07 的對照資料）
