# 基礎 08 · 第一次面對四萬顆細胞：Drop-seq 小鼠視網膜

**難度**：★ ｜ **預估時間**：1 個工作天 ｜ **對應集數**：B5–B10（B8、B9 重點）

## 背景與研究主題

Macosko 等人（2015）發明 Drop-seq 時拿小鼠視網膜當展示品：44,808 顆細胞，比你前面做過的任何一題大一個量級。這一題有兩個主題。第一是**規模**：dense txt 矩陣直接讀進來記憶體會很難看，你要學會「讀入即轉 sparse」的習慣，並親手量測差多少。第二是**解析度**：視網膜是分層組織，rod 佔絕對多數，bipolar/amacrine 底下藏著大量細胞數很少的亞型——nPC 夠不夠（B8）、resolution 切多細（B9），直接決定你看得到看不到它們。原論文的 cluster 註解檔公開在 McCarroll lab 網站，這題可以逐 barcode 對答案。你的研究主題：

1. 44,808 顆細胞的 dense txt，怎麼讀才不會把 16 GB 記憶體吃光？sparse 化省了幾倍？
2. 視網膜的主要細胞大類（感光、雙極、無軸突、節細胞、膠細胞…）分得開嗎？
3. nPC 與 resolution 怎麼影響「bipolar 亞型」這種少數群的解析度？你的亞型跟作者的 cluster 對得上嗎？

## 資料集

- **GSE63472**（Macosko；Drop-seq 平台）。小鼠 P14 視網膜，44,808 cells，7 個 replicates，merged DGE txt（檔名以 GEO 頁為準）。
- 下載：GEO accession 頁 Supplementary files；**cluster 註解檔**（barcode → cluster 對照）在 McCarroll lab 網站 `mccarrolllab.org/dropseq`（檔名以網站為準）。
- 格式注意：DGE 是 gene × cell 的 **dense txt**，建議用 `data.table::fread()` 讀、立刻轉 `dgCMatrix` 再建物件；本題建議 16 GB RAM，讀入階段不要同時開一堆大物件。

## 任務

### 階段 A：記憶體意識的讀入與 QC（對應 B5–B7）

1. 用 `fread()` 讀入 DGE，第一欄設為基因名，其餘轉 matrix 後 `as(x, "dgCMatrix")` 轉 sparse。用 `object.size()`（或 `lobstr::obj_size()`）記錄 dense 與 sparse 兩個版本的大小，算出倍率，轉完把 dense 版 `rm()` + `gc()` 掉——這段量測寫進筆記。
2. QC：小鼠 `^mt-` 前綴（抓不到就先 grep 查）；4 萬顆細胞的 nFeature/nCount 分布長尾比小資料更明顯，閾值三件事照舊，另外回答：Drop-seq 的每細胞深度跟你做過的 10x/C1 比起來如何？
3. Normalize→HVG→Scale。注意 `ScaleData` 預設只 scale HVG——這在大資料上是刻意的省記憶體設計，想想為什麼夠用。

### 階段 B：大類註解與對答案（對應 B8–B10）

4. PCA：4 萬顆細胞、十幾種型別的組織，ElbowPlot 的「拐點」比 PBMC 模糊。試 nPC = 10 與 30 各跑一次下游，看哪些群在 nPC 太少時消失——這就是 B8 說的「PC 數決定你帶多少結構進下游」。
5. 分群＋UMAP 後查 marker 註解大類（來源記進筆記）。起點：rod（Rho）、cone（Opn1sw/Arr3）、bipolar（Vsx2/Otx2）、amacrine（Tfap2a/Gad1）、horizontal（Lhx1）、節細胞 RGC（Rbpms）、Müller glia（Rlbp1/Glul）、微膠與血管細胞自查。rod 會是壓倒性大群，屬正常。
6. 對答案：讀入作者 cluster 註解檔，按 barcode join 到你的 metadata（對不上的 barcode 有多少？先想為什麼——作者的 QC 跟你的不同）。做你的 cluster × 作者 cluster 交叉表，檢討大類層級的一致性。

### 階段 C：resolution 與 bipolar 亞型（對應 B9 策略）

7. 全資料試至少 3 個 resolution（如 0.5/1/2），追蹤 bipolar 在各 resolution 下裂成幾群——用 `clustree` 或交叉表呈現。
8. `subset()` 出 bipolar 大群，重算 HVG/PCA 後以較高 resolution 重分群（B9 先大後小）。與作者註解檔中的 bipolar 亞型 cluster 對答案，寫 200–300 字筆記：全資料調高 resolution 與「子集重分群」兩條路，哪條把亞型拆得更乾淨？rod 這種巨無霸群在兩條路裡各扮演什麼角色（想想 HVG 被誰主導）？

## 繳交物

1. 可重跑的 R 專案（含讀入與 sparse 化程式碼、記憶體量測紀錄）。
2. 圖：QC 分布、大類註解 UMAP、marker DotPlot、resolution 比較圖、bipolar 子集重分群 UMAP（英文標籤）。
3. 分析筆記：參數三件事＋nPC 10 vs 30 的比較＋與作者 cluster 的交叉表檢討。
4. bipolar 亞型解析度筆記（200–300 字）。

## 自我檢核點

- [ ] dense→sparse 的大小差有實測數字，讀入後 dense 版有釋放
- [ ] nPC 的選擇有 10 vs 30 的實驗證據，不是「教學用 10 所以用 10」
- [ ] 大類註解每型都有小鼠 marker 證據與出處
- [ ] 與作者註解 join 時對不上的 barcode 數量有查、有解釋
- [ ] 能用自己的圖說明「全資料高 resolution」與「子集重分群」對亞型解析度的差別

## 提示（卡關再看）

<details><summary>提示 1：fread 與 sparse</summary>
`dt <- data.table::fread("DGE.txt.gz")`；`genes <- dt[[1]]; m <- as.matrix(dt[, -1]); rownames(m) <- genes; sm <- as(m, "dgCMatrix"); rm(dt, m); gc()`。若連 dense matrix 這一步都撐不住，改用分塊讀入或 `Matrix::readMM` 類路線；最重要的原則：同一時間記憶體裡只留一份大東西。
</details>

<details><summary>提示 2：rod 淹沒一切</summary>
rod 佔多數時，HVG 與前幾個 PC 都會被「rod vs 其他」主導，少數型別的結構被擠到後面的 PC——這正是 nPC 要拉高、以及子集重分群有效的原因。看 `DimHeatmap` 逐 PC 檢查各 PC 在分什麼。
</details>

<details><summary>提示 3：barcode 對不上</summary>
檢查兩邊 barcode 格式是否一致（有無 replicate 前綴、大小寫、`-1` 後綴）。修格式後仍對不上的細胞通常是「你留了作者丟掉的」或反過來——把這群細胞的 QC 指標畫出來看，往往正是兩套 QC 標準的差異所在。
</details>

## 進階挑戰

- 對 amacrine 大群重複階段 C 的子集重分群，比較 bipolar 與 amacrine 兩個家族誰更難拆、為什麼。
- 挑一個你拆出的 bipolar 亞型，查文獻找它的功能注釋（ON/OFF、rod/cone bipolar），驗證你的 marker 邏輯。
- 這裡練的「大矩陣讀入＋取捨記憶體」是進階07（肺腺癌 20 萬細胞大檔）與進階24（2 GB 壓縮矩陣）的生存技能——先在 4 萬顆細胞上把習慣養好。

## 參考文獻

- Macosko EZ, et al. Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. *Cell* (2015). GEO: GSE63472；cluster 註解檔見 mccarrolllab.org/dropseq。
