# 進階 59 · 從相關到因果：CROP-seq 的擾動分析範式

**難度**：★★★ ｜ **預估時間**：2 個工作天 ｜ **對應集數**：B12、B18、B19 ｜ **資料集**：GSE92872 Datlinger CROP-seq（本題獨用）

## 背景與研究主題

到目前為止，題庫裡所有結論都是觀察性的——你看到亞群 X 高表達基因 Y，但沒有人動過細胞，因果方向只能用文獻圓。Datlinger 等人（2017）的 CROP-seq 改變了遊戲：每顆細胞帶一支 CRISPR gRNA（可從轉錄本直接讀出指派），敲掉一個基因、再測整個轉錄組——**gRNA 是你親手放進去的處理變數**，這是題庫唯一有因果推論資格的資料集。他們在 Jurkat T 細胞株上敲 TCR 訊號通路的調控者，並設 TCR 刺激 vs 未刺激兩種條件：5,905 顆有唯一 gRNA 指派的細胞，單一 csv 僅 17.3 MB——題庫最小的檔案，裝著範式轉移最大的一題。分析思維要整個換檔：這裡沒有細胞型別多樣性（全是同一株細胞），**分析單位是 perturbation，不是型別**。你的研究主題：

1. 每個 target 的敲除在轉錄組上留下什麼足跡？哪些 target 的足跡大、哪些幾乎無感（gRNA 沒效？基因本來就不表達？）？
2. TCR 刺激改變了哪些 target 的效應？「只有刺激下才看得到的敲除效應」是條件依賴因果的直接展示。
3. 把 target 依效應足跡聚類，能不能重建 TCR 通路的功能結構（同路的 target 足跡相似）？

## 資料集

- **GSE92872**（CROP-seq，Drop-seq 讀出；Jurkat 人類 T 細胞株；5,905 cells，gRNA 唯一指派；TCR 刺激 vs 未刺激）。
- 下載：GEO Supplementary files 的**單一 csv.gz，僅 17.3 MB**——不用取子集、不用大記憶體，本題的門檻全在分析設計不在硬體。
- 檔案把細胞註解（條件、gRNA 指派與其目標基因）與表達矩陣裝在同一個 csv 裡——結構解剖是真實數據關卡第 1 項；target 清單與每 target 細胞數以檔案與原論文為準。

## 真實數據關卡

真實資料在你跑第一行分析之前就開始出題。這一節照做並記錄：

1. **格式解剖**：`readLines(n=10)` 先看清楚——哪幾列（或欄）是註解（condition、gRNA、目標基因）、表達值從哪開始、尺度是什麼（counts？以數值分布與 GEO 頁描述佐證）。把註解拆成 metadata、建物件；統計設計矩陣：條件 × target 的細胞數交叉表——每格幾顆細胞，直接決定你後面每個比較的檢定力。
2. **gRNA 指派欄的解析**：一支 gRNA 對一個 target，但一個 target 通常有多支 gRNA——gRNA 層級與 target 層級是兩層。找出 non-targeting／control gRNA 是哪些（命名以檔案為準）——它們是所有因果比較的基線。同 target 的不同 gRNA 效應一致嗎？這是內建的重複驗證，設計分析時就要留這一手。
3. **細胞株資料的長相**：Jurkat 是永生化細胞株——cell cycle 訊號極強、沒有型別多樣性。跑一次標準分群你會得到幾個「群」，但它們多半是週期相位不是生物亞群（`CellCycleScoring` 驗證）。這逼你放棄「分群 → 註解」的預設劇本：本題的分組變數是 gRNA/target 與條件，都在 metadata 裡，不用 UMAP 找。cell cycle 要不要迴歸掉？它可能也是敲除的下游效應——迴歸掉會不會把因果訊號洗掉？決定與理由進日誌。
4. **決策日誌**：開 `decisions.md`，凡是「資料逼你做的決定」（格式發現、設計矩陣、control 定義、cell cycle 處置）都記一筆：發現了什麼 → 選項 → 你選了什麼 → 理由。這是繳交物，也是論文 Methods 的草稿。

## 任務

### 階段 A：設計矩陣與效力檢查（B12）

1. 完成真實數據關卡，產出條件 × target 細胞數交叉表；細胞數太少的格子（自訂門檻並說理由）標記為低檢定力，後面解讀降級。
2. **敲除有效性檢查**：對每個 target，看它自身的表達在帶該 gRNA 的細胞裡有沒有掉（CRISPR 敲除不保證 mRNA 消失——想想為什麼，寫進筆記）。做一張「target 自身表達 vs control」總表，把「疑似沒打中」的 target 標出來。

### 階段 B：per-target 差異表達（B12）

3. 每個 target vs non-targeting control、**分條件各做一次** DE。方法選擇要過 B12 的腦：這裡的「重複」是什麼（細胞？gRNA？）——statistical 單位問題在擾動資料換了臉但沒消失。至少對一個 target 用「gRNA 當重複」的 pseudobulk 思路重跑，跟 cell-level 檢定比較 p 值膨脹程度。
4. 效應量地圖：每個 target 的 DE 基因數 × 條件的總覽圖。哪些 target 只在刺激下有足跡？挑 2–3 個講出通路故事（它是刺激誘導程式的正調控者還是負調控者——方向從 DE 的正負號直接讀）。
5. 同 target 多 gRNA 的一致性檢查：效應相關嗎？有沒有 outlier gRNA（off-target 或無效的候選）？

### 階段 C：通路重建與範式反思（發表導向；B18、B19）

6. 把每個 target 的效應向量（如 DE 的 log fold change 譜）互相關聯、聚類——功能相近的 target 應該聚在一起。重建出的結構跟已知 TCR 通路對得上嗎？有沒有意外的配對（兩個文獻上不同路的 target 足跡卻高度相似）？
7. 穩健性：換 DE 方法、換 control 定義（全部 non-targeting vs 隨機抽同數量）、丟掉低效力格子重跑，通路結構還在嗎？
8. 寫 300 字方法反思：「B12 的 DE 工具箱哪些直接可用、哪些假設在擾動資料上要重想？」＋寫「發表路徑」評估（見下節，300–500 字）。

## 繳交物

1. 可重跑的 R 專案（renv + set.seed(1234)）＋ `decisions.md` 決策日誌（含設計矩陣）。
2. 圖：條件 × target 細胞數表、敲除有效性總表、效應量地圖、2–3 個 target 的 DE 火山圖、target 聚類熱圖（英文標籤）。
3. 分析筆記：各階段回答＋每個參數三件事＋300 字方法反思。
4. 發表路徑評估（300–500 字）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- 你的發現屬於指南第二節的哪一類？擾動資料的優勢是訊號自帶因果方向——但「意外配對」仍要過穩健性與多 gRNA 一致性這兩關才算數。
- **驗證設計**：題庫內沒有第二個擾動資料集——這正是本題的功課：自查更大的 Perturb-seq 資源（Replogle 2022 全基因組 Perturb-seq 是起點；accession 與資料可得性**自行查證**），評估你的意外配對能否在獨立擾動資料中重現；通路層級的旁證可用觀察性資料（基礎03 的 IFN 刺激 PBMC 是「刺激誘導程式」的現成參照）。
- **novelty 定位**：GSE92872 是方法學經典，重分析者眾——查「CROP-seq reanalysis」與該通路文獻；你的機會在分析範式的比較（statistical 單位、control 設計對結論的影響）或特定意外配對的深挖。
- 若要成文：缺什麼（獨立擾動資料集驗證？效應的劑量關係？機制實驗）？方法比較型的合理落點是計算/組學期刊（指南第四節第二層）。

## 自我檢核點

- [ ] 格式解剖、設計矩陣、control 定義、cell cycle 處置全部進了 decisions.md
- [ ] 敲除有效性檢查做在 DE 之前，「沒打中」的 target 有標記且解讀降級
- [ ] DE 的統計單位想過並寫明，至少一個 target 做過 pseudobulk 對照
- [ ] 條件依賴效應有圖有量化，不是單看一個基因說故事
- [ ] 通路聚類換過方法與 control 定義，結構有重現性
- [ ] 發表路徑評估自查過 Perturb-seq 文獻與資料源，novelty 判斷有依據

## 提示（卡關再看）

<details><summary>提示 1：讀這個 csv</summary>
`read.csv` 全讀後先看前幾列與前幾欄的內容判斷方向（細胞在欄還是在列）。註解列/欄拆出來做 metadata（condition、gRNA、gene 三個關鍵欄），剩下的轉 sparse matrix 建 Seurat 物件。5,905 顆細胞很小，任何筆電都跑得動——把省下的算力花在多做幾組敏感度分析上。
</details>

<details><summary>提示 2：per-target DE 的骨架</summary>
最簡形式：`FindMarkers(obj, ident.1 = "targetX", ident.2 = "CTRL", group.by = "gene", subset = 條件)`。但先想 B12 教的：cell-level Wilcoxon 的 p 值把每顆細胞當獨立樣本——同一支 gRNA 感染的細胞獨立嗎？pseudobulk 版把「gRNA × 條件」當樣本聚合 counts 再上 DESeq2/edgeR；n 會很小，這是擾動資料的真實痛點，兩種做法的結果差異本身就是階段 C 反思的素材。
</details>

<details><summary>提示 3：效應向量聚類</summary>
對每個 target 取「vs control 的 log fold change」向量（限制在一組共同的可變基因上），target × 基因矩陣做相關聚類熱圖。只在刺激條件下做一版、未刺激一版——通路結構應該在刺激下才浮現（TCR 通路沒被啟動時，敲它的調控者理應無感——這個預期本身就是個內建陰性對照）。
</details>

## 進階挑戰

- 用線性模型把「條件 × target」交互作用項直接建進 DE（~ condition * target），跟分層分析的結果比——交互作用項顯著的基因才是「條件依賴效應」的正式統計版。
- 對「疑似沒打中」的 target 深挖：是 gRNA 效率、mRNA 不降解（NMD 逃逸）、還是表達本來就低？三種假說各找證據。
- 查 Replogle 2022（或更新的全基因組 Perturb-seq），寫 200 字：規模從數十 target 到全基因組，分析範式哪些直接放大、哪些要重新設計（多重檢定、效應量先驗、E 級細胞數的工程問題）？

## 參考文獻

- Datlinger P, et al. Pooled CRISPR screening with single-cell transcriptome readout. *Nature Methods* (2017). GEO: GSE92872.
- Replogle JM, et al. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell* (2022).（延伸驗證方向；資料可得性自行查證）
- Kang HM, et al.（基礎03 的刺激 PBMC——觀察性刺激程式參照）GEO: GSE96583.
