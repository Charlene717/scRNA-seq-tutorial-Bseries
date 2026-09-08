# 進階 25 · 換一種細胞、換一套規則：阿茲海默的細胞核

**難度**：★★★ ｜ **預估時間**：1.5–2 個工作天 ｜ **對應集數**：B5、B12、B19（snRNA 規則與五個決策點）｜ **資料集**：GSE138852（本題獨用）

## 背景與研究主題

腦組織拿不到完整的活細胞——神經元經不起解離，凍存檢體更沒得選。所以腦研究測的是**細胞核**：snRNA-seq。Grubman 等人（2019）定序了 12 位個體（阿茲海默症 AD 與對照）的 entorhinal cortex，共 13,214 顆細胞核。這一題一半是生物學（AD 改變了哪些細胞型別的什麼），另一半是方法學：**snRNA 與 scRNA 的規則差異**——細胞核裡幾乎沒有粒線體 RNA（percent.mt 這個 QC 支柱直接失效）、RNA 總量少（nFeature 整體偏低）、大量 intron reads（比對時算不算 pre-mRNA 影響你看到的一切）、但換來 scRNA 拿不到的神經元。B5 策略室與 B19 都預告過 snRNA 規則要重想——這裡就是實戰場。另外還有一個實驗設計陷阱：**8 個 10x libraries、每個 library 混了 2 位個體**——個體與批次糾纏在一起，統計單位問題比腫瘤題更刁。你的研究主題：

1. snRNA 的 QC 規則哪裡跟 scRNA 不同？percent.mt 在這裡量到的到底是什麼？
2. AD vs 對照：哪些型別變的是組成、哪些變的是狀態（astrocyte/microglia 的反應程式）？
3. 每 library 混 2 人——你的 n 到底是 12 還是 8？結論怎麼寫才誠實？

## 資料集

- **GSE138852**（10x snRNA-seq；13,214 nuclei；12 位個體，entorhinal cortex；**8 個 libraries、每 library 混 2 位個體**）
- 格式：**dense csv + covariates csv**。
- 規模小、單一腦區——本題不取子集；小資料的誠實限制留到發表路徑面對。

## 真實數據關卡

真實資料在你跑第一行 Seurat 之前就開始出題。照做並記錄：

1. **dense csv 讀入與尺度確認**：`data.table::fread` 讀、轉 sparse 再建物件。數值是整數 counts 還是已被處理過？讀 GEO 頁與原論文 Methods 確認，並用分布佐證——猜錯尺度，後面全錯。
2. **covariates csv 與個體拆分**：每 library 混 2 位個體，先解剖 covariates 檔：細胞層級有沒有個體指派欄？還是只有 library 層級的資訊（以檔案實況為準——先 `readLines` 看清楚再下結論）？如果個體只能追到 library，你的統計單位就是 **8 個 libraries 而不是 12 位個體**，而且「個體」與「library 批次」完全糾纏、拆不開——這不是你的錯，是實驗設計欠的債，但你的每一個組間結論都要背著它寫。把你對這份檔案的判讀與後果寫進 decisions.md。
3. **snRNA 的 QC 規則改寫**：細胞核幾乎沒有粒線體 RNA——所以 percent.mt 的邏輯**反轉**：在 scRNA 裡高 mt 是垂死細胞，在 snRNA 裡 mt 明顯不為零的「核」反而可疑（可能混入了胞質或 ambient RNA）。畫這份資料的 percent.mt 分布驗證這件事，然後重訂你的 QC 三件套：mt 閾值往哪個方向設？nFeature 下限比照 10x scRNA 的數字合理嗎（核內 RNA 本來就少）？每個閾值：數字、理由、紀錄。
4. **決策日誌**：開 `decisions.md`，凡是「資料逼你做的決定」都記一筆：發現了什麼 → 選項 → 你選了什麼 → 理由。它是繳交物，也是論文 Methods 的草稿。

## 任務

### 階段 A：snRNA 的第一次 QC 與註解（對應 B5）

1. 讀入、掛 covariates、用改寫後的規則 QC——並排畫「照抄 scRNA 規則會砍掉誰」vs「snRNA 規則留下誰」，量化兩者差多少細胞。
2. 分群註解腦的主要型別（神經元亞型、astrocyte、microglia、oligodendrocyte、OPC、內皮…）。注意 snRNA 的 marker 行為：核內轉錄本組成與全細胞不同，某些教科書 marker 會變弱——註解卡住時先懷疑這個。
3. 檢查 library 效應：UMAP 按 library 著色，8 個 libraries 有沒有誰自成一格？需要整合校正嗎（校正單位是 library）——決策與理由。

### 階段 B：AD vs 對照（對應 B12）

4. 組成比較：每個 library（或個體，依關卡 2 的判讀）各型別佔比，AD vs 對照。**統計單位用你在關卡 2 決定的那個**，點數誠實畫在圖上。
5. 狀態比較：astrocyte 與 microglia 各自子集重分群，找 AD 富集的反應狀態（marker 自查文獻：反應性 astrocyte 如 GFAP 上調程式、疾病相關 microglia 程式等），DE 用 pseudobulk——n 只有 8 個 libraries，檢定力有限是事實，效果量與方向一致性比 p 值誠實。
6. 神經元那邊呢？entorhinal cortex 是 AD 最早受損的腦區——神經元亞型的組成有沒有 AD 相關的缺失訊號？（組成變化在 snRNA 裡也受解離/取樣偏誤影響，解讀收斂一點。）

### 階段 C：跨疾病的髓系對話（發表導向）

7. 把你的 AD microglia 反應程式做成 signature，與**進階24 的 IPF 疾病巨噬程式**做基因層級比較（Jaccard、打分相關；沒做進階24 的人：用該題原論文的巨噬 marker 附表替代）——「疾病相關髓系狀態」跨器官通用的部分與腦特有的部分各是什麼？這是指南第一節「新組合」型的切入。
8. 穩健性：換 QC 閾值（snRNA 規則內合理範圍）、去掉一個 library、換 resolution，你的反應程式與跨題比較結論還在嗎？
9. 寫「發表路徑」評估（見下節）。

## 繳交物

1. 可重跑的 R 專案（renv + `set.seed(1234)`）＋ `decisions.md` 決策日誌。
2. 圖：scRNA vs snRNA QC 規則對照圖、註解 UMAP（含 library 著色）、AD vs 對照組成與狀態圖組、microglia 程式跨題比較圖（英文標籤）。
3. 分析筆記：參數三件事＋各階段回答＋關卡 2 的統計單位判讀。
4. 發表路徑評估（300–500 字）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- **訊號分類與穩健性**：你的 microglia/astrocyte 程式屬於指南第二節哪一類？以 library 為單位（n=8）畫出來的組間差異，有幾分把握不是某一兩個 library 撐起來的？去 library 檢查的結果如何？
- **驗證設計**：AD 的更大隊列如 Mathys 等人（2019）等——**存取有門檻（管控資料需申請），自行查證 accession 與申請條件**，把「哪個隊列、怎麼拿、拿到能驗什麼」寫進評估；查的過程本身是練習。跨疾病那半邊的比較場在題庫內：**進階24 的 IPF 巨噬狀態**（見階段 C）。
- **novelty 定位**：AD snRNA 文獻極厚（Grubman 與 Mathys 同年，之後爆發）——PubMed 查「Alzheimer snRNA-seq microglia astrocyte」；本資料集只有 13,214 nuclei、單一腦區、n=8 libraries，靠它單獨立論幾乎不可能。你的 novelty 更可能在「跨疾病髓系比較」或「snRNA QC 決策的系統性展示」這種組合/方法角度。
- **缺什麼＋目標期刊層級**：跨疾病比較若擴到 3+ 個疾病資料集（腦、肺、肝各一），是一個能往計算/組學期刊走的骨架；只有兩兩比較，落點是 Scientific Reports／Frontiers 層級或先掛 bioRxiv。任何 AD 生物學主張都缺大隊列驗證——沒拿到管控資料前，措辭停在「與文獻一致」。

## 自我檢核點

- [ ] 尺度與 covariates 的判讀（個體能不能拆到細胞層級）進了 decisions.md，統計單位的選擇有明寫
- [ ] QC 用 snRNA 規則重訂，且量化了「照抄 scRNA 規則」的後果
- [ ] percent.mt 的邏輯反轉能講出為什麼，不是背口訣
- [ ] 所有 AD vs 對照結論用 library/個體層級統計，n 誠實標在圖上
- [ ] microglia 程式與進階24（或其論文附表）做過跨疾病比較
- [ ] 發表路徑評估查證過更大 AD 隊列的存取條件，novelty 定位誠實

## 提示（卡關再看）

<details><summary>提示 1：dense csv 讀入</summary>
`fread` 後第一欄通常是基因名：`genes <- dt[[1]]; m <- as.matrix(dt[,-1]); rownames(m) <- genes` → `as(m, "dgCMatrix")`。covariates 用細胞名 `match()` 對齊，對不上的細胞數要報告。
</details>

<details><summary>提示 2：mt 幾乎為零時的 QC</summary>
`PercentageFeatureSet(pattern="^MT-")` 照算，但看分布：多數核應貼近 0。閾值改成「排除 mt 異常高的核」（例如 >5% 就很可疑，實際數字看分布訂）。QC 主力改壓在 nFeature/nCount 的雙尾與 doublet 偵測上。
</details>

<details><summary>提示 3：n=8 的組間比較</summary>
pseudobulk 到 library 層級後每組只有 4 個點——正式檢定選項有限（Wilcoxon 在 n=4+4 幾乎沒力）。誠實的做法：報效果量＋每點透明畫出＋方向一致性（幾個 library 同向），把「檢定力不足」明寫進筆記，而不是硬擠 p 值。
</details>

## 進階挑戰

- 找一個公開的人腦 scRNA（非核）資料集（自查 GEO），比較同型別在 sn vs sc 的可測性與 marker 差異——把「snRNA 換走了什麼、換來了什麼」做成量化附表。
- 對 doublet 加一層偵測（DoubletFinder 等）：每 library 混 2 人的設計下，跨個體 doublet 理論上可用基因型拆——查 demuxlet/souporcell 的思路，評估這份資料（無基因型檔時）能做到哪一步。

## 參考文獻

- Grubman A, et al. A single-cell atlas of entorhinal cortex from individuals with Alzheimer's disease reveals cell-type-specific gene expression regulation. *Nat Neurosci* (2019). GEO: GSE138852.
- Mathys H, et al. Single-cell transcriptomic analysis of Alzheimer's disease. *Nature* (2019).（更大 AD 隊列；存取條件自行查證）
- Adams TS, et al. *Sci Adv* (2020).（跨疾病髓系比較場，見進階24）GEO: GSE136831.
