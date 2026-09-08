# 進階 02 · 多病人 GBM：整合的兩難

**難度**：★★★ ｜ **預估時間**：1.5–2 個工作天 ｜ **對應集數**：B11（整合＋⑤+ 策略室）、B9、B16、B19 ｜ **資料集**：GSE84465（本題獨用）

## 背景與研究主題

Darmanis 等人（2017）對 4 位 GBM 病人的腫瘤核心（tumor core）與周邊（peripheral tissue）做了單細胞定序。這份資料是 B11 策略室那句話的最佳試煉場：**「惡性細胞按病人分群不是 batch effect，是生物學」**。你的研究主題：

1. 這 4 位病人的細胞，哪些群是「跨病人共享」的，哪些是「病人特有」的？
2. 對這種資料，「整合」與「不整合」各自回答什麼問題？各自掩蓋什麼？
3. **核心 vs 周邊的微環境差在哪裡？**——尤其是 myeloid 區塊：核心的巨噬細胞與周邊的 microglia 是同一群細胞的兩種狀態，還是兩群細胞？

## 資料集

- **GSE84465**（Smart-seq2，3,589 cells／4 病人，腫瘤核心＋周邊）
- 下載：GEO accession 頁 Supplementary files 的 count matrix（單一大檔；確切檔名以 GEO 頁為準）＋ series matrix（細胞的病人與部位 metadata 藏在這裡）。
- 注意：Smart-seq2 是全長、無 UMI 的資料——QC 的尺度與 10x 完全不同，詳見真實數據關卡。

## 真實數據關卡

真實資料在你跑第一行 Seurat 之前就開始出題。照做並記錄：

1. **格式解剖**：先 `readLines(n=5)` 看清楚再讀。這是 gene × cell 的 count matrix 文字檔——幾千個欄的大檔用 `read.csv` 會等到懷疑人生，改用 `data.table::fread` 再轉 sparse matrix。讀入後驗貨：`dim()` 對不對得上 3,589 cells？值是整數 counts 還是已被處理過（`min()`、是否有小數）？確認完才有資格建物件。
2. **metadata 藏在哪**：這個資料集沒有一個現成的「metadata.csv」。病人與部位資訊要嘛從 **series matrix** 檔解析（`GEOquery::getGEO()` 或手動讀 `!Sample_characteristics` 列），要嘛從**細胞名稱的編碼**還原——先把兩條路都看一眼，選一條，然後**用另一條交叉驗證**：兩邊解析出來的病人/部位標籤一致嗎？矩陣的欄名與 metadata 的細胞名對得起來嗎（大小寫、分隔符、前綴）？對不齊的細胞怎麼辦，記錄決定。
3. **無 UMI 的 QC 尺度**：Smart-seq2 的 nCount 是 read 數不是 UMI 數，動輒數十萬起跳；nFeature 也普遍比 10x 高。B5 教的 10x 閾值數字在這裡全部作廢——QC 看這份資料自己的分布抓相對離群，且核心與周邊、不同病人的分布可能系統性不同，閾值要不要分層訂？決定並記錄。
4. **決策日誌**：開一個 `decisions.md`，凡是「資料逼你做的決定」都記一筆：發現了什麼 → 選項 → 你選了什麼 → 理由。這份日誌是繳交物，也是之後論文 Methods 的草稿。

## 任務

### 階段 A：不整合，先看真相（對應 B5–B10）

1. 建 Seurat 物件，帶入病人與部位 metadata（經過交叉驗證的版本）。QC（每個閾值：數字、理由、紀錄）。
2. 標準流程到 UMAP＋分群（不做任何整合）。
3. 以病人上色、以部位上色、以初步註解上色，各出一張 UMAP。回答：哪些群按病人分開？哪些群跨病人混合？

### 階段 B：判定惡性，解釋分群（對應 B10、B16）

4. 用 marker（免疫 PTPRC、血管、寡樹突 MBP、神經元…）註解非惡性群。
5. 對剩餘疑似惡性的群跑 inferCNV（以免疫細胞當 reference——「跨病人混合良好」本身就是 CNV 正常的旁證）。回答：按病人分開的群，是不是正好就是 CNV 異常的群？
6. 產出一張「病人 × 部位 × 惡性狀態」的分群總表。

### 階段 C：只整合微環境，深挖核心 vs 周邊（發表導向）

7. **整合範圍決策**：做兩個版本——(a) 全部細胞整合（Harmony 或 Seurat anchors）；(b) **只整合非惡性細胞**。比較：整合全部時惡性細胞被揉在一起失去什麼？只整合微環境時你能回答哪些原本答不了的問題？這是 B19「整合範圍」決策點的實作。
8. **myeloid 深挖**：在 (b) 的微環境空間 subset 出 myeloid 重新分群，用 microglia（P2RY12/TMEM119）vs 血源性巨噬細胞（CD163 等）的 signature 打分。核心 vs 周邊的 myeloid 組成與狀態差異是什麼？每個結論檢查：是 4 位病人一致，還是被一位病人撐起來的（統計單位！只有 4 位病人，per-patient 的圖誠實畫出來）？
9. 穩健性：換整合方法（Harmony ↔ anchors）、去掉貢獻細胞最多的病人再跑一次，核心 vs 周邊的 myeloid 訊號還在嗎？寫「發表路徑」評估（見下節）。

## 繳交物

1. 可重跑的 R 專案（renv + `set.seed(1234)`、相對路徑）＋ `decisions.md` 決策日誌。
2. 圖：階段 A 三張 UMAP、inferCNV heatmap、(a)(b) 兩版整合對比、myeloid 深挖圖組（含 per-patient 呈現；英文標籤）。
3. 分析筆記：參數三件事＋各階段回答。
4. 發表路徑評估（300–500 字；涵蓋原「策略決定書」的整合設計論證）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- **訊號分類與穩健性**：你的核心 vs 周邊 myeloid 差異屬於指南第二節的哪一類？它通過「4 位病人一致」與換方法/去病人的穩健性檢查了嗎？只有 4 位病人——這個 n 撐得起什麼強度的主張，誠實寫。
- **驗證設計**：微環境的部位差異這條線，現成的驗證場：**GSE131907（進階07，肺腺癌原發 vs 轉移多部位）**驗「腫瘤內外 myeloid 狀態梯度」是否跨癌種成立；**GSE72056（進階03，19 顆黑色素瘤）**驗 myeloid 亞型結構的跨腫瘤一致性。你的 signature 拿過去打分就能初步驗。
- **novelty 定位**：GBM 的 microglia vs 血源性巨噬細胞之爭文獻很厚（PubMed 查「glioblastoma tumor-associated macrophage microglia single-cell」）——原論文與後續工作做過什麼、你的角度（核心/周邊配對＋只整合微環境的設計、或某個被帶過的亞群）還剩多少新？
- **缺什麼＋目標期刊層級**：4 病人的發現要成文，缺的第一塊幾乎一定是外部驗證（上面兩個資料集）＋更大的 n；補上後合理落點是計算/組學期刊或觀察型（Scientific Reports 層級）；若能接臨床變數（進階10 的 TCGA 反卷積接存活）才有機會往領域期刊（Neuro-Oncology 系）談。

## 自我檢核點

- [ ] metadata 用兩條路解析並交叉驗證過，對不齊的細胞有處理決定且進了 decisions.md
- [ ] QC 閾值是看這份 Smart-seq2 資料的分布訂的，能說出與 10x 尺度差異的原因
- [ ] 能明確指出哪些群「按病人分」、哪些「跨病人混」，且 inferCNV 結論與 marker 註解交叉驗證過
- [ ] 兩版整合的差異用「回答什麼問題」而不是「圖好不好看」來比較
- [ ] myeloid 結論有 per-patient 呈現，且做過去掉最大病人的穩健性檢查
- [ ] 發表路徑評估誠實面對 n=4，並指名了驗證資料集

## 提示（卡關再看）

<details><summary>提示 1：Smart-seq2 讀入與 metadata</summary>
`fread` 讀入後第一欄通常是基因名，設為 rownames 再 `as.sparse`。series matrix 用 `GEOquery::getGEO("GSE84465")` 取 `pData()`，或直接解壓文字檔找 `!Sample_characteristics_ch1` 列；細胞名編碼那條路，先 `head(colnames())` 觀察命名規則再寫解析。兩條路對不上的地方就是 decisions.md 的素材。
</details>

<details><summary>提示 2：只整合微環境</summary>
`subset()` 出非惡性細胞另建物件、**重跑 HVG/PCA** 再整合；惡性細胞的分析回到各病人自己的空間做。
</details>

<details><summary>提示 3：4 位病人的統計</summary>
n=4 做不了像樣的檢定——重點放在「方向一致性」：每位病人核心 vs 周邊的 myeloid 狀態比例各畫一條線，4 條線同方向就是你最強的證據；任何 cell-level p 值在這裡都是自欺。回想 B12 策略室的統計單位。
</details>

## 進階挑戰

- 對每位病人的惡性細胞單獨分群，對照 Neftel 四狀態（MES/AC/OPC/NPC-like）打分——與進階01 的單一檢體結果比：狀態組成像嗎？
- 核心 vs 周邊的 myeloid 做 pseudobulk DE（4 位病人、配對設計——設計矩陣怎麼寫？統計單位是誰？）。

## 參考文獻

- Darmanis S, et al. Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma. *Cell Reports* (2017). GEO: GSE84465.
- Neftel C, et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. *Cell* (2019).（進階挑戰用）
- Kim N, et al. *Nature Communications* (2020). GEO: GSE131907.（驗證資料集，見進階07）
