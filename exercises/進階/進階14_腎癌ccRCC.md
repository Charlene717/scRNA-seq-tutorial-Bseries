# 進階 14 · 腎癌 ccRCC：細胞起源、缺氧程式與一顆離群腫瘤

**難度**：★★★ ｜ **預估時間**：2 個工作天 ｜ **對應集數**：B10、B14、B16、B19 ｜ **資料集**：GSE159115（本題獨用）

## 背景與研究主題

亮細胞腎細胞癌（ccRCC）是「細胞起源」故事講得最完整的癌症之一：惡性細胞帶著近曲小管（proximal tubule）的身分證，加上 VHL 失活驅動的假性缺氧程式（HIF 目標基因全開）——不用缺氧也活得像缺氧。Zhang 等人（2021）的資料含約 30,000 顆細胞：8 顆腫瘤（**7 ccRCC + 1 chromophobe**）+ 6 份正常腎。正常腎給了你完整的腎小管上皮圖譜當參照，讓「惡性細胞最像哪一段小管」變成可以直接檢驗的問題；而那 1 顆 chromophobe（嫌色細胞癌，起源與生物學都不同）是題目送你的天然「離群腫瘤」教材——**留還是剔除，本身就是一個要寫進日誌的決策**，不是理所當然。你的研究主題：

1. ccRCC 惡性細胞的轉錄身分最接近正常腎的哪段小管？「cell of origin」的表達層證據能做到多強？
2. VHL/缺氧程式在 7 顆 ccRCC 間的一致性與病人間差異——哪些是共同核心、哪些是病人私有？
3. 那顆 chromophobe 在你的每一步分析裡是資產還是干擾？什麼問題該留它、什麼問題該剔除？

## 資料集

- **GSE159115**（10x 3' v2，~30,000 cells；8 腫瘤 = 7 ccRCC + 1 chromophobe，+ 6 正常腎）
- 下載：GEO Supplementary files。**per-sample H5 檔 + 官方細胞註解 csv**——每個樣本一個 H5，要批次讀入再合併；作者註解 csv 是現成的「對答案」材料。
- 樣本層級的臨床細節（分期、性別等）以 GEO 頁為準。

## 真實數據關卡

真實資料在你跑第一行 Seurat 之前就開始出題。照做並記錄：

1. **per-sample H5 批次讀入**：先 `list.files` 盤點——幾個 H5、檔名怎麼編碼樣本身分（GSM 編號？樣本名？腫瘤/正常標記在哪？）。寫一個迴圈 `Read10X_h5` → `CreateSeuratObject` → 加樣本 metadata → `merge`，**細胞 barcode 加樣本前綴防撞名**（同一 barcode 出現在兩個樣本是 10x 常態）。合併後硬檢查：總細胞數、每樣本細胞數對得上 GEO 頁嗎？
2. **官方註解 csv 對齊**：註解 csv 的細胞 ID 格式跟你合併後的細胞名多半**不一樣**（前綴規則、分隔符、後綴 -1 有無）——寫出對齊規則，`table(對上, 對不上)` 攤開來看，對不上的細胞數與處理方式記錄。作者註解此後只當「答案卷」，不進你的分析主線。
3. **chromophobe 決策（現在就做第一版）**：讀 GEO 頁確認哪個樣本是 chromophobe。先寫下預設方案（例：整合與註解階段保留、ccRCC 特異分析剔除、階段 B 拿它當離群對照），跑的過程中若證據要求改變方案，更新決策並記錄前後版本——決策可以改，但每一版都要留痕。
4. **決策日誌**：`decisions.md` 記下每一筆「資料逼你做的決定」（讀入迴圈、前綴規則、對齊損失、chromophobe 方案）：發現了什麼 → 選項 → 你選了什麼 → 理由。這是繳交物，也是論文 Methods 的草稿。

## 任務

### 階段 A：批次讀入、整合與版圖（B10）

1. 完成真實數據關卡。整合範圍決策：14 個樣本、腫瘤 + 正常——batch 單位是什麼？先寫預期再跑，檢查整合後正常腎小管各段是否仍分得開（過度整合會把小管段抹平——這裡是本題的整合品質試金石）。
2. 自己的 marker 流程註解：正常腎要拆出小管各段（近曲小管、髓袢、遠曲小管、集尿管等）與內皮/免疫/基質；與官方註解 csv 對答案，混淆矩陣進筆記。
3. 寫 100 字：對答案後你最不放心自己哪一類註解、為什麼？

### 階段 B：惡性判定、細胞起源與缺氧程式（B16、B19）

4. 惡性判定靠 CNV 不是 marker：inferCNV 以免疫細胞（與正常上皮）當 reference。ccRCC 有教科書級的染色體 3p 缺失——你的 CNV 圖看得到嗎？chromophobe 的 CNV 長相跟 ccRCC 一樣嗎（它是天然的陰性對照/離群對照，用起來）？
5. cell of origin：把 7 顆 ccRCC 的惡性細胞與正常腎各小管段做轉錄相似度分析（相關性、或以各段 signature 打分），惡性細胞是否一致地最靠近近曲小管？chromophobe 的惡性細胞靠近哪一段——與 ccRCC 相同嗎？
6. VHL/缺氧程式：HIF 目標與缺氧基因清單（如 VEGFA、CA9、NDUFA4L2 等，清單來源記錄）對惡性細胞打分，病人內與病人間的分布攤開；拆出「7 顆 ccRCC 共同的核心程式」與「病人私有程式」（病人內先做、再跨病人找共性——B19）。

### 階段 C：被忽略的角落（發表導向）

7. 從三個方向擇一深挖：(a) 腫瘤內皮——ccRCC 是高度血管化的癌，腫瘤 vs 正常腎內皮的程式差異與 B14 通訊分析（誰在餵 VEGF 訊號？）；(b) 腫瘤相關 myeloid 的亞型結構與跨病人一致性；(c) 正常腎近曲小管的異質性——起源段內部有沒有更細的結構能對應惡性程式的病人間差異？
8. 穩健性：換 resolution、去掉貢獻最大的病人、剔除/保留 chromophobe 重跑，訊號還在嗎？
9. 寫「發表路徑」評估（見下節，300–500 字）。

## 繳交物

1. 可重跑的 R 專案（renv + set.seed(1234)）＋ `decisions.md`（含 chromophobe 決策的版本紀錄）。
2. 圖：註解對答案混淆矩陣、inferCNV 圖（標出 3p 區域與 chromophobe 對照）、cell-of-origin 相似度圖、缺氧程式病人內/間分布圖、階段 C 圖組（英文標籤）。
3. 分析筆記：各階段回答＋參數三件事。
4. 發表路徑評估（300–500 字）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- 你的訊號屬於指南第二節的哪一類？7 顆 ccRCC 的 n 很小——結論寫成「7 顆中 X 顆一致」的誠實句式，而不是過度外推。
- **驗證設計**：本題的天然出口是 **TCGA-KIRC**（沿進階10 的反卷積 + 存活流程）——ccRCC 擁有 bulk 文獻裡最出名的免疫浸潤—預後關聯之一，你的微環境訊號（內皮程式、myeloid 亞型比例）接上 KIRC 存活是現成的驗證設計；單細胞端的獨立 ccRCC cohort 需自查（PubMed/GEO 搜尋，這是本題的文獻功課）。
- **novelty 定位**：ccRCC 的 cell of origin 與缺氧程式是原論文主線，不是你的發表點——你的機會在階段 C 的角落（內皮通訊細節、myeloid 結構、起源段內部異質性）＋「單細胞訊號 × KIRC 預後」的新組合。PubMed 查你的具體訊號。
- 若要成文：缺哪一塊（獨立單細胞 cohort／KIRC 方向一致性／Cox 多變數校正）？對照指南第四節評估合理層級。

## 自我檢核點

- [ ] H5 批次讀入有防撞名前綴與硬檢查，官方 csv 對齊損失有數字
- [ ] chromophobe 決策有第一版與（若有）修訂版，每一版理由留痕
- [ ] 惡性判定看得到 3p 缺失（或誠實報告看不到與可能原因），chromophobe 當對照用過
- [ ] cell-of-origin 與缺氧程式結論以病人為單位陳述，n=7 的限制寫明
- [ ] 階段 C 訊號通過至少三種穩健性檢查（含 chromophobe 剔除敏感度）
- [ ] 發表路徑評估指名 TCGA-KIRC 的具體驗證動作，novelty 查過文獻

## 提示（卡關再看）

<details><summary>提示 1：批次讀入骨架</summary>
`files <- list.files(dir, pattern="\\.h5$", full.names=TRUE)`，迴圈內 `m <- Read10X_h5(f); obj <- CreateSeuratObject(m, project=樣本名)`，合併用 `merge(x, y, add.cell.ids=樣本名向量)`——add.cell.ids 就是防撞名前綴。樣本名從檔名 parse，parse 規則寫進 decisions.md。
</details>

<details><summary>提示 2：小管段 signature</summary>
近曲小管看 LRP2/CUBN 類、集尿管看 AQP2 類——完整清單從人腎圖譜文獻自查並記錄出處。cell-of-origin 分析用「正常腎各段的 pseudobulk profile × 每顆腫瘤惡性細胞 pseudobulk」的相關性矩陣最穩，避免 cell-level 相似度的雜訊。
</details>

<details><summary>提示 3：離群腫瘤的用法</summary>
chromophobe 別急著丟——「一個已知不同起源的腫瘤在你的 cell-of-origin 分析裡落在哪」是你整套方法的 sanity check：它若也貼近近曲小管，你的相似度分析大概在量測別的東西（文庫、細胞週期、壓力程式）。陰性對照的價值正在於此。
</details>

## 進階挑戰

- 沿進階10 流程建 signature 反卷積 TCGA-KIRC：你的 myeloid/內皮亞型比例與存活的關聯，跟文獻中出名的 KIRC 免疫浸潤—預後方向一致嗎？
- 對 6 份正常腎做一次「純圖譜」分析（不看腫瘤），把近曲小管拆到你敢拆的最細，回頭檢驗起源段內部結構與 7 顆腫瘤程式差異的對應。

## 參考文獻

- Zhang Y, et al. Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response. *PNAS* (2021). GEO: GSE159115.
- Colaprico A, et al. TCGAbiolinks（TCGA-KIRC 驗證用，流程見進階10）. *Nucleic Acids Research* (2016).
