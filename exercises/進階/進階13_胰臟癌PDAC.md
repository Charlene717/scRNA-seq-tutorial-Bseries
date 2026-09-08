# 進階 13 · 胰臟癌 PDAC：免疫沙漠與一場檔案偵探戰

**難度**：★★★ ｜ **預估時間**：2–2.5 個工作天 ｜ **對應集數**：B10、B14、B16、B19 ｜ **資料集**：CRA001160／Zenodo 3969339（本題獨用）

## 背景與研究主題

胰管腺癌（PDAC）在單細胞視角下有兩張著名的臉：**免疫沙漠**（effector T 細胞稀少、抑制性 myeloid 當家）與**CAF 亞型分工**（myCAF/iCAF 的空間與功能分野）。Peng 等人（2019）的 57,530 顆細胞、24 位 PDAC + 11 份正常胰臟，是少數帶著正常器官對照的癌症大隊列——正常胰臟的 ductal cell 給了你「導管上皮基準線」，讓 PDAC 最難的問題（腫瘤導管細胞 vs 正常導管細胞怎麼分）有了著力點。但這題最大的特色在拿到資料之前：原始檔存放在中國 GSA、只有 FASTQ 且需登入，**主路線改走 Zenodo 的處理後 h5ad——而那個檔案的檔名誤標成 CRC（大腸癌）**。查證「我下載的到底是什麼」，本身就是本題第一課。你的研究主題：

1. 有正常胰臟當參照，PDAC 的 ductal cell 惡性判定能做到多乾淨？「導管化生程度」是連續譜還是二分？
2. PDAC 的免疫沙漠在這個 cohort 裡有多沙漠——跟正常胰臟比、跨病人比，怎麼量化？
3. myCAF/iCAF 亞型結構與它們對微環境的通訊角色（B14），有沒有被原論文輕描淡寫的細節？

## 資料集

- **CRA001160**（PRJCA001063；scRNA-seq 10x 型，57,530 cells；24 PDAC + 11 正常胰臟）
- **下載路線是本題最大特色**：GSA（CRA 開頭）僅提供 FASTQ 且需註冊登入——不走。**主路線：Zenodo record 3969339 的處理後 h5ad**（注意：檔名誤標 CRC，內容確為 PDAC——查證程序見真實數據關卡）。
- **替代路線**：GEO GSE155698（Steele 2020，PDAC，~124,575 cells QC 後）——Zenodo 路線不通時的備案，也是階段 C 的外部驗證集。
- h5ad 是 Python/AnnData 格式，進 R 需 zellkonverter 或 SeuratDisk 轉檔——轉檔本身是關卡之一。

## 真實數據關卡

真實資料在你跑第一行 Seurat 之前就開始出題。照做並記錄：

1. **來源查證（檔名 ≠ 內容）**：下載前先讀 Zenodo record 頁的描述與關聯論文；下載後用三個獨立證據確認內容是 PDAC 而非檔名寫的 CRC：(a) 細胞數與樣本數對得上原論文嗎？(b) metadata 的樣本名長相（有沒有 normal pancreas 類樣本）？(c) 胰臟 marker（如 PRSS1 腺泡、INS 內分泌）有訊號、大腸 marker 沒有？把這套「查證程序」寫成 decisions.md 的第一條——以後你再遇到可疑檔案，這就是 SOP。
2. **h5ad → Seurat 轉檔**：zellkonverter（經 SingleCellExperiment）或 SeuratDisk（h5ad→h5Seurat）擇一，寫下選擇理由。轉檔後硬檢查：細胞數、基因數、稀疏性保留了嗎？`X` 裝的是 counts 還是已正規化的值？AnnData 常把 raw counts 放在 `layers` 或 `raw` 而 `X` 是處理過的——搞錯這一點，下游全錯。
3. **metadata 考古**：obs 欄位有哪些？作者的型別註解、樣本/病人欄叫什麼、編碼規則？跟原論文圖對一對（型別種類、比例量級），確認註解欄可信再拿來對答案。
4. **非 GEO 平台的引用**：Data availability 要寫 GSA accession（CRA001160）、Zenodo record 與下載日期——處理後資料的出處鏈比 GEO 複雜，現在就記清楚。
5. **決策日誌**：`decisions.md` 記下每一筆「資料逼你做的決定」（查證結論、轉檔工具、X/layers 判定、metadata 對齊）：發現了什麼 → 選項 → 你選了什麼 → 理由。

## 任務

### 階段 A：查證、轉檔與版圖重建（B10）

1. 完成真實數據關卡：查證通過、轉檔成功、帶作者註解的 Seurat 物件就位。
2. 自己的 marker 流程做大類註解（ductal、acinar、endocrine、T/B/myeloid、fibroblast、endothelial…），與作者註解對答案，混淆矩陣進筆記。
3. 樣本盤點表：24 PDAC + 11 正常，每樣本細胞數、各大類比例——這張表是階段 B 免疫沙漠量化的底稿。

### 階段 B：惡性導管、免疫沙漠與 CAF 分工（B16、B19、B14）

4. 惡性判定靠 CNV 不是 marker：inferCNV 以正常胰臟的 ductal cell（加免疫細胞）當 reference，判定 PDAC 樣本裡的惡性導管細胞。正常樣本裡的 ductal cell 有沒有被誤判成惡性的？這是你的流程的假陽性率估計。
5. 免疫沙漠量化：以病人為單位算 T/NK、effector 程式分數、抑制性 myeloid 比例，PDAC vs 正常胰臟比較；病人間的「沙漠程度」排序出來——哪些病人其實不沙漠？
6. CAF 亞型：fibroblast 子分群，用 myCAF（ACTA2、TAGLN 類）/iCAF（IL6、CXCL 家族類）marker 定位亞型；接 B14 通訊分析，iCAF 對免疫區塊、myCAF 對惡性導管的訊號各長什麼樣？通訊結果照 B14 的批判性檢視原則過濾。

### 階段 C：被忽略的角落（發表導向）

7. 從三個方向擇一深挖：(a) 「不沙漠」的 PDAC 病人——他們的微環境哪裡不同，有沒有一致的程式？(b) 正常胰臟 ductal cell 的異質性——惡性判定的 reference 本身有結構嗎，影響判定嗎？(c) myCAF/iCAF 之外的 fibroblast 狀態——跨病人反覆出現嗎？
8. 穩健性：換 resolution、去掉貢獻最大的病人、換 marker 定義，訊號還在嗎？
9. 寫「發表路徑」評估（見下節，300–500 字）。

## 繳交物

1. 可重跑的 R 專案（renv + set.seed(1234)）＋ `decisions.md`（第一條必須是檔案查證程序）。
2. 圖：註解對答案混淆矩陣、inferCNV 圖（含正常 ductal 假陽性檢查）、免疫沙漠病人排序圖、CAF 亞型與通訊圖組（英文標籤）。
3. 分析筆記：各階段回答＋參數三件事＋非 GEO 資料的 Data availability 草稿。
4. 發表路徑評估（300–500 字）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- 你的訊號屬於指南第二節的哪一類？以病人為單位站得住嗎（免疫沙漠排序類的結論特別容易被單一病人撐起）？
- **驗證設計**：外部驗證首選 **GSE155698（Steele 2020）**——同癌種、獨立 cohort、規模更大，你的沙漠量化與 CAF 亞型流程可原樣重跑；bulk 端接 **TCGA-PAAD**（沿進階10 的反卷積 + 存活流程），看你的訊號與存活是否掛鉤。
- **novelty 定位**：PDAC 的免疫沙漠與 CAF 亞型是紅海——PubMed 查你的具體角度（不沙漠病人的共同程式？正常導管異質性對惡性判定的影響？）是否已被做過；方法學角度（reference 結構如何影響 CNV 判定）在指南第一節屬「新方法」切入點，別忽略它。
- 若要成文：缺哪一塊（Steele 重現／TCGA-PAAD 掛鉤／機制故事）？兩個獨立 cohort + bulk 驗證的計算發現，對照指南第四節評估合理層級。

## 自我檢核點

- [ ] 檔案查證有三個獨立證據，程序寫進 decisions.md 第一條
- [ ] 轉檔後做過硬檢查，X 與 raw counts 的判定有證據
- [ ] 惡性判定用了正常 ductal 當 reference，且報告了正常樣本的假陽性檢查
- [ ] 免疫沙漠與 CAF 結論以病人為統計單位，通訊結果經過批判性過濾
- [ ] 階段 C 訊號通過至少三種穩健性檢查
- [ ] 發表路徑評估指名 GSE155698/TCGA-PAAD 的具體驗證動作，Data availability 寫清楚出處鏈

## 提示（卡關再看）

<details><summary>提示 1：zellkonverter 路線</summary>
`sce <- zellkonverter::readH5AD("file.h5ad")` 得 SingleCellExperiment；`assayNames(sce)` 看有哪些矩陣、`counts` 在不在。`as.Seurat(sce)` 或手動 `CreateSeuratObject(counts = assay(sce, "X"))`——但先確認 X 的尺度再決定放 counts 還是 data slot，此決策進 decisions.md。SeuratDisk 路線（`Convert` + `LoadH5Seurat`）遇到版本相容問題時，換 zellkonverter 通常更穩。
</details>

<details><summary>提示 2：沙漠的量化</summary>
「免疫沙漠」不是一個數字就能定案：T/NK 佔比、CD8 effector 分數、Treg 與抑制性 myeloid 比例至少三個維度都算，病人排序在不同維度下一致才叫穩。排序圖建議畫成 per-patient 的多指標熱圖。
</details>

<details><summary>提示 3：myCAF/iCAF 不是二分</summary>
兩亞型最初來自 PDAC 文獻，但實際資料裡常是連續譜加中間態。先打分再看分布，別急著硬分兩群——「這個 cohort 裡 myCAF/iCAF 是離散還是連續」本身就是可以寫進筆記的觀察。
</details>

## 進階挑戰

- 走替代路線把 GSE155698 完整跑一遍，兩個 cohort 的免疫沙漠排序與 CAF 亞型結構做正式的一致性分析——這就是發表路徑的驗證章節初稿。
- 沿進階10 流程用你的型別註解建 signature，反卷積 TCGA-PAAD 接存活——PDAC 的預後極差，看看微環境組成能不能分出層次。

## 參考文獻

- Peng J, et al. Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma. *Cell Research* (2019). GSA: CRA001160（PRJCA001063）；處理後資料：Zenodo record 3969339.
- Steele NG, et al.（驗證資料集，GEO 替代路線）GEO: GSE155698.
