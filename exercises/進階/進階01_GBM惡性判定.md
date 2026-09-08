# 進階 01 · 這群是不是癌：單一 GBM 的惡性判定

**難度**：★★ ｜ **預估時間**：1–1.5 個工作天 ｜ **對應集數**：B10、B16 ｜ **資料集**：10x Human Glioblastoma 5k（本題獨用）

## 背景與研究主題

拿到一份腫瘤檢體的單細胞資料，第一個要回答的問題永遠是：**哪些細胞是惡性的？** 直覺的做法是找「癌 marker」——但 B16 已經警告過：GBM 的惡性細胞會模仿神經發育譜系，許多「神經幹性 marker」在正常膠質細胞也會亮。這一題用一份單一病人的 GBM 資料，讓你親身驗證那句鐵則：**惡性判定靠 CNV，不是 marker**。你的研究主題：

1. 這份 GBM 檢體裡，腫瘤微環境（免疫、血管）與疑似惡性細胞各佔多少？
2. 用 inferCNV 判定惡性時，reference 怎麼選才站得住腳？**換一組 reference，結論會不會跟著換？**
3. 惡性細胞內部的異質性：這顆腫瘤的惡性亞群，對得上 Neftel 等（2019）的 GBM 四狀態（MES/AC/OPC/NPC-like）嗎？

## 資料集

- **10x Human Glioblastoma 5k**（10x 3' v3），5,604 cells，單一 GBM 檢體。
- 下載：`https://www.10xgenomics.com/datasets` 搜尋「Human Glioblastoma Multiforme」，填 email 後下載 **Filtered feature-barcode matrix（HDF5 或 MEX tar.gz）**。
- 注意：單一病人、無配對正常組織——inferCNV 的 reference 必須從檢體內部找，這正是本題的考點之一。

## 真實數據關卡

真實資料在你跑第一行 Seurat 之前就開始出題。照做並記錄：

1. **格式解剖**：h5 與 MEX 三件套是同一份資料的兩種包裝。若抓 h5：`Read10X_h5()` 需要 `hdf5r` 套件，裝不起來時改抓 MEX（`Read10X()` 讀資料夾）。讀入後先驗貨：`dim()` 的 barcode 數對不對得上官方頁的 5,604？基因名有沒有重複或底線改名（`make.unique` 的痕跡）？「filtered」是 CellRanger 已做過 cell calling 的意思——它替你做了什麼、你還需不需要再做一層 QC，想清楚。
2. **沒有 metadata 可依賴**：GEO 題目有 series matrix、有註解檔；這裡什麼都沒有——單一檢體、單一批次，所有 metadata 都得由你產生（QC 指標、分群、註解、CNV 判定）。這不是輕鬆，是危險：沒有病人欄可以幫你發現批次混淆，也沒有作者答案可對。你的每個判定只能靠證據鏈自我支撐，這正是本題要練的。
3. **腫瘤檢體的 QC 陷阱**：GBM 手術檢體解離壓力大，percent.mt 分布普遍右移，ambient RNA（垂死細胞漏出的 mRNA）也比 PBMC 髒——你會看到免疫群「微量表達」膠質基因。閾值看這份資料自己的分布訂，別照抄基礎01 的 PBMC 數字；同時記錄：你放寬了什麼、憑什麼相信放進來的不是垃圾。
4. **決策日誌**：開一個 `decisions.md`，凡是「資料逼你做的決定」都記一筆：發現了什麼 → 選項 → 你選了什麼 → 理由。這份日誌是繳交物，也是之後論文 Methods 的草稿。

## 任務

### 階段 A：標準流程與微環境註解（對應 B5–B10）

1. 標準流程到 UMAP＋分群（參數三件事照規矩；QC 對腫瘤檢體放寬/收緊了什麼，記進 decisions.md）。
2. 先註解**有把握的非惡性群**：免疫細胞（PTPRC；myeloid/microglia 如 CD14、AIF1/P2RY12；淋巴球 CD3D）、血管內皮（PECAM1/CLDN5）、寡樹突細胞（MBP/PLP1/MOG）。
3. 剩下註解不動的群先標「疑似惡性/待判定」。出一張初步註解 UMAP。回答：此刻你對哪些群的身分最沒把握？為什麼？

### 階段 B：inferCNV 惡性判定（對應 B16）

4. 跑 inferCNV：reference 的選擇寫成一小段理由——你選了誰、為什麼相信它們 CNV 正常、如果選錯會怎樣。
5. 讀 heatmap：GBM 的教科書級變異是 chr7 增加與 chr10 缺失——你的疑似惡性群有沒有？寡樹突與免疫群乾不乾淨？
6. 把 CNV 判定寫回 metadata，出「惡性/非惡性」UMAP，與階段 A 註解交叉比對，產出最終註解總表（cluster × 身分 × 證據）。順手驗證一個 marker 誤導案例（GFAP、SOX2、OLIG2、EGFR 擇一畫 FeaturePlot 對照 CNV 判定）。

### 階段 C：單一腫瘤能問什麼（發表導向）

7. **惡性亞群 × Neftel 四狀態**：對惡性細胞單獨重跑 HVG/PCA/分群，用 Neftel 四狀態基因集 `AddModuleScore` 打分。你的亞群對得上四狀態嗎？有沒有對不上的亞群（例如高壓力/低品質假亞群，或狀態混合的細胞）？
8. **方法學比較：reference 敏感度**：換至少兩組 reference（免疫 only；免疫＋內皮；加或不加寡樹突）重跑 inferCNV，做判定結果的混淆比較。哪些細胞的惡性判定會隨 reference 翻盤？它們是誰？
9. 寫「發表路徑」評估（見下節）——重點是誠實回答：**單一樣本的發現，離可發表還缺什麼？**

## 繳交物

1. 可重跑的 R 專案（renv + `set.seed(1234)`；inferCNV 輸出目錄保留）＋ `decisions.md` 決策日誌。
2. 圖：初步註解 UMAP、inferCNV heatmap、惡性/非惡性 UMAP、Neftel 四狀態打分圖、reference 敏感度比較圖（英文標籤）。
3. 分析筆記：參數三件事＋reference 選擇理由＋最終註解總表。
4. 發表路徑評估（300–500 字）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- **訊號分類與穩健性**：你在階段 C 的發現（某個四狀態組成、或 reference 敏感的細胞群）屬於指南第二節的哪一類？先誠實面對：n=1 病人，「以病人為單位仍成立」這條你根本無法檢查——所以穩健性只能靠方法內部（換 reference、換 resolution、換打分方式）撐。
- **驗證設計**：單一樣本的觀察要成為主張，必須到多病人資料集重現。現成的驗證場就在題庫裡：**進階02（GSE84465，4 位 GBM 病人）**驗 GBM 內的推廣性；**進階03（GSE72056，19 顆黑色素瘤）**驗「reference 選擇影響 CNV 判定」這個方法學主張是否跨癌種成立。
- **novelty 定位**：Neftel 四狀態與 inferCNV 都是已發表的成熟工作——用 PubMed 查「inferCNV reference selection / benchmark」，看方法學比較這條線有沒有人系統做過。單一樣本的生物學發現幾乎不可能是 novelty；但「判定流程對 reference 的敏感度量化」有機會往方法評估型文章走。
- **缺什麼＋目標期刊層級**：誠實寫：這題本身最多是一篇論文的一小節或一個 supplementary 分析。若把 reference 敏感度分析擴到多資料集（進階02、進階03），才夠格投方法比較類（Briefings in Bioinformatics 層級）或輕量觀察類（Scientific Reports 層級）。

## 自我檢核點

- [ ] 讀檔驗貨（barcode 數、基因名、filtered 的含義）與 QC 決策全部進了 decisions.md
- [ ] 微環境註解每群都有 marker 證據，「待判定」群誠實標出而非硬塞名字
- [ ] inferCNV reference 的選擇有理由，且做了至少兩組 reference 的敏感度比較
- [ ] 能在 heatmap 上指出具體染色體事件（如 chr7/chr10）並對應到 cluster
- [ ] Neftel 打分的結論有考慮「假亞群」的可能（壓力/品質驅動的群）
- [ ] 發表路徑評估誠實面對 n=1 的侷限，並指名了驗證資料集

## 提示（卡關再看）

<details><summary>提示 1：inferCNV 輸入</summary>
需要三樣：raw counts 矩陣、細胞註解檔（cell → group，reference 群標明）、gene ordering 檔（官方 wiki 有現成下載）。細胞多時跑得慢，可對每群下採樣或開多執行緒；參數以官方預設起步，改動照三件事記錄。
</details>

<details><summary>提示 2：reference 站不站得住</summary>
好的內部 reference 要滿足：註解證據強（canonical marker 清楚）、與惡性譜系距離遠（免疫 > 內皮 > 寡樹突，越後者越可能被質疑——寡樹突與 GBM 同屬神經譜系，正是階段 C 敏感度比較值得納入它的原因）。
</details>

<details><summary>提示 3：Neftel 打分的品質陷阱</summary>
MES-like 狀態的基因集與缺氧/壓力反應高度重疊——先檢查你的「MES-like 亞群」是不是同時 percent.mt 高、nFeature 低。分不清時把壓力基因集（如 Fos/Jun 家族、熱休克蛋白）另外打分當對照，兩者相關性高就要在 decisions.md 誠實記錄。
</details>

## 進階挑戰

- 用 copykat 或 SCEVAN 重做惡性判定，與 inferCNV（多組 reference 版本）做三方混淆矩陣：三法不合的細胞落在哪裡？
- 把你的 reference 敏感度分析寫成可重用的函式，直接帶去進階02、進階03 重跑——這就是把練習升級成方法學專案的第一步。

## 參考文獻

- 10x Genomics Datasets：Human Glioblastoma Multiforme, 3' v3（官方資料集頁；報告引用時寫明資料集名稱與下載日期）。
- Tirosh I, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* (2016).（inferCNV 方法起源之一）
- Neftel C, et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. *Cell* (2019).
- Darmanis S, et al. *Cell Reports* (2017). GEO: GSE84465.（驗證資料集，見進階02）
