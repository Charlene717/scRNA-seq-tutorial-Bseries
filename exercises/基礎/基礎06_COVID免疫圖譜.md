# 基礎 06 · 疾病免疫圖譜：重症 COVID-19

**難度**：★★ ｜ **預估時間**：1–1.5 個工作天 ｜ **對應集數**：B12（組成）、B13

## 背景與研究主題

Wilk 等人（2020）在疫情最初幾個月定序了 7 位重症 COVID-19 病人與 6 位健康者的 PBMC，是「疾病 vs 健康免疫圖譜」的教科書級設計。這題的主戰場是 B12 反覆強調的那件事：**組成比較的統計單位是病人，不是細胞**——幾萬顆細胞會給你美到不真實的 p 值，但你的 n 其實是 13。你的研究主題：

1. 重症 COVID-19 病人的 PBMC 組成，跟健康者差在哪裡？哪些差異撐得起統計檢定？
2. 在組成之外，同一型細胞的「狀態」有什麼疾病相關變化（DE + 富集）？
3. 只有 13 個人、而且病人間差異很大——這對你每一條結論的強度有什麼影響？

## 資料集

- **GSE150728**（Seq-Well 產出、10x 型矩陣；7 重症 COVID + 6 健康 PBMC）
- 下載：GEO Supplementary files 的 per-sample matrices（或作者提供的整合物件；用整合物件可以省整合步驟，但階段 A 的 QC 與檢查仍要自己做）。
- 注意：這是 2020 年早期的資料，PBMC 裡出現一些「教科書說不該在 PBMC 的東西」（如發育中嗜中性球相關訊號）——看到怪群不要急著當 doublet 刪掉，先查證。**你自己資料裡看到什麼就寫什麼，不要照論文抄結論。**

## 任務

### 階段 A：整合與註解（對應 B5–B11）

1. 讀入樣本（metadata 帶病人 ID 與疾病狀態）、QC（三件事）、跨樣本整合（Harmony 或 anchors——13 個 PBMC 樣本該不該整合？跟進階02 的腫瘤情境有何不同？想清楚寫下來）。
2. 分群與註解：至少分出 T（CD4/CD8）、NK、B、漿細胞/漿母細胞、單核球（classical/non-classical）、DC、血小板。
3. 對「不太像教科書 PBMC」的群做偵查：marker、QC 指標、在哪些病人出現。給它一個你能辯護的標籤（哪怕是 "developing neutrophil-like, needs validation"）。

### 階段 B：組成分析——統計單位是病人（對應 B12 策略室）

4. 算每位「病人」的細胞型別比例，畫 boxplot（每個點是一個人，分健康/重症兩組）。
5. 檢定組成差異：用 `speckle::propeller`，或手動 per-patient 比例 + Wilcoxon/t 檢定（比例資料的 transformation 要不要做？查 propeller 的做法）。回答：細胞層級卡方檢定 vs 病人層級檢定，p 值差多少？哪個才誠實？
6. 報告至少 2 個顯著（或明顯但不顯著）的組成變化，並檢視是不是被單一極端病人拉出來的（把病人 ID 標在點上）。

### 階段 C：狀態差異與個體差異（對應 B13）

7. 挑 1–2 個細胞型別（例如 CD14 單核球）做**同型別的疾病 vs 健康 DE**：用 pseudobulk（per-patient 加總）做，回想 B13 為什麼細胞層級 DE 的 p 值不可信。
8. 對 DE 結果做富集分析（GO/GSEA，工具自選），描述重症的功能特徵（干擾素反應？HLA class II 下調？以你自己跑出來的為準）。
9. 寫 300–400 字短文：把你的每一條主要結論標上「在幾位病人身上成立」——13 人中若一個訊號只由 3 位病人貢獻，你的結論句要怎麼改寫才不超賣？

## 繳交物

1. 可重跑的 R 專案（renv + set.seed(1234)、相對路徑）。
2. 圖：註解 UMAP（分疾病狀態）、per-patient 組成 boxplot（點標病人）、pseudobulk DE 火山圖 + 富集圖（英文標籤）。
3. 分析筆記：每個參數三件事＋細胞層級 vs 病人層級檢定的對比結果。
4. 個體差異短文（300–400 字）。

## 自我檢核點

- [ ] 整合與否的決策有理由，且能說出這裡跟腫瘤資料整合兩難的差別
- [ ] 怪群沒有被無腦刪掉，偵查過程有紀錄
- [ ] 組成檢定以病人為單位，並實際展示了「細胞層級 p 值 vs 病人層級 p 值」的落差
- [ ] DE 用 pseudobulk，且能說出為什麼細胞當重複會假陽性爆炸
- [ ] 每條主要結論都標了「由幾位病人支撐」，句子的強度與之相稱

## 提示（卡關再看）

<details><summary>提示 1：propeller</summary>
`speckle` 套件（Bioconductor）：`propeller(clusters = ..., sample = patient_id, group = disease_status)`，它會做 logit/asin 轉換再以病人為單位檢定，輸出附 FDR。比較組別時記得 n = 7 vs 6，功效本來就低——「不顯著」不等於「沒差異」。
</details>

<details><summary>提示 2：pseudobulk DE</summary>
對每個「病人 × 細胞型別」把 counts 加總成一個樣本（`AggregateExpression`），對目標型別得到 13 欄的矩陣，丟 DESeq2/edgeR 做 7 vs 6 的比較。設計矩陣裡只有疾病狀態；發現某病人是離群值時，試著拿掉重跑看結論穩不穩。
</details>

<details><summary>提示 3：怪群偵查</summary>
先看它的 top marker 在 Human Protein Atlas / 文獻查是什麼譜系；再看它是不是集中在少數病人、QC 指標是否極端。重症感染時骨髓緊急造血會把未成熟細胞放進血液——「PBMC 不該有」的教科書前提在重症病人身上會失效。
</details>

## 進階挑戰

- 用作者 metadata 中的病人臨床資訊（如通氣狀態，以 GEO 頁提供者為準）把重症組再細分，看組成變化是否有梯度。
- 對漿母細胞群做 BCR 相關基因（IGHG/IGHA 比例）的病人間比較，討論「擴張」在轉錄組層面能推多遠。
- 本題的組成分析思路（per-patient 比例、propeller、pseudobulk DE）在進階04 與進階06 會原封不動地搬進腫瘤免疫場景。

## 參考文獻

- Wilk AJ, et al. A single-cell atlas of the peripheral immune response in patients with severe COVID-19. *Nature Medicine* (2020). GEO: GSE150728.
- Phipson B, et al. propeller: testing for differences in cell type proportions in single cell data. *Bioinformatics* (2022).（工具）
