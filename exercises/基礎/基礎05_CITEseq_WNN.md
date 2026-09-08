# 基礎 05 · RNA 之外：CITE-seq 與 WNN

**難度**：★★ ｜ **預估時間**：1 個工作天 ｜ **對應集數**：B4、B18

## 背景與研究主題

CITE-seq 在同一顆細胞上同時量 RNA 與表面蛋白（ADT）。流式細胞儀時代的免疫學註解本來就是建立在表面蛋白上的——很多 T 細胞亞群（CD4 vs CD8、naive vs memory）在 RNA 層面表達量低又吵，蛋白層面卻一刀兩斷。WNN（weighted nearest neighbor）讓兩個 modality 各自發揮：每顆細胞自己學出「該多信 RNA、多信蛋白」的權重。你的研究主題：

1. RNA-only、ADT-only、WNN 三種分群，各自看到（與看不到）什麼？
2. 哪些族群是「RNA 分不開、蛋白分得開」的？反過來的有沒有？
3. modality weight 在不同細胞型別間怎麼變化？這告訴你什麼？

## 資料集

- **SeuratData `bmcite`**（人類骨髓 CITE-seq，~30,672 cells，25 個 ADT）
- 下載：R 內 `install.packages("SeuratData", repos="https://seurat.nygenome.org")` 後 `InstallData("bmcite")`；`LoadData("bmcite")` 得到含 `RNA` 與 `ADT` 兩個 assay 的 Seurat 物件。
- 注意：ADT 只有 25 個蛋白，是「挑過的 panel」不是全轉錄組——它的正規化方式（CLR）與降維參數（nPC 不能太多）都跟 RNA 不同，這正是 B4 說「不同 modality 有不同統計性質」的實例。物件內附作者註解欄位，留到階段 C 對答案。

## 任務

### 階段 A：兩個 modality 各自跑（對應 B4、B5–B9）

1. RNA assay：標準流程（NormalizeData → HVG → ScaleData → PCA），紀錄 nPC 的三件事。
2. ADT assay：`NormalizeData(..., normalization.method = "CLR", margin = 2)` → ScaleData → PCA（25 個蛋白，nPC 上限是多少？想一下再訂）。
3. 各自跑 UMAP + 分群（resolution 的三件事），各出一張 UMAP，先用幾個熟悉的 marker（CD3、CD4、CD8、CD14、CD19 的 RNA 與 ADT 版本）粗略定位主要族群。

### 階段 B：WNN 聯合分群（對應 B18）

4. `FindMultiModalNeighbors`（兩個 modality 各給哪些 dims？寫理由）→ 基於 wsnn 的分群 → `RunUMAP(nn.name = "weighted.nn")`。
5. 三套分群（RNA-only、ADT-only、WNN）互相對照：用列聯表或 alluvial 圖看群與群的對應關係。找出至少一個「RNA 一團、WNN 拆成多群」的案例，用 ADT 的 FeaturePlot 證明拆分是真的（例如 T 細胞裡的 CD4/CD8/memory 亞群、CD8 naive vs effector）。
6. 檢查有沒有反向案例：ADT panel 只有 25 個蛋白，哪些族群是 RNA 分得開、蛋白分不開的（例如 progenitor 亞群、pDC vs cDC）？

### 階段 C：權重解讀與註解定稿（對應 B18、B10）

7. 畫 `RNA.weight` 的 VlnPlot（按 WNN 分群分組）。回答：哪些族群更信蛋白、哪些更信 RNA？跟你在 5、6 找到的案例對得起來嗎？
8. 給 WNN 分群做最終註解（RNA marker + ADT 一起用），與物件附的作者註解對答案；比較三種分群各自能穩定分辨的亞群數，寫 200–300 字總結：「多一個 modality」在這份資料上實際買到了什麼解析度、代價是什麼。

## 繳交物

1. 可重跑的 R 專案（renv + set.seed(1234)、相對路徑）。
2. 圖：RNA/ADT/WNN 三張 UMAP、分群對應圖（列聯表或 alluvial）、關鍵案例的 ADT FeaturePlot、RNA.weight VlnPlot（英文標籤）。
3. 分析筆記：每個參數三件事＋各階段回答。
4. 200–300 字的 modality 價值總結。

## 自我檢核點

- [ ] ADT 的正規化用 CLR 且能說出為什麼不能照 RNA 的 LogNormalize 做
- [ ] 兩個 modality 的 nPC/dims 各自有理由，不是同一組數字複製貼上
- [ ] 至少一個「蛋白拆開 RNA 拆不開」的案例有圖為證，且有檢查過反向案例
- [ ] 權重解讀與案例互相呼應（更信蛋白的族群正是靠蛋白拆開的那些）
- [ ] 總結誠實提到 25-plex panel 的限制（panel 沒放的蛋白就是看不到）

## 提示（卡關再看）

<details><summary>提示 1：CLR 與 margin</summary>
ADT 有很強的細胞間背景差異，CLR（centered log-ratio）是 Seurat 對 ADT 的建議正規化；`margin = 2` 表示按細胞做。ADT 的 PCA 用全部 25 個蛋白當 features（不做 HVG 篩選），dims 給到 18–25 之間自己看 ElbowPlot 訂。
</details>

<details><summary>提示 2：WNN 呼叫</summary>
`FindMultiModalNeighbors(bm, reduction.list = list("pca", "apca"), dims.list = list(1:30, 1:18))`（數字是示意，換成你自己有理由的版本）。分群用 `FindClusters(graph.name = "wsnn")`，UMAP 用 `RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap")`。
</details>

<details><summary>提示 3：找拆分案例</summary>
先鎖定 T 細胞大群：RNA-only 常把 CD4/CD8 memory 揉在一起；WNN 分群後看 adt_CD4、adt_CD8a、adt_CD45RA、adt_CD45RO 的 FeaturePlot，界線通常乾淨得多。列聯表用 `table(rna_clusters, wnn_clusters)` 就夠。
</details>

## 進階挑戰

- 對同一群細胞比較 RNA 的 CD8A 與 ADT 的 CD8a 訊號分布，量化「RNA dropout、蛋白清楚」的程度（例如兩者的雙峰分離度）。
- 把 resolution 掃一個範圍（0.4–2.0），看三種分群的亞群數成長曲線：WNN 是不是在同 resolution 下穩定給出更多「可註解」的群，而非只是更多群？
- 這裡練出的免疫亞群精細註解功力（尤其 T 細胞），在進階09（肝癌浸潤 T 細胞）會直接派上用場。

## 參考文獻

- Stuart T, Butler A, et al. Comprehensive Integration of Single-Cell Data. *Cell* (2019).（bmcite 資料來源）
- Hao Y, et al. Integrated analysis of multimodal single-cell data. *Cell* (2021).（WNN 方法）
