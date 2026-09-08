# 進階 30 · 圖譜時代的註解：把你的資料 map 到人類肺圖譜

**難度**：★★★ ｜ **預估時間**：2 個工作天 ｜ **對應集數**：B10、B11、B18、B19 ｜ **資料集**：HLCA core（CZ CELLxGENE；本題獨用）

## 背景與研究主題

前面二十九題你都在「自己分群、自己註解」。但單細胞領域已經進入圖譜時代：Sikkema 等人（2023）的人類肺圖譜（HLCA）core 用 scANVI 整合了 14 個 datasets、107 人、約 58 萬顆細胞，附帶專家共識的多層級註解。有了這種資源，註解新資料的主流做法翻轉了——**不重算整合，把你的查詢資料 map 到現成參照圖譜上，讓圖譜把標籤投影過來**（scArches／Azimuth／Seurat reference mapping）。這快、可比、標準化，但也埋著本題要你親手挖出的雷：**圖譜註解不是真理**——參照圖譜是健康為主的組成，你的查詢資料若帶著圖譜沒見過的狀態（疾病特有亞群），mapping 只會把它塞進「最不離譜的健康標籤」裡，而且不會告訴你它在說謊。你的研究主題：

1. reference mapping 給出的註解，跟你自己分群註解的結果差在哪？分歧集中在哪些型別、什麼解析度？
2. mapping 的不確定性指標（label transfer score／uncertainty）能不能事先指出「圖譜心虛」的細胞？
3. 把疾病細胞 map 到健康圖譜會發生什麼——找不到家的細胞被安置到哪裡去了？

## 資料集

- **HLCA core**：CZ CELLxGENE（Sikkema 2023 collection）下載 h5ad（**約 5.6 GB**）；584,444 cells（論文數字；portal 顯示 584.9k——引用時寫明出處）／107 人／14 datasets。collection 頁點 Download 取得永久連結，文件未要求登入。
- **只下載、不重算整合**（資料集總覽明訂）：圖譜的價值就在整合已經做完，你要用的是它的 latent space 與註解。
- **查詢集（擇一）**：(a) 從 HLCA core 內挑 1–2 個 dataset 抽出來當「假想查詢集」（標準答案就在圖譜註解裡，適合先驗證流程）；(b) 你在**進階24**（GSE136831 IPF）做過的子集當真查詢集（疾病 vs 健康圖譜的對決，階段 C 主戰場）。建議 a 練功、b 上場。
- h5ad 進 R 需 zellkonverter/SeuratDisk 轉檔；或全程 Python（scanpy + scArches）。磁碟預留 15–20 GB。

## 真實數據關卡

真實資料在你跑第一行程式之前就開始出題。這一節照做並記錄：

1. **5.6 GB 的 h5ad 不是用 read 全載的**：下載（`wget -c` 斷點續傳）後，先用 Python `anndata.read_h5ad(..., backed="r")` 以 backed mode 開檔——只載 metadata 不載矩陣。偵查 `obs` 的欄位：dataset 來源、供體、多層級註解欄（實際欄位名以檔案為準，逐欄看過並記下每欄是什麼）、還有 `obsm` 裡的整合 embedding。**先在 Python 端按 obs 條件取好子集、寫出小 h5ad，再考慮進 R**——在 16–32 GB 的機器上全量載入等於自殺。
2. **轉檔是有損壓縮**：h5ad → Seurat 的每條路（zellkonverter、SeuratDisk、自己拆 h5）都可能掉東西：多層 obsm 的 embedding、uns 裡的 colormap 與模型資訊、raw vs X 的尺度混淆（X 是 counts 還是 lognorm？轉完再驗一次數值分布）。轉檔前後各做一次盤點（細胞數、基因數、metadata 欄、embedding 有沒有跟過來），損耗清單進 decisions.md。全程 Python 可避開轉檔，但你的 B 系列工具箱在 R——這個取捨本身就是決策。
3. **查詢集的「乾淨切割」**：走路線 a 時，把查詢 dataset 從參照裡**真正拿掉**（參照不能包含查詢自己，否則是作弊的自我 mapping）；走路線 b 時，確認 HLCA core 與 GSE136831 的關係——查 HLCA 論文的 dataset 清單，確認你的查詢資料沒有被收進 core（若有重疊，處理方式與理由寫清楚）。這一步偷懶，整題的結論都是循環論證。
4. **決策日誌**：開 `decisions.md`，凡是「資料逼你做的決定」（backed mode 偵查發現、子集與轉檔路線、損耗清單、查詢/參照切割）都記一筆：發現了什麼 → 選項 → 你選了什麼 → 理由。這是繳交物，也是論文 Methods 的草稿。

## 任務

### 階段 A：圖譜偵查與查詢集準備（B18）

1. backed mode 偵查 HLCA core：畫出（或表列）dataset × 供體 × 細胞數的結構、註解層級各有幾類。挑好你的假想查詢 dataset 並完成乾淨切割，參照端與查詢端各存一份工作檔。
2. 對查詢集做一次「傳統路線」：自己 QC → 分群 → 用 marker 手動註解（B10 的功夫，不准偷看圖譜標籤）——這份「自力註解」是階段 B 分歧分析的對照組。
3. 選定 mapping 工具並寫理由：scArches（scANVI surgery，正宗 HLCA 路線、全 Python）、Azimuth（有現成肺 reference，自行查證其版本與 HLCA 的關係）、或 Seurat anchor-based label transfer（R 原生、但不是用 HLCA 的 latent space——它在你自己算的空間裡對齊，這個差異要寫進筆記）。

### 階段 B：mapping 與分歧解剖（B10、B11）

4. 執行 mapping，把圖譜標籤（至少兩個層級：粗與細）與信心分數投影到查詢細胞上。
5. 三方對照：圖譜投影標籤 vs 你的自力註解 vs （路線 a 才有的）原圖譜標準答案。列聯表＋分歧地圖：分歧集中在哪些型別？是解析度問題（圖譜分得細你分不出）還是真衝突（兩邊指向不同譜系）？挑 2–3 個分歧案例用 marker 表達當裁判，判給誰、為什麼。
6. 檢驗不確定性指標：mapping 信心低的細胞是不是正好落在分歧熱區／型別邊界？畫「信心分數 vs 是否分歧」的關係——如果信心指標抓不到分歧，它就不能當防線，這個結論本身很重要。

### 階段 C：圖譜的邊界——疾病細胞的審判（發表導向；B19）

7. 把**進階24** 的 IPF 子集（含你在該題註解過的疾病相關族群，特別是 aberrant basaloid 一類健康肺沒有的細胞）map 到 HLCA core。追蹤這些「圖譜裡沒有家」的細胞被安置到哪些健康標籤下、信心分數如何。
8. 寫出「圖譜註解不是真理」的證據報告：健康圖譜怎麼安置 IPF 特有狀態？如果一個不知情的使用者直接採信 mapping 標籤，會做出什麼錯誤結論？uncertainty 指標這次救得了他嗎？穩健性：換一個 mapping 工具（或換信心閾值）重跑第 7 步，被誤安置的模式一致嗎？
9. 寫「發表路徑」評估（見下節，300–500 字）。

## 繳交物

1. 可重跑的專案（R 部分 renv + set.seed(1234)；Python 部分附 environment.yml 與腳本）＋ `decisions.md` 決策日誌（含轉檔損耗清單、查詢/參照切割證明）。
2. 圖：HLCA 結構總覽圖、三方註解列聯表／alluvial、分歧案例 marker 裁判圖、信心分數 vs 分歧關係圖、IPF 細胞的安置去向圖（英文標籤）。
3. 分析筆記：各階段回答＋每個參數三件事。
4. 發表路徑評估（300–500 字）。

## 發表路徑

讀 `_從練習到投稿指南.md` 後回答：

- 你的發現屬於指南第二節的哪一類？本題天然長出兩種：mapping 方法的系統性行為（新方法/方法比較型）或「健康圖譜安置疾病細胞」的具體失效模式（對圖譜使用社群有用的警示型）。
- **驗證設計**：核心驗證場就在題庫內——**進階24**（GSE136831 IPF）：aberrant basaloid 在健康圖譜裡沒有家，是「mapping 失效」最好的正面教材；反向驗證用路線 a 的假想查詢集（標準答案在圖譜內，可量化 mapping 的基準錯誤率）。兩個方向都做才完整。
- **novelty 定位**：HLCA 原論文自己就示範過疾病資料的 mapping 與新狀態偵測——查它與 scArches 相關文獻，確認你的失效案例分析角度（特定亞群的安置去向？uncertainty 指標的盲區？）沒被做過；「重複圖譜團隊自己的示範」不是發表。
- 若要成文：這類工作合理落點是計算/組學期刊或資源評測型文章（指南第四節第二、三層）；缺的通常是多個查詢資料集的系統性重複與多工具比較——誠實列出。

## 自我檢核點

- [ ] backed mode 偵查、轉檔損耗清單、查詢/參照乾淨切割全部進了 decisions.md
- [ ] 自力註解在看 mapping 結果之前完成，三方對照誠實呈現（包括自己輸的格子）
- [ ] 分歧案例有 marker 層級的裁判，不是只報列聯表
- [ ] 信心/uncertainty 指標與分歧的關係有量化，結論不超賣指標的能力
- [ ] IPF 安置分析換過工具或閾值，失效模式的描述有重現性
- [ ] 發表路徑評估查過 HLCA 與 scArches 文獻，novelty 判斷有依據

## 提示（卡關再看）

<details><summary>提示 1：backed mode 與子集落地</summary>
`adata = anndata.read_h5ad(path, backed="r")` 後，`adata.obs` 隨便玩都不吃記憶體；取子集用 `adata[mask].to_memory()` 再 `.write_h5ad(...)`。先把 obs 的註解欄位名全部 `adata.obs.columns` 印出來抄進筆記——HLCA 的欄位很多，哪一欄是共識註解、哪一欄是原研究註解，搞錯了後面全錯。
</details>

<details><summary>提示 2：三條 mapping 路線的真實差異</summary>
scArches 是把查詢資料「手術」進 scANVI 模型的 latent space——用的是 HLCA 真正的整合空間；Seurat 的 `FindTransferAnchors`/`MapQuery` 則是在 anchor 空間重新對齊，沒有用到 HLCA 的模型。兩者標籤可以都合理但空間意義不同。你只需走通一條，但要在筆記寫清楚你那條路「用了圖譜的什麼、沒用到什麼」。
</details>

<details><summary>提示 3：追蹤沒有家的細胞</summary>
對 IPF 查詢集，把你在進階24 給的標籤當列、mapping 給的健康標籤當欄做列聯表——aberrant basaloid 那一列的分佈就是「安置去向」。再疊上信心分數：如果它被高信心地塞進某個健康上皮標籤，這就是本題最有力的一張警示圖。
</details>

## 進階挑戰

- 用 scArches 路線時，試 HLCA 論文示範的新狀態偵測思路（查詢細胞到參照鄰居的距離／uncertainty 分佈），看它能否自動標紅 aberrant basaloid。
- 把 mapping 註解與自力註解各自餵給同一個下游分析（如組成比較），看「註解來源」這個上游決策對下游結論的擾動有多大。
- 拿題庫其他組織的資料重演本題（如把基礎06 的 PBMC map 到 Azimuth PBMC reference），寫 200 字：組織不同、圖譜成熟度不同，mapping 的可信度怎麼變？

## 參考文獻

- Sikkema L, et al. An integrated cell atlas of the lung in health and disease. *Nature Medicine* (2023).（HLCA；CZ CELLxGENE collection）
- Lotfollahi M, et al. Mapping single-cell data to reference atlases by transfer learning (scArches). *Nature Biotechnology* (2022).（工具）
- Adams TS, et al. Single-cell RNA-seq reveals ectopic and aberrant lung-resident cell populations in idiopathic pulmonary fibrosis. *Science Advances* (2020). GEO: GSE136831.（進階24／本題查詢集）
