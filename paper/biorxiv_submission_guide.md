# bioRxiv 预印本提交指南

> **今日流程：** 先交 bioRxiv → **同一天** 公开推送 GitHub。  
> 一页行动清单见 [`PREPRINT_TODAY.md`](PREPRINT_TODAY.md)；粘贴字段见 [`biorxiv_paste/`](biorxiv_paste/)。

## 提交信息一览

| 项目 | 内容 |
|------|------|
| **标题** | ADA Predictor: An Integrated Sequence-Intrinsic Immunogenicity Risk Predictor for Protein Therapeutics with a 289-Drug Clinical Reference Library |
| **学科分类** | Bioinformatics |
| **稿件文件** | `ADA_manuscript_Bioinformatics.docx`（英文稿） |
| **摘要词数** | ~198 词（bioRxiv 上限 300 词） |
| **代码仓库** | https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor |
| **许可（推荐）** | CC BY 4.0 |
| **利益冲突** | None |
| **资助** | None |

## 作者信息

| # | 姓名 | 单位 |
|---|------|------|
| 1 | Kaitong Li | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 2 | Min Tang | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 3 | Feilin Wu | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 4 | Mei Yang | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 5 | Lijie Jiang | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 6 | Guoyou Liu | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 7 | Zonglei Zhou | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 8 | Peng Guo | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 9 | Jiaona Guang | JOINN Medical Testing Laboratories (Beijing) Co., Ltd. |
| 10 | Jingwen Dong | China Pharmaceutical University |
| 11 | Zhen Zhang | China Pharmaceutical University |
| 12 | Lin Pang | JOINN Pharmaceutical Quality Research and Testing (Beijing) Co., Ltd. |
| 13 | Chunming Rao | JOINN Pharmaceutical Quality Research and Testing (Beijing) Co., Ltd. |

**通讯作者：** Kaitong Li \<likaitong@joinn-lab.com\>

单位完整写法（含城市/国家）见 `biorxiv_paste/authors.csv`。

## 摘要（直接复制；亦见 `biorxiv_paste/abstract.txt`）

Immunogenicity — the development of anti-drug antibodies (ADA) against therapeutic proteins — drives reduced efficacy, altered pharmacokinetics, and adverse events. Existing tools address single axes (MHC binding or B-cell epitopes) in isolation, and none integrate sequence-derived T/B-cell epitope signals with antibody CDR/germline calibration and a large clinical ADA reference. We present ADA Predictor, a web platform on the immunopred engine that combines MHC-I/II epitope prediction, B-cell linear epitope consensus scoring, antibody CDR delimitation and germline V-gene alignment, multi-signal risk flags with capped uplift, and an independent NAb risk score; a 289-biologic reference library (IDC-DS-V1) provides clinical ADA/NAb rates for background comparison without participating in the sequence score. For infliximab (chimeric anti-TNFα), ADA Predictor reported a sequence-intrinsic risk of 58.2% (High) versus a clinical ADA of 34.23%; for adalimumab (fully human), 14.4% (Low) against a clinical ADA of 36.67% — illustrating that clinical ADA incidence is driven by non-sequence factors (dosing, aggregation, assay) and cannot be inferred from the molecule alone, while CDR/germline calibration distinguished chimeric, humanized, and fully human lineages. ADA Predictor offers transparent, reproducible, sequence-intrinsic risk assessment that complements clinical immunogenicity studies; decoupling molecular from clinical determinants enables honest early-stage candidate screening.

## 关键词

immunogenicity prediction; anti-drug antibody; MHC epitope; B-cell epitope; antibody CDR; germline; neutralizing antibody; protein therapeutics

## 提交步骤

1. 访问 https://www.biorxiv.org/ 并登录/注册  
2. 点击 "Submit a preprint"  
3. 上传：`ADA_manuscript_Bioinformatics.docx`  
4. 标题 / 摘要 / 关键词：从 `biorxiv_paste/` 复制  
5. 作者按上表顺序；通讯作者邮箱 `likaitong@joinn-lab.com`  
6. Subject：**Bioinformatics**  
7. Conflict of Interest：**None**；Funding：**None**  
8. License：推荐 **CC BY 4.0**  
9. 提交  
10. **立刻** 推送 GitHub（见 `PREPRINT_TODAY.md` §B），保证筛查期内链接可打开  

## 与 GitHub 的时间线（按你的计划）

1. **今天上午/下午**：提交 bioRxiv  
2. **提交后马上**：`git push` 公开仓库 +（可选）Release `v0.1.0`  
3. **预印本上线后（约 1–2 工作日）**：回填 DOI → `CITATION.cff` / `README.md` / 下方 Cover Letter  
4. **之后**：投 *Bioinformatics*，Cover Letter 写明预印本 DOI  

> 说明：bioRxiv 筛查期间会抽查 Data Availability。若 GitHub 仍为空或私有，可能延误发布——故「交完预印就推」比「隔天再推」更稳。

## 注意事项

- bioRxiv 审核通常 1–2 个工作日，通过后立即发布并获得 DOI  
- 发布后可在 Author Area 直接 "Submit to Journal" 选择 Bioinformatics  
- 预印本可随时更新版本（v2, v3…）  
- 中文稿仅供内部/单位存档，**不要**作为 bioRxiv 主文件上传  

## 后续投稿 Bioinformatics 用 Cover Letter 模板

```
Dear Editors,

We are pleased to submit our manuscript entitled "ADA Predictor: An Integrated 
Sequence-Intrinsic Immunogenicity Risk Predictor for Protein Therapeutics 
with a 289-Drug Clinical Reference Library" for consideration as an 
Original Paper in Bioinformatics.

This work presents ADA, an integrated immunogenicity prediction platform 
that addresses a critical gap in current tools: the fragmentation of 
T-cell, B-cell, and antibody-humanness analyses, and the methodological 
temptation to overfit sequence features directly to clinical ADA incidence. 
By explicitly decoupling sequence-intrinsic risk from clinical ADA rates, 
ADA provides honest, reproducible early-stage candidate screening.

We believe this work is well-suited to Bioinformatics because:
- It addresses the journal's core scope of sequence analysis and 
  immunoinformatics.
- It builds on and extends methods previously published in Bioinformatics 
  (PickPocket, POPI, EpiDope, Hu-mAb, BioPhi), providing an integrated 
  alternative to single-axis tools.
- All code and data are openly available at https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor.

A preprint of this manuscript has been deposited on bioRxiv 
[DOI: to be filled after preprint publication].

We declare no competing interests. No funding was received for this work.

Thank you for your consideration.

Sincerely,
Kaitong Li
On behalf of all co-authors
likaitong@joinn-lab.com
```
