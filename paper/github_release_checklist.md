# GitHub 仓库发布清单

**目标仓库：** https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor

手稿 Data Availability、Cover Letter 与 `CITATION.cff` 均已指向该地址。

**与预印本的顺序：** 先交 bioRxiv（见 `PREPRINT_TODAY.md`），**提交后同一天**完成本清单推送，确保筛查期内链接可打开。

## 一、首次推送前（本地）

- [x] `LICENSE` — MIT
- [x] `README.md` — 安装 / 结构 / 论文与数据说明
- [x] `CITATION.cff` — GitHub “Cite this repository”
- [x] `paper/` — 中英手稿（md + DOCX）+ 基准 CSV + bioRxiv 指南
- [x] `.gitignore` — 排除运行库、mhcflurry 运行时、Word 锁文件、临时稿
- [ ] 配置 remote（若尚未配置）：

```powershell
cd D:\Project\ADA2
git remote add origin https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor.git
git branch -M main
```

- [ ] 首次提交并推送（需你确认后再执行）：

```powershell
git add -A
git status   # 确认无密钥、无超大无关文件
git commit -m "Initial public release of ADA Predictor"
git push -u origin main
```

**建议纳入公开仓库：**

| 路径 | 说明 |
|------|------|
| `ada/` | 核心包（含 `immunopred`） |
| `data/ada_reference.csv` 等参照表 | 289 药临床参照 |
| `data/germlines/` | IMGT V 基因池 |
| `paper/*.md` / 正式 DOCX / `benchmark_*.csv` | 手稿与可复现基准 |
| `paper/_iedb_cache.json` | IEDB 复现缓存（~0.7 MB） |
| `tests/` `tools/`（除 get-pip / py310） | 测试与校准脚本 |

**不要推送：** `.trae/`、`data/ada.db`、`tools/py310/`、`paper/~$*.docx`、`*_tmp.docx`、旧版 `*-3.docx`（已在 `.gitignore`）。

## 二、GitHub 网页设置（登录后）

1. **About**
   - Description: `Sequence-intrinsic immunogenicity (ADA) risk predictor for protein therapeutics — MHC-I/II, B-cell epitopes, CDR/germline, 289-drug clinical reference`
   - Website: （可选）部署后的演示 URL
   - Topics: `immunogenicity` `anti-drug-antibody` `bioinformatics` `mhc` `epitope-prediction` `fastapi` `python`
2. **Settings → General**
   - Features: Issues ✅ / Discussions（可选）
   - Releases: 推送后打 `v0.1.0`
3. **Settings → Pages**（可选）：若日后做静态文档站再开
4. **Social preview**：可上传一张 Web UI 截图

## 三、与预印本联动

1. bioRxiv 发布获得 DOI 后：
   - 更新 `CITATION.cff` 中 `preferred-citation.doi`
   - 更新 `README.md` 徽章与 “Citation” 小节
   - 更新手稿 Data Availability（若需）
2. （推荐）[Zenodo GitHub 集成](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)：授权后 Create release → 自动得软件 DOI，可写入手稿与 bioRxiv。

## 四、Release 说明模板（v0.1.0）

```
## ADA Predictor v0.1.0

First public release accompanying the bioRxiv preprint.

### Highlights
- immunopred engine: MHC-I/II, B-cell consensus, CDR/germline, NAb score
- 289-drug clinical ADA/NAb reference (background only; not used in scoring)
- Bilingual (zh/en) web UI (FastAPI)
- Reproducible benchmarks under paper/

### Cite
See CITATION.cff (update with bioRxiv DOI after preprint is live).
```
