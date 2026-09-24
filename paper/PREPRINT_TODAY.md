# 今日操作清单：先 bioRxiv → 马上推 GitHub

按你的顺序：**先交预印本，同一天立刻推公开仓库**。bioRxiv 筛查约 1–2 个工作日，在此期间仓库必须已可访问（稿件 Data Availability 写了 GitHub 链接）。

---

## A. bioRxiv 提交（现在）

### 上传文件
| 用途 | 文件 |
|------|------|
| **主稿（必须）** | `paper/ADA_manuscript_Bioinformatics.docx` |
| 中文对照（不必上传） | `paper/ADA_manuscript_Bioinformatics_zh.docx` |

### 表单粘贴（已拆好，直接打开复制）
| 字段 | 文件 |
|------|------|
| Title | [`biorxiv_paste/title.txt`](biorxiv_paste/title.txt) |
| Abstract（~198 词，上限 300） | [`biorxiv_paste/abstract.txt`](biorxiv_paste/abstract.txt) |
| Keywords | [`biorxiv_paste/keywords.txt`](biorxiv_paste/keywords.txt) |
| Authors | [`biorxiv_paste/authors.csv`](biorxiv_paste/authors.csv) |

### 表单其它选项（照填）
| 项 | 填法 |
|----|------|
| Subject area | **Bioinformatics** |
| Corresponding author | Kaitong Li · likaitong@joinn-lab.com |
| Conflict of interest | **None** / No competing interests |
| Funding | **None** / No funding to declare |
| License | 推荐 **CC BY 4.0**（与数据许可一致；也可用 CC BY-NC） |
| New results? | Yes（工具/方法论文） |
| Human subjects / clinical trial | **No** |
| Dual use / pathogens | **No** |

### 提交步骤（简）
1. https://www.biorxiv.org/ → Submit a preprint  
2. 上传 `ADA_manuscript_Bioinformatics.docx`  
3. 粘贴 title / abstract / keywords；按 CSV 顺序录入 13 位作者与 3 个单位  
4. Subject = Bioinformatics；COI / Funding 如上  
5. 选 CC BY 4.0 → Submit  

更细的说明见 [`biorxiv_submission_guide.md`](biorxiv_submission_guide.md)。

---

## B. 提交后立刻推 GitHub（同一天）

仓库：https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor  

本地已配置 `origin`。在 Cursor 里说一声 **「推送 GitHub」**（或自己执行）：

```powershell
cd D:\Project\ADA2
git add -A
git status
git commit -m "Initial public release of ADA Predictor (v0.1)"
git push -u origin main
```

然后在 GitHub 网页：
1. About 填 Description + Topics（见 `github_release_checklist.md`）  
2. 打 Release **v0.1.0**（可选但推荐）  
3. 确认仓库为 **Public**，链接可打开  

---

## C. 预印本上线后（1–2 天后）

1. 记下 bioRxiv DOI（形如 `10.1101/2026.xx.xx.xxxxxx`）  
2. 回填：`CITATION.cff`、`README.md` Citation、本目录 Cover Letter 占位  
3. 可选：Zenodo ↔ GitHub 绑软件 DOI  
4. 再投 *Bioinformatics*（Cover Letter 模板在 `biorxiv_submission_guide.md`）  

---

## 自检（提交前 30 秒）

- [ ] 打开 DOCX：标题、作者、单位、通讯作者邮箱正确  
- [ ] Abstract / Keywords 与 paste 文件一致  
- [ ] Funding = None；COI = None  
- [ ] Data Availability 含 `https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor`  
- [ ] 计划：点完 Submit 后 **马上** push（勿隔夜）  
