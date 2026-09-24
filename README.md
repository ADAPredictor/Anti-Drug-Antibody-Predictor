# ADA Predictor — 蛋白药物免疫原性（抗药抗体 ADA）预测平台

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/)
[![GitHub](https://img.shields.io/badge/GitHub-ADAPredictor-181717?logo=github)](https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor)

**Repository:** https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor

对蛋白/多肽药物的 T 细胞表位（MHC-I / MHC-II）、B 细胞表位、聚集倾向、化学风险位点、给药模态先验与序列新颖性进行多因子加权，输出 0–100 的免疫原性风险分与 FDA 已批准药物的同模态对照。

> 仅供科研与早期候选物筛选，**不构成临床或注册决策依据**。ADA 发生率高度依赖检测方法、采样时间、适应症、剂量与合并用药，不可直接跨研究比较。

### English summary

ADA Predictor is an integrated, sequence-intrinsic immunogenicity risk platform for protein therapeutics. It combines MHC-I/II epitope prediction, B-cell linear epitope consensus scoring, antibody CDR delimitation and germline V-gene alignment, capped multi-signal risk uplift, and an independent neutralizing-antibody (NAb) score. A 289-biologic clinical reference library supplies ADA/NAb rates for **background comparison only** and does not enter the sequence score.

## 技术栈

| 模块 | 实现 |
|---|---|
| 核心运行时 | Python 3.7+ |
| Web | FastAPI + Uvicorn |
| MHC-I | 可选 [mhcflurry 2.2.x](https://github.com/openvax/mhcflurry)（Apache-2.0，pan-allele 神经网络，运行于隔离的 Python ≥3.10 子进程）；不可用时回退到 HLA 超型锚定基序启发式 |
| MHC-II | 9 肽核心 P1/P4/P6/P9 口袋疏水性 15 肽基线（预留 NetMHCIIPan 接口） |
| B 细胞 | Parker 亲水性 + Chou-Fasman 转角双条件，阈值用真实 RSA 校准 |
| 数据 | 内置 289 个 FDA 批准蛋白药物 ADA 参照表（IDC-DS-V1，CC BY 4.0），支持 openFDA / FDA SPL 实时抽取 |

## 快速开始

```bash
# 1. 安装核心依赖
pip install -r requirements.txt

# 2. 启动服务（首次运行自动建库并加载种子数据）
python run.py --host 127.0.0.1 --port 8000
```

浏览器打开 <http://127.0.0.1:8000>，粘贴序列即可预测。

## 启用高精度 MHC-I 后端（可选，推荐）

mhcflurry 需要 Python ≥3.10 与 PyTorch，与本项目的 Python 3.7 核心**互不污染**——以常驻子进程 + JSON-lines 协议通信，模型只加载一次。

Windows 一键脚本（会下载 embeddable Python 3.10、安装 mhcflurry==2.2.1、拉取 models_class1_presentation，约 1 GB）：

```powershell
$env:ADA_PIP_INDEX = "https://pypi.tuna.tsinghua.edu.cn/simple"   # 可选，国内加速
powershell -ExecutionPolicy Bypass -File tools\setup_mhcflurry.ps1
```

安装后 `api/meta` 的 `class1_backend` 会变为 `mhcflurry-2.2.x-subprocess`。可用环境变量 `ADA_DISABLE_MHCFLURRY=1` 强制关闭，`ADA_MHCFLURRY_PYTHON` 覆盖解释器路径。

## 项目结构

```
ada/            核心包（pipeline / scoring / epitopes / web / etl / db）
data/           运行期 SQLite 与结构文件（ada.db 由 seed 自动生成，不入库）
paper/          论文手稿（中英 markdown + DOCX）与基准数据
tests/          单元测试
tools/          mhcflurry 一键安装脚本、校准脚本、评估脚本
requirements.txt
run.py          服务入口
```

## 测试

```bash
python -m unittest discover -s tests
```

## 校准脚本（可复现，非 ADA 率拟合）

- `tools/calibrate_bcell_rsa.py` — B 细胞阈值基于真实 RSA 校准
- `tools/calibrate_mhc1_saturation.py` — MHC-I 强结合率饱和尺度（20%）基于 mhcflurry 在 10 个已批准蛋白上的实测结合背景
- `tools/eval_mhcflurry.py` — 10 蛋白双跑（启发式 vs mhcflurry）基准与流感 GILGFVFTL/A*02:01 阳性对照断言

## 论文与数据

本仓库随附论文手稿（`paper/` 目录）：

- `ADA_manuscript_en.md` / `ADA_manuscript_zh.md` — 中英文手稿源文件
- `ADA_manuscript_Bioinformatics.docx` / `ADA_manuscript_Bioinformatics_zh.docx` — Bioinformatics 格式排版
- `benchmark_composite.csv` — 9 药基准验证的权威数值源
- `build_docx.py` — 双语 DOCX 构建脚本
- `biorxiv_submission_guide.md` — bioRxiv 预印本提交信息与 Cover Letter
- `github_release_checklist.md` — GitHub 首次发布清单

**数据可用性**：289 药临床参照表（IDC-DS-V1，CC BY 4.0）随本仓库分发；IEDB 缓存（`paper/_iedb_cache.json`）供复现使用。

### Citation

预印本 DOI 发布前，请引用本仓库（GitHub → *Cite this repository*，或见 [`CITATION.cff`](CITATION.cff)）。获得 bioRxiv DOI 后将更新该文件与本小节。

```
Li K, Tang M, Wu F, et al. ADA Predictor: An Integrated Sequence-Intrinsic
Immunogenicity Risk Predictor for Protein Therapeutics with a 289-Drug
Clinical Reference Library. GitHub.
https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor
```

## 许可

MIT License — 见 [LICENSE](LICENSE)。临床参照表中来自 IDC-DS-V1 的内容遵循 CC BY 4.0。
