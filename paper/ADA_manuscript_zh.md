# ADA Predictor：集成式蛋白药物序列内在免疫原性风险预测平台（含 289 条上市药物临床参照库）

**作者：** 李开通¹，唐敏¹，吴飞林¹，杨梅¹，姜立杰¹，刘国友¹，周宗雷¹，郭鹏¹，光姣娜¹，董靖文²，张振²，庞琳³，饶春明³

**单位：**
¹ 北京昭衍新药检测实验室有限公司（JOINN Medical Testing Laboratories (Beijing) Co., Ltd.），北京
² 中国药科大学（China Pharmaceutical University），南京
³ 北京昭衍药物质量研究检测有限公司（JOINN Pharmaceutical Quality Research and Testing (Beijing) Co., Ltd.），北京

**通讯作者：** 李开通，likaitong@joinn-lab.com

---

## 摘要

免疫原性——即治疗性蛋白诱导患者产生抗药物抗体（ADA）——是疗效降低、药代动力学改变与不良事件的主要诱因。现有工具多仅覆盖单一维度（MHC 结合或 B 细胞表位），尚无工具能将序列层面的 T/B 细胞表位信号与抗体 CDR/种系校准及大规模临床 ADA 参照库整合于同一平台。我们提出 ADA Predictor——基于 *immunopred* 引擎的 Web 平台，集成 MHC-I/II 表位预测、B 细胞线性表位共识评分、抗体 CDR 划分与种系 V 基因比对、多信号风险标志（封顶上浮）与独立的 NAb 风险评分；289 条上市生物药参照库（整合自 IDC-DS-V1 [12]、FDA/EMA 标签、Thera-SAbDab 等公开来源）提供临床 ADA/NAb 率作背景对照，不参与序列评分。对英夫利西单抗（嵌合抗 TNFα），系统给出序列内在风险 58.2%（高），临床 ADA 参照 34.23%；对阿达木单抗（全人源），预测 14.4%（低），临床 ADA 36.67%——表明临床 ADA 发生率受非序列因素（给药、聚集、检测方法）驱动，无法仅凭分子推断，而 CDR/种系校准正确区分嵌合、人源化与全人源谱系。ADA Predictor 提供透明、可复现的序列内在免疫原性风险评估，可补充临床免疫原性研究；将分子风险与临床因素解耦，使早期候选物筛选的信号诚实、不过拟合。

**关键词：** 免疫原性预测；抗药物抗体；MHC 表位；B 细胞表位；抗体 CDR；种系；中和抗体；治疗性蛋白

---

## 1. 引言

治疗性蛋白（尤其是单克隆抗体）已深刻改变了肿瘤、自身免疫病与代谢性疾病的治疗格局。然而，由于这类分子可被患者免疫系统识别为外源物，相当比例的患者会产生抗药物抗体（ADA），导致药物中和、药代动力学改变或输注/超敏反应 [1, 14]。ADA 发生率差异巨大——部分全人源 mAb 低于 1%，而细菌来源酶类可超过 80%——并受药物分子来源、给药途径与频次、合并免疫抑制治疗以及患者 HLA 基因型等多因素影响 [14]。

近年来，大量研究致力于在药物研发早期通过计算方法预测免疫原性，以减少对昂贵、耗时的湿实验的依赖。主流范式包括：

- **MHC 结合预测**（如 NetMHCpan [9]、IEDB MHC [10]）：按预测的 HLA 结合亲和力对肽段排序，聚焦 T 细胞表位维度。
- **B 细胞表位预测**（如 Kolaskar–Tongaonkar 抗原性 [1]、Parker 亲水性 [2]、Emini 表面可及性 [3]）：识别线性表位区段。
- **抗体特异分析**：考察 CDR 异源性与种系邻近度 [5, 6]，认识到治疗性抗体的可变区是免疫原性的主要决定因素。

每一单一维度方法都有可测量的局限。MHC 结合对免疫原性是必要但不充分的：POPI 从 428 条人 MHC I 类结合肽中选出 23 项理化性质做分类，留一法准确率为 64.72%，仅略高于两种序列比对基线的 54.91% 与 53.23% [22]。线性 B 细胞预测即便引入深度学习仍然偏弱：EpiDope 报告的 ROC 曲线下面积为 0.67 ± 0.07 [23]。抗体一侧，人源度分数携带真实但有限的信号：Hu-mAb 基于抗体库的随机森林解释了报告 ADA 方差的 34%（R = −0.58），BioPhi 的 OASis 9 肽流行度分数解释了 28%（R = −0.53）；关键的是，若只统计 OASis 判定为非人源、且被 NetMHCIIpan 判为结合的 9 肽，已解释方差塌缩至 0.01（R = 0.44）[20]，表明把异源度与 MHC 结合简单相乘，非但不会增强、反而会抹掉序列与 ADA 之间的信号。

然而，这些工具彼此割裂：研发人员需分别运行各维度分析并人工整合结果。更关键的是，许多方法试图直接从序列特征预测*临床 ADA 发生率*，将分子内在表位负荷与外在临床因素（检测方法灵敏度、给药方案、人群 HLA 分布）混为一谈，而后者无法仅凭序列确定。

这种割裂并非偶然，而有其历史根源：T 细胞表位、B 细胞表位与抗体人源度分属不同方法学社区，各自在自有指标上迭代优化，缺乏跨维度整合的动机；而临床 ADA 发生率长期以零散、非机器可读的形式散见于药品说明书与文献，直至近年才被系统整理为可计算数据集 [12]，此前难以与序列工具对接。正因如此，将序列特征直接拟合为临床 ADA 发生率的回归目标在方法学上具有诱惑力——样本量小、维度众多、噪声结构复杂——却极易产生过拟合且不可解释的模型，并把来源迥异的方差折叠进同一目标。

我们提出 **ADA Predictor**——一个集成式 Web 平台，其贡献如下：

1. 将 MHC-I、MHC-II、B 细胞表位、抗体 CDR/种系及理化聚集信号整合为单一 0–100% 的**序列内在免疫原性风险**评分。
2. 严格解耦序列衍生风险与临床 ADA 发生率：289 条药物参照库仅提供临床 ADA/NAb 率作*背景对照*，不参与评分计算。
3. 提供透明、可复现的评分与中英双语注释，使计算与非计算背景用户均能理解各信号的贡献。

---

## 2. 系统与方法

### 2.1 系统架构

ADA Predictor 基于 FastAPI 实现，前端为单页应用。后端 *immunopred* 为 Python 包，含 23 个核心模块，负责序列解析、特征提取、校准与评分。前端以卡片式 UI 呈现结果，含可折叠技术解释、评分构成条形图与风险速览表。

参照库（`ada_reference.csv`）含 289 条上市生物药，字段包括：药物名称、INN、药物类别、ADA 率（%）、NAb 率（%）、给药途径、批准状态、序列（抗体含重链/轻链）及文献引用。数据整合自多个公开来源：IDC-DS-V1 临床免疫原性数据集（CC BY 4.0 [12]）、FDA/EMA 药品标签（公共领域监管数据）、Thera-SAbDab 抗体序列及 UniProt/PDB 记录。IDC-DS-V1 中所有治疗性分子里，96 个被纳入；其余因缺乏完整的信息被排除。

### 2.2 MHC-I 表位预测

以 9 肽滑窗扫描输入序列，对标 7 种常见 HLA-A/B 等位基因（HLA-A\*02:01、A\*24:02、A\*03:01、B\*07:02、B\*08:01、B\*35:01、B\*44:02）的口袋偏好矩阵。矩阵编码锚定位点（P2 与 C 端）的位置特异性氨基酸偏好。肽段按蛋白内百分位排序，top 10% 且矩阵分值 ≥0.15 者标记为强结合。

### 2.3 MHC-II 表位预测

采用 15 肽滑窗 + 9 肽核心扫描，对标 6 种 HLA-DRB1 等位基因（DRB1\*01:01、\*03:01、\*04:01、\*07:01、\*11:01、\*15:01）的位置特异性评分矩阵（PSSM）。PSSM 编码 P1、P4、P6、P9 口袋位点的偏好。强结合阈值 = 等位基因特异阈值 × 0.5 且百分位 ≥90。

### 2.4 B 细胞线性表位共识评分

四种经典算法以 z-score 均值方式组合：

| 组成 | 方法 | 文献 |
|---|---|---|
| 抗原性 | Kolaskar–Tongaonkar | [1] |
| 亲水性 | Parker 等 | [2] |
| 表面可及性 | Emini 等 | [3] |
| β-转角倾向 | Chou–Fasman | [4] |

连续 ≥6 氨基酸且共识 z-score ≥0.5 的区段判定为候选 B 细胞表位。需要指出，>90% 的 B 细胞表位为构象表位；线性评分仅用于候选区筛查。

### 2.5 抗体 CDR 划分与种系校准

CDR 由保守 Cys/Trp 残基与 WGXG/FGXG 锚点基序划分（类 Chothia [5]）。可变区序列通过 Biopython [8] 局部比对与内置 IMGT 人源 V 基因库 [6] 比对，获得：

- 框架区（FR）错配分率
- CDR 错配分率
- 种系 V 基因最佳命中（链、基因、一致性百分比）

谱系（嵌合/人源化/全人源/鼠源）由 INN 后缀（-ximab/-zumab/-umab [11, 13]）与种系 FR 一致度双重推断，模态不匹配时向用户告警。

### 2.6 多信号风险标志与封顶上浮

系统评估 9 项风险标志（表 1）。单一 MHC/B 细胞表位密度标志不触发上浮；需种系覆盖或框架距离的交叉印证。上浮封顶 +15%，防止无节制放大。

机制依据在于：单一表位密度信号假阳性偏高——任意蛋白均含若干高抗原性或高结合潜力肽段；而要求该表位同时落在种系错配残基（提示其位于外源区段）或框架远离低 ADA 人源库的背景上，方能提高其被递呈并触发 T 细胞识别的后验概率。封顶机制则避免任一强标志经多次放大而主导总分，使最终风险仍受原始表位负荷锚定。

| 标志 | 信号 |
|---|---|
| t_cdr | CDR 非人源种系表位负荷偏高 |
| b_epi | B 细胞线性表位覆盖偏高 |
| germline_hetero | 相对最佳命中 V 种系的异源残基偏高 |
| t_on_hetero | MHC-II 肽段覆盖种系 CDR 错配残基 |
| b_on_hetero | B 细胞表位覆盖种系 CDR 错配残基 |
| agg | 聚集相关风险偏高 |
| fr_distant | 框架远离低 ADA 人源库 |
| charge | 绝对电荷密度偏高 |
| hetero_cdr | CDR 异源性代理偏高（回退方案） |

### 2.7 评分模型

**第一步——原始表位信号加权：**

$$\text{全长表位负荷} = 0.45 \cdot T_{\text{密度}} + 0.20 \cdot T_{\text{分值}} + 0.20 \cdot B_{\text{覆盖}} + 0.08 \cdot M1_{\text{密度}} + 0.07 \cdot A_{\text{聚集}}$$

**第二步——路径校准：**

- *抗体路径：* 序列异源度（骨架 + FR/CDR 错配 + 格式暴露）× CDR 加权表位负担 × 谱系耐受系数（全人源 0.55、人源化 0.75、嵌合 1.15）→ 基线；叠加多信号上浮（封顶 +15%）。
- *蛋白路径：* 类别特异异源度 × 表位 × 类别耐受系数。
- *启发式路径：* 近内源 ×0.35；人源样抗体 ×0.40。

谱系耐受系数编码可变区来源的先验免疫识别预期：全人源抗体源自人源种系 V 基因、处于免疫耐受库内；人源化抗体保留小鼠 CDR 嫁接于人源框架；嵌合抗体保留完整小鼠可变结构域。该系数与上一步测得的异源度（FR/CDR 错配）相乘，使最终风险同时受"来源先验"与"实测异源度"约束，而非任一单独决定——从而避免仅凭来源归类或仅凭残基比对即定论的片面性。

在施加谱系系数之前，抗体路径将异源度与表位负担分解为带权重的命名分量。

**异源度（0–100）** 混合五个分量：种系骨架距离（权重 0.35）、异源残基比例（0.20）、相对最佳人 V 基因命中的 CDR 错配（0.20）、超出 0.78 种系阈值的超额 CDR 分歧（0.05）、格式暴露（0.20）。

**表位负担（0–100）** 混合 CDR 局部 MHC-II 密度（0.32）、非种系残基负担（0.33）、线性 B 细胞覆盖（0.20）和框架距离项（0.15）。其对最终分数的相关度由异源度门控：

$$\text{相关度} = 0.20 + 0.80 \cdot \frac{\text{异源度}}{100}$$

即可变区越接近人种系，表位负担的贡献越低。基线分值为：

$$\text{基线} = \text{截断}\left(0.45 \cdot \text{异源度} \cdot \lambda + 0.55 \cdot \text{表位} \cdot \text{相关度},\ 1,\ 80\right)$$

对非常接近人源的支架施加软上限：全人源且骨架 < 0.40 时不超过 22；人源化且骨架 < 0.45 时不超过 16；人源化且骨架 < 0.55 时不超过 19；嵌合基线不超过 58。这些上限是结构先验，不是对 289 个临床 ADA 率拟合得到的系数。

**第三步——需相互印证的上浮与机制情境：** 九个风险标志（表 1）贡献封顶上浮：单独的 MHC 或 B 细胞密度标志加 0 分；强 MHC-II 肽或 B 细胞区段与种系错配残基重叠时给出小幅上浮；两个标志同时出现加 7 分；三个及以上加 12 分；总和封顶 +15。全人源和人源化谱系还会把上浮再压低，除非框架距离已经很高。

随后按抗体类别施加单独的情境乘数：CD19/CD20 B 细胞耗竭抗体 ×0.12（上限 6）；双特异、ADC 或 T 细胞衔接器 ×0.48（上限 28）；肿瘤单抗 ×0.40（上限 24）；急性片段 ×0.40（上限 18）。这是序列分数中唯一引入类别水平临床格局的位置；它是手工规则，不是任一药物 ADA 率的回归残差。

**分级：** <20% 低；20–44.9% 中；≥45% 高。

### 2.8 NAb 风险评分

独立的三因子加权评分：

$$\text{NAb} = w_1 \cdot S_{\text{结构}} + w_2 \cdot A_{\text{ADA发生率}} + w_3 \cdot L_{\text{文献}}$$

- **结构/类别因子：** 基线 25 + 格式暴露加成 + 异源度加成 + 类别特异规则（酶替代 +35、细菌/毒素 +40、嵌合 +25、人源化 mAb −10、胰岛素/激素 −15）。
- **ADA 发生率因子：** 由最邻近参照药的临床 ADA 率推导。
- **文献因子：** 来自 IDC-DS-V1 NAb 注释与已发表报告。

分级：低（<25%）、中（25–55%）、高（≥55%）。

### 2.9 参照库与 kNN 匹配

序列匹配按链感知方式进行：
- **Fv（重链 + 轻链）相似度 ≥95%：** 高置信同药命中。
- **仅 VH 或仅 VL ≥85%：** 标注"仅单链框架相似"，并提示人源抗体框架彼此相似度常达 85% 以上。
- **k-mer 局部片段：** 仅作提示。

匹配药物的临床 ADA/NAb 值仅作背景对照展示并用于 NAb 发生率因子，绝不回灌序列风险评分。

### 2.10 与成熟预测器的基准对比

为定位 ADA Predictor 内置口袋/PSSM 启发式方法相对基于神经网络的机器学习预测器的水平，我们针对参照库中 9 个具有已验证 VH/VL 序列的上市 mAb（横跨三种谱系：嵌合——英夫利西单抗、利妥昔单抗、西妥昔单抗；人源化——曲妥珠单抗、贝伐珠单抗、奥马珠单抗；全人源——阿达木单抗、帕尼单抗、戈利木单抗），调用 IEDB Analysis Resource 公共 API（NetMHCpan 用于 MHC-I、NetMHCIIpan 用于 MHC-II [9, 16]）进行基准对比。对每条重链与轻链（共 18 条链、90 个链-等位基因组合、9,389 个评分窗口；MHC-I 等位基因 HLA-A\*02:01、A\*01:01、B\*07:02；MHC-II 等位基因 HLA-DRB1\*01:01、\*04:01），比较：(i) 逐肽排序一致性（ADA Predictor 内部分值与 IEDB 反向百分位秩之间的 Spearman ρ）；(ii) 各方法 top 10% 肽段集合的重叠度（Jaccard 指数）；(iii) 各自阈值下的强结合肽数目（ADA Predictor：蛋白内 top 10% 加引擎绝对阈值；NetMHCpan/NetMHCIIpan：百分位秩 ≤2）。完全可复现的脚本与原始 API 响应存档于 `paper/benchmark_iedb.py` 与 `paper/_iedb_cache.json`。

---

## 3. 结果

### 3.1 嵌合抗体案例：英夫利西单抗

英夫利西单抗（Remicade®/类克）为嵌合抗 TNFα IgG1 mAb，可变区源自小鼠、恒定区为人源。输入其 VH 与 VL 序列后系统输出：

| 指标 | 数值 |
|---|---|
| 序列内在风险 | 58.2%（高） |
| 临床 ADA 参照 | 34.23% |
| 临床 NAb 参照 | 15.0% |
| 预测 NAb 风险 | 23.2%（低） |
| 推断谱系 | 嵌合 |
| 评分路径 | 抗体 CDR/种系校准 |
| 序列异源度 | 49.5% |
| 基线分值 | 46.2% |
| 多信号上浮 | +12.0%（封顶 +15%） |
| 校准前全长表位负荷 | 68.9% |

高风险主要由 CDR 非人源种系表位负荷偏高（异源度 49.5%）与框架远离低 ADA 人源库驱动。5 项并发标志触发了多信号上浮（+12.0%）：CDR 非人源种系表位负荷偏高、种系 V 基因偏离度高、MHC-II 肽段覆盖种系异源 CDR 残基、B 细胞表位覆盖种系异源 CDR 残基、框架远离低 ADA 人源库。

### 3.2 全人源抗体案例：阿达木单抗

阿达木单抗（Humira®/修美乐）为全人源抗 TNFα mAb。系统预测：

| 指标 | 数值 |
|---|---|
| 序列内在风险 | 14.4%（低） |
| 临床 ADA 参照 | 36.67% |
| 临床 NAb 参照 | 5.0% |
| 预测 NAb 风险 | 17.0%（低） |
| 推断谱系 | 全人源 |
| 谱系耐受系数 | 0.55 |
| 序列异源度 | 20.7% |
| 多信号上浮 | 0.0%（未触发上浮） |
| 校准前全长表位负荷 | 62.8% |

尽管校准前表位负荷为 62.8%——与英夫利西单抗的 68.9% 相当——全人源谱系耐受系数（0.55）与低异源度（20.7%）将最终风险校准至 14.4%。值得注意的是，阿达木单抗在 IDC-DS-V1 数据集中的临床 ADA 为 36.67%（159 个队列的患者加权均值），*高于*其序列内在风险，这一偏离可归因于非序列因素（皮下给药、长期高剂量、聚集倾向及检测方法）而非分子免疫原性。该案例体现了核心设计原则：序列内在风险与临床 ADA 发生率度量的是不同维度，不应混为一谈。

### 3.3 参照库覆盖度

289 条参照库覆盖以下药物类别：

| 类别 | 数量 |
|---|---|
| 全人源 mAb | 26 |
| 人源化 mAb | 68 |
| 嵌合 mAb | 19 |
| 其他抗体/片段 | 3 |
| 酶（替代疗法/溶栓） | 37 |
| 激素/多肽 | 58 |
| 细胞因子/生长因子 | 30 |
| 凝血因子 | 22 |
| 抗体偶联药物（ADC） | 7 |
| Fc 融合蛋白 | 6 |
| 其他（毒素、血浆蛋白） | 13 |

*谱系（全人源/人源化/嵌合）由框架区异源性比例与人类胚系 V 基因池比对推断；非抗体类别按 `drug_class` 字段分组。*

### 3.4 与 NetMHCpan / NetMHCIIpan 的基准对比

在 9 个 mAb 的对比组中，ADA Predictor 内置的轻量级口袋/PSSM 启发式方法与神经网络预测器表现出中等程度的逐肽排序一致性：在 90 个链-等位基因组合（9,389 条肽段）中，MHC-I 的 Spearman ρ 平均为 0.485（中位数 0.473，标准差 0.108），MHC-II 平均为 0.385（中位数 0.399，标准差 0.135）；90 组比较中 89 组呈正相关，96%（MHC-I）/ 75%（MHC-II）超过 ρ = 0.3（表 2）。Top 10% 肽段集合重叠的 Jaccard 指数：MHC-I 为 0.334，MHC-II 为 0.250。一致性最高的是 HLA-A\*02:01（ρ = 0.602），即表征最完善的等位基因。

**表 2.** ADA Predictor 内部启发式与 IEDB 预测器在 9 个上市 mAb（18 条 VH/VL 链）上的逐肽一致性。

| 维度 | 比较组数 | 肽段总数 | ρ 均值（中位数 ± 标准差） | Top 10% Jaccard | ρ > 0.3 占比 |
|---|---|---|---|---|---|
| MHC-I vs NetMHCpan | 54 | 约 5,600 | 0.485（0.473 ± 0.108） | 0.334 | 96% |
| MHC-II vs NetMHCIIpan | 36 | 约 3,800 | 0.385（0.399 ± 0.135） | 0.250 | 75% |

各等位基因平均 ρ：HLA-A\*02:01 为 0.602，HLA-A\*01:01 为 0.447，HLA-B\*07:02 为 0.405；HLA-DRB1\*01:01 为 0.424，DRB1\*04:01 为 0.347。

ADA Predictor 标记的强结合肽数目约为 IEDB 百分位秩 ≤2 阈值的 3.2 倍（MHC-I 为 593 对 187；MHC-II 为 380 对 119），表明其敏感度更高而特异度更低——这对早期筛选启发式是合理的取向，但不能替代 NetMHCpan 用于等位基因特异的去免疫化决策。

在复合评分层面，九药对比组仅作示例而非验证集（n = 9）。ADA Predictor 序列内在风险保持了预期的分子排序——嵌合 29.1% > 人源化 15.6% > 全人源 11.2%（组均值）——而临床 ADA 发生率则不然（14.1%、2.1%、16.3%；序列风险与临床 ADA 之间的 Spearman ρ = −0.13）。两处最大的分歧在机制上均可解释：利妥昔单抗（嵌合，但为 CD20 B 细胞耗竭抗体；序列分被基于作用机制的情境耐受压低 12 倍，而临床 ADA 保持约 5%）与阿达木单抗（全人源序列，但该数据集中临床 ADA 达 36.67%，反映给药、聚集与序列外检测因素）。这一分歧正是 §2.9 所述设计原则的预期结果，并与已发表观察一致——即使生产级复合评分也只能解释临床 ADA 方差的一部分：一篇 2026 年行业综述报告诺和诺德 iPREDICT 的免疫原性风险分在 217 个 mAb 上与临床 ADA 相关系数 r = 0.35、反向人源度 r = 0.50（系 2024 年会议报告的二次转述 [19]），BioPhi OASis 人源度在其基准中仅解释临床 ADA 约 28% 的方差 [20]。

---

## 4. 讨论

### 4.1 与现有工具的比较

免疫原性预测工具的格局分散在三个层面（表 3）。**MHC/T 细胞工具**最为成熟：NetMHCpan/NetMHCIIpan [9, 16] 与 IEDB [10] 提供精确的逐等位基因结合与配体洗脱预测；PickPocket [25] 通过受体口袋相似性插值，扩展了对表征不足等位基因的覆盖，可作为神经网络预测器的补充。近期深度学习系统如 T-SCAPE [18] 甚至能在无显式 MHC 输入的情况下预测 T 细胞活化与 ADA 诱导潜力。**B 细胞工具**仍局限于单一维度，有线性的 BepiPred，也有基于构象结构的 DiscoTope-3.0 [17]；即使是后者也仅达到中等 F1（约 0.2–0.3）。**抗体层面工具**评估人源度或对整个分子分类：BioPhi OASis/Sapiens [20] 提供针对天然 OAS 库的可解释 9 肽人源度评估（人源度分类 ROC-AUC 0.966；对临床 ADA 方差解释约 28%）；Hu-mAb 与 TAP [24] 建立在 ANARCI 编号框架 [26] 之上，面向人源化与可开发性；AbImmPred [21] 基于 AntiBERTy 嵌入对抗体免疫原性分类（在 22 个 mAb 的独立测试集上准确率 72.7%）。制药企业内部管线（诺和诺德 iPREDICT、诺华 iSHAPe、EpiVax EpiMatrix）将 MHC 递呈与自体相似度及人群加权整合，但大多为专有系统 [19]。

**表 3.** 代表性免疫原性相关工具的功能覆盖。

| 工具 | MHC-I | MHC-II | B 细胞 | CDR/种系人源度 | 聚集 | NAb | 临床参照 | 集成 0–100 评分 | 临床 ADA 入评分 | 获取方式 |
|---|---|---|---|---|---|---|---|---|---|---|
| NetMHCpan 4.1 [9] | ✓ | – | – | – | – | – | – | – | – | 免费（IEDB） |
| NetMHCIIpan 4.x [16] | – | ✓ | – | – | – | – | – | – | – | 免费（IEDB） |
| IEDB 工具箱 [10] | ✓ | ✓ | ✓ | – | – | – | – | – | – | 免费 |
| BepiPred / DiscoTope | – | – | ✓ | – | – | – | – | – | – | 免费 Web |
| DiscoTope-3.0 [17] | – | – | ✓（3D） | – | – | – | – | – | – | 免费服务器 |
| BioPhi（OASis）[20] | – | – | – | ✓ | – | – | 部分 | – | – | 开源 |
| TAP / Hu-mAb（OPIG）[24] | – | – | – | ✓ | ✓ | – | – | – | – | 免费 Web |
| AbImmPred [21] | – | – | – | ✓ | – | – | 有（标签） | ✓（分类） | 是 | Web |
| iPREDICT / iSHAPe [19] | ✓ | ✓ | – | ✓ | – | – | 有 | ✓ | 是 | 专有 |
| T-SCAPE [18] | ✓ | ✓ | – | – | – | – | 有 | ✓ | 是 | 科研 |
| **ADA Predictor（本工作）** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | **289 药、解耦** | ✓ | **否（解耦）** | 开放 Web |

ADA Predictor 的独特定位不在于逐等位基因精度——§3.4 的基准测试证实神经预测器在该任务上仍占优势——而在于*覆盖度与可解释性*：单一、轻依赖、完全透明的平台，同时覆盖 T 细胞、B 细胞、抗体种系、聚集与独立的 NAb 维度，每个复合评分组件都以加权条形图直观呈现，每个阈值均有文档记录以保证可复现性。关键在于，ADA Predictor 是唯一在提供大规模临床 ADA/NAb 参照库的同时*从结构上排除*临床数据进入序列评分的平台，从而避免了小样本研究容易陷入的"序列特征预测临床结局"这一过度乐观且混杂的叙事框架。

### 4.2 序列风险与临床 ADA 的区分

核心设计原则是严格分离序列内在风险与临床 ADA 发生率。临床 ADA 率受非序列因素（给药途径、频次、合并免疫抑制、检测方法灵敏度、人群 HLA 分布）影响，无法仅凭分子序列确定。通过仅将临床率作背景对照展示——绝不回灌评分——ADA Predictor 避免了以片面分子信息预测临床结局的方法学错误。

需要辩证看待的是，已有工作确实证明序列特征对临床 ADA 携带*部分*预测信号：BioPhi OASis 人源度在其基准上解释了约 28% 的临床 ADA 方差 [20]，AbImmPred 基于 AntiBERTy 嵌入的分类在 22 个 mAb 上达 72.7% 准确率 [21]，一篇 2026 年行业综述报告 iPREDICT 的免疫原性风险分在 217 个 mAb 上与临床 ADA 相关系数 r = 0.35（系二次转述 [19]）。这些信号真实且有用，但 28% 的解释方差意味着 72% 归于非序列因素，n = 22 的分类结论尚属初步。因此 ADA Predictor 的立场并非"序列与临床 ADA 无关"，而是"序列是必要而非充分的信息"——将其与临床发生率分开展示，恰恰保留了两者各自的解释力，而非把弱相关性强行拟合为确定性预测。NAb 评分中的 ADA 发生率因子是唯一例外，其中最邻近药物的临床 ADA 率作为*邻近度加权先验*而非直接预测因子。

### 4.3 局限性

1. **MHC 结合启发式：** 口袋偏好矩阵对特定等位基因-肽段对的精度不及 NetMHCpan 神经网络模型。
2. **线性 B 细胞表位：** >90% 的 B 细胞表位为构象表位；线性评分的表面区分能力有限。
3. **CDR 划分：** 锚点启发式近似但不替代 IMGT/ANARCI 正式编号。
4. **无三维结构：** 格式暴露由序列架构推断，非晶体结构对接。
5. **临床数据异质性：** 参照库中的 ADA/NAb 率受不同检测方法与试验设计影响；跨药比较需谨慎。

### 4.4 未来方向

- 针对 §4.3 局限 1 与 §3.4 基准（MHC-I ρ = 0.485、特异度低于 NetMHCpan），集成 NetMHCpan-4.1 作为可选 MHC 后端，以在等位基因特异去免疫化场景下提供更高精度，同时保留内置启发式作为轻量默认。
- 针对 §4.3 局限 3，采用 ANARCI 替代锚点启发式进行正式 IMGT 编号，以满足申报文档的可追溯要求。
- 针对 §4.3 局限 2（>90% 的 B 细胞表位为构象表位），基于结构预测模型引入构象 B 细胞表位预测，补足线性评分的表面区分短板。
- 针对 §4.3 局限 5 与新模态涌现，扩充参照库以纳入生物类似药与双抗、CAR-T 构建体，既增大临床背景样本量，也为新模态提供解耦式参照。

---

## 致谢

感谢 ADA Predictor 开发团队在平台设计、实现与验证过程中的协作努力。感谢以下开源工具与公共资源使本工作得以实现：Biopython [8] 用于序列比对，NumPy 与 pandas 用于数值计算，FastAPI 用于 Web 服务，以及 IEDB [10]、IMGT [6]、PDB [7] 数据库提供的免疫学与结构参照数据。MHC 结合预测依赖 NetMHCpan/NetMHCIIpan [9, 16] 及可选的 mhcflurry 后端。289 条药物临床参照表整合自 IDC-DS-V1 [12]、FDA/EMA 标签及 Thera-SAbDab，以 CC BY 4.0 协议发布。

## 资助

无资助。

## 利益冲突

作者声明无利益冲突。

## 数据与代码可用性

*immunopred* 源代码见 https://github.com/ADAPredictor/Anti-Drug-Antibody-Predictor。本研究使用的 289 条药物临床参照表整合了 IDC-DS-V1 数据集（Immunogenicity Database Collaborative，CC BY 4.0 [12]，https://immunogenicitydb.org/）、FDA/EMA 药品标签（公共领域监管数据）、Thera-SAbDab 抗体序列及 UniProt/PDB 记录的临床 ADA/NAb 率。该参照表连同案例研究使用的基准复合数值（`paper/benchmark_composite.csv`）以 CC BY 4.0 协议发布于同一 GitHub 仓库。IEDB 百分位秩缓存源自公开 IEDB 数据。

---

## 参考文献

1. Kolaskar AS, Tongaonkar PC. A semi-empirical method for prediction of antigenic determinants on protein antigens. *FEBS Lett.* 1990;276(1-2):172-176.
2. Parker JM, Guo D, Hodges RS. New hydrophilicity scale derived from high-performance liquid chromatography peptide retention data. *Biochemistry.* 1986;25(19):5425-5432.
3. Emini EA, Hughes JV, Perlow DS, Boger J. Induction of hepatitis A virus-neutralizing antibody by a virus-specific synthetic peptide. *J Virol.* 1985;55(3):836-839.
4. Chou PY, Fasman GD. Prediction of the secondary structure of proteins from their amino acid sequence. *Adv Enzymol.* 1978;47:45-148.
5. Chothia C, Lesk AM. Canonical structures for the hypervariable regions of immunoglobulins. *J Mol Biol.* 1987;196(4):901-917.
6. Lefranc MP, Giudicelli V, Ginestoux C, et al. IMGT, the international ImMunoGeneTics information system. *Nucleic Acids Res.* 2009;37(suppl_1):D1006-D1012.
7. Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids Res.* 2000;28(1):235-242.
8. Cock PJ, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics.* 2009;25(11):1422-1423.
9. Nielsen M, Lundegaard C, Blicher T, et al. NetMHCpan, a method for quantitative prediction of peptide binding to any HLA-A and -B locus. *PLoS One.* 2008;3(3):e1797.
10. Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB): 2018 update. *Nucleic Acids Res.* 2019;47(D1):D339-D343.
11. World Health Organization. International Nonproprietary Names (INN) for biological and biotechnological substances. *INN Working Document 05.179.* 2005.
12. Immunogenicity Database Collaborative. IDC DS V1: Clinical immunogenicity dataset of marketed biologics. *CC BY 4.0.* https://immunogenicitydb.org/
13. Mayrhofer P, Kunert R. Nomenclature of humanized mAbs: early concepts, current challenges and future perspectives. *Hum Antibodies.* 2018;27(1):37-51.
14. Jawa V, Cousens LP, Awwad M, Wakshull E, Kropshofer H, De Groot AS. T-cell dependent immunogenicity of protein therapeutics: preclinical assessment and mitigation. *Clin Immunol.* 2013;149(3):534-555.
15. Rosenberg AS, Sauna ZE. Immunogenicity assessment during the development of protein therapeutics. *J Pharm Pharmacol.* 2018;70(5):584-594.
16. Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions across MHC alleles using refined data and analysis. *Nucleic Acids Res.* 2020;48(W1):W449-W454.
17. Sanchez-Trincado JL, Gomez-Perosanz M, Reche PA. Fundamentals and methods for T- and B-cell epitope prediction. *J Immunol Res.* 2017;2017:2680160.
18. Kim J, et al. T-SCAPE: a novel T-cell activation prediction method identifying immunogenic regions in protein therapeutics. *Sci Adv.* 2025;11(49):eadz8759.
19. Hickling TP, Nielsen M, Meysman P, Rose RH, Obrezanova O. Review: application and opportunities for machine learning and artificial intelligence in preclinical immunogenicity risk assessment. *Front Immunol.* 2026;17:1720928.
20. Prihoda D, Maamary J, Waight A, Juan V, Fayadat-Dilman L, Svozil D, Bitton DA. BioPhi: a platform for antibody design, humanization, and humanness evaluation based on natural antibody repertoires and deep learning. *MAbs.* 2022;14(1):2020203.
21. Wang H, Hao X, He Y, Fan L. AbImmPred: an immunogenicity prediction method for therapeutic antibodies using AntiBERTy-based sequence features. *PLoS One.* 2024;19(2):e0296737.
22. Tung CW, Ho SY. POPI: predicting immunogenicity of MHC class I binding peptides by mining informative physicochemical properties. *Bioinformatics.* 2007;23(8):942-949.
23. Collatz M, Mock F, Barth E, Hölzer M, Sachse K, Marz M. EpiDope: a deep neural network for linear B-cell epitope prediction. *Bioinformatics.* 2021;37(4):448-455.
24. Marks C, Hummer AM, Chin M, Deane CM. Humanization of antibodies using a machine learning approach on large-scale repertoire data. *Bioinformatics.* 2021;37(22):4041-4047.
25. Zhang H, Lund O, Nielsen M. The PickPocket method for predicting binding specificities for receptors based on receptor pocket similarities: application to MHC-peptide binding. *Bioinformatics.* 2009;25(10):1293-1299.
26. Dunbar J, Deane CM. ANARCI: antigen receptor numbering and receptor classification. *Bioinformatics.* 2016;32(2):298-300.

---

*稿件整理日期：2026-09-22。拟投稿至计算生物学/生物信息学领域同行评审期刊。*
