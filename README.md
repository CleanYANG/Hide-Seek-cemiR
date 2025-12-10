下面给你一份**完整、美观、专业级别的 README.md**，格式符合生物信息学工具的常用标准（类似 RNAhybrid、IntaRNA、TargetScan 的风格），也适用于你的 GitHub 项目首页。

---

Hide-Seek-cemiR

### *Structure-aware miRNA–RNA interaction prediction with seed rules, accessibility, ΔG, and shared-target analysis*

**Hide-Seek-cemiR** 是一个轻量但功能完整的 miRNA 靶点预测工具，它综合：

* **三类 Canonical Seed**（8mer, 7mer-m8, 7mer-A1）
* **RNAfold MFE 二级结构预测**（可及性：unpaired/partial/paired）
* **RNAduplex 计算局部结合自由能 ΔG**
* **多 RNA 聚合分析（per-miRNA-per-RNA best site）**
* **寻找共享靶点 miRNAs（shared-miRNA detection）**
* **PubMed 自动文献注释（可选）**

该工具最初为研究小鼠睾丸特异性 lncRNA **Stgart** 与 **Star 3′UTR** 之间的竞争性内源 RNA（ceRNA）关系而开发，但其方法完全通用于任何物种、任何 miRNA–RNA 体系。

---

# 🚀 **Features**

### 🔹 **1. Canonical miRNA seed scanning**

支持三种强结合 seed 类型：

| Seed Type   | Rule               | 描述       |
| ----------- | ------------------ | -------- |
| **8mer**    | A + rc(miRNA[2–8]) | 最强结合     |
| **7mer-m8** | rc(miRNA[2–8])     | 2–8 完全配对 |
| **7mer-A1** | A + rc(miRNA[2–7]) | 5′A 锚定   |

---

### 🔹 **2. RNAfold secondary-structure accessibility**

工具会对每条 RNA：

* 调用 **RNAfold** 预测结构
* 抽取所有连续 unpaired 区间
* 生成 `unpaired.region.csv`
* 构建 binary mask：

  * 1 = unpaired
  * 0 = paired

每个 seed 位点会分类成：

| access_cat | 解释                   |
| ---------- | -------------------- |
| **1**      | 完全暴露（fully open）     |
| **2**      | 部分暴露（partially open） |
| **3**      | 完全埋藏（paired region）  |

---

### 🔹 **3. RNAduplex ΔG calculation**

对每个 seed 命中位点：

* 截取一个可配置的局部窗口（默认 -10/+15 nt）
* 用 **RNAduplex** 计算 miRNA–RNA 双链形成自由能 ΔG
* ΔG 越负 → 结合越稳定

---

### 🔹 **4. Best-site aggregation (per-miRNA-per-RNA)**

对每个 (miRNA, RNA)：

* 根据 seed 强度、可及性、ΔG 排序
* 选择一个最强位点作为代表
* 输出 `per_miRNA_per_RNA_best_site.xlsx`

---

### 🔹 **5. Shared miRNA detection**

寻找同时强烈靶向 RNA-A 和 RNA-B 的 miRNAs。

例如：
✔ **Stgart**
✔ **3′UTR-Star**

过滤标准包括：

* access_cat ≤ 2
* ΔG ≤ -10
* seed 类型优先 8mer/7mer-m8
* 按 worst-accessibility + ΔG 总和排序

输出：

`shared_miRNAs_between_A_B.xlsx`

---

### 🔹 **6. PubMed annotation (可选)**

自动为 top shared miRNAs：

* 查询 PubMed
* 获取摘要、标题、PMID
* 输出：

`shared_miRNAs_with_pubmed.xlsx`

用于快速收集生物学证据。

---

# 🧪 **Installation**

以下操作基于 mamba / conda 环境：

```bash
git clone https://github.com/CleanYANG/Hide-Seek-cemiR.git
cd Hide-Seek-cemiR

mamba env create -f environment.yml
mamba activate hide-seek-cemir
```

environment.yml 安装：

* Python
* Biopython
* pandas
* requests
* **ViennaRNA**（RNAfold + RNAduplex）

---

# ⚙️ **Usage**

目前 pipeline 使用脚本内置参数（下一版将提供 CLI 参数接口）。

运行：

```bash
python hide_seek_cemir.py
```

默认输出包括：

```
*.unpaired.region.csv
*.all_sites.xlsx
per_miRNA_per_RNA_best_site.xlsx
shared_miRNAs_between_A_B.xlsx
shared_miRNAs_with_pubmed.xlsx
```

---

# 📂 **Input Format**

### MiRNA FASTA (e.g. mmu_mature.fa)

```
>mmu-miR-150-5p
UCUCCCAACCCUUGUACCAGUG
>mmu-let-7b-5p
UGAGGUAGUAGGUUGUGUGGUU
```

### RNA FASTA or plain TXT

#### FASTA (recommended)

```
>Stgart
ATGCGTAAAAAA...
```

#### TXT

文件内容仅为序列本体。

---

# 📤 **Output Files Overview**

| 文件                                 | 含义                          |
| ---------------------------------- | --------------------------- |
| `*.unpaired.region.csv`            | 每条 RNA 的所有 unpaired 区间      |
| `*.all_sites.xlsx`                 | 每个 seed 位点的 access_cat + ΔG |
| `per_miRNA_per_RNA_best_site.xlsx` | 每个 (miRNA, RNA) 最佳位点        |
| `shared_miRNAs_between_A_B.xlsx`   | 同时靶向 RNA-A & RNA-B 的 miRNA  |
| `shared_miRNAs_with_pubmed.xlsx`   | 上述结果 + PubMed 注释            |

---

# 📈 **Workflow Diagram**

```
miRNA FASTA           RNA sequences
     │                     │
     ▼                     ▼
 Seed motif generation   RNAfold structure prediction
     │                     │
     └─────► Accessibility mask ◄──────┘
                 ▼
          Seed site scanning
                 ▼
        Local ΔG (RNAduplex)
                 ▼
    Per-miRNA-per-RNA aggregation
                 ▼
       Shared miRNA detection
                 ▼
         PubMed annotation
                 ▼
              Results


---

*Citation

If you use *Hide-Seek-cemiR* in your research, please cite:

> Hong Yang, Hide-Seek-cemiR: Structure-aware prediction of miRNA–RNA interactions.
> GitHub repository: [https://github.com/CleanYANG/Hide-Seek-cemiR](https://github.com/CleanYANG/Hide-Seek-cemiR)

---
Author
Hong Yang, Hokkaido University


