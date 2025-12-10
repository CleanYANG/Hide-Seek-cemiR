#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Hide-Seek-cemiR — multi-RNA miRNA accessibility & energy scanner
---------------------------------------------------------------
当前版本功能：
  1. 读取成熟 miRNA（FASTA）
  2. 读取多条 RNA 序列（FASTA 或 txt）
  3. 使用 RNAfold 预测每条 RNA 的 MFE 结构
  4. 从 dot-bracket 中抽取所有连续 unpaired 区间
  5. 生成 unpaired.region.csv
  6. 为每条 RNA 构建 binary accessibility mask (1=unpaired, 0=paired)
  7. 对每条 miRNA，定义三种强 seed 类型：
        - 8mer
        - 7mer-m8
        - 7mer-A1
  8. 在所有 RNA 上扫描 8mer/7mer-m8/7mer-A1 motif
  9. 结合结构 mask 判断 seed 是否位于开放区域 (access_cat 1/2/3)
 10. 使用 RNAduplex 计算 miRNA 与局部目标序列的 ΔG
 11. 输出所有位点到 all_miRNA_sites_with_energy.xlsx

后续可以在此基础上：
  - 聚合到 (miRNA, RNA) 层面打分排序
  - 找出同时强烈作用于 RNA-A 和 RNA-B 的 miRNA
  - 调用 PubMed 注释相关文献
"""

import os
import re
import io
import time
import subprocess
import pandas as pd

from Bio import SeqIO
from Bio import Medline
import requests


# =========================
# 0. 文件名 & 外部命令配置
# =========================

# 成熟 miRNA FASTA
MIRNA_FASTA = "mmu_mature.fa"

# RNA 输入文件（可包含多条）
RNA_INPUT_FILES = [
    "Stgart.txt",
    "3UTR-Star.txt",
    # 以后可以继续添加其它 RNA 文件
]

# 物种过滤前缀（例如 "mmu-", "hsa-", "rno-" 等）
SPECIES_PREFIX = "mmu-"

# RNAfold & RNAduplex 命令（需在 PATH 中）
RNAFOLD_CMD = "RNAfold"
RNADUPLEX_CMD = "RNAduplex"

# 输出文件
UNPAIRED_CSV = "unpaired.region.csv"
OUTPUT_SITES_EXCEL = "all_miRNA_sites_with_energy.xlsx"


# =========================
# 1. 读取 miRNA 序列
# =========================

print("Loading miRNA fasta...")
mir_records = list(SeqIO.parse(MIRNA_FASTA, "fasta"))

if SPECIES_PREFIX:
    mir_records = [r for r in mir_records if r.id.startswith(SPECIES_PREFIX)]

mirna_dict = {r.id: str(r.seq) for r in mir_records}
print(f"Total miRNAs loaded ({SPECIES_PREFIX}): {len(mirna_dict)}")


# =========================
# 2. 读取 RNA 序列（多文件、多条）
# =========================

def load_rna_sequences(files):
    """
    从给定的文件列表中读取 RNA 序列。
    支持：
      - FASTA: 可包含多条序列，id 来自 fasta header
      - TXT  : 单条序列，id 来自文件名（去掉扩展名）
    返回:
      rna_seqs: dict[rna_id] = DNA 字母序列（大写，U->T）
    """
    rna_seqs = {}
    for path in files:
        if not os.path.exists(path):
            raise FileNotFoundError(f"RNA input file not found: {path}")

        ext = os.path.splitext(path)[1].lower()
        if ext in (".fa", ".fasta", ".fna"):
            for rec in SeqIO.parse(path, "fasta"):
                sid = rec.id
                seq = str(rec.seq).strip().upper().replace("U", "T")
                if not seq:
                    continue
                if sid in rna_seqs:
                    print(f"[Warning] Duplicate RNA ID {sid}, overwriting.")
                rna_seqs[sid] = seq
        else:
            # 简单 txt：单条序列
            with open(path) as f:
                seq = f.read().strip().upper().replace("U", "T")
            if not seq:
                continue
            sid = os.path.splitext(os.path.basename(path))[0]
            if sid in rna_seqs:
                print(f"[Warning] Duplicate RNA ID {sid}, overwriting.")
            rna_seqs[sid] = seq

    return rna_seqs


print("Loading RNA sequences...")
rna_seqs = load_rna_sequences(RNA_INPUT_FILES)
print(f"Total RNA sequences loaded: {len(rna_seqs)}")
for rid, seq in rna_seqs.items():
    print(f"  - {rid}: length {len(seq)}")


# =========================
# 3. RNAfold：预测 MFE 结构
# =========================

def run_rnafold(seq, rna_id=None, cmd=RNAFOLD_CMD):
    """
    使用 RNAfold 预测单条 RNA 的 MFE 结构。
    输入:
      seq   : DNA 风格序列（A/T/G/C），内部会自动把 T→U
      rna_id: 仅用于报错/打印
    返回:
      dotbracket: 结构字符串，如 "..((...))..."
      energy    : 浮点数，自由能（kcal/mol），若解析失败则为 None
    """
    rna_seq = seq.replace("T", "U")

    try:
        p = subprocess.run(
            [cmd, "--noPS"],
            input=(rna_seq + "\n").encode("utf-8"),
            capture_output=True,
            check=True,
        )
        out = p.stdout.decode("utf-8").strip().splitlines()
        if len(out) < 2:
            raise RuntimeError(f"Unexpected RNAfold output for {rna_id}: {out}")

        structure_line = out[1]
        # 例如: "..((...)).."  (-12.30)
        m = re.match(r"^([\.\(\)]+)\s+\(([-0-9\.]+)\)", structure_line)
        if not m:
            struct = structure_line.split()[0]
            energy = None
        else:
            struct = m.group(1)
            try:
                energy = float(m.group(2))
            except ValueError:
                energy = None

        return struct, energy

    except subprocess.CalledProcessError as e:
        print(f"[Error] RNAfold failed for {rna_id}: {e}")
        print("stdout:", e.stdout.decode("utf-8", errors="ignore"))
        print("stderr:", e.stderr.decode("utf-8", errors="ignore"))
        return None, None


# =========================
# 4. 从 dot-bracket 抽取 unpaired 区间 & 构建 mask
# =========================

def extract_unpaired_regions(dotbracket):
    """
    从 dot-bracket 字符串中抽取所有连续的 '.' 区间。
    返回:
      ranges: list[(start, end)]，1-based 闭区间
    """
    ranges = []
    in_region = False
    start = None

    for i, ch in enumerate(dotbracket, start=1):  # 1-based
        if ch == ".":
            if not in_region:
                in_region = True
                start = i
        else:
            if in_region:
                ranges.append((start, i - 1))
                in_region = False
                start = None

    if in_region and start is not None:
        ranges.append((start, len(dotbracket)))

    return ranges


def build_mask_from_dotbracket(dotbracket):
    """
    根据 dot-bracket 构建 binary accessibility mask。
    返回:
      mask: list[int]，长度 = len(dotbracket) + 1
            mask[i] = 1 表示第 i 个 nt 为 '.', 0 表示 '(' 或 ')'
            mask[0] 闲置不用
    """
    n = len(dotbracket)
    mask = [0] * (n + 1)
    for i, ch in enumerate(dotbracket, start=1):
        mask[i] = 1 if ch == "." else 0
    return mask


print("Running RNAfold & extracting unpaired regions...")

rna_structures = {}    # rna_id -> dotbracket
rna_energies = {}      # rna_id -> MFE energy
rna_masks = {}         # rna_id -> mask[0..len]

unpaired_records = []  # 用于写 unpaired.region.csv

for rna_id, seq in rna_seqs.items():
    struct, energy = run_rnafold(seq, rna_id=rna_id)
    if struct is None:
        print(f"[Warning] Skip {rna_id} because RNAfold failed.")
        continue

    if len(struct) != len(seq.replace("U", "T")):
        print(f"[Warning] Length mismatch for {rna_id}: seq={len(seq)}, struct={len(struct)}")

    rna_structures[rna_id] = struct
    rna_energies[rna_id] = energy

    mask = build_mask_from_dotbracket(struct)
    rna_masks[rna_id] = mask

    ranges = extract_unpaired_regions(struct)
    for idx, (s, e) in enumerate(ranges, start=1):
        unpaired_records.append({
            "rna_id": rna_id,
            "region_id": idx,
            "start": s,
            "end": e,
            "length": e - s + 1,
        })

    print(f"  - {rna_id}: length={len(seq)}, unpaired_bases={sum(mask)}, regions={len(ranges)}, MFE={energy}")


if unpaired_records:
    df_unp = pd.DataFrame(unpaired_records)
    df_unp.to_csv(UNPAIRED_CSV, index=False)
    print(f"Saved unpaired regions to: {UNPAIRED_CSV}")
else:
    print("[Warning] No unpaired regions found; CSV not written.")


# =========================
# =========================
# 5. 计算 miRNA–target duplex ΔG (RNAduplex) — 正确输入格式版
# =========================

def compute_duplex_energy(mirna_seq, target_dna_subseq, cmd=RNADUPLEX_CMD):
    """
    使用 RNAduplex 计算 miRNA 与目标局部 RNA 序列的结合能 ΔG。
    输入:
      mirna_seq        : miRNA 序列（A/U/G/C；若含 T 会自动转 U）
      target_dna_subseq: 目标 RNA 的局部 DNA 风格子序列（A/T/G/C）
    返回:
      deltaG (float)   : kcal/mol；解析失败时返回 None
    说明:
      RNAduplex 期望的是“两条 RNA 序列”，分行或 FASTA 格式，不是 "seq1&seq2" 一行。
    """
    # 转成 RNA
    mir_rna = mirna_seq.upper().replace("T", "U")
    tgt_rna = target_dna_subseq.upper().replace("T", "U")

    # ✅ 正确的输入格式：两条序列，FASTA 风格
    # 也可以不用 >header，只写两行序列，但 FASTA 更保险
    duplex_input = f">target\n{tgt_rna}\n>mirna\n{mir_rna}\n"

    try:
        p = subprocess.run(
            [cmd],
            input=duplex_input.encode("utf-8"),
            capture_output=True,
            check=True,
        )
        out_text = p.stdout.decode("utf-8", errors="ignore").strip()
        if not out_text:
            print("[Warning] RNAduplex output is empty. Input was:")
            print(duplex_input)
            return None

        lines = out_text.splitlines()
        last = lines[-1]

        # 如有需要，可以暂时打开调试输出看看格式
        # print("RNAduplex raw output:\n", out_text)

        # 1) 优先匹配括号里的数字: (...)，例如 (-21.90)
        m = re.search(r"\(([-+]?\d+\.\d+)\)", last)
        if m:
            return float(m.group(1))

        # 2) 如果没匹配上，尝试在最后一行中找一个带小数点的数字
        m2 = re.search(r"([-+]?\d+\.\d+)", last)
        if m2:
            return float(m2.group(1))

        # 3) 再向前看看前几行（某些版本可能把能量放在别的行）
        for line in reversed(lines[:-1]):
            m3 = re.search(r"\(([-+]?\d+\.\d+)\)", line)
            if m3:
                return float(m3.group(1))
            m4 = re.search(r"([-+]?\d+\.\d+)", line)
            if m4:
                return float(m4.group(1))

        print("[Warning] Cannot parse ΔG from RNAduplex output. Last line was:")
        print(last)
        return None

    except FileNotFoundError:
        print(f"[ERROR] RNAduplex command '{cmd}' not found. 请确认已安装 ViennaRNA 并在 PATH 中。")
        return None
    except subprocess.CalledProcessError as e:
        print(f"[Warning] RNAduplex failed (CalledProcessError).")
        print("stdout:", e.stdout.decode("utf-8", errors="ignore"))
        print("stderr:", e.stderr.decode("utf-8", errors="ignore"))
        return None
    except Exception as e:
        print(f"[Warning] RNAduplex failed with unexpected error: {e}")
        return None


# =========================
# 6. 定义 seed 类型：8mer / 7mer-m8 / 7mer-A1
# =========================

def get_seed_motifs(mirna_seq):
    """
    根据 TargetScan 规范，为一条 miRNA 生成三种强结合 seed motif：
      - 8mer     : A + rc(miRNA[2–8])
      - 7mer-m8  :     rc(miRNA[2–8])
      - 7mer-A1  : A + rc(miRNA[2–7])

    返回:
      motifs: dict[type] = motif_seq (DNA 风格: A/T/G/C)
              type ∈ {"8mer", "7mer-m8", "7mer-A1"}

    miRNA 序列默认视为 RNA（A/U/G/C），内部会自动 T->U。
    """
    s = mirna_seq.upper().replace("T", "U")
    if len(s) < 8:
        return {}

    seed_2_8 = s[1:8]  # 2–8
    seed_2_7 = s[1:7]  # 2–7

    comp_table = str.maketrans("AUGC", "UACG")
    rc_2_8 = seed_2_8.translate(comp_table)[::-1]
    rc_2_7 = seed_2_7.translate(comp_table)[::-1]

    rc_2_8_dna = rc_2_8.replace("U", "T")
    rc_2_7_dna = rc_2_7.replace("U", "T")

    motifs = {
        "7mer-m8": rc_2_8_dna,
        "8mer":   "A" + rc_2_8_dna,
        "7mer-A1": "A" + rc_2_7_dna,
    }
    return motifs


# =========================
# 7. accessibility 分类函数
# =========================

def classify_site(start, length, mask):
    """
    对一个 seed site 按可及性分类。
    输入:
      start : 1-based 起始位置
      length: motif 长度
      mask  : 对应 RNA 的 mask[0..n], 1=unpaired, 0=paired

    返回:
      1 = 全不配对（窗口内全是 1）
      2 = 部分不配对（至少一个 1，但不全是 1）
      3 = 全配对（全是 0 或越界视为 0）
    """
    positions = range(start, start + length)
    vals = []
    for i in positions:
        if 0 <= i < len(mask):
            vals.append(mask[i])
        else:
            vals.append(0)

    if all(v == 1 for v in vals):
        return 1
    elif any(v == 1 for v in vals):
        return 2
    else:
        return 3


# =========================
# 8. 多 RNA 通用扫描 + ΔG
# =========================

def scan_motifs_in_all_rnas(rna_seqs, rna_masks, mirna_dict,
                            flank_left=10, flank_right=15):
    """
    在所有 RNA 上扫描所有 miRNA 的 8mer / 7mer-m8 / 7mer-A1 位点，
    并计算：
      - seed 是否位于 unpaired/partial/paired 区
      - miRNA 与目标局部序列的 ΔG（RNAduplex）

    输入:
      rna_seqs   : dict[rna_id] = DNA 序列 (A/T/G/C)
      rna_masks  : dict[rna_id] = mask[0..n], 1=unpaired, 0=paired
      mirna_dict : dict[mir_id] = miRNA 序列 (mature)

    返回:
      sites: list of dict，每个 dict 是一条结合位点记录：
        {
          "miRNA"       : mir_id,
          "miRNA_seq"   : mir_seq,
          "rna_id"      : rna_id,
          "rna_len"     : len(seq),
          "seed_type"   : "8mer"/"7mer-m8"/"7mer-A1",
          "motif"       : motif_seq,
          "start"       : start_1based,
          "end"         : end_1based,
          "access_cat"  : 1/2/3,
          "window_start": win_start,
          "window_end"  : win_end,
          "deltaG"      : ΔG (float 或 None),
        }
    """
    all_sites = []

    for mir_id, mir_seq in mirna_dict.items():
        motifs = get_seed_motifs(mir_seq)
        if not motifs:
            continue

        for rna_id, seq in rna_seqs.items():
            mask = rna_masks.get(rna_id)
            if mask is None:
                continue

            for seed_type, motif in motifs.items():
                mlen = len(motif)
                if mlen == 0:
                    continue

                idx = seq.find(motif)
                while idx != -1:
                    start = idx + 1  # 1-based
                    end = start + mlen - 1

                    access_cat = classify_site(start, mlen, mask)

                    win_start = max(1, start - flank_left)
                    win_end = min(len(seq), end + flank_right)
                    sub_dna = seq[win_start - 1: win_end]  # 0-based slice

                    deltaG = compute_duplex_energy(mir_seq, sub_dna)

                    site = {
                        "miRNA": mir_id,
                        "miRNA_seq": mir_seq,
                        "rna_id": rna_id,
                        "rna_len": len(seq),
                        "seed_type": seed_type,
                        "motif": motif,
                        "start": start,
                        "end": end,
                        "access_cat": access_cat,
                        "window_start": win_start,
                        "window_end": win_end,
                        "deltaG": deltaG,
                    }
                    all_sites.append(site)

                    idx = seq.find(motif, idx + 1)

    return all_sites


print("Scanning all RNAs for 8mer / 7mer-m8 / 7mer-A1 sites and computing ΔG...")

sites = scan_motifs_in_all_rnas(rna_seqs, rna_masks, mirna_dict)
print(f"Total binding sites found: {len(sites)}")

df_sites = pd.DataFrame(sites)
print(df_sites.head())

df_sites.to_excel(OUTPUT_SITES_EXCEL, index=False)
print(f"Saved per-site results to: {OUTPUT_SITES_EXCEL}")

# =========================
# 9. 按 (miRNA, RNA) 聚合，挑选“最佳位点”
# =========================

if df_sites.empty:
    print("[Info] No binding sites found. Skip aggregation.")
else:
    # 为 seed 类型定义一个“强度等级”，方便排序：8mer 最强
    seed_rank_map = {"8mer": 0, "7mer-m8": 1, "7mer-A1": 2}
    df_sites["seed_rank"] = df_sites["seed_type"].map(seed_rank_map).fillna(3)

    # 有些 deltaG 可能为 None，用一个中性值进行排序占位（0 大概介于弱和强之间）
    df_sites["deltaG_for_sort"] = df_sites["deltaG"].fillna(0.0)

    # 按 miRNA, RNA, 可及性类别, seed 强度, ΔG 排序
    df_sites_sorted = df_sites.sort_values(
        ["miRNA", "rna_id", "access_cat", "seed_rank", "deltaG_for_sort"]
    )

    # 对每个 (miRNA, RNA) 组合，选排序后的第一条作为“最佳位点”
    df_best = df_sites_sorted.groupby(["miRNA", "rna_id"], as_index=False).first()

    # 统计每个 (miRNA, RNA) 有多少个 seed 命中
    sizes = df_sites.groupby(["miRNA", "rna_id"]).size().reset_index(name="n_sites")

    # 合并回 df_best
    df_best = df_best.merge(sizes, on=["miRNA", "rna_id"], how="left")

    # 精简一下列，保持关键信息
    cols_keep = [
        "miRNA", "miRNA_seq", "rna_id", "rna_len",
        "seed_type", "access_cat", "deltaG",
        "start", "end", "window_start", "window_end",
        "n_sites",
    ]
    df_best = df_best[cols_keep]

    df_best.to_excel("per_miRNA_per_RNA_best_site.xlsx", index=False)
    print("Saved per-miRNA-per-RNA best sites to per_miRNA_per_RNA_best_site.xlsx")

    # =========================
    # 10. 寻找“同时强烈 target RNA-A 和 RNA-B 的 miRNA”
    # =========================

    # 默认把前两条 RNA 视为 A 和 B（你也可以手动改成具体名字）
    rna_ids = list(rna_seqs.keys())
    if len(rna_ids) < 2:
        print("[Info] Less than 2 RNAs. Skip shared miRNA analysis.")
        df_shared_sorted = pd.DataFrame()
    else:
        RNA_A_ID = rna_ids[0]
        RNA_B_ID = rna_ids[1]
        print(f"[Info] Treating '{RNA_A_ID}' as RNA-A, '{RNA_B_ID}' as RNA-B for shared-miRNA search.")

        df_A = df_best[df_best["rna_id"] == RNA_A_ID].copy()
        df_B = df_best[df_best["rna_id"] == RNA_B_ID].copy()

        # 左右拼成 A/B 两列
        df_shared = df_A.merge(df_B, on="miRNA", suffixes=("_A", "_B"))

        if df_shared.empty:
            print("[Info] No miRNAs that bind both RNA-A and RNA-B (before filtering).")
            df_shared_sorted = pd.DataFrame()
        else:
            # 过滤条件：你可以根据实际调整阈值
            MAX_ACCESS_CAT = 2   # 1=全开放, 2=部分开放, 3=全配对
            MIN_DELTA_G   = -10  # ΔG <= -10 视为比较靠谱的结合

            df_shared_filt = df_shared[
                (df_shared["access_cat_A"] <= MAX_ACCESS_CAT) &
                (df_shared["access_cat_B"] <= MAX_ACCESS_CAT) &
                (df_shared["deltaG_A"] <= MIN_DELTA_G) &
                (df_shared["deltaG_B"] <= MIN_DELTA_G)
            ].copy()

            if df_shared_filt.empty:
                print("[Info] No shared miRNAs pass the accessibility/ΔG filters.")
                df_shared_sorted = pd.DataFrame()
            else:
                # 计算 worst_access（两边里较差的一侧）和 ΔG 总和（越负越好）
                df_shared_filt["worst_access"] = df_shared_filt[
                    ["access_cat_A", "access_cat_B"]
                ].max(axis=1)
                df_shared_filt["sum_deltaG"] = df_shared_filt["deltaG_A"] + df_shared_filt["deltaG_B"]

                # 排序：先按 worst_access（1→2→3），再按 sum_deltaG（越负越好）
                df_shared_sorted = df_shared_filt.sort_values(
                    ["worst_access", "sum_deltaG", "miRNA"]
                )

                df_shared_sorted.to_excel("shared_miRNAs_between_A_B.xlsx", index=False)
                print("Saved shared miRNAs between RNA-A and RNA-B to shared_miRNAs_between_A_B.xlsx")


# =========================
# 9. PubMed 注释工具（暂不自动调用）
# =========================

SEARCH_TERMS = [
    "steroidogenesis",
    "Leydig",
    "StAR",
    "Cyp11a1",
    "Hsd3b1",
    "testis",
    "spermatogenesis",
    "ovary",
    "gonad",
    "reproduction"
]

PUBMED_SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_FETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def pubmed_search_and_fetch(miRNA_name, max_papers=2):
    """
    对一条 miRNA：
      1) 用 miRNA 名 + 生殖/甾体关键词 在 PubMed 检索
      2) 拿到总 hits 数 + 前 max_papers 篇文献（标题+摘要+PMID）

    返回:
      total_hits: int，总命中数
      papers    : list[{"pmid", "title", "abstract"}]
    """
    query = f'{miRNA_name} AND (' + " OR ".join(SEARCH_TERMS) + ")"

    try:
        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
        }
        r = requests.get(PUBMED_SEARCH_URL, params=params, timeout=10)
        r.raise_for_status()
        js = r.json()
        total_hits = int(js["esearchresult"]["count"])
        idlist = js["esearchresult"].get("idlist", [])

        if not idlist:
            return total_hits, []

        pmids = idlist[:max_papers]

        params_fetch = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "text",
            "rettype": "medline",
        }
        r2 = requests.get(PUBMED_FETCH_URL, params=params_fetch, timeout=15)
        r2.raise_for_status()

        handle = io.StringIO(r2.text)
        records = Medline.parse(handle)

        papers = []
        for rec in records:
            pmid = rec.get("PMID", "")
            title = rec.get("TI", "")
            abstract = rec.get("AB", "")
            papers.append({
                "pmid": pmid,
                "title": title,
                "abstract": abstract
            })

        return total_hits, papers

    except Exception as e:
        print(f"[Warning] PubMed request failed for {miRNA_name}: {e}")
        return 0, []


def annotate_miRNAs(df, max_papers=2):
    """
    在给定的 DataFrame 上增加 4 列 PubMed 注释：
      - pubmed_hits   : 相关文献数量
      - top_titles    : 前 max_papers 篇文献的标题（用 ' || ' 连接）
      - top_pmids     : PMID（用 '; ' 连接）
      - top_abstracts : 多个摘要，用 '---' 拼接

    注意：这里 df 应该至少包含 'miRNA' 一列。
    """
    hits_list = []
    titles_list = []
    pmids_list = []
    abstracts_list = []

    for mir in df["miRNA"]:
        total_hits, papers = pubmed_search_and_fetch(mir, max_papers=max_papers)
        hits_list.append(total_hits)

        if papers:
            titles = " || ".join(p["title"] for p in papers if p["title"])
            pmids = "; ".join(p["pmid"] for p in papers if p["pmid"])
            abstracts = "\n---\n".join(p["abstract"] for p in papers if p["abstract"])
        else:
            titles = ""
            pmids = ""
            abstracts = ""

        titles_list.append(titles)
        pmids_list.append(pmids)
        abstracts_list.append(abstracts)

        time.sleep(0.34)

    df["pubmed_hits"]   = hits_list
    df["top_titles"]    = titles_list
    df["top_pmids"]     = pmids_list
    df["top_abstracts"] = abstracts_list

    return df

# =========================
# 11. 对 shared miRNAs 进行 PubMed 注释
# =========================

try:
    if "df_shared_sorted" in globals() and not df_shared_sorted.empty:
        # 只对排名前 N 个做文献注释，避免 PubMed 请求太多
        TOP_N = 50
        df_top = df_shared_sorted.head(TOP_N).copy()
        print(f"[Info] Annotating top {len(df_top)} shared miRNAs with PubMed info...")

        df_top_annot = annotate_miRNAs(df_top, max_papers=2)
        df_top_annot.to_excel("shared_miRNAs_with_pubmed.xlsx", index=False)
        print("Saved PubMed-annotated shared miRNAs to shared_miRNAs_with_pubmed.xlsx")
    else:
        print("[Info] No shared miRNAs to annotate with PubMed.")
except Exception as e:
    print(f"[Warning] PubMed annotation step failed: {e}")

