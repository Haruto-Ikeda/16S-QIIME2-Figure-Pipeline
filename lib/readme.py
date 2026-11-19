#!/usr/bin/env python3
"""readme.py — 出力の各フォルダに README.md (そのフォルダの中身と図の読み方) を書く。

- 書く場所: output/16S/・Qiime2/・figure/ver<N>/ と、その版の中にある図のフォルダ (01_read-depth/ 〜 12_ancombc-summary/)。
  図のフォルダは、あるものにだけ書く (README のためにフォルダは作らない。図のフォルダの有無で「描いた図」を判定しているため)。
- 文面は lib/readme_templates/<名前>.md。「この実行の設定」の表などの {{名前}} を、実行ごとの値に置き換える。
- 値は Step 10 と同じもの (config・群・06 の深さと除外) と、図のフォルダの CSV から取る (既存の図を飛ばした実行でも同じ値になるように)。
- 実行のたびに書き直す。config の pipeline.run.readme が false なら書かない (前の実行の README は消さない)。
- 図の中身は変えないので、図の版の判定 (VERSION_KEYS) には入れない。
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

TEMPLATES = Path(__file__).resolve().parent / "readme_templates"
NAME = "README.md"
RANK_JA = {"domain": "ドメイン", "kingdom": "界", "phylum": "門", "class": "綱", "order": "目", "family": "科",
           "genus": "属", "species": "種"}
METHOD_JA = {"tukey": "Tukey HSD (全ペア)", "welch": "Welch の t 検定 (2 群)"}


def _render(name: str, values: dict) -> str:
    """雛形の {{名前}} を値に置き換える。置き換え残しや、値に無い名前があれば止める (中身の欠けた README を書かない)。"""
    text = (TEMPLATES / f"{name}.md").read_text(encoding="utf-8")
    missing = sorted(set(re.findall(r"\{\{(\w+)\}\}", text)) - set(values))
    if missing:
        raise KeyError(f"README の雛形 {name}.md の値がありません: {missing}")
    return re.sub(r"\{\{(\w+)\}\}", lambda m: str(values[m.group(1)]), text)


def _pipeline_version() -> str:
    try:
        line = (Path(__file__).resolve().parent.parent / "VERSION").read_text(encoding="utf-8").splitlines()[0]
        m = re.search(r"\d+(?:\.\d+)+", line)
        return m.group(0) if m else line.strip()
    except (OSError, UnicodeDecodeError, IndexError):
        return "不明"


def _qiime_version() -> str:
    import importlib.metadata as md
    try:
        return md.version("qiime2")
    except md.PackageNotFoundError:
        return "不明"


def _ranks(ctx: dict) -> list[str]:
    return [str(r).strip().lower() for r in ctx["cfg"]["plotting"]["ranks"]]


def _top_n_text(v) -> str:
    return "all (有意なタクソンをすべて)" if isinstance(v, str) and v.strip().lower() == "all" else str(v)


def _taxon_bar_count(d: Path, ranks: list[str]) -> str:
    """#10 の枚数 (階級ごと)。csv/10_taxon-bar_<階級>.files.csv と stats.csv から数える。"""
    parts = []
    for r in ranks:
        f = d / "10_taxon-bar" / "csv" / f"10_taxon-bar_{r}.files.csv"
        if not f.exists():
            continue
        n = len(pd.read_csv(f))
        text = f"{n:,}"
        s = d / "10_taxon-bar" / "csv" / f"10_taxon-bar_{r}.stats.csv"
        if s.exists():
            st = pd.read_csv(s, dtype=str, keep_default_na=False)
            z = st.loc[st["reason"] == "全検体で 0", "taxon"].nunique() if "reason" in st else 0
            if z:
                text += f" (すべての検体で 0 だった {z:,} タクソンは描いていません)"
        parts.append(text if len(ranks) == 1 else f"{r}: {text}")
    return "・".join(parts) if parts else "不明"


def _ancom(d: Path, ranks: list[str], G) -> tuple[str, str]:
    """#11・#12 の「比べた群」と「計算の方法」。csv の all.csv (実際に計算した結果) から読む。"""
    for folder, stem in (("11_ancombc", "11_ancombc"), ("12_ancombc-summary", "12_ancombc-summary")):
        for r in ranks:
            f = d / folder / "csv" / f"{stem}_{r}.all.csv"
            if not f.exists():
                continue
            t = pd.read_csv(f, dtype=str, keep_default_na=False)
            if t.empty:
                continue
            got = set(t["group"])
            comps = [g for g in G.order if g in got]
            how = "対照群と各群の比較 (Dunnett 型)" if "Dunnett" in t["method"].iloc[0] else "対照群と 1 つの群の比較 (2 群)"
            return "・".join(comps), f"ANCOM-BC2 (R の ANCOMBC {t['ancombc_version'].iloc[0]})。{how}"
    return "不明", "ANCOM-BC2 (R の ANCOMBC)"


def _values(ctx: dict, G, d: Path) -> dict:
    from . import figures as FG
    from . import plotcfg as PC
    from .steps import runs_of
    cfg, out = ctx["cfg"], ctx["out"]
    ranks = _ranks(ctx)
    depth = int((out / "06_rarefaction" / "sampling_depth.txt").read_text().strip())
    floor = int(cfg["diversity"]["sampling_depth_floor"])
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["sample-id"])
    n_all = len(G.label)
    n_excl = int(G.label.index.isin(list(excl)).sum())
    n_used = n_all - n_excl
    left = G.label[~G.label.index.isin(list(excl))]
    n_groups = len([g for g in G.order if (left == g).any()])

    if G.method == "dunnett":
        method = f"Dunnett (対照群 {G.control} と各群)"
    else:
        method = METHOD_JA.get(G.method, str(G.method))
    heat = method + ("。マスの中の文字で示します" if G.method == "tukey" else "。対照群と比べて有意な群のマスに、星で示します")
    has_venn = any((d / "09_venn-upset").glob("09_venn_*.pdf"))
    comps, ancom_method = _ancom(d, ranks, G)
    plot = cfg["plotting"]
    if G.use_subgroup:
        unit = "群 × サブグループ (Group_name と Subgroup_name の組)"
        color = "サブグループ (群 × サブグループの組) ごと" if G.color_by_subgroup else "もとの群 (Group_name)"
    else:
        unit, color = "群 (Group_name)", "群"
    if G.no_color:
        color = "なし (群の色を付けない設定)"
    return {
        "ver": d.name,
        "pipeline_version": _pipeline_version(),
        "qiime_version": _qiime_version(),
        "run_mode": ("絞り込み (解析を終えた別の実験の結果から、一部の検体を選んだもの)" if not _is_none(cfg.get("parent"))
                     else "全工程 (生データから)"),
        "runs": "・".join(r["m"] for r in runs_of(cfg)),
        "classifier": Path(str(cfg["classifier"]["path"])).name,
        "stats_method": method,
        "heatmap_test": heat,
        "min_n": G.min_n,
        "depth": f"{depth:,}",
        "floor": f"{floor:,}",
        "n_samples": f"{n_all:,}",
        "n_excluded": f"{n_excl:,}",
        "samples_used": (f"{n_used:,} (全 {n_all:,} 検体のうち、深さの下限に届かなかった {n_excl:,} 検体を除く)" if n_excl
                         else f"{n_all:,} (全検体。深さの下限に届かなかった検体はありません)"),
        "samples_fig01": (f"{n_all:,} (うち、下限に届かず以降の解析から外れた {n_excl:,} 検体は白抜き)" if n_excl
                          else f"{n_all:,} (下限に届かなかった検体はありません)"),
        "samples_total": (f"{n_all:,} (うち、深さの下限に届かず 02 以降から外れた検体 {n_excl:,})" if n_excl
                          else f"{n_all:,} (深さの下限に届かなかった検体はありません)"),
        "ranks": "・".join(f"{r} ({RANK_JA[r]})" if r in RANK_JA else r for r in ranks),
        "top_n": plot["top_n_bacteria"],
        "ancombc_top_n": _top_n_text(plot["ancombc_top_n"]),
        "ancombc_summary_top_n": _top_n_text(plot["ancombc_summary_top_n"]),
        "group_unit": unit,
        "group_order": "・".join(G.order),
        "control": PC._control_written(cfg) or "なし",   # G.control は tukey では None (検定に使わないため)。#11・#12 の対照群は config の値
        "comps": comps,
        "ancombc_method": ancom_method,
        "pcoa_color": color,
        "pcoa_marker": "サブグループ" if G.show_subgroup else "すべて丸 (サブグループを記号で分けない設定)",
        "nperm": f"{FG.NPERM:,}",
        "perm_seed": FG.PERM_SEED,
        "venn": f"{'あり' if has_venn else 'なし'} ({n_groups} 群)",
        "taxon_bar_count": _taxon_bar_count(d, ranks),
    }


def _is_none(v) -> bool:
    return v is None or (isinstance(v, str) and v.strip().lower() == "none")


def write_all(ctx: dict, G, d: Path) -> list[Path]:
    """README.md を書き、書いたファイルの一覧を返す。d = 図の版のフォルダ (output/16S/figure/ver<N>)。"""
    from .figures import FIG_DIRS
    out = ctx["out"]                                   # output/16S/Qiime2
    values = _values(ctx, G, d)
    targets = [("_output-16S", out.parent), ("_Qiime2", out), ("_figure-ver", d)]
    targets += [(name, d / name) for name in FIG_DIRS.values() if (d / name).is_dir()]
    written = []
    for name, folder in targets:
        text = _render(name, values)                   # 先に全部を組み立ててから書く (途中で止まって一部だけ新しくならないように)
        written.append((folder / NAME, text))
    for path, text in written:
        path.write_text(text, encoding="utf-8")
    return [p for p, _ in written]
