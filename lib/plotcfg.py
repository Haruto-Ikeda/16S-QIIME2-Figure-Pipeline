#!/usr/bin/env python3
"""plotcfg.py — 作図の config (群構成・対照群・検定・図ごとのスイッチ)を読み、検証する。

設計の正本は ../README.md「作図パイプラインの設計」節。
群の扱いはメタボローム (Scripts/figure/metabolome/lib/load.py)と同じ規則にしてある:
  - use_subgroup=false → 群 = Group_name
  - use_subgroup=true  → 群 = "Group_name Subgroup_name" (半角スペースでつなぐ)。並び = group_order × subgroups の定義順
  - 色は常に Group_name で決める
**16S には群・検体の除外が無い** (260917 決定。除外するなら QIIME 2 解析からやり直す)ので、
データにある群が group_order に無ければ停止する。不整合は黙って直さず止める。
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

# 実装済みの図のスイッチ（pipeline.run のキー）。図を実装したらここに足す。
# 14 個のスイッチ（図の番号は #1〜#12。#4 と #7 は 2 種類ずつ）がそろっていなければ停止する。
IMPLEMENTED = ["read_retention", "rarefaction", "alpha", "pcoa", "pcoa_matrix", "permanova", "dendrogram", "stacked_group", "stacked_sample", "heatmap", "venn", "taxon_barplot", "ancombc", "ancombc_summary"]
ALL_FIGURES = ["read_retention", "rarefaction", "alpha", "pcoa", "pcoa_matrix", "permanova", "dendrogram",
               "stacked_group", "stacked_sample", "heatmap", "venn", "taxon_barplot", "ancombc", "ancombc_summary"]
METHODS = ("tukey", "dunnett", "welch")   # welch は群が 2 つのときだけ（261009 本人決定）
# Tol discrete rainbow（29 色）と、色数 n = 1〜23 のときに使う番号（公式 Fig. 25。tol-colors 2.2.0 から転記・公式の 29 色と一致を確認 260924）
TOL_DISCRETE_RAINBOW = ['#E8ECFB', '#D9CCE3', '#D1BBD7', '#CAACCB', '#BA8DB4', '#AE76A3', '#AA6F9E', '#994F88', '#882E72', '#1965B0', '#437DBF', '#5289C7', '#6195CF', '#7BAFDE', '#4EB265', '#90C987', '#CAE0AB', '#F7F056', '#F7CB45', '#F6C141', '#F4A736', '#F1932D', '#EE8026', '#E8601C', '#E65518', '#DC050C', '#A5170E', '#72190E', '#42150A']
TOL_DISCRETE_RAINBOW_INDEX = [[9], [9, 25], [9, 17, 25], [9, 14, 17, 25], [9, 13, 14, 17, 25], [9, 13, 14, 16, 17, 25], [8, 9, 13, 14, 16, 17, 25], [8, 9, 13, 14, 16, 17, 22, 25], [8, 9, 13, 14, 16, 17, 22, 25, 27], [8, 9, 13, 14, 16, 17, 20, 23, 25, 27], [8, 9, 11, 13, 14, 16, 17, 20, 23, 25, 27], [2, 5, 8, 9, 11, 13, 14, 16, 17, 20, 23, 25], [2, 5, 8, 9, 11, 13, 14, 15, 16, 17, 20, 23, 25], [2, 5, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25], [2, 5, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 27], [2, 4, 6, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 27], [2, 4, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 27], [2, 4, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 26, 27], [1, 3, 4, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 26, 27], [1, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 19, 21, 23, 25, 26, 27], [1, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 25, 26, 27], [1, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 25, 26, 27, 28], [0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 25, 26, 27, 28]]
# Tol smooth rainbow（公式 Fig. 23 の 34 色。つないで連続色にし、必要な数だけ等間隔に取り出す＝何色でも作れる）
TOL_SMOOTH_RAINBOW = ['#E8ECFB', '#DDD8EF', '#D1C1E1', '#C3A8D1', '#B58FC2', '#A778B4', '#9B62A7', '#8C4E99', '#6F4C9B', '#6059A9', '#5568B8', '#4E79C5', '#4D8AC6', '#4E96BC', '#549EB3', '#59A5A9', '#60AB9E', '#69B190', '#77B77D', '#8CBC68', '#A6BE54', '#BEBC48', '#D1B541', '#DDAA3C', '#E49C39', '#E78C35', '#E67932', '#E4632D', '#DF4828', '#DA2222', '#B8221E', '#95211B', '#721E17', '#521A13']
PALETTE_NAMES = list(("bright", "high_contrast", "vibrant", "muted", "medium_contrast", "dark", "light", "pale")) + ["discrete_rainbow", "smooth_rainbow"]
NO_COLOR = "#222222"   # group_subgroup_colors: none（群の色なし）のときの全群の色 = 前景・軸と同じ（261013 本人決定）


def palette_colors(name: str, n: int) -> list[str]:
    """配色名から n 色を返す。色数が足りなければ止める (260924 本人決定)。
    24 色以上が要るときは smooth_rainbow を勧める (Tol の配色で 24 色以上を出せるのは smooth_rainbow だけ)。"""
    import matplotlib.colors as mc
    import numpy as np
    if name in TOL_QUALITATIVE and n <= len(TOL_QUALITATIVE[name]):
        return TOL_QUALITATIVE[name][:n]
    if name == "discrete_rainbow" and 1 <= n <= len(TOL_DISCRETE_RAINBOW_INDEX):
        return [TOL_DISCRETE_RAINBOW[i] for i in TOL_DISCRETE_RAINBOW_INDEX[n - 1]]
    if name != "smooth_rainbow":
        size = len(TOL_QUALITATIVE[name]) if name in TOL_QUALITATIVE else len(TOL_DISCRETE_RAINBOW_INDEX)
        hint = ("24 色以上が要るので smooth_rainbow にしてください (隣どうしの色は近くなります)" if n >= 24
                else "色の多い配色 (例 discrete_rainbow は 23 色まで)にするか、個別に指定してください")
        raise Stop(f"配色 {name} は {size} 色しかありません (必要 {n} 色)。{hint}")
    cm = mc.LinearSegmentedColormap.from_list("tol_smooth_rainbow", TOL_SMOOTH_RAINBOW)
    # 白に近い端 10% を外す (261016 本人決定・E6a。1 色目 #E8ECFB が白い背景でほぼ見えなかった。点検 F120・F125)
    return [mc.to_hex(cm(v)).upper() for v in np.linspace(0.1, 1, n)]
TOL_QUALITATIVE = {   # Paul Tol のカテゴリ用配色。公式ページの「Colours in default order」をそのまま転記（灰・欠損用は除く）
    # 出典: https://sronpersonalpages.nl/~pault/（260924 取得）。tol-colors 2.2.0 とは dark・pale の色と high/medium-contrast の順が違ったので公式に合わせた
    "bright": ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"],
    "high_contrast": ["#004488", "#DDAA33", "#BB5566"],
    "vibrant": ["#EE7733", "#0077BB", "#33BBEE", "#EE3377", "#CC3311", "#009988"],
    "muted": ["#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", "#999933", "#AA4499"],
    "medium_contrast": ["#6699CC", "#004488", "#EECC66", "#994455", "#997700", "#EE99AA"],
    "dark": ["#4477BB", "#CC4466", "#228833", "#2255AA", "#993366", "#006644", "#117788", "#882288", "#BB4488", "#775500"],
    "light": ["#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00"],
    "pale": ["#AACCEE", "#CCEEFF", "#BBDDBB", "#EEEEBB", "#FFBBCC", "#EEBBDD"],
}


class Stop(Exception):
    """config やデータの不整合。メッセージを出して止める。"""


def load_samples(out: Path) -> pd.DataFrame:
    """01_metadata/sample-metadata.txt から検体表を返す (列: sample-id, newID, Group_name, Subgroup_name)。"""
    md = pd.read_csv(out / "01_metadata" / "sample-metadata.txt", sep="\t", dtype=str)
    need = ["#SampleID", "newID", "Group_name", "Subgroup_name"]
    miss = [c for c in need if c not in md.columns]
    if miss:
        raise Stop(f"01_metadata/sample-metadata.txt に列がありません: {miss}")
    md = md[need].rename(columns={"#SampleID": "sample-id"})
    for c in md.columns:
        md[c] = md[c].str.strip()
    return md


def check_config(cfg: dict) -> None:
    """config を読むだけで分かる検査 (01 の前に呼ぶ。261021 本人決定・案 B)。
    データと照らし合わせる検査 (群名・対照群・色の数) は、手順 10 の Groups で行う。"""
    run = switches(cfg)
    g = _check_grouping(cfg)
    method, _ = _check_stats(cfg)
    written = _control_written(cfg)
    if method in ("dunnett", "welch") and not written:
        raise Stop(f"stats.method: {method} のときは comparison.control_group に対照群を書いてください (none は不可)")
    if (run.get("ancombc") or run.get("ancombc_summary")) and not written:   # 黙って図を飛ばさない (261021 本人決定・C3)
        raise Stop("#11・#12 ANCOM-BC を描くときは comparison.control_group に対照群を書いてください"
                   " (描かないなら pipeline.run の ancombc・ancombc_summary を false にしてください)")
    check_plotting(cfg)


def _is_none(v) -> bool:
    return isinstance(v, str) and v.strip().lower() == "none"


def _control_written(cfg: dict) -> str | None:
    """対照群 (comparison.control_group)。必須で、使わないときは none (261021 本人決定)。
    値はそのときの群の名前で書く: use_subgroup: false なら "Control"、true なら "Control SubA" (control_group_composite は廃止)。"""
    comp = cfg.get("comparison")
    if not isinstance(comp, dict) or "control_group" not in comp:
        raise Stop("config の comparison.control_group がありません (対照群の名前。使わないときは none)")
    v = comp["control_group"]
    if _is_none(v):
        return None
    if not isinstance(v, str) or not v.strip():
        raise Stop(f"comparison.control_group は群の名前か none で書いてください (現在: {v!r})")
    return v.strip()


def _pos_int_or_all(p: dict, k: str) -> None:
    v = p[k]
    if isinstance(v, str) and v.strip().lower() == "all":
        return
    if isinstance(v, bool) or not isinstance(v, int) or v < 1:
        raise Stop(f"plotting.{k} は正の整数か all で書いてください (現在: {v!r})")


def check_plotting(cfg: dict) -> None:
    """plotting の値 (261021 本人決定・C4。すべて必須)。"""
    p = cfg.get("plotting")
    if not isinstance(p, dict):
        raise Stop("config に plotting がありません")
    for k in ("ranks", "top_n_bacteria", "ancombc_top_n", "ancombc_summary_top_n"):
        if k not in p:
            raise Stop(f"config の plotting.{k} がありません")
    ranks = p["ranks"]
    if not isinstance(ranks, list) or not ranks:
        raise Stop(f"plotting.ranks は階級名のリスト (1 つ以上) で書いてください (現在: {ranks!r})")
    have = [str(x).strip().lower() for x in ((cfg.get("taxonomy") or {}).get("collapse_ranks") or [])]
    bad = [r for r in ranks if str(r).strip().lower() not in have]
    if bad:
        raise Stop(f"plotting.ranks に collapse していない階級があります: {bad} (taxonomy.collapse_ranks: {have})")
    n = p["top_n_bacteria"]
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise Stop(f"plotting.top_n_bacteria は正の整数で書いてください (現在: {n!r})")
    _pos_int_or_all(p, "ancombc_top_n")
    _pos_int_or_all(p, "ancombc_summary_top_n")


def _check_grouping(cfg: dict) -> dict:
    g = cfg.get("grouping")
    if not isinstance(g, dict):
        raise Stop("config に grouping がありません (group_order・group_subgroup_colors・use_subgroup)")
    for k in ("use_subgroup", "group_order", "group_subgroup_colors", "show_group_subgroup"):
        if k not in g:
            raise Stop(f"config の grouping.{k} がありません")
    for k in ("use_subgroup", "show_group_subgroup"):        # "false" などの文字は true 扱いになるので型で確かめる (261021 本人決定・C1)
        if not isinstance(g[k], bool):
            raise Stop(f"grouping.{k} は true / false で書いてください (現在: {g[k]!r})")
    # group_order は群名のリストか false。同じ群を 2 回書かない (261021 本人決定・C2)
    go = g["group_order"]
    if go is not False:
        if not isinstance(go, list) or not go:
            raise Stop(f"grouping.group_order は群名のリストか false で書いてください (現在: {go!r})")
        names = [str(x) for x in go]
        dup_g = sorted({x for x in names if names.count(x) > 1})
        if dup_g:
            raise Stop(f"grouping.group_order に同じ群が 2 回以上あります: {', '.join(dup_g)}")
    # 色は #RRGGBB の辞書か、配色名か、none だけ (色の名前 red などは受け付けない。261021 本人決定・C2)
    spec = g["group_subgroup_colors"]
    if isinstance(spec, str):
        name = spec.strip().lower().replace("-", "_")
        if name != "none" and name not in PALETTE_NAMES:
            raise Stop(f"grouping.group_subgroup_colors の配色名が不正です: {spec} (使える名前: {', '.join(PALETTE_NAMES)}、または none)")
    elif isinstance(spec, dict):
        bad = {str(k): v for k, v in spec.items() if not (isinstance(v, str) and re.fullmatch(r"#[0-9A-Fa-f]{6}", v.strip()))}
        if bad:
            raise Stop("grouping.group_subgroup_colors の色は #RRGGBB (例 \"#228833\") で書いてください: "
                       + ", ".join(f"{k} = {v!r}" for k, v in bad.items()))
    else:
        raise Stop("grouping.group_subgroup_colors は、全部の色の辞書か、配色名 (例 muted)か、none (色なし)で書いてください")
    # subgroups・subgroup_open_note は必須。使わないときは none (261021 本人決定)。
    # 2 つのスイッチが false でも辞書を書いておいてよい (あとで true にするときのため)
    for k in ("subgroups", "subgroup_open_note"):
        if k not in g:
            raise Stop(f"config の grouping.{k} がありません (使わないときは none)")
    if _is_none(g["subgroups"]):
        if g["use_subgroup"] or g["show_group_subgroup"]:
            raise Stop("use_subgroup または show_group_subgroup が true のときは grouping.subgroups にサブグループを書いてください (none は不可)")
        subs = {}
    elif isinstance(g["subgroups"], dict) and g["subgroups"]:
        subs = {str(k): v for k, v in g["subgroups"].items()}
    else:
        raise Stop(f"grouping.subgroups はサブグループの辞書か none で書いてください (現在: {g['subgroups']!r})")
    note = g["subgroup_open_note"]
    if not isinstance(note, str) or not note.strip():
        raise Stop(f"grouping.subgroup_open_note は白抜きの意味の文か none で書いてください (現在: {note!r})")
    # サブグループごとの marker・open・scale は必須 (marker は auto = 自動で割り当て。261021 本人決定)
    for name, spec in subs.items():
        if not isinstance(spec, dict):
            raise Stop(f"grouping.subgroups.{name} は {{marker: …, open: …, scale: …}} で書いてください (現在: {spec!r})")
        miss = [k for k in ("marker", "open", "scale") if k not in spec]
        if miss:
            raise Stop(f"grouping.subgroups.{name} に {', '.join(miss)} がありません (marker は記号か auto・open は true / false・scale は正の数)")
        if not isinstance(spec["marker"], str) or not spec["marker"].strip():
            raise Stop(f"grouping.subgroups.{name}.marker は記号 (例 \"o\") か auto で書いてください (現在: {spec['marker']!r})")
        if not isinstance(spec["open"], bool):
            raise Stop(f"grouping.subgroups.{name}.open は true / false で書いてください (現在: {spec['open']!r})")
        sc = spec["scale"]
        if isinstance(sc, bool) or not isinstance(sc, (int, float)) or sc <= 0:
            raise Stop(f"grouping.subgroups.{name}.scale は正の数で書いてください (現在: {sc!r})")
    marks = {}
    for name, spec in subs.items():
        if spec["marker"].strip().lower() != "auto":
            marks.setdefault(spec["marker"], []).append(name)
    dup_m = {m: v for m, v in marks.items() if len(v) > 1}
    if dup_m:
        raise Stop("grouping.subgroups で同じ記号 (marker) が複数のサブグループに使われています: "
                   + "; ".join(f"{m} = {', '.join(v)}" for m, v in dup_m.items()))
    return g


def _check_stats(cfg: dict) -> tuple[str, int]:
    s = cfg.get("stats")
    if not isinstance(s, dict) or "method" not in s:
        raise Stop("config に stats.method がありません (tukey | dunnett | welch)")
    method = str(s["method"]).strip().lower()
    if method not in METHODS:
        raise Stop(f"stats.method が不正です: {s['method']!r} (tukey | dunnett | welch)")
    # 検定に必要な最小の n（261009 本人決定・D2a）。これ未満の群は検定から外し、図に「n < 3」と書く。必須 (261021)
    if "min_n" not in s:
        raise Stop("config の stats.min_n がありません (2 以上の整数)")
    raw = s["min_n"]
    if isinstance(raw, bool) or not isinstance(raw, int) or raw < 2:
        raise Stop(f"stats.min_n は 2 以上の整数で書いてください (現在: {raw!r})")
    return method, raw


class Groups:
    """検体 → 群ラベル・群の並び・色・対照群・検定法をまとめて持つ。"""

    def __init__(self, cfg: dict, samples: pd.DataFrame):
        g = _check_grouping(cfg)
        self.use_subgroup = g["use_subgroup"]
        self.show_subgroup = g["show_group_subgroup"]      # 中の名前は show_subgroup のまま
        data_groups = sorted(set(samples["Group_name"]))
        if g["group_order"] is False:
            # false = 検体情報 CSV の行の順（No. の順）で最初に出てきた順（260928 本人決定）。
            # ID (M####P####N####) の P → N の順 (CSV が複数なら P 番号の順 → No. の順。261022 に新しい ID へ)
            pn = samples["sample-id"].str.extract(r"^M\d{4}P(\d{4})N(\d{4})$").astype(int)
            order = list(dict.fromkeys(samples.assign(_p=pn[0].values, _n=pn[1].values).sort_values(["_p", "_n"])["Group_name"]))
        else:
            order = [str(x) for x in (g["group_order"] or [])]
        unlisted = [x for x in data_groups if x not in order]
        if unlisted:
            raise Stop(f"grouping.group_order に無い群がデータにあります: {', '.join(unlisted)} (データにある群はすべて書いてください)")
        extra = [x for x in order if x not in data_groups]
        if extra:
            raise Stop(f"grouping.group_order にデータに無い群があります (綴りの違いの可能性): {', '.join(extra)}")

        self.subgroups = {} if _is_none(g["subgroups"]) else {str(k): dict(v) for k, v in g["subgroups"].items()}
        self.open_note = "" if _is_none(g["subgroup_open_note"]) else g["subgroup_open_note"].strip()
        if self.use_subgroup or self.show_subgroup:
            subs = list(self.subgroups)
            data_subs = sorted(set(samples["Subgroup_name"]))
            unl = [x for x in data_subs if x not in subs]
            if unl:
                raise Stop(f"grouping.subgroups に無い Subgroup がデータにあります (綴りの違いの可能性): {unl}"
                           f" (16S のメタデータの綴りで書いてください。例: SubA・SubA-SubB)")
        if self.use_subgroup:
            subs = list(self.subgroups)
            label = samples["Group_name"] + " " + samples["Subgroup_name"]
            full = [f"{a} {b}" for a in order for b in subs]
            self.order = [x for x in full if x in set(label)]
        else:
            label = samples["Group_name"]
            self.order = order
        self.label = pd.Series(label.values, index=samples["sample-id"].values)          # sample-id → 群ラベル
        self.label_by_newid = pd.Series(label.values, index=samples["newID"].values)     # newID → 群ラベル
        self.group_by_newid = pd.Series(samples["Group_name"].values, index=samples["newID"].values)
        self.subgroup_by_newid = pd.Series(samples["Subgroup_name"].values, index=samples["newID"].values)
        self.group_order = order                       # Group_name の並び

        # 色は grouping.group_subgroup_colors の 1 つだけ（260924 本人決定）。何に色を付けるかは設定で決まる:
        #   use_subgroup: true かつ show_group_subgroup: false → 処理×サブグループごとに色
        #   それ以外 → 元の群（Group_name）ごとに色（サブグループは元の群の色を受け継ぐ）
        # 書き方は「全部を個別に指定（辞書）」か「配色名（文字列）」。配色の色数が足りなければ止める（24 色以上は smooth_rainbow を勧める）
        self.color_by_subgroup = self.use_subgroup and not self.show_subgroup
        units = self.order if self.color_by_subgroup else order
        what = "サブグループ" if self.color_by_subgroup else "群"
        spec = g["group_subgroup_colors"]
        self.color_note = ""
        # none = 群の色なし（261013 本人決定・案 A）: 全群を前景色 #222222 で描き、色だけで群を見分ける図と群の凡例は描かない。
        # キーを書かない（null を含む）場合は従来どおり止める（config は明示的に書かないと止める方針）
        self.no_color = isinstance(spec, str) and spec.strip().lower() == "none"
        if self.no_color:
            colors = {x: NO_COLOR for x in units}
            self.color_note = (f"[注記] group_subgroup_colors: none のため、全{what}を {NO_COLOR} で描きます。色だけで群を見分ける図"
                               " (#2・#4・#6 検体ごと・#7 検体ごと・#9 ベン図) と群の凡例は描きません")
        elif isinstance(spec, str):
            pal = spec.strip().lower().replace("-", "_")
            if pal not in PALETTE_NAMES:
                raise Stop(f"grouping.group_subgroup_colors の配色名が不正です: {spec} (使える名前: {', '.join(PALETTE_NAMES)})")
            colors = dict(zip(units, palette_colors(pal, len(units))))
            self.color_note = f"[注記] {what}の色は配色 {pal} から並び順に自動で割り当てました ({len(units)} 色)"
        elif isinstance(spec, dict):
            colors = {str(k): str(v) for k, v in spec.items()}
            miss = [x for x in units if x not in colors]
            if miss:
                raise Stop(f"grouping.group_subgroup_colors に色の無い{what}があります: {miss}"
                           f" (全{what}を個別に指定するか、配色名を書いてください)")
        else:
            raise Stop("grouping.group_subgroup_colors は、全部の色の辞書か、配色名 (例 muted)か、none (色なし)で書いてください")
        dup = {}
        for x in units:
            dup.setdefault(colors[x].strip().upper(), []).append(x)
        dup = {c: v for c, v in dup.items() if len(v) > 1}
        if dup and not self.no_color:                                            # 同じ色が 2 つ以上にあれば止める（260924 本人指示）
            raise Stop(f"grouping.group_subgroup_colors で同じ色が複数の{what}に使われています: "
                       + "; ".join(f"{c} = {', '.join(v)}" for c, v in dup.items()))
        if self.color_by_subgroup:
            self.color = {x: colors[x] for x in self.order}
            self.group_color = {}                          # 元の群の色は使わない
        else:
            self.group_color = {x: colors[x] for x in order}
            base = samples.set_index("sample-id")["Group_name"]
            first_group = {lab: base[sid] for sid, lab in self.label.items()}
            self.color = {lab: colors[first_group[lab]] for lab in self.order}

        self.method, self.min_n = _check_stats(cfg)
        # 群が 2 つなら welch だけ、welch は群が 2 つのときだけ（261009 本人決定。tukey は 2 群で Student の t と同じになり、
        # dunnett は 2 群で Welch に切り替わっていたため、設定しだいで検定が変わっていた）
        if len(self.order) == 2 and self.method != "welch":
            raise Stop(f"群が 2 つ ({', '.join(self.order)})のときは stats.method: welch を指定してください (現在: {self.method})")
        if self.method == "welch" and len(self.order) != 2:
            raise Stop(f"stats.method: welch は群が 2 つのときだけ使えます (現在 {len(self.order)} 群)。tukey か dunnett を指定してください")

        # 対照群は、書いてあれば方法によらずデータの群にあるか確かめる (tukey でも #11・#12 で使う。261021 本人決定・C3)
        written = _control_written(cfg)
        if written and written not in self.order:
            hint = " (use_subgroup: true のときは「群 サブグループ」の形で書きます。例 Control SubA)" if self.use_subgroup else ""
            raise Stop(f"対照群 '{written}' (comparison.control_group) がデータの群にありません。存在する群: {self.order}{hint}")
        self.control = None
        if self.method in ("dunnett", "welch"):     # welch も #8 の星を付ける相手として対照群が要る（261009 本人決定）
            if not written:
                raise Stop(f"stats.method: {self.method} のときは comparison.control_group に対照群を書いてください (none は不可)")
            self.control = written                  # tukey では検定に対照群を使わない (#8 は文字で示す)

    def check_control(self, excluded: set[str]) -> None:
        """dunnett の対照群の n を、rarefaction で除外した検体を引いた数で確かめる (261009 本人決定・D2c・時点 B)。
        足りなければ設定の誤りとみなして止める。excluded = 06_rarefaction/excluded_samples.tsv の sample-id。"""
        if self.method not in ("dunnett", "welch"):
            return
        n_ctrl = int(((self.label == self.control) & ~self.label.index.isin(list(excluded))).sum())
        if n_ctrl < self.min_n:
            raise Stop(f"対照群 '{self.control}' の検体数 {n_ctrl} (rarefaction の除外後)が stats.min_n ({self.min_n})"
                       "未満です。別の対照群にするか、min_n を下げてください")

    def describe(self) -> list[str]:
        n = self.label.value_counts()
        lines = ([self.color_note] if self.color_note else []) + [f"群の単位 : {'Group × Subgroup' if self.use_subgroup else 'Group_name'} ({len(self.order)} 群)",
                 "群       : " + ", ".join(f"{x}(n={int(n.get(x, 0))})" for x in self.order),
                 f"検定     : {self.method}" + (f" (対照群 = {self.control})" if self.control else " (全ペア)")]
        return lines


def switches(cfg: dict) -> dict[str, bool]:
    """pipeline.run を読む。実装済みの図のキーが欠けていたら止める。
    未実装の図のキーが書かれていても読まない (ログで知らせる)。"""
    run = (cfg.get("pipeline") or {}).get("run")
    if not isinstance(run, dict):
        raise Stop("config に pipeline.run がありません")
    unknown = [k for k in run if k not in ALL_FIGURES]
    if unknown:
        raise Stop(f"pipeline.run に知らないキーがあります: {unknown} (使える名前: {', '.join(ALL_FIGURES)})")
    missing = [k for k in IMPLEMENTED if k not in run]
    if missing:
        raise Stop(f"pipeline.run に実装済みの図のキーがありません: {missing}")
    for k in IMPLEMENTED:
        if not isinstance(run[k], bool):
            raise Stop(f"pipeline.run.{k} は true / false で書いてください (現在: {run[k]!r})")
    return {k: run[k] for k in IMPLEMENTED}
