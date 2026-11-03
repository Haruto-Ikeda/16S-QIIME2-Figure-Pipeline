#!/usr/bin/env python3
"""figures.py — Step 10 作図 (QIIME 2 解析の出力から図を描く)。

設計の正本は ../README.md「作図パイプラインの設計」節。
- 出力先は output/16S/figure/ver<N>/ (固定。260928 に 10_figures → figure)。ファイル名は「作図一覧の番号_図の名前」。
- **版は設定から自動で決める** (260917 判断・案 B): 図の中身を変える設定 (grouping・comparison・stats・plotting と上流の条件 upstream。261014 E1b)が
  既存の版の settings.yaml と同じならその版に描き足し (図のフォルダがある番号は飛ばす)、どれとも違えば新しい版を作る。実行の途中で尋ねない。
- 図のフォルダ (例 04_pcoa/)が版の中に既にあれば、その番号の図は飛ばす (260929 決定・案 A)。
  描き直すときは、その図のフォルダか版のフォルダをゴミ箱へ移して実行する (--redraw・--force は 261015 に廃止)。描くときは検定も config の群設定から計算し直す (図と検定を一致させるため)。
- 図 1 枚につき PDF・図の根拠の値 (.values.csv)・検定結果 (.stats.csv)・楕円 (.ellipses.csv)を出す。
  **凡例・settings.yaml・名前の対応表 (taxa_names_<階級>.csv。261021) は ver<N>/00_info/、図は番号ごとのフォルダ (例 04_pcoa/)の直下に PDF・csv/ に CSV** (260928)。
"""
from __future__ import annotations

import itertools
import math
import re
import sys
from pathlib import Path

import matplotlib
import numpy as np
import yaml
import pandas as pd
import scipy

from . import groupstats as GS
from . import plotcfg as PC
from . import plotstyle as PS
from .steps import merged, read_denoising_stats


RULE = "-" * 65
FILE_PATTERNS = {                     # 同じ図を指標・距離・階級ごとに描くときの見出しのファイル名（260925 ログ整理）
    "#2 レアファクションカーブ": "02_rarefaction_<指標>.pdf",
    "#3 α多様性": "03_alpha_<指標>.pdf",
    "#4 PCoA PCo1×PCo2": "04_pcoa_<距離>.pdf",
    "#4 PCoA 交差形式 PCo1〜4": "04_pcoa-matrix_<距離>.pdf",
    "#5 PERMANOVA / PERMDISP 行列": "05_permanova_<距離>.pdf",
    "#6 デンドログラム": "06_dendrogram_<距離>.pdf",
    "#6 デンドログラム・群の平均": "06_dendrogram-group_<距離>.pdf",
    "#7 積み上げ棒・群平均": "07_stacked-group_<階級>.pdf",
    "#7 積み上げ棒・検体ごと": "07_stacked-sample_<階級>.pdf",
    "#8 ヒートマップ": "08_heatmap_<階級>.pdf",
    "#9 ベン図 / UpSet": "09_venn_<階級>.pdf・09_upset_<階級>.pdf",
    "#10 タクソン別の棒": "10_taxon-bar/10_taxon-bar_<階級>_<タクソン>.pdf",
    "#11 ANCOM-BC": "11_ancombc_<階級>.pdf",
    "#12 ANCOM-BC まとめ": "12_ancombc-summary_<階級>.pdf",   # 261019・F045・F088
}


class FigLog:
    """Step 10 の記録を 09 以前と同じ形に整える (260925 決定)。
    1 図 = 見出し 1 行 (指標・距離・階級ごとの繰り返しは 1 行にまとめる)＋ 要点 1 行 (sum() で集めて次の見出しの前に出す)。
    細かい数値は CSV に任せて出さない。[注意]・[エラー]・[注記] の行だけはそのまま出す。 (版のフォルダの log.txt は 260928 に廃止。記録は logs/ にある)"""

    def __init__(self, base):
        self.base, self.lines, self.parts, self.n_warn = base, [], [], 0

    def _out(self, msg: str) -> None:
        self.base(msg)
        self.lines.append(msg)

    def __call__(self, msg: str = "") -> None:
        if "[エラー]" in msg:
            self._out(msg)
        elif "[注意]" in msg or "[注記]" in msg:            # 警告・[注記] は最後にまとめて出す（Log が集める。260927 決定）
            self.base(msg)

    def show(self, msg: str = "") -> None:        # 見出しのブロックなど、必ず出す行
        self._out(msg)

    def sum(self, text: str) -> None:             # 要点は出さない（260925 指示: 細かい数値は CSV に任せる）
        pass

    def flush(self) -> None:
        if self.parts:
            self._out("      " + " / ".join(self.parts))
            self.parts = []

    def step(self, step: str, n: int, total: int, done: bool, name: str, detail: str = "") -> None:
        m = re.fullmatch(r"(.+) \((\d+)/(\d+)\)", name)
        if m:
            if m.group(2) != "1":
                return
            name = m.group(1)
            detail = f"{FILE_PATTERNS.get(name, detail)} ({m.group(3)} 種)"
        self.flush()
        mark = "スキップ" if done else "　実行　"
        self._out(f"{step} {n:>2}/{total}|{mark}|{name}" + (f": {detail}" if detail else ""))

    def head(self, title: str) -> None:
        self._out("")
        self._out(title)


def step10_figures(ctx: dict) -> None:
    import logging
    logging.getLogger("fontTools").setLevel(logging.ERROR)   # PDF のフォント処理の無害なメッセージを出さない（260925）
    log = FigLog(ctx["log"])
    ctx = {**ctx, "log": log}
    out, cfg = ctx["out"], ctx["cfg"]
    log.head("10. Figures (作図)")
    root = out.parent / "figure"                    # output/16S/figure/（260928 決定。out = output/16S/Qiime2）
    root.mkdir(parents=True, exist_ok=True)       # ステップ冒頭で明示（アーカイブ準拠）
    try:
        run = PC.switches(cfg)
        G = PC.Groups(cfg, PC.load_samples(out))
        G.check_control(set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["sample-id"]))
    except PC.Stop as e:
        log(f"  [エラー] {e}")
        raise SystemExit(1)
    log.show(RULE)
    d, ver_text = resolve_version(root, cfg, log, _upstream(ctx))
    info = d / "00_info"                                 # 凡例・settings.yaml（260928 決定）・名前の対応表 taxa_names_<階級>.csv（261021）
    nn = G.label.value_counts()
    excl_ids = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["sample-id"])
    na = G.label[~G.label.index.isin(list(excl_ids))].value_counts()   # 解析 (#1 以外) で使う n (261019・F039・F046)
    log.show(f"版          {ver_text}")
    log.show(f"群          {'Group × Subgroup' if G.use_subgroup else 'Group_name'} {len(G.order)} 群: "
             + ", ".join(f"{x} {int(nn.get(x, 0))}" + (f" (解析 {int(na.get(x, 0))})" if na.get(x, 0) != nn.get(x, 0) else "")
                         for x in G.order))
    log.show(f"表示        use_subgroup {str(G.use_subgroup).lower()} / show_group_subgroup {str(G.show_subgroup).lower()}"
             f" (色 = {'サブグループ' if G.color_by_subgroup else '群'})")
    log.show(f"検定        {G.method}" + (f" (対照群 = {G.control})" if G.control else " (全ペア)"))
    if G.color_note:
        log(f"  {G.color_note}")

    if G.no_color:                                                   # 群の色なしでは群の凡例を出さない（261013 決定）
        pass
    else:
        PS.save_legend(info / "00_legend_group.pdf", G, log=log)              # 凡例は別ファイル（全図共通の群の凡例）
    if G.show_subgroup and not G.no_color:
        subs = [s for s in G.subgroups if s in set(G.subgroup_by_newid)]
        PS.save_legend(info / "00_legend_group-subgroup.pdf", G, PS.subgroup_markers(G.subgroups), subs, log=log)
    log.show("凡例        " + ("なし (group_subgroup_colors: none)" if G.no_color else
                              "00_legend_group.pdf" + ("・00_legend_group-subgroup.pdf" if G.show_subgroup else "")))
    _write_taxa_names(ctx, info)
    log.show("名前の対応  00_info/taxa_names_<階級>.csv (図の表示名 ↔ 元の SILVA 系統)")
    n_si = _write_sample_information(ctx, info)
    log.show(f"検体情報    00_info/sample_information.csv (解析した {n_si} 検体・リード数つき)")
    log.show(RULE)

    total = 12                          # 図の番号は #1〜#12（#12 は 261007 追加）（#4 と #7 はそれぞれ 2 種類の図。260924 指示で番号を詰めた）
    # 版のフォルダに図のフォルダが既にあれば、その番号の図は飛ばす（260929 決定・案 A）。描き直すときはフォルダをゴミ箱へ移して実行する（261015）。
    # 判定は描き始める前のフォルダの有無で決める（#4・#6・#7 は 1 つのフォルダに 2 種類の図が入るため）
    existing = {k for k, v in FIG_DIRS.items() if (d / v).is_dir()}
    # 1 つのフォルダに 2 種類の図が入る番号 (#4・#7) は、その図の PDF があるかで判定する (261016 決定・E9c・点検 F271・F287)。
    # 片方だけ描いたあとに pipeline.run でもう片方を ON にしても描かれるように。判定は描き始める前に決める
    pdf_exists = {pat: any((d / FIG_DIRS[pat[:2]]).glob(pat)) for pat in SHARED_FIG_PATTERNS}

    def no_color(num: int, title: str) -> bool:
        """群の色なしで描かない図 (色だけで群を見分ける図。261013 決定)。"""
        if G.no_color:
            log.step("10", num, total, True, title, "group_subgroup_colors: none")
            log(f"  [注記] {title}: group_subgroup_colors: none のため描きません (色だけで群を見分ける図のため)")
            return True
        return False

    def done(num: int, title: str, pattern: str = "") -> bool:
        k = f"{num:02d}"
        if (pdf_exists[pattern] if pattern else k in existing):
            log.step("10", num, total, True, title, f"既存) {FIG_DIRS[k]}/")
            return True
        return False

    if not run["read_retention"]:
        log.step("10", 1, total, True, "#1 最終リード数と sampling depth", "config の pipeline.run で OFF")
    elif not done(1, "#1 最終リード数と sampling depth"):
        fig01_read_depth(ctx, G, d, 1, total)
    if run["rarefaction"]:
        if not no_color(2, "#2 レアファクションカーブ") and not done(2, "#2 レアファクションカーブ"):
            fig02_rarefaction(ctx, G, d, 2, total)
    else:
        log.step("10", 2, total, True, "#2 レアファクションカーブ", "config の pipeline.run で OFF")
    if not run["alpha"]:
        log.step("10", 3, total, True, "#3 α多様性", "config の pipeline.run で OFF")
    elif not done(3, "#3 α多様性"):
        fig03_alpha(ctx, G, d, 3, total)
    if run["pcoa"]:
        if not no_color(4, "#4 PCoA PCo1×PCo2") and not done(4, "#4 PCoA PCo1×PCo2", "04_pcoa_*.pdf"):
            fig04_pcoa(ctx, G, d, 4, total)
    else:
        log.step("10", 4, total, True, "#4 PCoA PCo1×PCo2", "config の pipeline.run で OFF")
    if run["pcoa_matrix"]:
        if not no_color(4, "#4 PCoA 交差形式 PCo1〜4") and not done(4, "#4 PCoA 交差形式 PCo1〜4", "04_pcoa-matrix_*.pdf"):
            fig05_pcoa_matrix(ctx, G, d, 4, total)
    else:
        log.step("10", 4, total, True, "#4 PCoA 交差形式 PCo1〜4", "config の pipeline.run で OFF")
    if not run["permanova"]:
        log.step("10", 5, total, True, "#5 PERMANOVA / PERMDISP 行列", "config の pipeline.run で OFF (計算も省く)")
    elif not done(5, "#5 PERMANOVA / PERMDISP 行列"):   # Subgroup を群にしても描く（260921 指示）
        fig06_permanova(ctx, G, d, 5, total)
    if run["dendrogram"]:
        if not done(6, "#6 デンドログラム"):
            if not no_color(6, "#6 デンドログラム・検体ごと"):
                fig07_dendrogram(ctx, G, d, 6, total)
            fig07_dendrogram_group(ctx, G, d, 6, total)     # 群の平均の版（260921 追加。同じスイッチで出す）
        else:                                              # 飛ばしたときも 2 種類の図をそれぞれ 1 行出す (261019・F256)
            done(6, "#6 デンドログラム・群の平均")
    else:
        log.step("10", 6, total, True, "#6 デンドログラム", "config の pipeline.run で OFF")
        log.step("10", 6, total, True, "#6 デンドログラム・群の平均", "config の pipeline.run で OFF")
    for key, num, title, fn in (("stacked_group", 7, "#7 積み上げ棒・群平均", fig08_stacked_group),
                                ("stacked_sample", 7, "#7 積み上げ棒・検体ごと", fig09_stacked_sample),
                                ("heatmap", 8, "#8 ヒートマップ", fig10_heatmap),
                                ("venn", 9, "#9 ベン図 / UpSet", fig11_venn_upset),
                                ("taxon_barplot", 10, "#10 タクソン別の棒", fig12_taxon_bar),
                                ("ancombc", 11, "#11 ANCOM-BC", fig13_ancombc),
                                ("ancombc_summary", 12, "#12 ANCOM-BC まとめ", fig14_ancombc_summary)):
        if not run[key]:
            log.step("10", num, total, True, title, "config の pipeline.run で OFF")
            continue
        if done(num, title, {"stacked_group": "07_stacked-group_*.pdf", "stacked_sample": "07_stacked-sample_*.pdf"}.get(key, "")):
            continue
        if key == "stacked_sample" and no_color(num, title):
            continue
        try:
            fn(ctx, G, d, num, total)
        except PC.Stop as e:
            log(f"  [エラー] {e}")
            raise SystemExit(1)
    log.flush()


VERSION_KEYS = ["grouping", "comparison", "stats", "plotting", "parent", "select", "upstream"]   # 図の中身を変える設定（pipeline.run は含めない）
# parent・select = 絞り込みの形 (イ) の親と選ぶ検体 (261026 決定・7-5)。通常の形ではどちらも none
# upstream = 上流の条件（261014 決定・E1b 案 ウ）。config からではなく 01〜09 の条件と出力から作る（_upstream）
# 261021: config の設定はすべて必須にした (書かないと 01 の前に止まる) ので、版を比べるときの既定値 (VERSION_DEFAULTS) は廃止した


def _upstream(ctx: dict) -> dict:
    """図の中身を変える上流の条件 (261014 決定・E1b 案 ウ)。違えば新しい版にする。"""
    from .steps import runs_of
    cfg, out = ctx["cfg"], ctx["out"]
    runs = runs_of(cfg)                              # ランごとの quality_control と exclude (261026 決定)
    dv = cfg.get("diversity") or {}
    excl = pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)
    return {"quality_control": {r["name"]: dict(r["qc"]) for r in runs},
            "exclude": {r["name"]: list(r["exclude"]) for r in runs},
            "sampling_depth_floor": dv.get("sampling_depth_floor"),
            "sampling_depth": int((out / "06_rarefaction" / "sampling_depth.txt").read_text().strip()),
            "excluded_samples": len(excl),
            "random_seed": dv.get("random_seed"),
            "classifier": Path(str(cfg["classifier"]["path"])).name,
            "taxonomy_confidence": (cfg.get("taxonomy") or {}).get("confidence")}


def _settings(cfg: dict, upstream: dict | None = None) -> dict:
    st = {k: cfg[k] for k in VERSION_KEYS if k in cfg}
    # 辞書の比較は並びを無視するので、subgroups の並び (凡例・#6・#7 の検体の並び・use_subgroup の群の並びを変える) を別に残す (261021・C8)
    subs = (cfg.get("grouping") or {}).get("subgroups")
    st["subgroup_order"] = list(subs) if isinstance(subs, dict) else subs
    if upstream is not None:
        st["upstream"] = upstream
    return st


def _upstream_diff(last: dict, cur: dict) -> str:
    """版の行に出す上流の条件の違い (変わった値だけ)。"""
    a, b = (last or {}).get("upstream") or {}, (cur or {}).get("upstream") or {}
    if not a:
        return "上流の条件 (前の版には記録が無い)"

    def flat(dct, pre=""):                          # 入れ子 (quality_control.<ラン>.trim_left_f) も 1 段ずつ開く (261026)
        res = {}
        for k, v in dct.items():
            if isinstance(v, dict):
                res.update(flat(v, f"{pre}{k}."))
            else:
                res[f"{pre}{k}"] = v
        return res
    def fmt(v):                                       # リストは [ ] や ' ' を出さずに並べる (261023 決定に合わせる)
        return (", ".join(map(str, v)) or "なし") if isinstance(v, list) else v
    a, b = flat(a), flat(b)
    parts = []
    for k in list(b) + [k for k in a if k not in b]:
        if a.get(k) != b.get(k):
            parts.append(f"{k}: {fmt(a.get(k, '記録なし'))} → {fmt(b.get(k, '記録なし'))}")   # 記録なし = その版の記録にキーが無い (空のリストの「なし」と区別)
    return "上流の条件 (" + ", ".join(parts) + ")"


def resolve_version(root: Path, cfg: dict, log, upstream: dict | None = None) -> tuple[Path, str]:
    """設定から版を決め、その版のフォルダを返す (260917 判断・案 B)。
    既存の ver<N>/00_info/settings.yaml と設定が同じ → その版に描き足す / どれとも違う → 新しい版。"""
    cur = _settings(cfg, upstream)
    versions = {}
    for p in root.glob("ver*/00_info/settings.yaml"):
        m = re.fullmatch(r"ver(\d+)", p.parent.parent.name)
        if m:
            try:
                versions[int(m.group(1))] = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
            except yaml.YAMLError as e:
                log(f"  [エラー] {p} を読めません: {e}")
                raise SystemExit(1)
    same = [v for v, s in versions.items() if s == cur]
    if same:
        n = same[0]
        text = f"ver{n} (既存の版と設定が同じ → この版に描き足す。図のフォルダがある番号は飛ばす)"   # 261019・F257
    else:
        n = max(versions, default=0) + 1
        if versions:
            last = versions[max(versions)]
            diff = [(_upstream_diff(last, cur) if k == "upstream" else k)
                    for k in VERSION_KEYS if (last or {}).get(k) != cur.get(k)]
            if "grouping" not in diff and (last or {}).get("subgroup_order") != cur.get("subgroup_order"):
                diff.append("grouping (subgroups の並び)")
            text = f"ver{n} (新しい版。最新の版 ver{max(versions)} と違う設定: {', '.join(diff)})"
        else:
            text = f"ver{n} (最初の版)"
    d = root / f"ver{n}"
    (d / "00_info").mkdir(parents=True, exist_ok=True)
    header = ("# この版の作図に使った設定 (Step 10 が自動で書く。編集しない)\n"
              "# 設定が同じなら同じ版に描き足し、違えば新しい版を作る (README「作図パイプラインの設計」)\n")
    (d / "00_info" / "settings.yaml").write_text(header + yaml.safe_dump(cur, allow_unicode=True, sort_keys=False),
                                     encoding="utf-8")
    return d, text


def _shown_labels(G, order: list | None = None) -> list[str]:
    """横軸の群名。use_subgroup: true・show_group_subgroup: true ならサブグループ名だけ (処理は色で分かる。他の図と同じ)。
    order = 並べる群 (None なら G.order)。"""
    order = list(G.order) if order is None else list(order)
    return [g.split(" ", 1)[1] for g in order] if (G.use_subgroup and G.show_subgroup) else order


def _warn_labels(log, labels: list[str]) -> None:
    """45° の群名が下余白 20 mm に収まらない可能性があれば警告する (260921 指示。図は見切れうる)。"""
    over = PS.labels_overflow(labels)
    if over:
        log(f"  [注意] 群名が長く、下余白 {PS.MARGIN_B} mm で見切れている可能性があります (群名を短くしてください): "
            + ", ".join(f"{s} (約 {h:.1f} mm)" for s, h in over))


FIG_DIRS = {"01": "01_read-depth", "02": "02_rarefaction", "03": "03_alpha", "04": "04_pcoa",
            "05": "05_permanova", "06": "06_dendrogram", "07": "07_stacked-bar", "08": "08_heatmap",
            "09": "09_venn-upset", "10": "10_taxon-bar", "11": "11_ancombc",
            "12": "12_ancombc-summary"}


SHARED_FIG_PATTERNS = ("04_pcoa_*.pdf", "04_pcoa-matrix_*.pdf", "07_stacked-group_*.pdf", "07_stacked-sample_*.pdf")


def _fig_dir(d: Path, name: str) -> Path:
    """図の番号ごとのフォルダ (260928 決定: 同じ番号の図は同じフォルダ。PDF はこの直下、CSV は csv/)。"""
    sub = d / FIG_DIRS[name[:2]]
    sub.mkdir(exist_ok=True)
    return sub


def _csv_dir(d: Path, name: str) -> Path:
    """CSV は図のフォルダの中の csv/ (260928 決定)。"""
    sub = _fig_dir(d, name) / "csv"
    sub.mkdir(exist_ok=True)
    return sub


def _write_stats(res: dict, path: Path, extra: dict) -> pd.DataFrame:
    df = pd.DataFrame(res["rows"], columns=["groupA", "groupB", "n_A", "n_B", "p", "test", "reason"])
    df["stars"] = df["p"].map(PS.p_stars)
    for k, v in extra.items():
        df[k] = v
    df.to_csv(path, index=False)
    return df


def _tops(data: list) -> list[float]:
    """群ごとの上端 = max(最大値, 平均 + SD) (SD のエラーバーが点より上に出ることがあるため)。"""
    out = []
    for v in data:
        v = np.asarray(v, float)
        out.append(max(v.max(), v.mean() + (v.std(ddof=1) if len(v) > 1 else 0)) if len(v) else np.nan)
    return out


def _ylim_for(G, res: dict, data: list) -> tuple[float, float, float]:
    """データ範囲 (平均 + SD を含む)と有意差の数から y 軸を決める (ブラケットは段数、対照群との星は 1 段ぶん確保)。"""
    allv = np.concatenate([x for x in data if len(x)] + [np.array([t for t in _tops(data) if np.isfinite(t)])])
    levels = len(_sig_items(G, res))                 # どの方法もブラケットを段に積む（261009 決定）。「n < 3」の印も 1 段
    return PS.choose_ylim(float(allv.min()), float(allv.max()), levels)


def _draw_significance(ax, G, res: dict, xpos: dict, data: list, lo: float, hi: float, order: list | None = None) -> float:
    """order = data の並び (None なら G.order)。#10 は検体の残った群だけで並べるので、その並びを渡す (261019・F183・F234)。"""
    order = list(G.order) if order is None else list(order)
    idx = {g: i for i, g in enumerate(order)}
    tops = {g: (t if np.isfinite(t) else lo) for g, t in zip(order, _tops(data))}
    data_top = max(tops.values())
    return PS.draw_brackets(ax, xpos, [(idx[a], idx[b], s, line) for a, b, s, line in _sig_items(G, res)],
                            data_top, lo, hi)


def _sig_items(G, res: dict) -> list:
    """ブラケットに描くもの: (群 A, 群 B, 文字, 線を引くか)。
    tukey = 有意な全ペア、dunnett = 対照群と各群、welch = 2 群の 1 組。どれもブラケットで結ぶ
    (261009 決定。それまで dunnett は群の上に星だけだった)。
    検定から外した群 (n < stats.min_n)は線を引かず、どの方法でもその群の真上に「n < 3」とだけ書く
    (ブラケットの最下段と同じ高さ。261009 決定)。"""
    items = [(r["groupA"], r["groupB"], PS.p_stars(r["p"]), True) for r in res["rows"] if PS.p_stars(r["p"])]
    label = f"n < {res.get('min_n', 2)}"
    return items + [(g, g, label, False) for g in res["skipped"]]


# ============================================================ #1
def fig01_read_depth(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#1 最終リード数と sampling depth (260917 決定)。
    縦軸 = キメラ除去後に残ったリード数。depth (破線)と下限 floor (点線)を引き、
    rarefaction で除外された検体は白抜き。検定は stats.method に従い、除外された検体も含める (260917 決定)。"""
    log, out, cfg = ctx["log"], ctx["out"], ctx["cfg"]
    name = "01_read-depth"
    log.step("10", n, total, False, "#1 最終リード数と sampling depth", f"{name}.pdf")

    t = read_denoising_stats(merged(ctx)["stats"])
    depth = int((out / "06_rarefaction" / "sampling_depth.txt").read_text().strip())
    floor = int(cfg["diversity"]["sampling_depth_floor"])
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["sample-id"])

    missing = sorted(set(t["sample-id"]) - set(G.label.index))
    if missing:
        log(f"  [エラー] denoising-stats にあってメタデータに無い検体があります: {missing}")
        raise SystemExit(1)
    t["group"] = t["sample-id"].map(G.label)
    t["excluded"] = t["sample-id"].isin(excl)
    if len(t) != len(G.label):
        log(f"  [エラー] 検体数が合いません: denoising-stats {len(t)} / メタデータ {len(G.label)}")
        raise SystemExit(1)

    col = "non-chimeric"
    data = [t.loc[t.group == g, col].to_numpy(float) for g in G.order]
    hollow = [t.loc[t.group == g, "excluded"].to_numpy() for g in G.order]
    labels = _shown_labels(G)                     # 群名に n 数は併記しない（260917 指示・全図共通。n は values.csv とログに残す）

    _warn_labels(log, labels)
    fig, ax, xpos = PS.box_panel(labels, "Non-chimeric reads")
    PS.dots_mean_sd(ax, data, [G.color[g] for g in G.order], hollow)
    ax.axhline(depth, color="#999999", lw=0.7, ls=(0, (4, 2)), zorder=0)   # 破線と注記は #999999（260921 指示）
    ax.axhline(floor, color="#999999", lw=0.6, ls=(0, (1, 2)), zorder=0)

    res = GS.compare(dict(zip(G.order, data)), G.order, G.method, G.control, G.min_n)
    lo, hi, step = _ylim_for(G, res, data)
    ax.set_ylim(lo, hi)
    hi = _draw_significance(ax, G, res, xpos, data, lo, hi)
    lo, hi = PS.fixed_ticks(ax, lo, hi)               # 目盛 5 個（260924 決定）
    PS.keep_brackets(ax)                         # 伸ばした軸でも、ブラケットは描いたときの見た目のまま (261020・F122)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1000:g}k" if v else "0"))
    # 注記は depth だけ、破線の下・左端に書く（260921 指示。floor の数値は書かない＝点線だけ残す）
    ax.annotate(f"depth {depth:,}", xy=(-0.45, depth), xytext=(0, -2.5), textcoords="offset points",
                fontsize=PS.FS_NOTE, va="top", ha="left", color="#999999")   # 破線から 2.5 pt 下げる（260921 指示）
    fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
    matplotlib.pyplot.close(fig)

    t[["sample-id", "group", col, "excluded"]].to_csv(_csv_dir(d, name) / f"{name}.values.csv", index=False)
    st = _write_stats(res, _csv_dir(d, name) / f"{name}.stats.csv",
                      {"metric": col, "dunnett_seed": GS.DUNNETT_SEED if G.method == "dunnett" else ""})

    n_plot = int(sum(len(v) for v in data))
    if n_plot != len(t):                            # 網羅の不変条件
        log(f"  [エラー] 図に載った検体数が denoising-stats の検体数と合いません: {n_plot} ≠ {len(t)}")
        raise SystemExit(1)
    log.sum(f"検体 {n_plot} (不変条件 OK)")
    log(f"  検体 {n_plot} = denoising-stats {len(t)} (不変条件 OK)/ depth {depth:,} / floor {floor:,}"
        f" / 白抜き (除外検体){int(t.excluded.sum())} 件 (検定には含める)")
    _log_skipped(log, res)
    sig = st[st["stars"] != ""]
    log.sum(f"有意 {len(sig)}/{int(st['p'].notna().sum())}")
    for _, r in sig.iterrows():                   # ログには有意なものだけ（use_subgroup だと数百比較になるため）
        log(f"    {r.groupA} vs {r.groupB}: p = {r.p:.3g} {r.stars}")


# ============================================================ #2
# 指標の表示名の単一正本（260921 決定）。ファイル名は QIIME 2 の名前のまま、図に出す名前だけをここで決める。
# 人名・手法名（Shannon・Faith・Pielou・Jaccard・Bray–Curtis・UniFrac）は大文字、それ以外は文頭だけ大文字。
# Bray–Curtis はエンダッシュ（2 人の人名）。
METRIC_LABEL = {
    "observed_features": "Observed ASVs", "shannon": "Shannon", "faith_pd": "Faith's PD",
    "evenness": "Pielou's evenness",
    "jaccard": "Jaccard", "bray_curtis": "Bray–Curtis",
    "unweighted_unifrac": "Unweighted UniFrac", "weighted_unifrac": "Weighted UniFrac",
}
RAREFACTION_METRICS = [(m, METRIC_LABEL[m]) for m in ("observed_features", "shannon", "faith_pd")]


def fig02_rarefaction(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#2 レアファクションカーブ (260917 決定・案 A)。
    群の平均 ± SD を 1 指標 1 枚 (observed_features・shannon・faith_pd)。検定は付けない。
    値は Step 06 の alpha-rarefaction (Qiime2/HTML_06_alpha-rarefaction.qzv/<指標>.csv)。
    各深さの 10 反復を検体ごとに平均してから、群で平均・SD を取る。
    [注意] カーブは --p-max-depth = sampling depth で止まる (depth より先は描けない)。"""
    log, out = ctx["log"], ctx["out"]
    src = out / "HTML_06_alpha-rarefaction.qzv"
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["newID"])
    for k, (metric, ylabel) in enumerate(RAREFACTION_METRICS, 1):
        name = f"02_rarefaction_{metric}"
        log.step("10", n, total, False, f"#2 レアファクションカーブ ({k}/{len(RAREFACTION_METRICS)})", f"{name}.pdf")
        t = pd.read_csv(src / f"{metric}.csv", index_col=0)
        t = t[[c for c in t.columns if c.startswith("depth-")]]
        depths = sorted({int(re.search(r"depth-(\d+)_", c).group(1)) for c in t.columns})
        per = pd.DataFrame({dp: t[[c for c in t.columns if c.startswith(f"depth-{dp}_")]].mean(axis=1)
                            for dp in depths})                       # 検体ごとに反復を平均
        per = per.drop(index=list(excl), errors="ignore")       # 除外した検体を明示的に外す (261019・F181・F215。リードが残っていても)
        per = per.dropna(how="all")
        unknown = sorted(set(per.index) - set(G.label_by_newid.index))
        if unknown:
            log(f"  [エラー] {metric}.csv にあってメタデータに無い検体があります: {unknown}")
            raise SystemExit(1)
        expect = len(G.label_by_newid) - len(excl)
        if len(per) != expect:                                        # 網羅の不変条件
            log(f"  [エラー] 曲線の検体数が合いません: {len(per)} ≠ 全検体 {len(G.label_by_newid)} − 除外 {len(excl)} = {expect}")
            raise SystemExit(1)
        grp = per.index.map(G.label_by_newid)

        labels = list(G.order)                    # 群名に n 数は併記しない（260917 指示・全図共通）
        units = int(np.ceil(depths[-1] / PS.READS_PER_UNIT))          # 横幅 = 5,000 リード 16 mm × 区間の数（260929）
        xmax = units * PS.READS_PER_UNIT
        fig, ax = PS.fixed_axes(units * PS.MM_PER_UNIT)
        rows = []
        for g, lab in zip(G.order, labels):
            s = per[grp == g]
            m, sd = s.mean(), s.std(ddof=1)
            ax.fill_between(depths, m - sd, m + sd, color=G.color[g], alpha=0.2, lw=0)
            ax.plot(depths, m, color=G.color[g], lw=1.2, marker="o", ms=2.5, label=lab)
            for dp in depths:
                rows.append(dict(group=g, depth=dp, n=len(s), mean=m[dp], sd=sd[dp]))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_xlim(0, xmax)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(PS.READS_PER_UNIT))
        top = pd.DataFrame(rows).eval("mean + sd").max()
        PS.fixed_ticks(ax, 0, float(top))                              # 縦軸は目盛 5 個（規約・260929 に #2 にも適用）
        ax.set_xlabel("Sequencing depth (reads)")
        ax.set_ylabel(ylabel)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1000:g}k" if v else "0"))
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)
        pd.DataFrame(rows).to_csv(_csv_dir(d, name) / f"{name}.values.csv", index=False)

        last = (per[depths[-1]] - per[depths[-2]]) / per[depths[-2]] * 100
        log.sum(f"検体 {len(per)} (不変条件 OK)")
        log(f"  検体 {len(per)} = 全 {len(G.label_by_newid)} − 除外 {len(excl)} (不変条件 OK)/ 深さ {len(depths)} 段"
            f" (最大 {depths[-1]:,} = sampling depth)/ 最後の 2 段の増加率: 中央値 {last.median():.2f}%"
            f"・最大 {last.max():.2f}% ({last.idxmax()})")


# ============================================================ #3
ALPHA_METRICS = [(m, METRIC_LABEL[m]) for m in ("observed_features", "shannon", "faith_pd", "evenness")]


def _log_skipped(log, res: dict) -> None:
    """検定から外した群を注記に回す ([注記] は注記に集める印。最後の注記では印を外して出る。[注意] は使わない。261009 決定)。"""
    if res["skipped"]:
        n = {r["groupA"]: r["n_A"] for r in res["rows"]} | {r["groupB"]: r["n_B"] for r in res["rows"]}
        log(f"  [注記] n < {res['min_n']} のため検定から外した群があります: "
            + ", ".join(f"{g} (n = {n.get(g, 0)})" for g in res["skipped"]))
    if res.get("fallback"):
        log("  [注記] 群を外した結果 2 群だけ残ったので、Tukey の代わりに Welch の t 検定にしました")


def _log_stats(log, res: dict, st: pd.DataFrame, name: str) -> None:
    _log_skipped(log, res)
    sig = st[st["stars"] != ""]
    log.sum(f"有意 {name.split('_', 2)[2]} {len(sig)}/{int(st['p'].notna().sum())}")
    for _, r in sig.iterrows():                   # ログには有意なものだけ
        log(f"    {r.groupA} vs {r.groupB}: p = {r.p:.3g} {r.stars}")


def fig03_alpha(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#3 α多様性 (箱ひげ＋点・1 指標 1 枚・検定は stats.method)。
    値 = Step 07 の core-metrics (rarefied table・depth で揃えた後)の *_vector.tsv。
    rarefaction で除外された検体は、この時点で既に入っていない。"""
    log, out = ctx["log"], ctx["out"]
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["newID"])
    for k, (metric, ylabel) in enumerate(ALPHA_METRICS, 1):
        name = f"03_alpha_{metric}"
        log.step("10", n, total, False, f"#3 α多様性 ({k}/{len(ALPHA_METRICS)})", f"{name}.pdf")
        v = pd.read_csv(out / "07_diversity" / f"{metric}_vector.tsv", sep="\t", index_col=0).iloc[:, 0].dropna()
        unknown = sorted(set(v.index) - set(G.label_by_newid.index))
        if unknown:
            log(f"  [エラー] {metric}_vector.tsv にあってメタデータに無い検体があります: {unknown}")
            raise SystemExit(1)
        expect = len(G.label_by_newid) - len(excl)
        if len(v) != expect:                                          # 網羅の不変条件
            log(f"  [エラー] 検体数が合いません: {len(v)} ≠ 全検体 {len(G.label_by_newid)} − 除外 {len(excl)} = {expect}")
            raise SystemExit(1)
        grp = v.index.map(G.label_by_newid)
        order = [g for g in G.order if (grp == g).any()]   # 全検体が除外された群 (n = 0)は外す。1 ≤ n < 3 は残して「n < 3」(261020 決定)
        data = [v[grp == g].to_numpy(float) for g in order]
        labels = _shown_labels(G, order)          # 群名に n 数は併記しない（260917 指示・全図共通）

        if k == 1:
            _warn_labels(log, labels)
        fig, ax, xpos = PS.box_panel(labels, ylabel)
        PS.dots_mean_sd(ax, data, [G.color[g] for g in order])
        res = GS.compare(dict(zip(order, data)), order, G.method, G.control, G.min_n)
        lo, hi, step = _ylim_for(G, res, data)
        ax.set_ylim(lo, hi)
        hi = _draw_significance(ax, G, res, xpos, data, lo, hi, order)
        lo, hi = PS.fixed_ticks(ax, lo, hi)           # 目盛 5 個（260924 決定）
        PS.keep_brackets(ax)                         # 伸ばした軸でも、ブラケットは描いたときの見た目のまま (261020・F122)
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)

        pd.DataFrame({"sample": v.index, "group": grp, metric: v.values}).to_csv(_csv_dir(d, name) / f"{name}.values.csv", index=False)
        st = _write_stats(res, _csv_dir(d, name) / f"{name}.stats.csv",
                          {"metric": metric, "dunnett_seed": GS.DUNNETT_SEED if G.method == "dunnett" else ""})
        log.sum(f"検体 {len(v)} (不変条件 OK)")
        log(f"  検体 {len(v)} = 全 {len(G.label_by_newid)} − 除外 {len(excl)} (不変条件 OK)/ n: "
            + ", ".join(f"{g}={len(x)}" for g, x in zip(order, data)))
        _log_stats(log, res, st, name)


# ============================================================ #4
DISTANCES = [(m, METRIC_LABEL[m]) for m in ("jaccard", "bray_curtis", "unweighted_unifrac", "weighted_unifrac")]


def _read_pcoa(out: Path, dist: str):
    from skbio import OrdinationResults
    o = OrdinationResults.read(str(out / "07_diversity" / f"{dist}_pcoa_results.tsv"))
    return o.samples, (o.proportion_explained * 100).to_numpy()


def _scatter_panel(ax, G, x, y, grp, sub, markers, present_subs) -> list[dict]:
    """散布図 1 枚ぶん (原点の線・群ごとの楕円・点)。楕円の数値を返す (CSV とログに残す)。"""
    ax.axhline(0, **PS.ORIGIN_LINE)
    ax.axvline(0, **PS.ORIGIN_LINE)
    # 描く単位 = config の群（use_subgroup: true なら 処理×サブグループ）。色は G.color（260924 決定の 4 パターン）:
    #   show_group_subgroup: true  → 元の群の色＋サブグループの記号 / false → use_subgroup: true ならサブグループの色・丸、false なら群の色・丸
    ells = []
    for g in G.order:
        m = np.asarray(grp == g)
        if not m.any():                                # 全検体が除外された群 (n = 0)は外す (261020 決定)
            continue
        e = PS.data_ellipse(x[m], y[m])
        if e is not None:
            PS.draw_ellipse(ax, e, G.color[g])
            ells.append(dict(group=g, **e))
        else:
            ells.append(dict(group=g, n=int(m.sum())))
    for g in G.order:
        if G.show_subgroup:
            for s in present_subs:
                m = np.asarray((grp == g) & (sub == s))
                if m.any():
                    PS.scatter_subgroup(ax, x[m], y[m], G.color[g], G.subgroups[s], markers[s],
                                        s=PS.SCATTER_S, zorder=2)
        else:
            m = np.asarray(grp == g)
            ax.scatter(x[m], y[m], s=PS.SCATTER_S, color=G.color[g], lw=0.3, edgecolors="white", zorder=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return ells


def _load_scores(log, G, out: Path, dist: str, excl: set):
    sc, pct = _read_pcoa(out, dist)
    unknown = sorted(set(sc.index) - set(G.label_by_newid.index))
    expect = len(G.label_by_newid) - len(excl)
    if unknown or len(sc) != expect:                                  # 網羅の不変条件
        log(f"  [エラー] {dist}_pcoa_results.tsv の検体が合いません: {len(sc)} (期待 {expect})/ メタデータに無い: {unknown}")
        raise SystemExit(1)
    return sc, pct, sc.index.map(G.label_by_newid), sc.index.map(G.subgroup_by_newid)   # 群 = config の群（use_subgroup で変わる）


def _log_ellipses(log, ells: list[dict]) -> None:
    skipped = [f"{e['group']}(n={e['n']})" for e in ells if "width" not in e]
    if skipped:
        log(f"  [注記] n < 3 のため楕円を描かなかった群があります: {', '.join(skipped)}")


def fig04_pcoa(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#4 PCoA PCo1×PCo2 (4 距離 × 1 枚)。散布図の規約 (260921 確定)に従う:
    軸 48 × 48 mm・縦横の縮尺をそろえ**原点を中心**にする (範囲は点だけで決め、楕円ははみ出してよい)・目盛 9 pt / 軸ラベル 10 pt・原点の線 #DDDDDD 0.5 pt・点 14・
    **群ごとに 95% データ楕円** (χ² 近似・半径 2.45 SD で一定)・上余白に指標名のタイトル・凡例は別ファイル。
    色 = Group_name・記号 = サブグループ (show_subgroup: true)。検定は付けない。"""
    log, out = ctx["log"], ctx["out"]
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["newID"])
    markers = PS.subgroup_markers(G.subgroups) if G.show_subgroup else {}
    for k, (dist, title) in enumerate(DISTANCES, 1):
        name = f"04_pcoa_{dist}"
        log.step("10", n, total, False, f"#4 PCoA PCo1×PCo2 ({k}/{len(DISTANCES)})", f"{name}.pdf")
        sc, pct, grp, sub = _load_scores(log, G, out, dist, excl)
        present_subs = [s for s in G.subgroups if s in set(sub)] if G.show_subgroup else []
        x, y = sc.iloc[:, 0].to_numpy(), sc.iloc[:, 1].to_numpy()

        fig, ax = PS.fixed_axes(PS.SCATTER_AX_MM, PS.SCATTER_AX_MM)
        ells = _scatter_panel(ax, G, x, y, grp, sub, markers, present_subs)
        lim = max(np.abs(x).max(), np.abs(y).max()) * 1.08   # 原点中心・範囲は点だけで決める（楕円ははみ出してよい。260921 案 a）
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ticks = PS.symmetric_ticks(lim)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel(f"PCo1 ({pct[0]:.1f}%)")
        ax.set_ylabel(f"PCo2 ({pct[1]:.1f}%)")
        ax.set_title(title, fontsize=PS.TITLE_FS, loc="left", pad=2)
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)

        pd.DataFrame({"sample": sc.index, "group": grp, "subgroup": sub, "PCo1": x, "PCo2": y}
                     ).to_csv(_csv_dir(d, name) / f"{name}.values.csv", index=False)
        pd.DataFrame(ells).assign(method="data_ellipse_chi2", level=PS.ELLIPSE_LEVEL).to_csv(_csv_dir(d, name) / f"{name}.ellipses.csv",
                                                                                         index=False)
        log.sum(f"検体 {len(sc)} (不変条件 OK)")
        log(f"  検体 {len(sc)} (不変条件 OK)/ {title}: PCo1 {pct[0]:.1f}%・PCo2 {pct[1]:.1f}% / 範囲 ±{lim:.2f}"
            f" / 楕円 データ楕円 χ² {PS.ELLIPSE_LEVEL:.0%} ({len(ells)} 群ごと)")
        _log_ellipses(log, ells)


# ============================================================ #4
def fig05_pcoa_matrix(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#4 PCoA 交差形式 PCo1〜4 (4 距離 × 1 枚)。
    列 = 横軸 PCo1〜PCo3・行 = 縦軸 PCo2〜PCo4 の 3 × 3 の三角形。軸ラベル・目盛は左端の列と最下段のみ。
    散布図の規約は #4 と同じ (1 マス 48 mm・文字・点・原点の線・楕円)。**全マス同じ範囲・原点中心**
    (PCo1〜4 の座標の絶対値の最大値から毎回決める。楕円は範囲に含めない)。凡例は別ファイル。"""
    log, out = ctx["log"], ctx["out"]
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["newID"])
    markers = PS.subgroup_markers(G.subgroups) if G.show_subgroup else {}
    K, cell, gap = 4, float(PS.SCATTER_AX_MM), 4.0
    M = K - 1
    for k, (dist, title) in enumerate(DISTANCES, 1):
        name = f"04_pcoa-matrix_{dist}"
        log.step("10", n, total, False, f"#4 PCoA 交差形式 PCo1〜4 ({k}/{len(DISTANCES)})", f"{name}.pdf")
        sc, pct, grp, sub = _load_scores(log, G, out, dist, excl)
        present_subs = [s for s in G.subgroups if s in set(sub)] if G.show_subgroup else []

        grid_w = cell * M + gap * (M - 1)
        fig_w = PS.MARGIN_L + grid_w + PS.MARGIN_R
        fig_h = PS.MARGIN_B + grid_w + PS.MARGIN_T
        fig = matplotlib.pyplot.figure(figsize=(fig_w * PS.MM, fig_h * PS.MM))
        axes, all_ells, extents = [], [], [float(np.abs(sc.iloc[:, :K].to_numpy()).max())]
        for r in range(M):                        # 行 = PCo2..PCo4（縦軸）
            for c in range(r + 1):                # 列 = PCo1..PCo(r+1)（横軸）
                xi, yi = c, r + 1
                x0 = PS.MARGIN_L + c * (cell + gap)
                y0 = PS.MARGIN_B + (M - 1 - r) * (cell + gap)
                ax = fig.add_axes([x0 / fig_w, y0 / fig_h, cell / fig_w, cell / fig_h])
                ells = _scatter_panel(ax, G, sc.iloc[:, xi].to_numpy(), sc.iloc[:, yi].to_numpy(), grp, sub,
                                      markers, present_subs)
                all_ells += [dict(pair=f"PCo{xi + 1}xPCo{yi + 1}", **e) for e in ells]
                axes.append((ax, r, c, xi, yi))
        lim = max(extents) * 1.08                 # 範囲は点だけで決める（楕円ははみ出してよい。260921 案 a）
        ticks = PS.symmetric_ticks(lim)
        for ax, r, c, xi, yi in axes:
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            if c == 0:
                ax.set_ylabel(f"PCo{yi + 1} ({pct[yi]:.1f}%)")
            else:
                ax.set_yticklabels([])
            if r == M - 1:
                ax.set_xlabel(f"PCo{xi + 1} ({pct[xi]:.1f}%)")
            else:
                ax.set_xticklabels([])
            if r == 0 and c == 0:
                ax.set_title(title, fontsize=PS.TITLE_FS, loc="left", pad=2)
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)

        cols = {f"PCo{i + 1}": sc.iloc[:, i].values for i in range(K)}
        pd.DataFrame({"sample": sc.index, "group": grp, "subgroup": sub, **cols}).to_csv(_csv_dir(d, name) / f"{name}.values.csv",
                                                                                         index=False)
        pd.DataFrame(all_ells).assign(method="data_ellipse_chi2", level=PS.ELLIPSE_LEVEL).to_csv(
            _csv_dir(d, name) / f"{name}.ellipses.csv", index=False)
        log.sum(f"検体 {len(sc)} (不変条件 OK)")
        log(f"  検体 {len(sc)} (不変条件 OK)/ {title}: 寄与率 PCo1〜4 = "
            + " / ".join(f"{v:.1f}%" for v in pct[:K]) + f" (計 {pct[:K].sum():.1f}%)/ 共通範囲 ±{lim:.2f}")
        _log_ellipses(log, [e for e in all_ells if e["pair"] == "PCo1xPCo2"])


# ============================================================ #5
NPERM = 999
PERM_SEED = 0          # 置換検定の乱数の種（固定。同じ入力から必ず同じ q が出るように。260917 決定）
CELL_MM = PS.SCATTER_AX_MM   # 行列の 1 マス（交差形式の PCoA #4 と同じ 48 mm。260921 決定・260923 に 46 → 48）
CELL_GAP_MM = 4        # マスのあいだ（#4 と同じ 4 mm）


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def _qfmt(q: float) -> str:
    return "q ≤ 0.001" if q <= 1 / (NPERM + 1) + 1e-12 else f"q = {q:.3f}"


def fig06_permanova(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#5 PERMANOVA / PERMDISP 行列 (4 距離 × 1 枚。260917 決定・260921 に新しいレイアウト規則を適用)。
    ペアごとに PERMANOVA (位置の差)と PERMDISP (ばらつきの差)を scikit-bio で計算する (案 P)。
    置換 999 回・乱数の種を固定・BH 補正は検定ごとに全ペアの中で行う。R² = F / (F + (N − 2))。
    図: 左下の三角形 (群数 − 1) × (群数 − 1)。色 = PERMANOVA の R² (viridis の 0〜95%・0〜1 固定)、マスの中 =
    PERMANOVA の q と PERMDISP の q。星は使わない。図の中に注記は書かない (CSV とログに残す)。"""
    import itertools
    import time
    import warnings
    import skbio
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle
    from skbio import DistanceMatrix
    from skbio.stats.distance import permanova, permdisp
    log, out = ctx["log"], ctx["out"]
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["newID"])
    left = set(G.label_by_newid[[s for s in G.label_by_newid.index if s not in excl]])
    order = [g for g in G.order if g in left]          # 全検体が除外された群は組を作らない (261019・F234)
    Mg = len(order) - 1
    if Mg < 1:
        log.step("10", n, total, True, "#5 PERMANOVA / PERMDISP 行列", "群が 1 つしかないため比較できない")
        return
    cmap, norm = PS.viridis(), Normalize(0, 1)            # viridis の 0〜95%（全図共通・261009）
    idx = {g: i for i, g in enumerate(order)}
    for k, (dist, title) in enumerate(DISTANCES, 1):
        name = f"05_permanova_{dist}"
        log.step("10", n, total, False, f"#5 PERMANOVA / PERMDISP 行列 ({k}/{len(DISTANCES)})", f"{name}.pdf")
        if k == 1 and G.use_subgroup:                  # 見出しの後に出す（注記が #5 の下にまとまるように。260928）
            npair = len(G.order) * (len(G.order) - 1) // 2
            log(f"  [注記] use_subgroup: true ({len(G.order)} 群・{npair} ペア) のため、計算に時間がかかり、図も大きくなります")
        t0 = time.time()
        D = pd.read_csv(out / "07_diversity" / f"{dist}_distance_matrix.tsv", sep="\t", index_col=0)
        expect = len(G.label_by_newid) - len(excl)
        if len(D) != expect or set(D.index) - set(G.label_by_newid.index) or list(D.index) != list(D.columns):
            log(f"  [エラー] {dist}_distance_matrix.tsv の検体が合いません: {len(D)} (期待 {expect})")
            raise SystemExit(1)
        dm = DistanceMatrix(np.ascontiguousarray(D.to_numpy(dtype=float)), list(D.index))   # scikit-bio は C 連続の配列が必要
        grp = D.index.map(G.label_by_newid)
        rows = []
        for a, b in itertools.combinations(order, 2):
            keep = [s for s, g in zip(D.index, grp) if g in (a, b)]
            labs = [g for g in grp if g in (a, b)]
            na, nb = labs.count(a), labs.count(b)
            # 検定しない組（261009 決定・D2d）: n < stats.min_n の群を含む組は計算しない。
            # 並べ替えで取りうる最小の p（= 1 / 区別できる分け方の数。両群が同じ数なら分け方は半分）が 0.05 を超える組は、
            # どれだけ差があっても有意になれないので「not testable」として q の計算から外す（R² は色として描く）
            if na < G.min_n or nb < G.min_n:
                rows.append(dict(groupA=a, groupB=b, n_A=na, n_B=nb, permanova_F=np.nan, R2=np.nan, permanova_p=np.nan,
                                 permdisp_F=np.nan, permdisp_p=np.nan, min_p=np.nan, testable=False,
                                 reason=f"n < {G.min_n}"))
                continue
            ways = math.comb(na + nb, na) // (2 if na == nb else 1)
            min_p = 1.0 / ways
            sub = dm.filter(keep)
            with warnings.catch_warnings():        # scikit-bio の内部 PCoA の警告 (負の固有値・全次元の計算) は出さない。以前と同じ (261021 決定・F048 案 b)
                warnings.simplefilter("ignore")
                pa = permanova(sub, labs, permutations=NPERM, seed=PERM_SEED)
                pdp = permdisp(sub, labs, permutations=NPERM, seed=PERM_SEED)
            F = float(pa["test statistic"])
            ok = min_p <= 0.05
            rows.append(dict(groupA=a, groupB=b, n_A=na, n_B=nb,
                             permanova_F=F, R2=F / (F + (len(keep) - 2)), permanova_p=float(pa["p-value"]),
                             permdisp_F=float(pdp["test statistic"]), permdisp_p=float(pdp["p-value"]),
                             min_p=min_p, testable=ok, reason="" if ok else f"min p = {min_p:.2f} > 0.05"))
        T = pd.DataFrame(rows)
        T["permanova_q"] = np.nan
        T["permdisp_q"] = np.nan
        if T.testable.any():                           # BH は検定できる組の中だけで行う
            T.loc[T.testable, "permanova_q"] = _bh(T.loc[T.testable, "permanova_p"])
            T.loc[T.testable, "permdisp_q"] = _bh(T.loc[T.testable, "permdisp_p"])
        T = T.assign(distance=dist, permutations=NPERM, seed=PERM_SEED, p_adjust="BH (検定ごとに全ペアで)",
                     permdisp_test="median (scikit-bio 既定)", scikit_bio=skbio.__version__)

        # ---- 図: (群数 − 1) × (群数 − 1) の三角形。群名は左（行）と下（列）、どちらも 45° ----
        # 座標は mm（1 マス CELL_MM・マスのあいだ CELL_GAP_MM。交差形式の PCoA #4 と同じ並べ方）
        step_mm = CELL_MM + CELL_GAP_MM
        grid = CELL_MM * Mg + CELL_GAP_MM * (Mg - 1)
        # カラーバーは一番上のマスの右隣（1 マス分）。2 群（マス 1 つ）のときは右隣がマスの外なので図を広げる（260921）
        extra_w = step_mm if Mg == 1 else 0.0
        fig_w = PS.MARGIN_L + grid + extra_w + PS.MARGIN_R
        fig_h = PS.MARGIN_B + grid + PS.MARGIN_T
        fig = matplotlib.pyplot.figure(figsize=(fig_w * PS.MM, fig_h * PS.MM))
        ax = fig.add_axes([PS.MARGIN_L / fig_w, PS.MARGIN_B / fig_h, grid / fig_w, grid / fig_h])
        for _, r in T.iterrows():
            c, rr = idx[r.groupA], idx[r.groupB] - 1          # 列 = 前の群・行 = 後の群
            if np.isnan(r.R2):                          # n < min_n の群を含む組: 欠損の灰色（規約 #DDDDDD）に「n < 3」
                face, tc, txt = "#DDDDDD", PS.GREY, r.reason
            else:
                v = norm(r.R2)
                face, tc = cmap(v), "white"                # マスの中の文字は白（全図共通・261009）
                txt = (f"PERMANOVA\n{_qfmt(r.permanova_q)}\n\nPERMDISP\n{_qfmt(r.permdisp_q)}" if r.testable
                       else f"not testable\n(min p = {r.min_p:.2f})")      # 図の中は英語（261009 決定・案 a）
            ax.add_patch(Rectangle((c * step_mm, rr * step_mm), CELL_MM, CELL_MM, facecolor=face, lw=0))
            ax.text(c * step_mm + CELL_MM / 2, rr * step_mm + CELL_MM / 2, txt,
                    ha="center", va="center", fontsize=PS.FS_CELL, color=tc, linespacing=1.15)
        ax.set_xlim(0, grid)
        ax.set_ylim(grid, 0)
        ax.set_aspect("equal")
        for sp in ax.spines.values():
            sp.set_visible(False)
        centers = [k * step_mm + CELL_MM / 2 for k in range(Mg)]
        ax.set_xticks(centers)
        ax.set_xticklabels(order[:-1], fontsize=PS.XLABEL_FS, rotation=PS.XLABEL_ROT, ha="right",
                           rotation_mode="anchor")
        ax.set_yticks(centers)
        ax.set_yticklabels(order[1:], fontsize=PS.XLABEL_FS, rotation=PS.XLABEL_ROT, ha="right", va="center",
                           rotation_mode="anchor")
        ax.tick_params(length=0)
        _warn_labels(log, list(order[:-1]))
        over = [(g, h) for g in order[1:] if (h := PS.rotated_label_mm(g) + 1.5) > PS.MARGIN_L]
        if over:
            log("  [注意] 左の群名が長く、左余白 15 mm で見切れている可能性があります (群名を短くしてください): "
                + ", ".join(f"{g} (約 {h:.1f} mm)" for g, h in over))
        ax.set_title(title, fontsize=PS.TITLE_FS, loc="left", pad=2)
        # カラーバーは一番上のマスの右隣に、横向き・1 マス分の幅で置く（縦はそのマスの中央。260921 決定）
        cb_w = CELL_MM
        cb_x = PS.MARGIN_L + step_mm
        cb_y = PS.MARGIN_B + grid - CELL_MM / 2 - 1.5
        cax = fig.add_axes([cb_x / fig_w, cb_y / fig_h, cb_w / fig_w, 3.0 / fig_h])
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
        cax.xaxis.set_label_position("top")
        cb.outline.set_visible(False)
        cb.set_label("PERMANOVA R²", fontsize=PS.FS_CBAR)
        cb.ax.tick_params(labelsize=PS.FS_CBAR, length=2)
        cb.set_ticks([0, 0.5, 1])
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)
        T.to_csv(_csv_dir(d, name) / f"{name}.stats.csv", index=False)

        n_both = int(((T.permanova_q < 0.05) & (T.permdisp_q >= 0.05)).sum())
        n_small, n_untest = int(T.R2.isna().sum()), int((~T.testable & T.R2.notna()).sum())
        if k == 1 and n_small:
            small = [g for g in order if (grp == g).sum() < G.min_n]
            log(f"  [注記] n < {G.min_n} のため検定から外した群があります: " + ", ".join(f"{g} (n = {int((grp == g).sum())})" for g in small))
        if k == 1 and n_untest:
            log(f"  [注記] 並べ替えで取りうる最小の p が 0.05 を超えるため検定できない組が {n_untest} ペアあります (図では not testable と表示します)")
        log.sum(f"{title} q < 0.05 {int((T.permanova_q < 0.05).sum())}/{int(T.testable.sum())}")
        log(f"  {title}: {len(T)} ペア (検定できた {int(T.testable.sum())}・n 不足 {n_small}・not testable {n_untest})"
            f" / R² {T.R2.min():.2f}〜{T.R2.max():.2f} / PERMANOVA q < 0.05: "
            f"{int((T.permanova_q < 0.05).sum())} / PERMDISP q < 0.05: {int((T.permdisp_q < 0.05).sum())}"
            f" (PERMANOVA 有意で PERMDISP では差を検出せず = ばらつきが等しいことの証明ではない: {n_both})/ {time.time() - t0:.0f} 秒")
    log(f"  置換 {NPERM} 回・種 {PERM_SEED}・BH (検定ごと)・scikit-bio {skbio.__version__} (permanova / permdisp)")


LEAF_MM = 4            # 1 検体 = 4 mm（#6 検体版・#7 共通。260922 決定: 3 mm → 4 mm。群 8 mm = 検体 2 本分・記号に余裕 (260922 時点は番号も)。X01 で約 520 mm）
BAND_MM, BAND_GAP_MM = 3.0, 0.5   # 葉の下: 群の色の帯 → サブグループの記号 → 検体名 (45°)
# 葉の下の文字は検体名 (newID) を 45° (261022 決定。それまでは群の中の番号を 0°)。下の余白は広げず、長い名前は見切れてよい


def fig07_dendrogram(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#6 デンドログラム (4 距離 × 1 枚。260921 決定)。
    UPGMA (average)で距離行列をそのまま木にする (Ward はユークリッド距離が前提なので使わない)。
    縦向き・葉は下・葉 1 つ 4 mm・縦軸は距離の目盛を出す。cophenetic 相関は出さない (決定)。
    葉の下: 群の色の帯 → サブグループの記号 (show_subgroup: true のとき。#4 と同じ記号) → 検体名 45° (261022。見切れてよい)。
    凡例は既存の 00_legend_group-subgroup.pdf を使う。"""
    from matplotlib.patches import Rectangle
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform
    log, out = ctx["log"], ctx["out"]
    markers = PS.subgroup_markers(G.subgroups) if G.show_subgroup else {}
    for k, (dist, title) in enumerate(DISTANCES, 1):
        name = f"06_dendrogram_{dist}"
        log.step("10", n, total, False, f"#6 デンドログラム ({k}/{len(DISTANCES)})", f"{name}.pdf")
        D = pd.read_csv(out / "07_diversity" / f"{dist}_distance_matrix.tsv", sep="\t", index_col=0)
        Z = linkage(squareform(D.to_numpy(dtype=float), checks=False), method="average")
        leaves = list(D.index[dendrogram(Z, no_plot=True)["leaves"]])
        N = len(leaves)

        ax_w = N * LEAF_MM
        fig_w = PS.MARGIN_L + ax_w + PS.MARGIN_R
        fig_h = PS.MARGIN_B + PS.PLOT_H_MM + PS.MARGIN_T
        fig = matplotlib.pyplot.figure(figsize=(fig_w * PS.MM, fig_h * PS.MM))
        ax = fig.add_axes([PS.MARGIN_L / fig_w, PS.MARGIN_B / fig_h, ax_w / fig_w, PS.PLOT_H_MM / fig_h])
        dendrogram(Z, ax=ax, no_labels=True, color_threshold=0, above_threshold_color=PS.GREY,
                   link_color_func=lambda _: PS.GREY)
        for ln in ax.collections:
            ln.set_linewidth(0.6)
        for sp in ("top", "right", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.set_xlim(0, N * 10)
        PS.fixed_ticks(ax, 0, float(Z[:, 2].max()))   # 目盛 5 個（260924 決定）
        ax.set_xticks([])
        ax.set_ylabel("Distance (UPGMA)")
        ax.set_title(title, fontsize=PS.TITLE_FS, loc="left", pad=2)

        # 下余白 20 mm の中に mm 座標の軸を重ねて、帯・記号・検体名を置く
        bx = fig.add_axes([PS.MARGIN_L / fig_w, 0, ax_w / fig_w, PS.MARGIN_B / fig_h])
        bx.set_xlim(0, ax_w)
        bx.set_ylim(-PS.MARGIN_B, 0)
        bx.axis("off")
        y = -BAND_GAP_MM - BAND_MM
        for i, s in enumerate(leaves):
            bx.add_patch(Rectangle((i * LEAF_MM, y), LEAF_MM, BAND_MM, facecolor=G.color[G.label_by_newid[s]], lw=0))
        if markers:
            y -= 2.0
            for i, s in enumerate(leaves):
                sub = G.subgroup_by_newid[s]
                PS.scatter_subgroup(bx, [i * LEAF_MM + LEAF_MM / 2], [y], "#666666", G.subgroups.get(sub, {}),
                                    markers[sub], s=PS.SCATTER_S, clip_on=False)
            y -= 2.0
        else:
            y -= 1.0
        for i, s in enumerate(leaves):
            bx.text(i * LEAF_MM + LEAF_MM / 2, y, s, fontsize=PS.XLABEL_FS, rotation=PS.XLABEL_ROT, rotation_mode="anchor",
                    ha="right", va="top")   # 検体名 45° (261022 決定。見切れてよい)
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)

        cd = _csv_dir(d, name)
        pd.DataFrame({"leaf_order": range(1, N + 1), "sample": leaves,
                      "group": [G.group_by_newid[s] for s in leaves],
                      "subgroup": [G.subgroup_by_newid[s] for s in leaves],
                      "label": [G.label_by_newid[s] for s in leaves]}).to_csv(cd / f"{name}.leaves.csv", index=False)   # label = 図の色帯の単位 (261020・F252)
        pd.DataFrame(Z, columns=["child1", "child2", "height", "n_samples"]).assign(
            method="average (UPGMA)").to_csv(cd / f"{name}.linkage.csv", index=False)
        log.sum(f"検体 {N}")
    log("  葉の下 = 検体名 (45°)。検体の並びと群は dendrogram/*.leaves.csv")


def fig07_dendrogram_group(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#6 デンドログラム (群の平均・4 距離 × 1 枚。260921 追加)。検体版とは別に出す。
    群間の距離 = 群の中心どうしの距離 (PCoA の空間での重心間の距離。PERMANOVA と同じ考え方) → UPGMA (261009 決定・D7)。
      d_AB² = mean(d²_ij, i∈A, j∈B) − Σ_{i,i'∈A} d²_ii' / (2 n_A²) − Σ_{j,j'∈B} d²_jj' / (2 n_B²)。負になれば 0 (非ユークリッド距離)。
      それまでの「群 A の検体 × 群 B の検体の距離の平均」は群内のばらつきが足し込まれ、ばらつきの小さい群ほど近く見えた (点検 F064)。
    群の単位は config に従う (use_subgroup: false → Group_name、true → 処理 × サブグループ)。
    見た目: 1 群 8 mm・葉の真下に群の色の帯・その下に名前 45° (黒。use_subgroup: true ならサブグループ名だけ)。260921 指示。"""
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform
    log, out = ctx["log"], ctx["out"]
    for k, (dist, title) in enumerate(DISTANCES, 1):
        name = f"06_dendrogram-group_{dist}"
        log.step("10", n, total, False, f"#6 デンドログラム・群の平均 ({k}/{len(DISTANCES)})", f"{name}.pdf")
        D = pd.read_csv(out / "07_diversity" / f"{dist}_distance_matrix.tsv", sep="\t", index_col=0)
        grp = D.index.map(G.label_by_newid)
        order = [g for g in G.order if g in set(grp)]
        if len(order) < 3:
            log(f"  [注記] 群が {len(order)} つしかないため描きません (木にならないため)")
            return
        X = D.to_numpy(dtype=float)
        g_arr = np.asarray(grp)
        M = pd.DataFrame(0.0, index=order, columns=order)
        n_neg = 0
        for a in order:
            for b in order:
                if a != b:
                    ia, ib = g_arr == a, g_arr == b
                    d2 = ((X[np.ix_(ia, ib)] ** 2).mean() - (X[np.ix_(ia, ia)] ** 2).sum() / (2 * ia.sum() ** 2)
                          - (X[np.ix_(ib, ib)] ** 2).sum() / (2 * ib.sum() ** 2))
                    n_neg += d2 < 0
                    M.loc[a, b] = np.sqrt(max(d2, 0.0))
        within = {g: (float(X[np.ix_(g_arr == g, g_arr == g)][np.triu_indices(int((g_arr == g).sum()), 1)].mean())
                      if (g_arr == g).sum() > 1 else np.nan) for g in order}
        if n_neg:
            log(f"  [注記] {title}: 中心どうしの距離の 2 乗が負になった組 {n_neg // 2} 個を 0 としました (非ユークリッド距離のため)")
        Z = linkage(squareform(M.to_numpy(), checks=False), method="average")
        leaves = [order[i] for i in dendrogram(Z, no_plot=True)["leaves"]]
        N = len(leaves)

        ax_w = N * PS.MM_PER_GROUP
        fig_w = PS.MARGIN_L + ax_w + PS.MARGIN_R
        fig_h = PS.MARGIN_B + PS.PLOT_H_MM + PS.MARGIN_T
        fig = matplotlib.pyplot.figure(figsize=(fig_w * PS.MM, fig_h * PS.MM))
        ax = fig.add_axes([PS.MARGIN_L / fig_w, PS.MARGIN_B / fig_h, ax_w / fig_w, PS.PLOT_H_MM / fig_h])
        dendrogram(Z, ax=ax, no_labels=True, color_threshold=0, above_threshold_color=PS.GREY,
                   link_color_func=lambda _: PS.GREY)
        for ln in ax.collections:
            ln.set_linewidth(0.6)
        for sp in ("top", "right", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.set_xlim(0, N * 10)
        PS.fixed_ticks(ax, 0, float(Z[:, 2].max()))   # 目盛 5 個（260924 決定）
        ax.set_xticks([])
        ax.set_ylabel("Centroid distance (UPGMA)")
        ax.set_title(title, fontsize=PS.TITLE_FS, loc="left", pad=2)
        # 葉の真下に群の色の帯、その下に名前 45°（黒）。use_subgroup: true のときは名前 = サブグループ名だけ（260921 指示）
        from matplotlib.patches import Rectangle
        shown = [g.split(" ", 1)[1] for g in leaves] if (G.use_subgroup and G.show_subgroup) else leaves
        W = PS.MM_PER_GROUP
        bx = fig.add_axes([PS.MARGIN_L / fig_w, 0, ax_w / fig_w, PS.MARGIN_B / fig_h])
        bx.set_xlim(0, ax_w)
        bx.set_ylim(-PS.MARGIN_B, 0)
        bx.axis("off")
        y = -BAND_GAP_MM - BAND_MM
        for i, g in enumerate(leaves if not G.no_color else []):     # 群の色なしでは帯を描かない（261013）
            bx.add_patch(Rectangle((i * W, y), W, BAND_MM, facecolor=G.color[g], lw=0))
        for i, lab in enumerate(shown):
            bx.text(i * W + W / 2 + 1.0, y - 1.0, lab, fontsize=PS.XLABEL_FS, rotation=PS.XLABEL_ROT,
                    ha="right", va="top", rotation_mode="anchor")
        over = [(x, h) for x in shown if (h := PS.rotated_label_mm(x) + BAND_GAP_MM + BAND_MM + 1.0) > PS.MARGIN_B]
        if over:
            log(f"  [注意] 群名が長く、下余白 {PS.MARGIN_B} mm で見切れている可能性があります (群名を短くしてください): "
                + ", ".join(f"{x} (約 {h:.1f} mm)" for x, h in over))
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)

        cd = _csv_dir(d, name)
        n_by = grp.value_counts()
        M.to_csv(cd / f"{name}.group-distance.csv")
        pd.DataFrame({"leaf_order": range(1, N + 1), "group": leaves,
                      "n": [int(n_by[g]) for g in leaves],
                      "mean_within_group_distance": [within[g] for g in leaves]}).to_csv(cd / f"{name}.leaves.csv", index=False)
        pd.DataFrame(Z, columns=["child1", "child2", "height", "n_groups"]).assign(
            method="average (UPGMA)", group_distance="distance between group centroids (PCoA space)").to_csv(
            cd / f"{name}.linkage.csv", index=False)
        log.sum(f"群 {N}")


# ---------------- #7 積み上げ棒（260921 決定） ----------------
# 色 = Paul Tol の discrete rainbow の「載せるタクソン数 N 色用の組」を逆順（下の層ほど赤）＋ Others 灰
# （261016 決定・E6b。それまでは 20 色用の組の先頭 N 色で、N = 12 では暗い赤が続き下の層が見分けにくかった。点検 F123・F279）。
# discrete rainbow は 23 色までなので、24 以上は smooth_rainbow から等間隔に取る（261021 決定・案 ア。隣どうしの色は近くなる）
TAXA_MAX = 23
TITLE_MAX_LINES = 3     # #10 のタイトルは 3 行まで (上余白 13 mm。260923)


def _taxa_colors(n: int) -> list[str]:
    """#7 の属の色 (下の層ほど赤)＋ Others の灰。"""
    return PC.palette_colors("discrete_rainbow" if n <= TAXA_MAX else "smooth_rainbow", n)[::-1] + [OTHERS_COLOR]

OTHERS_COLOR = "#BBBBBB"
BAR_FRAC = 0.8          # 棒の幅（1 群 / 1 検体の幅に対する割合。棒のあいだに隙間を空ける）
BAR_ALPHA = 0.6         # #10 の棒の塗りの濃さ（群の色を薄くして、棒の中の点を見せる。260924 に 0.3・260928 決定で 0.6。規約 Skills/figure-plotting と同じ）
# #10 の縦軸の上端は決まった段から選ぶ（260924 決定）。1・2・4・6・8 の階段 (桁共通) で下限 0.4%（261020 決定: 1・2・4 から 6・8 を足し、4 → 10 などの隙間を埋めた）。
# 4 等分した目盛がどの段でも読みやすい数になる（1・2・5 の階段では 0.5 が 0.125 刻みになった）。X01 の属で軸は 13 種類 (261020)
TAXON_BAR_TOPS = (0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100)
# アーカイブ plot_02 の名前の整形（normalize_taxon_header）をそのまま移植。SILVA 144 の k__ にも対応させた
UNCULTURED_TERMS = {"uncultured", "uncultured bacterium", "metagenome", "gut metagenome", "soil metagenome",
                    "marine metagenome", "incertae_sedis"}


def _extract_label(part: str) -> str:
    import re
    p = (part or "").strip()
    if not p or p in {"__", "_", "-", "NA", "N/A"}:
        return ""
    if ":" in p:
        p = p.split(":", 1)[1].strip()
    m = re.match(r"(?:[Dd]_\d+|[dkpcfogsDKPCFOGS])__?(.*)$", p)      # アーカイブに k を足した（SILVA 144）
    if m:
        return m.group(1).strip()
    return "Unassigned" if p.lower() == "unassigned" else p


def _uncultured(v: str) -> bool:
    v = v.strip().lower().replace("-", "_")          # SILVA 144 のハイフン綴り Incertae-Sedis も印の語にする (261019・F012・F031)
    return bool(v) and (v in UNCULTURED_TERMS or v.endswith(" metagenome") or v == "metagenome")


def taxon_name_archive(s: str, max_level: int) -> str:
    """アーカイブと同じ規則で名前を付ける: 末端が有効ならその名前 / 空・uncultured・親と同名なら
    「有効な祖先;__子;__子」/ 先頭が Unassigned なら「Unassigned;__…」。"""
    parts = [p.strip() for p in s.split(";")]
    vals = [_extract_label(parts[i]) if i < len(parts) else "" for i in range(max_level + 1)]
    if vals and vals[0].lower() == "unassigned":
        return "Unassigned;" + ";".join(["__"] * max_level) if max_level > 0 else "Unassigned"
    for k in range(1, max_level + 1):
        if vals[k] == "":                            # 途中が空: 上へさかのぼり、印の語も飛ばして有効な祖先を探す (261019・F025)
            j = k - 1
            while j > 0 and (vals[j] == "" or _uncultured(vals[j])):
                j -= 1
            return f"{vals[j]};" + ";".join(["__"] * (max_level - j))
    target, parent = vals[max_level], (vals[max_level - 1] if max_level > 0 else "")
    if target and not _uncultured(target) and target != parent:
        return target
    j = max_level - 1
    while j >= 0 and (vals[j] == "" or _uncultured(vals[j]) or vals[j] == vals[j + 1]):
        j -= 1
    anchor = vals[j] if j >= 0 else (vals[0] if vals else "")
    # 印の語（uncultured・Incertae_Sedis など）は名前として残さず "__" にする（260923 決定・案 A。
    # アーカイブは小文字にそろえずに比べていたため "__Incertae_Sedis" と表示されていた）
    tail = ["__" if (not vals[t] or _uncultured(vals[t])) else "__" + vals[t] for t in range(j + 1, max_level + 1)]
    return anchor + ";" + ";".join(tail) if anchor else ";".join(tail)


def _unique_names(raw: list[str], names: list[str]) -> list[str]:
    """同じ名前になった行は足さず、名前で区別できるようにする (260923 決定)。
    区別のしかた: 名前の先頭 (基準の階級) より上の階級を、近い順に 1 つずつ前に足す (例 Lactobacillales;Enterococcus--other)。
    空の階級と印の語は足さない (261019・F010・F013・F032。以前は基準より下の階級や Incertae_Sedis を足していた)。
    それでも同じなら末尾に (2)(3) を付ける (例 属名が空の行と Incertae_Sedis の行。番号は表の行の順)。"""
    import collections
    out = list(names)
    dup = {n for n, c in collections.Counter(names).items() if c > 1}
    for n in dup:
        idx = [i for i, x in enumerate(names) if x == n]
        head = n.split(";")[0]
        ancs = []
        for i in idx:
            v = [_extract_label(p) for p in raw[i].split(";")]
            j = max((k for k, x in enumerate(v) if x == head), default=len(v))
            ancs.append([x for x in reversed(v[:j]) if x and not _uncultured(x)])
        for up in range(max(len(a) for a in ancs)):              # 上の階級を 1 つずつ足す
            cand = [f"{a[up]};{n}" if len(a) > up else n for a in ancs]
            if len(set(cand)) == len(cand):
                for i, c in zip(idx, cand):
                    out[i] = c
                break
        else:
            for k, i in enumerate(idx, 1):                     # 系統が同じなら番号で分ける
                out[i] = n if k == 1 else f"{n} ({k})"
    return out


def _display_names(out: Path, rank: str):
    """collapsed_table_<rank>.tsv と表示名 (アーカイブと同じ整形・同名は上の階級や (2) で区別)。
    返り値: (表, 整形しただけの名前, 区別した表示名)。#7〜#12 と対応表 (00_info/taxa_names_<rank>.csv) で共通。"""
    t = pd.read_csv(out / "09_taxonomy" / f"collapsed_table_{rank}.tsv", sep="\t", skiprows=1, index_col=0)
    max_level = max(len(s.split(";")) for s in t.index) - 1
    names = [taxon_name_archive(s, max_level) for s in t.index]
    return t, names, _unique_names(list(t.index), names)


def _write_sample_information(ctx: dict, info: Path) -> int:
    """01 の sample_information.csv (検体情報 CSV の全列) に、最終のリード数 (Reads = キメラ除去後。#1 の縦軸と同じ) を
    Subgroup_name のすぐ後ろへ足して 00_info/ に書く (261024 決定)。CSV は読み直さない (01 と食い違わないように)。"""
    out = ctx["out"]
    t = pd.read_csv(out / "01_metadata" / "sample_information.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)
    reads = read_denoising_stats(merged(ctx)["stats"]).set_index("sample-id")["non-chimeric"]
    t.insert(t.columns.get_loc("Subgroup_name") + 1, "Reads", t["#SampleID"].map(reads).astype(int))
    t.to_csv(info / "sample_information.csv", index=False, encoding="utf-8-sig")
    return len(t)


def _write_taxa_names(ctx: dict, info: Path) -> None:
    """表示名と元の SILVA 系統の対応表を階級ごとに 1 つ出す (261021 決定・F033/F086 案 a。(2) の行がどれか分かるように)。"""
    _, ranks = _plot_cfg(ctx)
    for rank in ranks:
        t, _, uniq = _display_names(ctx["out"], rank)
        pd.DataFrame({"display_name": uniq, "silva_lineage": list(t.index)}).to_csv(info / f"taxa_names_{rank}.csv", index=False)


def _taxa_relative(ctx: dict, rank: str):
    """collapsed_table_<rank>.tsv → 検体ごとの相対存在量 (%)。名前はアーカイブと同じ整形 (260921 決定)。
    同じ名前になった行は足さず、上の階級を足して区別する。
    返り値: (相対存在量 [タクソン × 検体], 除外検体の集合, 同名で足した行数)"""
    out = ctx["out"]
    t, names, uniq = _display_names(out, rank)      # 行は足さず、同じ名前になったら上の階級を足して区別する
    n_fixed = sum(a != b for a, b in zip(names, uniq))
    t = t.set_axis(pd.Index(uniq))
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["newID"])
    keep = [c for c in t.columns if c not in excl]
    ra = t[keep] / t[keep].sum() * 100
    return ra, excl, n_fixed


def _top_taxa(ra: pd.DataFrame, G, top_n: int) -> tuple[pd.DataFrame, list[str]]:
    """上位 N = 群ごとの平均を出してから、その合計の大きい順 (260922 決定・アーカイブ plot_03 と同じ)。
    群を同じ重みで扱う (検体数の多い群に引っ張られない)。群 = config の群。#7・#8 で共通。"""
    grp = ra.columns.map(G.label_by_newid)
    gm = pd.DataFrame({g: ra.loc[:, grp == g].mean(axis=1) for g in G.order if (grp == g).any()})
    top = gm.sum(axis=1).sort_values(ascending=False).head(top_n).index.tolist()
    st = ra.loc[top].copy()
    st.loc["Others"] = ra.drop(index=top).sum()
    return st, top


def _stack(ax, x, st: pd.DataFrame, width: float, colors: list[str]) -> None:
    """1 位を一番下、Others を一番上に積む。"""
    bottom = np.zeros(st.shape[1])
    for (name, row), c in zip(st.iterrows(), colors):
        ax.bar(x, row.to_numpy(), width, bottom=bottom, color=c, lw=0)
        bottom += row.to_numpy()


def _taxa_legend(path: Path, names: list[str], colors: list[str], log) -> None:
    """凡例は別ファイル (高さ 81 mm・横幅は実測)。上から Others → 下位 → 1 位 (棒の積み順と同じ向き)。
    1 列で高さに収まらなければ 2 列にする。"""
    from matplotlib.patches import Patch
    fig_h = PS.MARGIN_B + PS.PLOT_H_MM + PS.MARGIN_T
    handles = [Patch(facecolor=c, label=n) for n, c in zip(names, colors)][::-1]
    for ncol in (1, 2, 3):
        fig = matplotlib.pyplot.figure(figsize=(300 * PS.MM, fig_h * PS.MM))
        lg = fig.legend(handles=handles, fontsize=PS.FS_LEGEND, frameon=False, loc="upper left", ncol=ncol,
                        bbox_to_anchor=(2 / 300, 1 - PS.MARGIN_T / fig_h), borderaxespad=0, borderpad=0,
                        handletextpad=0.4, labelspacing=0.4, columnspacing=1.0, handlelength=1.2)
        fig.canvas.draw()
        bb = lg.get_window_extent()
        h_mm, w_mm = bb.height / fig.dpi * 25.4, bb.width / fig.dpi * 25.4
        if PS.MARGIN_T + h_mm <= fig_h - 2 or ncol == 3:
            break
        matplotlib.pyplot.close(fig)
    fig_w = 2 + w_mm + 2
    fig.set_size_inches(fig_w * PS.MM, fig_h * PS.MM)
    lg.set_bbox_to_anchor((2 / fig_w, 1 - PS.MARGIN_T / fig_h), transform=fig.transFigure)
    fig.savefig(path)
    matplotlib.pyplot.close(fig)
    log(f"  凡例     : {path.name} ({ncol} 列・幅 {fig_w:.0f} mm)")


def _plot_cfg(ctx: dict) -> tuple[int, list[str]]:
    """plotting の値 (検査は PC.check_plotting。01 の前にも同じ検査を通している)。"""
    PC.check_plotting(ctx["cfg"])
    p = ctx["cfg"]["plotting"]
    return p["top_n_bacteria"], [str(r).strip().lower() for r in p["ranks"]]


def _top_n_or_all(ctx: dict, k: str) -> float:
    """ancombc_top_n・ancombc_summary_top_n。all は上限なし (math.inf。261021 決定で #11 も all を可に)。"""
    v = ctx["cfg"]["plotting"][k]
    return math.inf if isinstance(v, str) and v.strip().lower() == "all" else v


def _bottom_band(fig, fig_w, fig_h, ax_w, pitch, colors_per_bar, draw: bool = True):
    """軸の下余白に mm 座標の軸を重ね、群の色の帯を描く。帯の下端の y (mm) を返す。
    draw=False (群の色なし)なら帯は描かず、位置だけ同じにする (群名・「n < 3」の位置は変えない。261013 決定)。"""
    from matplotlib.patches import Rectangle
    bx = fig.add_axes([PS.MARGIN_L / fig_w, 0, ax_w / fig_w, PS.MARGIN_B / fig_h])
    bx.set_xlim(0, ax_w)
    bx.set_ylim(-PS.MARGIN_B, 0)
    bx.axis("off")
    y = -BAND_GAP_MM - BAND_MM
    for i, c in enumerate(colors_per_bar if draw else []):
        bx.add_patch(Rectangle((i * pitch, y), pitch, BAND_MM, facecolor=c, lw=0))
    return bx, y


def _band_note(bx, y: float, pitch: float, j: int, color, text: str) -> None:
    """群の色の帯の中に短い文字を書く (例「n < 3」= 検定・計算から外した群。261009 決定)。
    文字の色は帯の色の明るさで白か前景色を選ぶ。"""
    r, g, b = matplotlib.colors.to_rgb(color)
    tc = PS.GREY if 0.299 * r + 0.587 * g + 0.114 * b > 0.55 else "white"
    bx.text(j * pitch + pitch / 2, y + BAND_MM / 2, text, ha="center", va="center", fontsize=PS.FS_NOTE, color=tc)


def fig08_stacked_group(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#7 積み上げ棒・群平均 (階級ごとに 1 枚。260921 決定)。
    検体ごとに相対存在量 (%) → 群で平均。除外検体は平均に入れない。上位 N = 群平均の合計順 (群を同じ重みで扱う)。
    1 群 8 mm・1 位が一番下・Others が一番上。葉の下に群の色の帯・その下に名前 45° (use_subgroup: true ならサブグループ名だけ)。"""
    log = ctx["log"]
    top_n, ranks = _plot_cfg(ctx)
    for k, rank in enumerate(ranks, 1):
        name = f"07_stacked-group_{rank}"
        log.step("10", n, total, False, f"#7 積み上げ棒・群平均 ({k}/{len(ranks)})", f"{name}.pdf")
        ra, excl, merged = _taxa_relative(ctx, rank)
        st, top = _top_taxa(ra, G, top_n)
        grp = ra.columns.map(G.label_by_newid)
        order = [g for g in G.order if g in set(grp)]
        gm = pd.DataFrame({g: st.loc[:, grp == g].mean(axis=1) for g in order})
        colors = _taxa_colors(len(top))

        W = PS.MM_PER_GROUP
        fig, ax = PS.fixed_axes(W * len(order))
        _stack(ax, np.arange(len(order)), gm, BAR_FRAC, colors)
        ax.set_xlim(-0.5, len(order) - 0.5)
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 25, 50, 75, 100])          # 目盛 5 個（260924 決定）
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
        ax.set_ylabel("Relative abundance (%)")
        ax.set_title(rank.capitalize(), fontsize=PS.TITLE_FS, loc="left", pad=2)
        fig_w, fig_h = fig.get_size_inches() / PS.MM
        bx, y = _bottom_band(fig, fig_w, fig_h, W * len(order), W, [G.color[g] for g in order], not G.no_color)
        shown = [g.split(" ", 1)[1] for g in order] if (G.use_subgroup and G.show_subgroup) else order
        for i, lab in enumerate(shown):
            bx.text(i * W + W / 2 + 1.0, y - 1.0, lab, fontsize=PS.XLABEL_FS, rotation=PS.XLABEL_ROT,
                    ha="right", va="top", rotation_mode="anchor")
        over = [(x, h) for x in shown if (h := PS.rotated_label_mm(x) + BAND_GAP_MM + BAND_MM + 1.0) > PS.MARGIN_B]
        if over:                                           # 名前が下余白に収まらないとき（260924 指示）
            log(f"  [注意] 群名が長く、下余白 {PS.MARGIN_B} mm で見切れている可能性があります (群名を短くしてください): "
                + ", ".join(f"{x} (約 {h:.1f} mm)" for x, h in over))
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)
        _taxa_legend(_fig_dir(d, name) / f"{name}_legend.pdf", top + ["Others"], colors, log)

        cd = _csv_dir(d, name)
        n_by = pd.Series(grp).value_counts()
        gm.T.assign(n=[int(n_by[g]) for g in order]).rename_axis("group").to_csv(cd / f"{name}.values.csv")
        log.sum(f"{rank}: 上位 {len(top)}")
        log(f"  {rank}: {len(ra)} タクソン (名前が重なって上の階級で区別した行 {merged})/ 上位 {len(top)} (群平均の合計順)の全検体平均の合計 "
            f"{st.loc[top].sum().mean():.1f}% / 除外検体 {len(excl)} は平均に入れない")


def fig09_stacked_sample(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#7 積み上げ棒・検体ごと (階級ごとに 1 枚。260921 決定)。
    1 検体 4 mm (#6 検体版と同じ)。並び = 群の順 → サブグループの定義順 → newID の末尾の番号。
    除外検体は枠だけ描き、中に「< 10,000 reads」 (下限 floor の値。260922 決定)。葉の下は #6 検体版と同じ (群の色の帯 → サブグループの記号 → 検体名 45°)。"""
    import re
    from matplotlib.patches import Rectangle
    log = ctx["log"]
    top_n, ranks = _plot_cfg(ctx)
    markers = PS.subgroup_markers(G.subgroups) if G.show_subgroup else {}
    subs = list(G.subgroups)
    srank = {x: i for i, x in enumerate(subs)}
    gidx = {g: i for i, g in enumerate(G.order)}

    def rep(s):
        m = re.search(r"_(\d+)$", s)
        return int(m.group(1)) if m else 0

    for k, rank in enumerate(ranks, 1):
        name = f"07_stacked-sample_{rank}"
        log.step("10", n, total, False, f"#7 積み上げ棒・検体ごと ({k}/{len(ranks)})", f"{name}.pdf")
        ra, excl, _ = _taxa_relative(ctx, rank)
        floor = int(ctx["cfg"]["diversity"]["sampling_depth_floor"])   # 除外検体の枠に書く（N.D. は「未検出」と読まれるため。260922 決定）
        st, top = _top_taxa(ra, G, top_n)
        colors = _taxa_colors(len(top))
        samples = sorted(list(ra.columns) + sorted(excl & set(G.label_by_newid.index)),
                         key=lambda s: (gidx[G.label_by_newid[s]], srank.get(G.subgroup_by_newid[s], len(subs)), rep(s), s))
        st_all = st.reindex(columns=samples)
        # 除外検体は st_all で NaN（棒を描かない）
        N = len(samples)

        fig, ax = PS.fixed_axes(N * LEAF_MM)
        x = np.arange(N)
        _stack(ax, x, st_all.fillna(0), BAR_FRAC, colors)
        for i, s in enumerate(samples):
            if s in excl:
                ax.add_patch(Rectangle((i - BAR_FRAC / 2, 0), BAR_FRAC, 100, fill=False, ec=PS.GREY, lw=0.5))
                ax.text(i, 50, f"< {floor:,} reads", rotation=90, rotation_mode="anchor", ha="center", va="center_baseline",
                        fontsize=PS.FS_NOTE)   # 文字の下の余白（ディセンダ）の分だけ左に寄るので、ベースラインの中央でそろえる（260922）
        ax.set_xlim(-0.5, N - 0.5)
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 25, 50, 75, 100])          # 目盛 5 個（260924 決定）
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
        ax.set_ylabel("Relative abundance (%)")
        ax.set_title(rank.capitalize(), fontsize=PS.TITLE_FS, loc="left", pad=2)
        fig_w, fig_h = fig.get_size_inches() / PS.MM
        bx, y = _bottom_band(fig, fig_w, fig_h, N * LEAF_MM, LEAF_MM, [G.color[G.label_by_newid[s]] for s in samples])
        if markers:
            y -= 2.0
            for i, s in enumerate(samples):
                sub = G.subgroup_by_newid[s]
                PS.scatter_subgroup(bx, [i * LEAF_MM + LEAF_MM / 2], [y], "#666666", G.subgroups.get(sub, {}),
                                    markers[sub], s=PS.SCATTER_S, clip_on=False)
            y -= 2.0
        else:
            y -= 1.0
        for i, s in enumerate(samples):
            bx.text(i * LEAF_MM + LEAF_MM / 2, y, s, fontsize=PS.XLABEL_FS, rotation=PS.XLABEL_ROT, rotation_mode="anchor",
                    ha="right", va="top")   # 検体名 45° (261022 決定。見切れてよい)
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)
        _taxa_legend(_fig_dir(d, name) / f"{name}_legend.pdf", top + ["Others"], colors, log)

        cd = _csv_dir(d, name)
        st_all.T.assign(group=[G.label_by_newid[s] for s in samples],
                        excluded=[s in excl for s in samples]
                        ).rename_axis("sample").to_csv(cd / f"{name}.values.csv")
        log.sum(f"{rank}: 検体 {N}")
        log(f"  {rank}: {N} 検体 (除外 {len(excl & set(samples))} は枠と「< {floor:,} reads」)/ 図の幅 {fig_w:.0f} mm")


LETTER_MAX = 4         # #8 の有意差の文字の上限（8 pt で 4 文字 6.1 mm・マス 8 mm。261009 決定）
HEAT_ROW_MM = 4        # ヒートマップの 1 行（1 タクソン）の高さ。8 pt の名前が重ならない最小の目安
HEAT_VMIN, HEAT_VMAX, HEAT_PSEUDO = -3, 0, 1e-6   # log10(相対存在量 + 1e-6)・範囲 −3〜0 (0.1〜100%)
# 261016 決定 (E6d): 範囲を 0.1〜100% (log10 の割合で −3〜0) にし、カラーバーの目盛は % で書く。それまでは −5〜0 (0.001%〜) で目盛は log10。
# 0.1% 未満のマスはいちばん暗い色になる (0 も含む)。#8・#11・#12 で共通。


def _pct_ticklabels() -> list[str]:
    """存在量のカラーバーの目盛 (割合の log10 → %)。例 −3〜0 → 0.1・1・10・100。"""
    return [f"{10 ** (t + 2):g}" for t in range(HEAT_VMIN, HEAT_VMAX + 1)]


def fig10_heatmap(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#8 ヒートマップ・群平均 (階級ごとに 1 枚。260922 決定)。
    行 = 上位 N タクソン (#7 と同じ選び方・同じ並び。1 位が一番上)・名前は右。列 = 群 (1 群 8 mm)・下に群の色の帯と名前 45°。
    値 = 群平均の相対存在量 (割合)を log10(x + 1e-6)・範囲 −3〜0 (0.1〜100%。261016)・viridis の 0〜95% (261009 に 0〜85% から変更。全図共通)。
    有意差 (stats.method。#1・#3・#10 と同じ検定。261009 決定):
      tukey   = マスに compact letter display (同じ文字を共有しない群どうしが p < 0.05)。1 マスの文字が LETTER_MAX を
                超える属が 1 つでもあれば、その図の文字はすべて出さない (注記に残す)
      dunnett・welch = 対照群との比較で有意な群のマスに星 (welch も対照群が必要。無ければ作図の最初で止まる)
    検定から外した群 (n < stats.min_n)は文字・星を付けず、帯の中に「n < 3」。"""
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    log, cfg = ctx["log"], ctx["cfg"]
    top_n, ranks = _plot_cfg(ctx)
    control = G.control                                # dunnett・welch の対照群（tukey は使わない）
    cmap = PS.viridis()                                # viridis の 0〜95%（全図共通・261009。それまで 0〜85%）
    norm = Normalize(HEAT_VMIN, HEAT_VMAX)
    for k, rank in enumerate(ranks, 1):
        name = f"08_heatmap_{rank}"
        log.step("10", n, total, False, f"#8 ヒートマップ ({k}/{len(ranks)})", f"{name}.pdf")
        ra, excl, _ = _taxa_relative(ctx, rank)
        st, top = _top_taxa(ra, G, top_n)
        grp = ra.columns.map(G.label_by_newid)
        order = [g for g in G.order if (grp == g).any()]
        gm = pd.DataFrame({g: ra.loc[top, grp == g].mean(axis=1) for g in order}) / 100    # 割合
        L = np.log10(gm + HEAT_PSEUDO)

        W, H = PS.MM_PER_GROUP, HEAT_ROW_MM
        ax_w, ax_h = W * len(order), H * len(top)
        name_w = max(PS.text_mm(t, PS.XLABEL_FS) for t in top) + 2
        # カラーバーは行列の左。バーの左端を他の図の Y 軸の位置（左余白 15 mm）に合わせ、目盛と見出しは左余白の中に置く（260923 決定）
        cb_w, cb_gap = 3, 3
        left = PS.MARGIN_L + cb_w + cb_gap
        fig_w = left + ax_w + name_w + PS.MARGIN_R
        fig_h = PS.MARGIN_B + ax_h + PS.MARGIN_T
        fig = matplotlib.pyplot.figure(figsize=(fig_w * PS.MM, fig_h * PS.MM))
        ax = fig.add_axes([left / fig_w, PS.MARGIN_B / fig_h, ax_w / fig_w, ax_h / fig_h])
        ax.imshow(L.to_numpy(), cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top, fontsize=PS.XLABEL_FS)
        ax.yaxis.tick_right()
        ax.tick_params(axis="y", length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(rank.capitalize(), fontsize=PS.TITLE_FS, loc="left", pad=2)

        rows, letter_rows, res = [], [], None
        for t_i, t in enumerate(top):
            vals = {g: ra.loc[t, grp == g].to_numpy() for g in order}
            res = GS.compare(vals, order, G.method, control, G.min_n)
            if G.method == "tukey":
                for r in res["rows"]:
                    rows.append(dict(taxon=t, **r))
                usable = [g for g in order if g not in res["skipped"]]
                lt = GS.letters(usable, res["rows"]) if len(usable) >= 2 else {}
                letter_rows.append(dict(taxon=t, **{g: lt.get(g, "") for g in order}))
                continue
            for r in res["rows"]:                          # dunnett・welch: 対照群との比較に星
                other = r["groupB"] if r["groupA"] == control else r["groupA"]
                rows.append(dict(taxon=t, control=control, group=other, n_control=len(vals[control]),
                                 n_group=len(vals[other]), p=r["p"], test=r["test"], reason=r.get("reason", "")))
                star = PS.p_stars(r["p"])
                if star:
                    ax.text(order.index(other), t_i + 0.25, star, ha="center", va="center",   # 星の字形は上寄りなので下げる
                            fontsize=PS.FS_STAR, color="white")   # 太字は使わない（規約: 太字は表の見出しだけ。260923）
        letters_shown = False
        if letter_rows:
            longest = max(len(v) for r in letter_rows for g, v in r.items() if g != "taxon")
            if longest <= LETTER_MAX:
                for t_i, r in enumerate(letter_rows):
                    for j, g in enumerate(order):
                        if r[g]:
                            ax.text(j, t_i, r[g], ha="center", va="center", fontsize=PS.FS_CELL, color="white")
                letters_shown = True
            else:
                log(f"  [注記] {rank}: 有意差の文字が最長 {longest} 文字で上限 {LETTER_MAX} を超えるため、#8 の文字は描きません"
                    " (文字は letters.csv に残しています)")
        skipped_cols = res["skipped"] if res else []    # 群の n は属によらない
        if skipped_cols:
            _log_skipped(log, res)
        bx, y = _bottom_band(fig, fig_w, fig_h, ax_w, W, [G.color[g] for g in order], not G.no_color)
        bx.set_position([left / fig_w, 0, ax_w / fig_w, PS.MARGIN_B / fig_h])
        shown = [g.split(" ", 1)[1] for g in order] if (G.use_subgroup and G.show_subgroup) else order
        for g in skipped_cols:                         # 検定から外した群は帯の中に「n < 3」（261009 決定）
            _band_note(bx, y, W, order.index(g), "white" if G.no_color else G.color[g], f"n < {G.min_n}")
        for i, lab in enumerate(shown):
            bx.text(i * W + W / 2 + 1.0, y - 1.0, lab, fontsize=PS.XLABEL_FS, rotation=PS.XLABEL_ROT,
                    ha="right", va="top", rotation_mode="anchor")
        over = [(x, h) for x in shown if (h := PS.rotated_label_mm(x) + BAND_GAP_MM + BAND_MM + 1.0) > PS.MARGIN_B]
        if over:                                           # 名前が下余白に収まらないとき（260924 指示）
            log(f"  [注意] 群名が長く、下余白 {PS.MARGIN_B} mm で見切れている可能性があります (群名を短くしてください): "
                + ", ".join(f"{x} (約 {h:.1f} mm)" for x, h in over))
        cb_len = min(PS.PLOT_H_MM, ax_h)          # 長さは軸高 48 mm に固定・上端を行列の上端にそろえる（規約 260929。低い行列ではその高さ）
        ext = cb_len * 0.05                           # 下端の三角 (既定の長さ = 本体の 5%) は本体の下に足す。本体は 48 mm で行列の上下の端にそろえる (261019 決定・案 b)
        cax = fig.add_axes([PS.MARGIN_L / fig_w, (PS.MARGIN_B + ax_h - cb_len - ext) / fig_h, cb_w / fig_w, (cb_len + ext) / fig_h])
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, extend="min")   # 下端の三角 = 0.1% 未満 (0 を含む) もこの色 (261016 決定)
        cax.yaxis.set_ticks_position("left")
        cax.yaxis.set_label_position("left")
        cb.outline.set_visible(False)
        cb.set_label("Relative abundance (%)", fontsize=PS.FS_CBAR)   # 目盛は % (261016 決定・E6d)
        cb.ax.tick_params(labelsize=PS.FS_CBAR, length=2)
        cb.set_ticks(range(HEAT_VMIN, HEAT_VMAX + 1))
        cb.set_ticklabels(_pct_ticklabels())
        fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
        matplotlib.pyplot.close(fig)

        cd = _csv_dir(d, name)
        (gm * 100).rename_axis("taxon").to_csv(cd / f"{name}.values.csv")   # CSV は % (#7・#10 とそろえる。色の計算は割合のまま。261021・F1)
        if rows:
            pd.DataFrame(rows).to_csv(cd / f"{name}.stats.csv", index=False)
        if letter_rows:                                # 文字は図に出したかどうかによらず残す
            pd.DataFrame(letter_rows).to_csv(cd / f"{name}.letters.csv", index=False)
        nsig = sum(r["p"] < 0.05 for r in rows)
        n_tested = sum(r["test"] != "not_tested" for r in rows)
        log.sum(f"{rank}: {len(top)} タクソン × {len(order)} 群")
        log(f"  {rank}: {len(top)} タクソン × {len(order)} 群 / 図 {fig_w:.0f} × {fig_h:.0f} mm / "
            + (f"文字 = Tukey HSD の全ペア (p < 0.05)・有意 {nsig} / {n_tested}・図に{'出した' if letters_shown else '出さない'}"
               if G.method == "tukey" else f"星 = 対照群 {control} との比較 ({G.method})・p < 0.05 {nsig} / {n_tested}"))


# ---------------- #9 ベン図 / UpSet（260922 決定） ----------------
VENN_AX_MM = PS.SCATTER_AX_MM   # ベン図の軸（散布図と同じ 48 × 48 mm）
UPSET_CELL_MM = 8      # UpSet の点の行列: 格子 8 × 8 mm（1 群 8 mm にそろえる。260923 決定・旧 4 mm）
UPSET_DOT_S, UPSET_LINE_LW = 60, 1.2   # 点の面積・つなぐ線の太さ（9 パターンから選択。260923）
UPSET_GAP_MM, UPSET_NUM_MM = 2, 4      # 行列と棒のあいだ・棒の先の数の分
# 4 群のベン図は楕円（中心 x, y・幅・高さ・角度。0〜1 の座標）。2・3 群は円
VENN_SHAPES = {
    2: [((0.37, 0.5), 0.5, 0.5, 0), ((0.63, 0.5), 0.5, 0.5, 0)],
    3: [((0.37, 0.58), 0.5, 0.5, 0), ((0.63, 0.58), 0.5, 0.5, 0), ((0.5, 0.36), 0.5, 0.5, 0)],
    4: [((0.350, 0.450), 0.72, 0.45, 140), ((0.450, 0.550), 0.72, 0.45, 140),   # 軸の中で上下の中央にそろえた
        ((0.544, 0.550), 0.72, 0.45, 40), ((0.644, 0.450), 0.72, 0.45, 40)],
    5: [((0.428, 0.449), 0.87, 0.50, 155), ((0.469, 0.543), 0.87, 0.50, 82), ((0.558, 0.523), 0.87, 0.50, 10),
        ((0.578, 0.432), 0.87, 0.50, 118), ((0.489, 0.383), 0.87, 0.50, 46)],   # pyvenn の 5 群の配置
}


def _presence(ra: pd.DataFrame, G) -> pd.DataFrame:
    """タクソン × 群の有無。「その群にいる」= 群の**全検体**で 1 リード以上 (260922 決定)。"""
    grp = ra.columns.map(G.label_by_newid)
    order = [g for g in G.order if (grp == g).any()]
    return pd.DataFrame({g: (ra.loc[:, grp == g] > 0).all(axis=1) for g in order})


def _in_shape(xx, yy, shape) -> np.ndarray:
    (cx, cy), w, h, ang = shape
    t = np.deg2rad(ang)
    dx, dy = xx - cx, yy - cy
    u = dx * np.cos(t) + dy * np.sin(t)
    v = -dx * np.sin(t) + dy * np.cos(t)
    return (u / (w / 2)) ** 2 + (v / (h / 2)) ** 2 <= 1


def _draw_venn(path: Path, P: pd.DataFrame, G, title: str) -> None:
    """2〜5 群のベン図。数は各領域の内側で境界から一番遠い点に書く。群名は書かない (凡例は別ファイル)。"""
    from matplotlib.patches import Ellipse
    from scipy.ndimage import distance_transform_edt
    order = list(P.columns)
    shapes = VENN_SHAPES[len(order)]
    fig, ax = PS.fixed_axes(VENN_AX_MM, VENN_AX_MM)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    for g, sh in zip(order, shapes):
        ax.add_patch(Ellipse(sh[0], sh[1], sh[2], angle=sh[3], facecolor=G.color[g], alpha=0.25, lw=0))
        ax.add_patch(Ellipse(sh[0], sh[1], sh[2], angle=sh[3], fill=False, ec=G.color[g], lw=0.8))
    n = 400
    xx, yy = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    inside = np.stack([_in_shape(xx, yy, sh) for sh in shapes])
    key = P.apply(lambda r: tuple(r.to_numpy()), axis=1)
    import itertools
    for combo in itertools.product([False, True], repeat=len(order)):
        if not any(combo):
            continue
        mask = np.all([inside[i] == combo[i] for i in range(len(order))], axis=0)
        if not mask.any():
            continue
        dist = distance_transform_edt(mask)
        iy, ix = np.unravel_index(np.argmax(dist), dist.shape)
        cnt = int((key == combo).sum())
        ax.text(xx[iy, ix], yy[iy, ix], str(cnt), ha="center", va="center", fontsize=PS.FS_CELL)
    ax.set_title(title, fontsize=PS.TITLE_FS, loc="left", pad=2)
    fig.savefig(path)
    matplotlib.pyplot.close(fig)


def _draw_upset(path: Path, P: pd.DataFrame, G, title: str, shown: list[str]) -> list[tuple]:
    """UpSet (260923 決定): 左 = 点の行列 (列 = 群・行 = 重なり・格子 8 × 8 mm)、右 = 重なりごとの数の横棒 (長さ 46 mm)。
    行列の下に横軸を置き、群の位置に目盛・その下に群名 45° (群の図と同じ)。棒には縦軸の線を入れる。
    点 = 面積 60・つなぐ線 = 1.2 pt (9 パターンを見比べて決定)。重なりは数の多い順にすべて出す。"""
    order = list(P.columns)
    key = P.apply(lambda r: tuple(r.to_numpy()), axis=1)
    counts = key[key.map(any)].value_counts()
    combos = list(counts.index)
    C, R, cell = len(combos), len(order), UPSET_CELL_MM
    mat_w, mat_h = R * cell, C * cell
    bar_w = PS.PLOT_H_MM
    fig_w = PS.MARGIN_L + mat_w + UPSET_GAP_MM + bar_w + UPSET_NUM_MM + PS.MARGIN_R
    fig_h = PS.MARGIN_B + mat_h + PS.MARGIN_T
    fig = matplotlib.pyplot.figure(figsize=(fig_w * PS.MM, fig_h * PS.MM))
    mx = fig.add_axes([PS.MARGIN_L / fig_w, PS.MARGIN_B / fig_h, mat_w / fig_w, mat_h / fig_h])
    mx.set_xlim(-0.5, R - 0.5)
    mx.set_ylim(C - 0.5, -0.5)
    for j, cb in enumerate(combos):
        on = [i for i in range(R) if cb[i]]
        mx.scatter(range(R), [j] * R, s=UPSET_DOT_S, color="#DDDDDD", lw=0, zorder=1)
        if len(on) > 1:
            mx.plot([min(on), max(on)], [j, j], color=PS.GREY, lw=UPSET_LINE_LW, zorder=2)
        mx.scatter(on, [j] * len(on), s=UPSET_DOT_S, c=[G.color[order[i]] for i in on], lw=0, zorder=3)
    mx.set_yticks([])
    mx.set_xticks(range(R))
    mx.set_xticklabels(shown, rotation=PS.XLABEL_ROT, ha="right", rotation_mode="anchor", fontsize=PS.XLABEL_FS)
    for sp in ("top", "right", "left"):
        mx.spines[sp].set_visible(False)
    mx.set_title(title, fontsize=PS.TITLE_FS, loc="left", pad=2)
    bx = fig.add_axes([(PS.MARGIN_L + mat_w + UPSET_GAP_MM) / fig_w, PS.MARGIN_B / fig_h, bar_w / fig_w, mat_h / fig_h])
    bx.barh(range(C), counts.to_numpy(), 0.6, color=PS.GREY, lw=0)
    for j, v in enumerate(counts.to_numpy()):
        bx.text(v, j, f" {v}", va="center", fontsize=PS.FS_NOTE)
    PS.fixed_ticks(bx, 0, float(counts.max()), axis="x")   # 目盛 5 個（260924 決定）
    bx.set_ylim(C - 0.5, -0.5)
    bx.set_yticks([])
    for sp in ("top", "right"):
        bx.spines[sp].set_visible(False)
    bx.set_xlabel("Taxa")
    fig.savefig(path)
    matplotlib.pyplot.close(fig)
    return [(cb, int(counts[cb])) for cb in combos]


def fig11_venn_upset(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#9 ベン図・UpSet (階級ごと。260922 決定)。
    「その群にいる」= 群の全検体で 1 リード以上。除外検体は入れない。名前はアーカイブと同じ整形。
    ベン図は 2〜5 群のときだけ (6 群以上は描かない。VENN_SHAPES の形がある群の数まで)。UpSet は群の数によらず描く。"""
    import matplotlib.ticker  # noqa: F401
    log = ctx["log"]
    _, ranks = _plot_cfg(ctx)
    for k, rank in enumerate(ranks, 1):
        log.step("10", n, total, False, f"#9 ベン図 / UpSet ({k}/{len(ranks)})", f"09_venn_{rank}.pdf・09_upset_{rank}.pdf")
        ra, excl, _ = _taxa_relative(ctx, rank)
        P = _presence(ra, G)
        order = list(P.columns)
        shown = [g.split(" ", 1)[1] for g in order] if (G.use_subgroup and G.show_subgroup) else order
        title = rank.capitalize()
        if G.no_color:                                  # 群の色なしではベン図を描かない（円の色だけで群を見分けるため。261013）
            log(f"  [注記] {rank}: group_subgroup_colors: none のため、ベン図は描きません (UpSet は描きます)")
        elif 2 <= len(order) <= max(VENN_SHAPES):
            _draw_venn(_fig_dir(d, "09") / f"09_venn_{rank}.pdf", P, G, title)
        else:
            log(f"  [注記] ベン図は 2〜{max(VENN_SHAPES)} 群のときだけ描くため、描きません ({len(order)} 群)")
        if P.any(axis=1).any():
            combos = _draw_upset(_fig_dir(d, "09") / f"09_upset_{rank}.pdf", P, G, title, shown)
        else:                                          # どの群の全検体にもいる属が無い (261019・F217)
            combos = []
            log(f"  [注記] {rank}: どの群の全検体にもいるタクソンが無いため、UpSet は描きません")
        _warn_labels(log, shown)                           # UpSet の群名が下余白に収まらないとき（260925 指示）
        cd = _csv_dir(d, f"09_venn-upset_{rank}")
        key = P.apply(lambda r: tuple(r.to_numpy()), axis=1)
        rows = [dict(groups=" & ".join(g for g, b in zip(order, cb) if b), n_taxa=c,
                     taxa="; ".join(sorted(key.index[key == cb]))) for cb, c in combos]
        pd.DataFrame(rows).to_csv(cd / f"09_venn-upset_{rank}.intersections.csv", index=False)
        P.rename_axis("taxon").to_csv(cd / f"09_venn-upset_{rank}.presence.csv")
        n_none = int((~P.any(axis=1)).sum())
        log.sum(f"{rank}: 重なり {len(combos)} 通り")
        log(f"  {rank}: {len(P)} タクソン / どこかの群の全検体にいる {len(P) - n_none} (どの群にもいない {n_none} は図に出さない)"
            f" / 重なり {len(combos)} 通り / 群ごと: " + ", ".join(f"{g} {int(P[g].sum())}" for g in order))


def _safe_name(s: str) -> str:
    """タクソン名をファイル名にする (記号を置き換える)。対応は CSV に残す。"""
    return re.sub(r"[^0-9A-Za-z._-]+", "-", s).strip("-")[:60]


def _taxon_bar_axis(ax, need: float) -> float:
    """#10 の縦軸。上端 = need 以上の決まった段 (TAXON_BAR_TOPS)・目盛 5 個。
    100 を超えるとき (平均 + SD やブラケット) は 25 刻みの次の段 (125・150 …) まで伸ばし、目盛は 0〜100 だけ
    (100 を超える数値は出さない。261016 E9a・261020 決定)。返り値: 上端。"""
    if need > 100:
        top = 100 + 25 * math.ceil((need - 100) / 25)
        ticks = [100 * k / PS.TICK_N for k in range(PS.TICK_N + 1)]
    else:
        top = next(v for v in TAXON_BAR_TOPS if v >= need)
        ticks = [top * k / PS.TICK_N for k in range(PS.TICK_N + 1)]   # 目盛 5 個
    ax.set_ylim(0, top)
    ax.set_yticks(ticks)
    return top


def fig12_taxon_bar(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#10 タクソン別の相対存在量の棒 (タクソンごとに 1 枚・表の全列。260923 決定)。
    棒 = 群平均 (群の色)・SD のエラーバー・検体の点 (#1・#3 と同じ描き方の棒版)。
    縦軸 = Relative abundance (%)。検定は stats.method (タクソンごと。タクソンのあいだは補正しない)。
    use_subgroup: true でも描く (指示)。"""
    log = ctx["log"]
    _, ranks = _plot_cfg(ctx)
    for k, rank in enumerate(ranks, 1):
        name = f"10_taxon-bar_{rank}"
        ra, excl, _ = _taxa_relative(ctx, rank)
        grp = ra.columns.map(G.label_by_newid)
        order = [g for g in G.order if (grp == g).any()]
        zero = set(ra.index[(ra[[c for c, g in zip(ra.columns, grp) if g in order]] == 0).all(axis=1)])
        taxa = [t for t in ra.index if t not in zero]   # 残った検体ですべて 0 の属は描かない (261020 決定・b)
        fname, seen = {}, {}                           # ファイル名が重なったら -2, -3 … を付ける (Mac は大文字小文字を区別しない。261019・F218)
        for t in taxa:
            base = f"{name}_{_safe_name(t)}"
            k2 = seen.get(base.casefold(), 0) + 1
            seen[base.casefold()] = k2
            fname[t] = f"{base}.pdf" if k2 == 1 else f"{base}-{k2}.pdf"
        log.step("10", n, total, False, f"#10 タクソン別の棒 ({k}/{len(ranks)})",
                 f"{name}_<タクソン>.pdf ({len(taxa)} 枚)")
        if zero:
            log(f"  [注記] {rank}: 残った検体ですべて 0 の属 {len(zero)} 個は描きませんでした (stats.csv には not_tested・全検体で 0 と書きます)")
        _warn_labels(log, [g.split(" ", 1)[1] for g in order] if (G.use_subgroup and G.show_subgroup) else order)   # 階級ごとに 1 回（260925 指示）
        sub = _csv_dir(d, name)                          # CSV は他の図と同じ「taxon-bar」
        pdfs = _fig_dir(d, name)                        # PDF（X01 で 1,027 枚）は 10_taxon-bar/ の直下
        rows, stat_rows, n_sig = [], [], 0
        for t in taxa:
            data = [ra.loc[t, grp == g].to_numpy(dtype=float) for g in order]
            res = GS.compare({g: v for g, v in zip(order, data)}, order, G.method, G.control, G.min_n)
            fig, ax = PS.fixed_axes(PS.MM_PER_GROUP * len(order))
            x = np.arange(len(order))
            means = np.array([v.mean() if len(v) else np.nan for v in data])
            sds = np.array([v.std(ddof=1) if len(v) > 1 else np.nan for v in data])   # n < 2 は SD を出さない (261020・F011・F259。#1・#3 と同じ)
            # 前後は 棒（一番後ろ）→ 点 → エラーバー（一番前）。棒は群の色を薄く、点は群の色（260924 決定・案 A）
            ax.bar(x, means, BAR_FRAC, color=[(*matplotlib.colors.to_rgb(G.color[g]), BAR_ALPHA) for g in order],
                   lw=0, zorder=1)
            rng = np.random.default_rng(PS.JITTER_SEED)
            for i, v in enumerate(data):
                ax.scatter(i + rng.uniform(-PS.DOT_JITTER, PS.DOT_JITTER, len(v)), v, s=PS.DOT_S,
                           color=G.color[order[i]], lw=0.3, edgecolors="white", zorder=2)
            ok = np.isfinite(sds)
            ax.errorbar(x[ok], means[ok], yerr=sds[ok], fmt="none", ecolor=PS.GREY, elinewidth=PS.MEAN_LW,
                        capsize=2, capthick=PS.MEAN_LW, zorder=3)
            lo, hi, _step = _ylim_for(G, res, data)
            hi = _taxon_bar_axis(ax, hi)
            ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
            hb = _draw_significance(ax, G, res, {i: i for i in range(len(order))}, data, 0, hi, order)
            if hb > hi:                                # ブラケットが入りきらなければ次の段へ上げる (261020 決定・F242/F184)
                _taxon_bar_axis(ax, hb)
                PS.keep_brackets(ax)                   # ブラケットは描いたときの見た目のまま (B2-3 と同じ)
            ax.set_xlim(-0.5, len(order) - 0.5)
            ax.set_xticks(x)
            shown = [g.split(" ", 1)[1] for g in order] if (G.use_subgroup and G.show_subgroup) else order
            ax.set_xticklabels(shown, rotation=PS.XLABEL_ROT, ha="right", rotation_mode="anchor",
                               fontsize=PS.XLABEL_FS)
            ax.set_ylabel("Relative abundance (%)")
            title = PS.wrap_title(t, PS.MM_PER_GROUP * len(order))
            if title.count("\n") >= TITLE_MAX_LINES:     # 上余白 13 mm に入るのは 3 行まで (261021・F4)
                log(f"  [注意] {rank}: {t} のタイトルが {title.count(chr(10)) + 1} 行になり、上で切れている可能性があります")
            ax.set_title(title, fontsize=PS.TITLE_FS, loc="left", pad=2)
            fig.savefig(pdfs / fname[t])
            matplotlib.pyplot.close(fig)
            for g, m, sd, v in zip(order, means, sds, data):
                rows.append(dict(taxon=t, group=g, n=len(v), mean=m, sd=sd))
            vals = dict(zip(order, data))
            for r in res["rows"]:
                if r["test"] != "not_tested" and np.isnan(r["p"]) and not (vals[r["groupA"]].any() or vals[r["groupB"]].any()):
                    r.update(p=np.nan, test="not_tested", reason="両群とも全検体で 0")   # 0 ÷ 0 で p が出ない組 (261020)
                stat_rows.append(dict(taxon=t, **r))
                n_sig += r["p"] < 0.05
        pairs = (list(itertools.combinations(order, 2)) if G.method == "tukey"
                 else [(G.control, g) for g in order if g != G.control] if G.method == "dunnett" else [tuple(order[:2])])
        for t in [t for t in ra.index if t in zero]:   # 検定もしない (0 ÷ 0 で p が NaN になるため)
            for a, b in pairs:
                stat_rows.append(dict(taxon=t, groupA=a, groupB=b, n_A=int((grp == a).sum()), n_B=int((grp == b).sum()),
                                      p=np.nan, test="not_tested", reason="全検体で 0"))
        pd.DataFrame(rows).to_csv(sub / f"{name}.values.csv", index=False)
        st = pd.DataFrame(stat_rows).assign(note="タクソンのあいだの補正なし")   # test 列は実際に使った検定のまま
        if taxa:
            _log_skipped(log, res)                   # 群の n はタクソンによらないので階級ごとに 1 回
        st.to_csv(sub / f"{name}.stats.csv", index=False)
        pd.DataFrame({"taxon": taxa, "file": [fname[t] for t in taxa]}).to_csv(
            sub / f"{name}.files.csv", index=False)
        n_test = int(st["p"].notna().sum()) if len(st) else 0
        log.sum(f"{rank}: {len(taxa)} 枚・有意 {n_sig}/{n_test} (タクソンのあいだの補正なし)")
        log(f"  {rank}: {len(taxa)} タクソン × {len(order)} 群 / 検定 {n_test} 件・p < 0.05 {n_sig} 件"
            f" ({G.method}・タクソンのあいだの補正なし) / 除外検体 {len(excl)} は入れない")


# ---------------- #11 ANCOM-BC バブルプロット（260924 決定・README #11） ----------------
# 261018 決定: 計算を ANCOM-BC (scikit-bio) から ANCOM-BC2 (R ANCOMBC・Dunnett 型) へ移した (BC2 実装手順書・lib/ancombc2.R)。
# 設定は lib/ancombc2.R に固定 (pseudo 0・struc_zero TRUE / neg_lb FALSE・BH・mdFDR holm / B 1000・prv_cut 0)
ANCOM_PRV = 0.0        # 出現率で外さない (261018 決定・既定 0.10 から変更。群の大きさが違う設計で小さい群に偏るため)
ANCOM_SEED = 0         # mdFDR の乱数の種 (パイプラインのほかの種とそろえる。261018 決定)
ANCOM_B = 1000         # mdFDR の乱数の回数 (lib/ancombc2.R と同じ値。記録用)
ANCOM_SIZES = [(0.001, "q ≤ 0.001", 70), (0.01, "0.001 < q ≤ 0.01", 38), (0.05, "0.01 < q ≤ 0.05", 16)]
ANCOM_NS_S = 4         # 有意でない点（灰色の小さな点）
BURD = ["#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"]   # Tol BuRd


def _ancom_control(ctx: dict, G) -> str | None:
    c = PC._control_written(ctx["cfg"])                   # control_group だけ (composite は廃止。261021 決定)
    if c and c not in G.order:
        raise PC.Stop(f"対照群 '{c}' がデータの群にありません。存在する群: {G.order}")
    return c


def _ancom_run(ctx: dict, G, rank: str, control: str):
    """ANCOM-BC2 の計算 (#11 と #12 で共通。1 回の実行の中では階級ごとに 1 回だけ計算して使い回す)。
    R ANCOMBC の ancombc2 (Dunnett 型) を lib/ancombc2.R で呼ぶ (261018 決定)。
    返り値: (結果の表 [taxon, group, ln_fc, SE, W, pvalue, qvalue, diff, passed_ss, sig, struc_zero, zero_in_control,
            zero_in_group, not_tested], 比較する群, 全タクソン数, 検定したタクソン数, n < stats.min_n で外した群 {群: n})。
    sig = diff_robust (疑似カウントの感度解析をすべて通った有意。図の有意はこれだけで決める)。
    構造的ゼロ (片方の群で全検体 0) の組は検定しない。どれか 1 群でも構造的ゼロの属は、ほかの比較も検定されない
    (ANCOM-BC2 の設計どおり。U7 決定 B1)。その組は not_tested = True (図では空欄)。
    対照群が n < min_n なら None (261009 決定・D2e)"""
    import subprocess
    import tempfile
    cache = ctx.setdefault("_ancom_cache", {})
    if (rank, control) in cache:
        return cache[(rank, control)]
    out, log = ctx["out"], ctx["log"]
    t, _, uniq = _display_names(out, rank)
    t.index = uniq
    excl = set(pd.read_csv(out / "06_rarefaction" / "excluded_samples.tsv", sep="\t", dtype=str)["newID"])
    X = t[[c for c in t.columns if c not in excl]].T
    n_all = X.shape[1]
    counts = X.index.map(G.label_by_newid).value_counts()
    # 1 ≤ n < min_n の群は計算から外し、列は残して「n < 3」。n = 0 の群は列も出さない (261020 決定)
    small = {g: int(counts.get(g, 0)) for g in G.order if 0 < counts.get(g, 0) < G.min_n}
    if counts.get(control, 0) < G.min_n:               # 対照群が足りなければ計算しない
        cache[(rank, control)] = None
        return None
    X = X[~X.index.map(G.label_by_newid).isin(list(small))]   # n < min_n の群はモデルから外す（列を出さない）
    nz = (X > 0).sum()
    X = X.loc[:, (nz > 0) & (nz >= ANCOM_PRV * len(X))]      # 全検体 0 の属は外す (prv_cut 0 でも)
    groups = [g for g in G.order if g in set(X.index.map(G.label_by_newid))]
    # 各群ちょうど 1 検体ずつだけに出る属は、0 でない検体数 = 群の数で自由度 0 になり、R が p = 1 (差なし)にする (261020 決定・b)。
    # そのとき R が出す NaN・NA・max の警告は、この注記に置き換える
    gi = X.index.map(G.label_by_newid)
    nz_g = (X > 0).groupby(gi).sum()
    dof0 = [tx for tx in X.columns if (nz_g[tx] == 1).all()]
    code = {g: f"g{i}" for i, g in enumerate(groups, 1)}      # 群名は空白などを含みうるので R にはコードで渡す
    meta = pd.DataFrame({"grp": X.index.map(G.label_by_newid).map(code)}, index=X.index)
    rscript = Path(sys.prefix) / "bin" / "Rscript"            # QIIME 2 の env の R (check_env で env の中にいることは確認済み)
    with tempfile.TemporaryDirectory() as wd:
        X.T.to_csv(Path(wd) / "counts.csv")
        meta.to_csv(Path(wd) / "meta.csv")
        p = subprocess.run([str(rscript), str(Path(__file__).resolve().parent / "ancombc2.R"), wd, code[control],
                            ",".join(code[g] for g in groups), str(ANCOM_SEED)], capture_output=True, text=True)
        if p.returncode != 0:
            raise PC.Stop(f"ANCOM-BC2 (R) が失敗しました ({rank}):\n{p.stderr[-2000:]}")
        dn = pd.read_csv(Path(wd) / "res_dunn.csv").set_index("taxon")
        zi = pd.read_csv(Path(wd) / "zero_ind.csv").set_index("taxon")
        warns = [w for w in (Path(wd) / "warnings.txt").read_text().splitlines() if w.strip()]
        version = (Path(wd) / "session.txt").read_text().splitlines()[0].split()[-1]
    if dof0:
        log(f"  [注記] {rank}: 各群 1 検体ずつしか出ないため、自由度 0 で p = 1 (差なし) として扱われた属があります: {', '.join(dof0)}")
        noise = ("NaN", "NA が生成", "NAs produced", "max の引数に有限な値", "no non-missing arguments to max")
        warns = [w for w in warns if not any(k in w for k in noise)]
    if len(groups) == 2:                               # 2 群では ANCOMBC が必ずこの 2 行を出す (主解析を使うので想定どおり。261020)
        warns = [w for w in warns if "< 3 categories" not in w and "will be deactivated" not in w]
    for w in warns:                                    # 収束しないなどの R の警告は残す (261018 決定・O7)
        log(f"  [注意] {rank}: ANCOM-BC2 (R) の警告: {w}")
    zero = {g: zi[f"structural_zero (grp = {code[g]})"].astype(bool) for g in groups}
    comps = [g for g in groups if g != control]
    rows = []
    for tx in zi.index:
        for g in comps:
            zc, zg = bool(zero[control][tx]), bool(zero[g][tx])
            r = dict(taxon=tx, group=g, n_control=int(counts[control]), n_group=int(counts[g]),   # 計算に使った検体数 (除外後。261021・F2)
                     zero_in_control=zc, zero_in_group=zg, struc_zero=zc or zg)
            c = f"grp{code[g]}"
            if tx in dn.index:
                r.update(ln_fc=dn.at[tx, f"lfc_{c}"], SE=dn.at[tx, f"se_{c}"], W=dn.at[tx, f"W_{c}"],
                         pvalue=dn.at[tx, f"p_{c}"], qvalue=dn.at[tx, f"q_{c}"], diff=bool(dn.at[tx, f"diff_{c}"]),
                         passed_ss=bool(dn.at[tx, f"passed_ss_{c}"]), sig=bool(dn.at[tx, f"diff_robust_{c}"]),
                         not_tested=False)
            else:                                      # 構造的ゼロのある属は検定されない (この組が構造的ゼロでなければ not_tested)
                r.update(ln_fc=np.nan, SE=np.nan, W=np.nan, pvalue=np.nan, qvalue=np.nan, diff=False, passed_ss=False,
                         sig=False, not_tested=not r["struc_zero"])
            rows.append(r)
    res = pd.DataFrame(rows)
    res.attrs["ancombc_version"] = version
    cache[(rank, control)] = (res, comps, n_all, len(dn), small)
    return cache[(rank, control)]


def _ancom_draw(G, d: Path, name: str, rank: str, res, comps: list, small: dict, control: str, rows: list,
                base_mean, log):
    """#11・#12 共通の描画 (261009 決定: #11 にも対照群の平均相対存在量の列を入れる。#12 と同じ形)。
    列 = 群の並びどおり (対照群の列 = 存在量、他 = バブル)。返り値: (図の幅, 図の高さ, 載せた行の結果)"""
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    cmap_fc = LinearSegmentedColormap.from_list("BuRd", BURD)
    cmap_ab = PS.viridis()                             # #8 と同じ（viridis の 0〜95%）
    norm_ab = Normalize(HEAT_VMIN, HEAT_VMAX)
    sub = res[res.taxon.isin(rows)]
    L = sub.pivot(index="taxon", columns="group", values="ln_fc").reindex(index=rows, columns=comps)
    Q = sub.pivot(index="taxon", columns="group", values="qvalue").reindex(index=rows, columns=comps)
    A = np.log10(base_mean.reindex(rows).fillna(0).to_numpy() + HEAT_PSEUDO)

    W, H = PS.MM_PER_GROUP, HEAT_ROW_MM
    cols = [g for g in G.order if g == control or g in comps or g in small]   # 計算から外した群も列は残す（261009）   # 対照群の列は群の並び (group_order) どおりの位置（261007 決定・案 B）
    jc = cols.index(control)
    ax_w, ax_h = W * len(cols), H * len(rows)
    name_w = max(PS.text_mm(tx, PS.XLABEL_FS) for tx in rows) + 2
    cb_w, cb_gap, cb_lab = 3, 3, 15                 # カラーバーは 2 本。左が存在量、右が ln fold change（目盛と見出しは各バーの左）
    cb2_x = PS.MARGIN_L + cb_w + cb_lab
    left = cb2_x + cb_w + cb_gap
    fig_w = left + ax_w + name_w + PS.MARGIN_R
    fig_h = PS.MARGIN_B + ax_h + PS.MARGIN_T
    fig = matplotlib.pyplot.figure(figsize=(fig_w * PS.MM, fig_h * PS.MM))
    ax = fig.add_axes([left / fig_w, PS.MARGIN_B / fig_h, ax_w / fig_w, ax_h / fig_h])
    vmax = float(np.nanmax(np.abs(L.to_numpy())))
    step = PS.tick_step(vmax, 0, 2)
    vmax = step * 2
    norm_fc = Normalize(-vmax, vmax)
    from matplotlib.patches import Rectangle
    SZ = res.set_index(["taxon", "group"])
    vals = base_mean.reindex(rows).fillna(0).to_numpy() * 100
    for i, tx in enumerate(rows):
        ax.add_patch(Rectangle((jc - 0.5, i - 0.5), 1, 1, fc=cmap_ab(norm_ab(A[i])), lw=0, zorder=1))   # 対照群の存在量
        v = vals[i]
        ax.text(jc, i, f"{v:.1f}" if v >= 0.1 else ("<0.1" if v > 0 else "0"), ha="center", va="center",
                fontsize=PS.FS_CELL, color="white", zorder=3)     # 値 (%)・文字は白（261007 指示）
        for g in comps:
            j = cols.index(g)
            q, v = Q.loc[tx, g], L.loc[tx, g]
            r = SZ.loc[(tx, g)]
            if r["struc_zero"]:                       # 構造的ゼロ → 白抜き (BC2 では q が無いので先に判定する。261018)
                _struc_mark(ax, j, i, r)
                continue
            if r["not_tested"] or pd.isna(q):         # ほかの群の構造的ゼロで検定されなかった組 → 空欄 (U7 B1)
                continue
            if r["sig"]:                              # 有意 = diff_robust (感度解析を通ったもの。261014 D3)
                sz = next(z for th, _, z in ANCOM_SIZES if q <= th)
                ax.scatter(j, i, s=sz, color=cmap_fc(norm_fc(v)), edgecolors=PS.GREY, linewidths=0.3, zorder=2)
            else:
                ax.scatter(j, i, s=ANCOM_NS_S, color="#BBBBBB", lw=0, zorder=2)
    ax.set_xlim(-0.5, len(cols) - 0.5)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=PS.XLABEL_FS)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", length=0, pad=2)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(rank.capitalize(), fontsize=PS.TITLE_FS, loc="left", pad=2)
    bx, y = _bottom_band(fig, fig_w, fig_h, ax_w, W, [G.color[g] for g in cols], not G.no_color)
    bx.set_position([left / fig_w, 0, ax_w / fig_w, PS.MARGIN_B / fig_h])
    for jj, g in enumerate(cols):
        if g in small:
            _band_note(bx, y, W, jj, "white" if G.no_color else G.color[g], f"n < {G.min_n}")
    shown = [g.split(" ", 1)[1] for g in cols] if (G.use_subgroup and G.show_subgroup) else cols
    for j, lab in enumerate(shown):
        bx.text(j * W + W / 2 + 1.0, y - 1.0, lab, fontsize=PS.XLABEL_FS, rotation=PS.XLABEL_ROT,
                ha="right", va="top", rotation_mode="anchor")
    over = [(x, h) for x in shown if (h := PS.rotated_label_mm(x) + BAND_GAP_MM + BAND_MM + 1.0) > PS.MARGIN_B]
    if over:
        log(f"  [注意] 群名が長く、下余白 {PS.MARGIN_B} mm で見切れている可能性があります (群名を短くしてください): "
            + ", ".join(f"{x} (約 {h:.1f} mm)" for x, h in over))
    cb_len = min(PS.PLOT_H_MM, ax_h)          # #8・#11 と同じ: 長さ 48 mm・上端を行列の上端にそろえる
    for x0, norm, cmap, label, ticks in ((PS.MARGIN_L, norm_ab, cmap_ab, "Relative abundance (%)",
                                          range(HEAT_VMIN, HEAT_VMAX + 1)),
                                         (cb2_x, norm_fc, cmap_fc, "Log fold change", [-vmax, -step, 0, step, vmax])):
        ext = cb_len * 0.05 if x0 == PS.MARGIN_L else 0.0   # 存在量のバーだけ、三角を本体の下に足す (261019 決定・案 b)
        cax = fig.add_axes([x0 / fig_w, (PS.MARGIN_B + ax_h - cb_len - ext) / fig_h, cb_w / fig_w, (cb_len + ext) / fig_h])
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                          extend="min" if x0 == PS.MARGIN_L else "neither")   # 存在量だけ下端に三角 (261016 決定)
        cax.yaxis.set_ticks_position("left")
        cax.yaxis.set_label_position("left")
        cb.outline.set_visible(False)
        cb.set_label(label, fontsize=PS.FS_CBAR)
        cb.ax.tick_params(labelsize=PS.FS_CBAR, length=2)
        cb.set_ticks(list(ticks))
        if x0 == PS.MARGIN_L:                     # 存在量の目盛は % (261016 決定・E6d)
            cb.set_ticklabels(_pct_ticklabels())
    fig.savefig(_fig_dir(d, name) / f"{name}.pdf")
    matplotlib.pyplot.close(fig)
    _ancom_legend(_fig_dir(d, name) / f"{name}_legend.pdf")
    return fig_w, fig_h, sub

def fig13_ancombc(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#11 ANCOM-BC バブルプロット (階級ごとに 1 枚)。計算は R ANCOMBC の ancombc2 (Dunnett 型・lib/ancombc2.R。261018 に scikit-bio の ancombc から移した)。
    列 = 比較 (対照群 vs 各群・群名だけ)、行 = タクソン、色 = ln fold change (Tol BuRd・0 が白)、
    円の大きさ = q の 3 段階 (mdFDR 補正後)。有意 = diff_robust (感度解析を通ったもの)。有意でないものは灰色の小さな点。
    載せるタクソン = 比較から交互に 1 つずつ、有意 (diff_robust) の中から |ln fold change| の大きい順に、全体で ancombc_top_n 個。
    行の並び = 選ばれた比較ごとにまとめ、その中は増加 (大きい順) → 減少 (その比較の値で並べる)。"""
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    log, cfg = ctx["log"], ctx["cfg"]
    control = _ancom_control(ctx, G)
    if not control:
        log.step("10", n, total, True, "#11 ANCOM-BC",
                 "comparison.control_group が none なので描かない")
        return
    _, ranks = _plot_cfg(ctx)
    top_n = _top_n_or_all(ctx, "ancombc_top_n")
    cmap = LinearSegmentedColormap.from_list("BuRd", BURD)
    for k, rank in enumerate(ranks, 1):
        name = f"11_ancombc_{rank}"
        log.step("10", n, total, False, f"#11 ANCOM-BC ({k}/{len(ranks)})", f"{name}.pdf")
        run_ = _ancom_run(ctx, G, rank, control)
        if run_ is None:                               # 対照群が n < stats.min_n（261009 決定・D2e）
            log(f"  [注意] {rank}: 対照群 {control} の検体数が stats.min_n ({G.min_n}) 未満のため、#11 は描きません")
            continue
        res, comps, n_all, n_test, small = run_
        # 計算から外した群（n < stats.min_n）も列としては残し、帯の中に「n < 3」と書く（261009 決定）
        pcomps = [g for g in G.order if g != control and (g in comps or g in small)]
        if small and k == 1:
            log(f"  [注記] n < {G.min_n} のため ANCOM-BC の計算から外した群があります (列は残し、帯に n < {G.min_n} と書きます): "
                + ", ".join(f"{g} (n = {m})" for g, m in small.items()))
        # 載せるタクソン = 比較から交互に 1 つずつ、|ln fold change| の大きい順に、全体で top_n 個（260924 決定・案 B）。
        # 行数を決まった数（X01 は 12 = 軸高 48 mm ÷ 4 mm）にして、図の高さを他の図とそろえる。all なら有意の属を全部 (261021)
        queue = {}
        for g in comps:
            r = res[(res.group == g) & res.sig & ~res.struc_zero]
            queue[g] = list(r.reindex(r.ln_fc.abs().sort_values(ascending=False).index)["taxon"])
        picked, chosen_by = [], {}
        while len(picked) < top_n and any(queue.values()):
            for g in comps:
                while queue[g] and queue[g][0] in picked:
                    queue[g].pop(0)
                if queue[g] and len(picked) < top_n:
                    picked.append(queue[g].pop(0))
                    chosen_by[picked[-1]] = g              # どの比較で選ばれたか（行の並びに使う）
        pick = set(picked)
        sub = res[res.taxon.isin(pick)]
        L = sub.pivot(index="taxon", columns="group", values="ln_fc").reindex(columns=pcomps)
        Q = sub.pivot(index="taxon", columns="group", values="qvalue").reindex(columns=pcomps)
        rows = []                                          # 選ばれた比較ごとにまとめ、その中は増加 → 減少（260924 指示）
        for g in comps:
            mem = [tx for tx in L.index if chosen_by[tx] == g]
            up = sorted([tx for tx in mem if L.loc[tx, g] > 0], key=lambda tx: -L.loc[tx, g])
            dn = sorted([tx for tx in mem if L.loc[tx, g] <= 0], key=lambda tx: L.loc[tx, g])
            rows += up + dn
        cd = _csv_dir(d, name)
        _ancom_all_csv(res, control, cd / f"{name}.all.csv")
        if not rows:
            log(f"  [注記] {rank}: 有意 (diff_robust) のタクソンが無いため、図は描きません")   # 261019・F014・F050・F189
            continue
        ra, _, _ = _taxa_relative(ctx, rank)
        base_mean = ra.loc[:, ra.columns.map(G.label_by_newid) == control].mean(axis=1) / 100
        fig_w, fig_h, sub = _ancom_draw(G, d, name, rank, res, comps, small, control, rows, base_mean, log)
        sub.merge((base_mean * 100).rename(f"mean_relative_abundance_{control}").rename_axis("taxon").reset_index(), on="taxon") \
            .to_csv(cd / f"{name}.plotted.csv", index=False)
        nsig = {g: int(((res.group == g) & res.sig).sum()) for g in comps}
        ndiff = {g: int(((res.group == g) & res["diff"]).sum()) for g in comps}
        log.sum(f"{rank}: 有意 (diff_robust): " + ", ".join(f"{g} {c}" for g, c in nsig.items()))
        log(f"  {rank}: 検定 {n_test} / {n_all} タクソン / 基準 {control} / 有意 (diff_robust): "
            + ", ".join(f"{g} {c}" for g, c in nsig.items()) + " (感度解析の前: " + ", ".join(f"{g} {c}" for g, c in ndiff.items())
            + f")/ 図に {len(rows)} 行 (比較から交互に選んで全体で {cfg['plotting']['ancombc_top_n']})/ 図 {fig_w:.0f} × {fig_h:.0f} mm")
    log(f"  {_ancom_method_text(ctx)}")


def fig14_ancombc_summary(ctx: dict, G, d: Path, n: int, total: int) -> None:
    """#12 ANCOM-BC まとめ (階級ごとに 1 枚。261007 追加)。#11 と同じ ANCOM-BC2 の結果に、対照群での存在量を 1 列のヒートマップで足す。
    載せるタクソン = どれかの比較で有意 (diff_robust) のもの。上から対照群の平均相対存在量の多い順に、ancombc_summary_top_n 個。
    左から: 存在量のカラーバー / ln fold change のカラーバー / バブル (列 = 対照群 vs 各群。#11 と同じ描き方)/
    対照群の平均相対存在量 1 列 (log10(x + 1e-6)・範囲 −3〜0・viridis の 0〜95%。#8 と同じ目盛り)/ タクソン名。"""
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    log, cfg = ctx["log"], ctx["cfg"]
    control = _ancom_control(ctx, G)
    if not control:
        log.step("10", n, total, True, "#12 ANCOM-BC まとめ",
                 "comparison.control_group が none なので描かない")
        return
    _, ranks = _plot_cfg(ctx)
    raw_n = cfg["plotting"]["ancombc_summary_top_n"]
    top_n = _top_n_or_all(ctx, "ancombc_summary_top_n")   # all = 有意 (diff_robust) のタクソンを全部載せる（261007 追加）
    top_n = None if top_n == math.inf else top_n
    cmap_fc = LinearSegmentedColormap.from_list("BuRd", BURD)
    cmap_ab = PS.viridis()                             # #8 と同じ（viridis の 0〜95%）
    norm_ab = Normalize(HEAT_VMIN, HEAT_VMAX)
    for k, rank in enumerate(ranks, 1):
        name = f"12_ancombc-summary_{rank}"
        log.step("10", n, total, False, f"#12 ANCOM-BC まとめ ({k}/{len(ranks)})", f"{name}.pdf")
        run_ = _ancom_run(ctx, G, rank, control)
        if run_ is None:                               # 対照群が n < stats.min_n（261009 決定・D2e）
            log(f"  [注意] {rank}: 対照群 {control} の検体数が stats.min_n ({G.min_n}) 未満のため、#12 は描きません")
            continue
        res, comps, n_all, n_test, small = run_
        if small and k == 1:
            log(f"  [注記] n < {G.min_n} のため ANCOM-BC の計算から外した群があります (列は残し、帯に n < {G.min_n} と書きます): "
                + ", ".join(f"{g} (n = {m})" for g, m in small.items()))
        ra, _, _ = _taxa_relative(ctx, rank)
        grp = ra.columns.map(G.label_by_newid)
        base_mean = ra.loc[:, grp == control].mean(axis=1) / 100            # 対照群の平均相対存在量（割合）
        sig = sorted(set(res.loc[res.sig & ~res.struc_zero, "taxon"]))
        cd = _csv_dir(d, name)
        _ancom_all_csv(res, control, cd / f"{name}.all.csv")    # #11 を描かない設定でも計算条件と検体数を残す (261021・F2)
        if not sig:
            log(f"  [注記] {rank}: 有意 (diff_robust) のタクソンが無いため、図は描きません")
            continue
        rows = sorted(sig, key=lambda tx: -base_mean.get(tx, 0.0))[:top_n]  # 対照群で多い順（261007 決定）。None = all
        fig_w, fig_h, sub = _ancom_draw(G, d, name, rank, res, comps, small, control, rows, base_mean, log)
        sub.merge((base_mean * 100).rename(f"mean_relative_abundance_{control}").rename_axis("taxon").reset_index(), on="taxon") \
            .to_csv(cd / f"{name}.plotted.csv", index=False)
        log(f"  {rank}: 有意 (diff_robust) のタクソン {len(sig)} 個 (検定 {n_test} / {n_all})のうち、{control} で多い順に "
            f"{len(rows)} 行 (ancombc_summary_top_n = {raw_n})/ 図 {fig_w:.0f} × {fig_h:.0f} mm")
        # 増加・減少の内訳（261009 決定・D9 a）: 対照群で多い属から選ぶので、絞ると主に減少が載る
        hit = res[res.sig & ~res.struc_zero]
        shown = hit[hit.taxon.isin(rows)]
        cnt = {"scope": ["all_significant", "plotted"], "taxa": [len(sig), len(rows)],
               "increase": [int((hit.ln_fc > 0).sum()), int((shown.ln_fc > 0).sum())],
               "decrease": [int((hit.ln_fc < 0).sum()), int((shown.ln_fc < 0).sum())]}
        pd.DataFrame(cnt).to_csv(cd / f"{name}.direction.csv", index=False)
        log(f"  [注記] {rank}: 有意 (diff_robust) の組 (タクソン × 比較) は全体で増加 {cnt['increase'][0]}・減少 {cnt['decrease'][0]} で、"
            f"図に載せたのは増加 {cnt['increase'][1]}・減少 {cnt['decrease'][1]} です ({control} で多い属から選ぶため、絞ると主に減少が載ります)")
    log(f"  {_ancom_method_text(ctx)}")


def _ancom_method_text(ctx: dict) -> str:
    """ログに残す計算条件 (261018・BC2)。"""
    v = next((r[0].attrs.get("ancombc_version", "?") for r in ctx.get("_ancom_cache", {}).values() if r), "?")
    return (f"ANCOM-BC2: R ANCOMBC {v}・Dunnett 型 (対照群 vs 各群。2 群のときは主解析で mdFDR なし)・疑似カウントなし (感度解析 0.1・0.5・1 を通った有意 = diff_robust)・"
            f"構造的ゼロは全検体 0 の群・BH・mdFDR holm (B = {ANCOM_B})・出現率で外さない・乱数の種 {ANCOM_SEED}・ln fold change (自然対数)")


def _ancom_all_csv(res, control: str, path: Path) -> None:
    """全タクソン × 比較の結果 (261018: BC2 の列に変更)。"""
    v = res.attrs.get("ancombc_version", "")
    two = res.group.nunique() == 1                      # 2 群のときは主解析 (mdFDR なし。261020 決定)
    res.assign(control=control, log_base="natural (ln)",
               method="ANCOM-BC2 primary (R ANCOMBC)" if two else "ANCOM-BC2 Dunnett (R ANCOMBC)", ancombc_version=v,
               pseudo=0, sensitivity="0.1, 0.5, 1", p_adjust="BH", mdfdr="none" if two else f"holm (B = {ANCOM_B})", prv_cut=ANCOM_PRV,
               neg_lb=False, seed=ANCOM_SEED)[
        ["taxon", "control", "group", "n_control", "n_group", "ln_fc", "log_base", "SE", "W", "pvalue", "qvalue", "diff", "passed_ss", "sig",
         "not_tested", "struc_zero", "zero_in_control", "zero_in_group", "method", "ancombc_version", "pseudo", "sensitivity",
         "p_adjust", "mdfdr", "prv_cut", "neg_lb", "seed"]].to_csv(path, index=False)


def _struc_mark(ax, x, y, r) -> None:
    """構造的ゼロの組の印 (261009・D5a。検定はせず、片方の群で検出されないことを示す)。
    白抜きの円: 比べる群で未検出 = 減少の色 (BuRd の青端)、対照群で未検出 = 増加の色 (赤端)、両方 = 灰色 (261009 決定・案 ②)。
    構造的ゼロの属は #11・#12 の行の候補にしない (261009 決定・案 a。q が擬似カウントで決まるため)。"""
    c = "#BBBBBB" if (r.zero_in_group and r.zero_in_control) else (BURD[0] if r.zero_in_group else BURD[-1])
    ax.scatter(x, y, s=ANCOM_SIZES[1][2], facecolors="none", edgecolors=c, linewidths=1.0, zorder=2)


def _ancom_legend(path: Path) -> None:
    """#11 の円の大きさの凡例 (別ファイル・高さは他の図と同じ)。"""
    from matplotlib.lines import Line2D
    fig_h = PS.MARGIN_B + PS.PLOT_H_MM + PS.MARGIN_T
    h = [Line2D([], [], ls="", marker="o", ms=np.sqrt(sz), markerfacecolor="#BBBBBB", markeredgecolor=PS.GREY,
                markeredgewidth=0.3, label=lab) for _, lab, sz in ANCOM_SIZES]
    h.append(Line2D([], [], ls="", marker="o", ms=np.sqrt(ANCOM_NS_S), color="#BBBBBB", label="n.s."))
    # 構造的ゼロの白抜きの見本は載せない (261019 決定・案 a。ANCOM-BC2 ではどれか 1 群でも構造的ゼロの属は検定されず
    # 行に選ばれないので、図に白抜きは出ない。描く処理 _struc_mark は念のため残す)
    fig = matplotlib.pyplot.figure(figsize=(60 * PS.MM, fig_h * PS.MM))
    lg = fig.legend(handles=h, fontsize=PS.FS_LEGEND, frameon=False, loc="upper left", title="q (mdFDR)",
                    title_fontsize=PS.FS_LEGEND, alignment="left", borderaxespad=0, borderpad=0,
                    bbox_to_anchor=(2 / 60, 1 - PS.MARGIN_T / fig_h), labelspacing=0.6)
    fig.canvas.draw()
    w = lg.get_window_extent().width / fig.dpi * 25.4 + 4
    fig.set_size_inches(w * PS.MM, fig_h * PS.MM)
    lg.set_bbox_to_anchor((2 / w, 1 - PS.MARGIN_T / fig_h), transform=fig.transFigure)
    fig.savefig(path)
    matplotlib.pyplot.close(fig)
