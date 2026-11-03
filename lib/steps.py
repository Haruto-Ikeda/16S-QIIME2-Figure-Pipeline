#!/usr/bin/env python3
"""steps.py — 各ステップ本体 (01〜09)。番号・出力名はアーカイブ準拠。

各関数は (ctx) を受け取る。ctx は run_16s.py が組み立てる dict:
  exp, base, out, cfg, log
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import yaml

from . import metadata as MD
from . import conditions as C
from . import qiime as Q


def threads_arg(cfg: dict, style: str) -> str:
    """config の `threads` (0 = 全コア)を各コマンドの流儀に変換する。
    style="zero": 0 が全コア (dada2 / classify-sklearn の n_jobs)
    style="auto": 'auto' が全コア (phylogeny / core-metrics)"""
    t = cfg["threads"]                     # 必須 (check_config で確かめ済み。261021)
    if t <= 0:
        return "0" if style == "zero" else "auto"
    return str(t)


RUN_DIR_RE = re.compile(r"^\d{6}_(M\d{4})$")          # ランのフォルダ名 <YYMMDD>_M<4 桁> (261021 本人決定・1-2)
ID_RE = re.compile(r"^M(\d{4})P(\d{4})N(\d{4})$")      # 検体の ID (261021 本人決定・②)
N_RE = re.compile(r"^N(\d{4})$")                        # absent_ok の書き方 (その CSV の中の No.。261022 本人決定で各 CSV の中へ)
CSV_NAME_RE = re.compile(r"^\d{6}_(\d{4})_.+\.csv$")   # 検体情報 CSV のファイル名 <YYMMDD>_<P 4 桁>_…csv (3-4)
RUN_KEYS = ("threads", "runs", "sample_info", "diversity", "classifier", "taxonomy")   # 解析し直しの形 (ア) だけで使うキー (絞り込みでは none)


def _die(msg: str) -> None:
    raise SystemExit(f"  [エラー] {msg}")


def abs_path(raw, key: str) -> Path:
    """config のパス。~ で始まるパスか絶対パスだけ受け付け、相対パスは止める (261021 本人決定・ウ 案 A)。"""
    if not isinstance(raw, str) or not raw.strip():
        _die(f"config の {key} はパスで書いてください (現在: {raw!r})")
    p = Path(raw.strip()).expanduser()
    if not p.is_absolute():
        _die(f"config の {key} は ~ で始まるパスか絶対パスで書いてください (現在: {raw!r})")
    return p


def is_none(v) -> bool:
    return isinstance(v, str) and v.strip().lower() == "none"


def mode(cfg: dict) -> str:
    """config の形 (261021 本人決定・7-1。261022 に雛形を 1 つにまとめた)。
    runs: と parent: はどちらも必須で、使わない方は none と書く。runs: = 解析し直し (ア。01〜10)、parent: = 絞り込み (イ。10 だけ)。"""
    for k in ("runs", "parent", "select"):
        if k not in cfg:
            _die(f"config の {k} がありません (使わないときは none)")
    use_r, use_p = not is_none(cfg["runs"]), not is_none(cfg["parent"])
    if use_r == use_p:
        _die("config の runs: (FASTQ から解析する) と parent: (親の実験の結果を絞って描く) は、どちらか一方だけを書き、もう一方は none にしてください"
             + (" (両方に値があります)" if use_r else " (両方 none です)"))
    return "runs" if use_r else "parent"


def classifier_path(cfg: dict) -> Path:
    """分類器のパス (261015 本人決定・E4a)。~ で始まるパスか絶対パスだけ受け付け、相対パスは止める。"""
    p = abs_path((cfg.get("classifier") or {}).get("path"), "classifier.path")
    if not p.exists():                     # 01 の前に止める (261021・C7。以前は 01〜07 のあとの 08 で初めて確かめていた)
        guide = Path(__file__).resolve().parent.parent / "MANUAL.md"   # パイプラインの置き場所から計算 (261019 本人決定・移植性。261027 に手引きの §2 へ)
        raise SystemExit(f"  [エラー] 分類器がありません: {p}\n  ※ {guide} の「2. 分類器」を見て、置き場所と config の classifier.path を確かめてください")
    return p


def rarefaction_seed(cfg: dict) -> int:
    """希薄化 (07 core-metrics)の乱数の種 (261009 本人決定・D11a)。config に明示が無ければ止める。"""
    seed = (cfg.get("diversity") or {}).get("random_seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise SystemExit(f"  [エラー] config の diversity.random_seed を 0 以上の整数で書いてください (現在: {seed!r})")
    return seed


def _none_or_list(v, key: str, pat: re.Pattern, form: str) -> list[str]:
    if isinstance(v, str) and v.strip().lower() == "none":
        return []
    if not isinstance(v, list) or not v:
        _die(f"config の {key} は {form} のリストか none で書いてください (現在: {v!r})")
    bad = [x for x in v if not (isinstance(x, str) and pat.fullmatch(x.strip()))]
    if bad:
        _die(f"config の {key} に形の違う値があります: {', '.join(map(str, bad))} ({form} の形で書いてください)")
    dup = sorted({x for x in v if v.count(x) > 1})
    if dup:
        _die(f"config の {key} に同じ値が 2 回以上あります: {', '.join(dup)}")
    return [x.strip() for x in v]


def runs_of(cfg: dict) -> list[dict]:
    """config の runs: を確かめ、[{name, m, path, fastq, sheet, qc, exclude}] を返す (261021 本人決定・①⑥)。"""
    runs = cfg.get("runs")
    if not isinstance(runs, list) or not runs:
        _die(f"config の runs: は使うランのリスト (1 つ以上) で書いてください (現在: {runs!r})")
    out, seen = [], {}
    for i, r in enumerate(runs, 1):
        key = f"runs[{i}]"
        if not isinstance(r, dict):
            _die(f"config の {key} は path・quality_control・exclude を書いてください (現在: {r!r})")
        for k in ("path", "quality_control", "exclude"):
            if k not in r:
                _die(f"config の {key} に {k} がありません")
        p = abs_path(r["path"], f"{key}.path")
        m = RUN_DIR_RE.fullmatch(p.name)
        if not m:
            _die(f"ランのフォルダ名は <YYMMDD>_M<4 桁> の形にしてください (例 250101_M0001。現在: {p.name})")
        if not (p / "fastq").is_dir() or not (p / "SampleSheet.txt").is_file():
            _die(f"ランのフォルダに fastq/ と SampleSheet.txt が必要です: {p}")
        if m.group(1) in seen:
            _die(f"同じラン {m.group(1)} が runs: に 2 回あります: {seen[m.group(1)]} と {p}")
        seen[m.group(1)] = p
        qc = r["quality_control"]
        if not isinstance(qc, dict):
            _die(f"config の {key}.quality_control は trim_left_f・trim_left_r・trunc_len_f・trunc_len_r を書いてください")
        for k in ("trim_left_f", "trim_left_r", "trunc_len_f", "trunc_len_r"):
            v = qc.get(k)
            if k not in qc or isinstance(v, bool) or not isinstance(v, int) or v < 0:
                _die(f"config の {key}.quality_control.{k} は 0 以上の整数で書いてください (現在: {v!r})")
        exc = _none_or_list(r["exclude"], f"{key}.exclude", ID_RE, "検体の ID (M0001P0001N0001)")
        other = [x for x in exc if ID_RE.fullmatch(x).group(1) != m.group(1)[1:]]
        if other:
            _die(f"config の {key}.exclude にほかのランの ID があります: {', '.join(other)} (このランは {m.group(1)})")
        out.append(dict(name=p.name, m=m.group(1), path=p, fastq=p / "fastq", sheet=p / "SampleSheet.txt",
                        qc={k: qc[k] for k in ("trim_left_f", "trim_left_r", "trunc_len_f", "trunc_len_r")}, exclude=exc))
    return out


def sample_info_of(cfg: dict) -> dict[int, dict]:
    """config の sample_info: (CSV のリスト) を確かめ、{P 番号: {path, absent_ok (No. の集合)}} を返す (261021 本人決定・③)。
    各 CSV に path と absent_ok (その CSV の検体のうち FASTQ が無くてよい No.。N0127 の形か none) を書く (261022 本人決定)。
    中身の照合 (2 行目の Project No.・No. の重複など) は 01 で行う。"""
    v = cfg.get("sample_info")
    if isinstance(v, dict):
        _die("config の sample_info: はリストにして、各 CSV に path と absent_ok を書いてください (261022 に形を変更。1 つでもリスト)")
    if not isinstance(v, list) or not v:
        _die(f"config の sample_info: は検体情報 CSV のリスト (1 つ以上) で書いてください (現在: {v!r})")
    out = {}
    for i, e in enumerate(v, 1):
        key = f"sample_info[{i}]"
        if not isinstance(e, dict) or "path" not in e or "absent_ok" not in e:
            _die(f"config の {key} は path と absent_ok を書いてください (absent_ok は FASTQ が無くてよい No. のリストか none。現在: {e!r})")
        p = abs_path(e["path"], f"{key}.path")
        if not p.is_file():
            _die(f"検体情報 CSV がありません: {p}")
        m = CSV_NAME_RE.fullmatch(p.name)
        if not m:
            _die(f"検体情報 CSV のファイル名は <YYMMDD>_<P 番号 4 桁>_…csv の形にしてください (例 260709_0503_Sample information document.csv。現在: {p.name})")
        pn = int(m.group(1))
        if pn in out:
            _die(f"P 番号 {m.group(1)} の検体情報 CSV が 2 つあります: {out[pn]['path'].name} と {p.name}")
        absent = _none_or_list(e["absent_ok"], f"{key}.absent_ok", N_RE, "No. (N0127)")
        out[pn] = dict(path=p, absent_ok={int(N_RE.fullmatch(x).group(1)) for x in absent})
    return out


def parent_of(cfg: dict) -> Path:
    """絞り込み (イ) の親の実験 (261021 本人決定・⑦)。親は 01〜09 を行った実験だけ (入れ子にしない)。"""
    p = abs_path(cfg.get("parent"), "parent")
    q = p / "output/16S/Qiime2"
    if not (p / "input/16S/config.yaml").is_file() or not (q / "07_diversity").is_dir() or not (q / "09_taxonomy").is_dir():
        _die(f"親の実験に 01〜09 の結果がありません: {p} (先に親の実験を解析してください)")
    pcfg = yaml.safe_load((p / "input/16S/config.yaml").read_text(encoding="utf-8")) or {}
    if not is_none(pcfg.get("parent", "none")):
        _die(f"親の実験も絞り込みの形 (parent:) です: {p} (親は runs: で解析した実験だけにしてください)")
    return p


def with_parent(cfg: dict) -> dict:
    """絞り込み (イ) では、作図の検査に親の taxonomy (collapse_ranks) を使う (書かないキーを親から借りるだけ。ファイルは変えない)。"""
    if is_none(cfg.get("parent", "none")):
        return cfg
    pcfg = yaml.safe_load((parent_of(cfg) / "input/16S/config.yaml").read_text(encoding="utf-8")) or {}
    return {**cfg, "taxonomy": pcfg.get("taxonomy")}


def check_config(cfg: dict) -> str:
    """01〜09 の config を 01 の前に確かめる。すべて必須で、書いていなければ止める (261021 本人決定・C4/C6)。返り値: 形 (runs / parent)。"""
    if "quality_control" in cfg:
        _die("config の quality_control は、runs: の各ランの中へ移しました (261021。ランごとに書きます)")
    if "absent_ok" in cfg:
        _die("config の absent_ok は、sample_info: の各 CSV の中へ移しました (261022。その CSV の No. を N0127 の形で書きます)")
    md = mode(cfg)
    if md == "parent":
        miss = [k for k in RUN_KEYS if k != "runs" and k not in cfg]
        if miss:
            _die(f"config の {', '.join(miss)} がありません (絞り込みの形では none と書きます)")
        extra = [k for k in RUN_KEYS if k != "runs" and not is_none(cfg[k])]
        if extra:
            _die(f"絞り込みの形 (parent:) では {', '.join(extra)} は none にしてください (親の実験の設定を使います)")
        parent_of(cfg)
        sel = cfg["select"]
        if not isinstance(sel, dict):
            _die("config の select: に Group_name・Subgroup_name・exclude を書いてください")
        for k in ("Group_name", "Subgroup_name"):
            v = sel.get(k)
            if k not in sel or not ((isinstance(v, str) and v.strip().lower() == "all") or (isinstance(v, list) and v)):
                _die(f"config の select.{k} は使う値のリストか all で書いてください (現在: {v!r})")
        if "exclude" not in sel:
            _die("config の select.exclude がありません (外す検体の ID のリストか none)")
        _none_or_list(sel["exclude"], "select.exclude", ID_RE, "検体の ID (M0001P0001N0001)")
        return md

    if not is_none(cfg["select"]):
        _die("config の select: は、FASTQ から解析する形 (runs:) では none にしてください (絞り込みのときだけ使います)")

    def need(sec: str, key: str):
        d = cfg.get(sec) if sec else cfg
        if not isinstance(d, dict) or key not in d:
            raise SystemExit(f"  [エラー] config の {sec + '.' if sec else ''}{key} がありません")
        return d[key]

    def is_int(v, lo: int) -> bool:
        return isinstance(v, int) and not isinstance(v, bool) and v >= lo

    v = need("", "threads")
    if not is_int(v, 0):
        raise SystemExit(f"  [エラー] config の threads は 0 以上の整数で書いてください (0 = 全コア。現在: {v!r})")
    runs_of(cfg)
    sample_info_of(cfg)
    v = need("diversity", "sampling_depth_floor")
    if not is_int(v, 1):
        raise SystemExit(f"  [エラー] config の diversity.sampling_depth_floor は正の整数で書いてください (現在: {v!r})")
    rarefaction_seed(cfg)
    classifier_path(cfg)
    v = need("classifier", "n_jobs")
    if not is_int(v, 1):
        raise SystemExit(f"  [エラー] config の classifier.n_jobs は正の整数で書いてください (現在: {v!r})")
    v = need("taxonomy", "confidence")
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not 0 <= v <= 1:
        raise SystemExit(f"  [エラー] config の taxonomy.confidence は 0〜1 の数で書いてください (現在: {v!r})")
    v = need("taxonomy", "collapse_ranks")
    if not isinstance(v, list) or not v:
        raise SystemExit(f"  [エラー] config の taxonomy.collapse_ranks は階級名のリスト (1 つ以上) で書いてください (現在: {v!r})")
    bad = [r for r in v if str(r).strip().lower() not in RANK_PREFIX]
    if bad:
        raise SystemExit(f"  [エラー] config の taxonomy.collapse_ranks に知らない階級名があります: {bad} (使える名前: {', '.join(RANK_PREFIX)})")
    return md


# ============================================================ 02 Import
def _runs_sorted(cfg: dict) -> list[dict]:
    """ランを M の順に (6-9: 各ランの表は ID の順に並ぶので、M の順につなげば全体が ID の順になる)。"""
    return sorted(runs_of(cfg), key=lambda r: r["m"])


def step02_import(ctx: dict) -> None:
    """ランごとに FASTQ を読み込む (261022・6-6: manifest 形式。exclude した検体は読み込まない)。
    FASTQ の名前・SampleSheet・検体情報の照合は 01 で済ませてある。"""
    log, out, cfg = ctx["log"], ctx["out"], ctx["cfg"]
    log.head("02. Import (データの読み込み)")
    C.check(ctx, "02", 1, 4, "manifest の生成", "<ラン>/manifest.tsv, demux.qza")   # 条件の照合（261014 E1a）
    md = pd.read_csv(out / "01_metadata" / "sample-metadata.txt", sep="\t", dtype=str)
    for r in _runs_sorted(cfg):
        log(f"[{r['name']}]")
        d = out / "02_import" / r["m"]
        d.mkdir(parents=True, exist_ok=True)   # ステップ冒頭で明示（アーカイブ準拠）
        manifest, demux_qza, demux_qzv = d / "manifest.tsv", d / "demux.qza", d / "demux.qzv"

        # 1/4 manifest (01 の表に残った検体 = exclude を除いた検体)
        if manifest.exists():
            log.step("02", 1, 4, True, "manifest の生成", "既存) manifest.tsv")
        else:
            log.step("02", 1, 4, False, "manifest の生成", "manifest.tsv")
            fq = MD.fastq_ids(r["fastq"])
            ids = md.loc[md["Run"] == r["m"], "#SampleID"].tolist()
            manifest.write_text("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n"
                                + "".join(f"{i}\t{fq[i]['R1']}\t{fq[i]['R2']}\n" for i in ids), encoding="utf-8")
            log(f"  {len(ids)} 検体" + (f" (exclude で外した {len(r['exclude'])} 検体を除く)" if r["exclude"] else ""))

        # 2/4 import（公式 how-to-import の PairedEndFastqManifestPhred33V2）
        Q.run(log, "02", 2, 4, "インポートファイルの生成", [demux_qza], [
            "tools", "import",
            "--type", "SampleData[PairedEndSequencesWithQuality]",
            "--input-path", str(manifest),
            "--input-format", "PairedEndFastqManifestPhred33V2",
            "--output-path", str(demux_qza),
        ])

        # 3/4 demux.qzv
        Q.run(log, "02", 3, 4, "可視化インポートファイルの生成", [demux_qzv], [
            "demux", "summarize",
            "--i-data", str(demux_qza),
            "--o-visualization", str(demux_qzv),
        ])

        # 4/4 HTML 化 (ランごとに名前を分ける)
        Q.export_qzv(log, "02", 4, 4, "インポートファイルの HTML 化", demux_qzv, out / f"HTML_02_demux_{r['m']}.qzv")


# ============================================================ 03 DADA2
def step03_dada2(ctx: dict) -> None:
    """QC・マージ (denoise) をランごとに (6-1)。ランが 2 つ以上なら 03_merge/ でつなげる。出力名はアーカイブ準拠。
    2026.7 の変更 3 点を反映:
    ① `--o-base-transition-stats` が必須出力に加わった (公式＝マニュアル)
    ② `feature-table summarize` は `--o-summary` ＋ 頻度 2 出力 ＋ `--m-metadata-file` (公式＝マニュアル)
    ③ `tabulate-seqs` に `--m-metadata-file feature-frequencies.qza` を渡すと ASV ごとの本数が並ぶ (公式のみ)
    """
    log, out, cfg = ctx["log"], ctx["out"], ctx["cfg"]
    log.head("03. DADA2 (データの品質管理)")
    C.check(ctx, "03", 1, 7, "DADA2 実行結果ファイルの生成", "<ラン>/table.qza, rep-seqs.qza, denoising-stats.qza, base-transition-stats.qza")   # 条件の照合（261014 E1a）
    threads = threads_arg(cfg, "zero")
    md = out / "01_metadata" / "sample-metadata.txt"
    runs = _runs_sorted(cfg)
    for r in runs:
        log(f"[{r['name']}]")
        qc, m = r["qc"], r["m"]
        d = out / "03_dada2" / m
        d.mkdir(parents=True, exist_ok=True)   # ステップ冒頭で明示（アーカイブ準拠）
        demux = out / "02_import" / m / "demux.qza"
        table, seqs = d / "table.qza", d / "rep-seqs.qza"
        stats, trans = d / "denoising-stats.qza", d / "base-transition-stats.qza"

        # 1/7 denoise（公式 gut-to-soil ＝ マニュアル §2。出力 4 つ）
        Q.run(log, "03", 1, 7, "DADA2 実行結果ファイルの生成", [table, seqs, stats, trans], [
            "dada2", "denoise-paired",
            "--i-demultiplexed-seqs", str(demux),
            "--p-trim-left-f", str(qc["trim_left_f"]),
            "--p-trim-left-r", str(qc["trim_left_r"]),
            "--p-trunc-len-f", str(qc["trunc_len_f"]),
            "--p-trunc-len-r", str(qc["trunc_len_r"]),
            "--p-n-threads", threads,
            "--o-table", str(table),
            "--o-representative-sequences", str(seqs),
            "--o-denoising-stats", str(stats),
            "--o-base-transition-stats", str(trans),
        ])

        # 2/7 table.qzv（2026.7 形式: --o-summary ＋ 頻度 2 出力）
        table_qzv = d / "table.qzv"
        sfreq, ffreq = d / "sample-frequencies.qza", d / "feature-frequencies.qza"
        Q.run(log, "03", 2, 7, "リード数ファイルの可視化", [table_qzv, sfreq, ffreq], [
            "feature-table", "summarize",
            "--i-table", str(table),
            "--m-metadata-file", str(md),
            "--o-summary", str(table_qzv),
            "--o-sample-frequencies", str(sfreq),
            "--o-feature-frequencies", str(ffreq),
        ])

        # 3/7 HTML 化
        Q.export_qzv(log, "03", 3, 7, "リード数ファイルの HTML 化", table_qzv, out / f"HTML_03_table_{m}.qzv")

        # 4/7 rep-seqs.qzv（公式: feature-frequencies を流し込むと配列ごとの本数が並ぶ）
        seqs_qzv = d / "rep-seqs.qzv"
        Q.run(log, "03", 4, 7, "代表配列ファイルの可視化", [seqs_qzv], [
            "feature-table", "tabulate-seqs",
            "--i-data", str(seqs),
            "--m-metadata-file", str(ffreq),
            "--o-visualization", str(seqs_qzv),
        ])

        # 5/7 HTML 化
        Q.export_qzv(log, "03", 5, 7, "代表配列ファイルの HTML 化", seqs_qzv, out / f"HTML_03_rep-seqs_{m}.qzv")

        # 6/7 denoising-stats.qzv
        stats_qzv = d / "denoising-stats.qzv"
        Q.run(log, "03", 6, 7, "品質統計ファイルの可視化", [stats_qzv], [
            "metadata", "tabulate",
            "--m-input-file", str(stats),
            "--o-visualization", str(stats_qzv),
        ])

        # 7/7 denoising-stats.tsv（export は stats.tsv を吐くので改名する。アーカイブ準拠）
        stats_tsv = d / "denoising-stats.tsv"
        if not stats_tsv.exists():
            log.step("03", 7, 7, False, "品質統計ファイルの .tsv 化", stats_tsv.name)   # 改名後の名前を出す (261019・F258)
            Q.run(log, "03", 7, 7, "品質統計ファイルの .tsv 化", [d / "stats.tsv"], [
                "tools", "export", "--input-path", str(stats), "--output-path", str(d),
            ], force=True, quiet_head=True)
            (d / "stats.tsv").replace(stats_tsv)
        else:
            log.step("03", 7, 7, True, "品質統計ファイルの .tsv 化", "既存) denoising-stats.tsv")

    if len(runs) == 1:
        log("")
        log(f"ランが 1 つ ({runs[0]['m']}) のため、結合は行いません (04 以降は 03_dada2/{runs[0]['m']}/ を使います)")
    else:
        _merge_runs(ctx, runs)


def merged(ctx: dict) -> dict[str, Path]:
    """04 以降が読む 03 の出力 (table・rep-seqs・denoising-stats.tsv)。ランが 1 つならそのラン、2 つ以上なら 03_merge/ (6-4)。"""
    runs = runs_of(ctx["cfg"])
    d = ctx["out"] / ("03_merge" if len(runs) > 1 else f"03_dada2/{runs[0]['m']}")
    return {"table": d / "table.qza", "rep-seqs": d / "rep-seqs.qza", "stats": d / "denoising-stats.tsv"}


def _merge_runs(ctx: dict, runs: list[dict]) -> None:
    """ランの表と代表配列をつなげる (6-1・6-9)。trim_left がランどうしでそろっているかは確かめない
    (261026 本人決定: 結合前の位置の検査はやめた。そろえるのは config を書く人。雛形のコメントにタグ分けのランは 23・22 と書いてある)。"""
    log, out = ctx["log"], ctx["out"]
    log("")
    log(f"[結合: {' + '.join(r['m'] for r in runs)}]")
    d = out / "03_merge"
    d.mkdir(parents=True, exist_ok=True)
    table, seqs, stats = d / "table.qza", d / "rep-seqs.qza", d / "denoising-stats.tsv"
    src = [out / "03_dada2" / r["m"] for r in runs]   # M の順 = ID の順 (6-9)

    # 1/3 table（重なる検体は 01 で止めているので、既定の overlap-method = error のまま）
    Q.run(log, "03", 1, 3, "リード数ファイルの結合", [table],
          ["feature-table", "merge"] + sum((["--i-tables", str(s / "table.qza")] for s in src), []) + ["--o-merged-table", str(table)])

    # 2/3 rep-seqs
    Q.run(log, "03", 2, 3, "代表配列ファイルの結合", [seqs],
          ["feature-table", "merge-seqs"] + sum((["--i-data", str(s / "rep-seqs.qza")] for s in src), []) + ["--o-merged-data", str(seqs)])

    # 3/3 denoising-stats.tsv（ランごとの tsv を ID の順につなげる。2 行目の #q2:types は 1 回だけ）
    if stats.exists():
        log.step("03", 3, 3, True, "品質統計ファイルの結合", "既存) denoising-stats.tsv")
    else:
        log.step("03", 3, 3, False, "品質統計ファイルの結合", "denoising-stats.tsv")
        parts = [(s / "denoising-stats.tsv").read_text(encoding="utf-8").rstrip("\n").split("\n") for s in src]
        stats.write_text("\n".join(parts[0][:2] + sum((p[2:] for p in parts), [])) + "\n", encoding="utf-8")


def read_denoising_stats(path: Path) -> pd.DataFrame:
    """denoising-stats.tsv を読む。2 行目の `#q2:types` を落とし、数値列を int/float にする。"""
    t = pd.read_csv(path, sep="\t", dtype=str)
    t = t[~t.iloc[:, 0].astype(str).str.startswith("#")].copy()
    for c in t.columns[1:]:
        t[c] = pd.to_numeric(t[c], errors="coerce")
    return t.rename(columns={t.columns[0]: "sample-id"})


# ============================================================ 04 Relabel
def step04_relabel(ctx: dict) -> None:
    """サンプル名を newID に付け替える (マニュアル §4・アーカイブ準拠)。

    [注意] `feature-table group --p-mode sum` は同じ newID を持つ検体を**合算**する。
    newID の重複は 01 で検証済みだが、念のため改名後の検体数が一致することを確かめる。
    """
    log, out = ctx["log"], ctx["out"]
    log.head("04. Relabel (サンプル名の変更)")
    C.check(ctx, "04", 1, 2, "リード数ファイルのサンプル名変更", "table_cn.qza")   # 条件の照合（261014 E1a）
    d = out / "04_relabel"
    d.mkdir(parents=True, exist_ok=True)   # ステップ冒頭で明示（アーカイブ準拠）
    table = merged(ctx)["table"]   # ランが 1 つならそのラン・2 つ以上なら 03_merge/ (261026)
    md, md_cn = out / "01_metadata" / "sample-metadata.txt", out / "01_metadata" / "sample-metadata_cn.txt"
    table_cn = d / "table_cn.qza"

    Q.run(log, "04", 1, 2, "リード数ファイルのサンプル名変更", [table_cn], [
        "feature-table", "group",
        "--i-table", str(table),
        "--p-axis", "sample",
        "--m-metadata-file", str(md),
        "--m-metadata-column", "newID",
        "--p-mode", "sum",
        "--o-grouped-table", str(table_cn),
    ])

    table_cn_qzv = d / "table_cn.qzv"
    sfreq, ffreq = d / "sample-frequencies_cn.qza", d / "feature-frequencies_cn.qza"
    Q.run(log, "04", 2, 2, "サンプル名変更済みリード数ファイルの可視化", [table_cn_qzv, sfreq, ffreq], [
        "feature-table", "summarize",
        "--i-table", str(table_cn),
        "--m-metadata-file", str(md_cn),
        "--o-summary", str(table_cn_qzv),
        "--o-sample-frequencies", str(sfreq),
        "--o-feature-frequencies", str(ffreq),
    ])

    # 検証: 改名前後で検体数が一致するか（合算事故の検出）
    n_before = len(pd.read_csv(md, sep="\t"))
    n_after = _count_samples(log, sfreq, d / "_sample-frequencies_cn")
    mark = "[OK]" if n_before == n_after else "[エラー] 検体数が減っています (newID の重複による合算の疑い)"
    log("-" * 65)
    log(f"検体数: 改名前 {n_before} → 改名後 {n_after}" + ("" if n_before == n_after else f" {mark}"))
    log("-" * 65)
    if n_before != n_after:
        raise SystemExit(1)


def _count_samples(log, freq_qza: Path, tmp_dir: Path) -> int:
    """sample-frequencies.qza (ImmutableMetadata)を export して行数を数える。"""
    Q.run_plain(log, ["qiime", "tools", "export", "--input-path", str(freq_qza), "--output-path", str(tmp_dir)])
    tsv = next(tmp_dir.glob("*.tsv"))
    t = pd.read_csv(tsv, sep="\t", dtype=str)
    import shutil
    shutil.rmtree(tmp_dir)                 # 数え終わったら一時フォルダを消す (07・09 の書き出しと同じ。261021・F5)
    return int((~t.iloc[:, 0].astype(str).str.startswith("#")).sum())


# ============================================================ 05 Phylogeny
def step05_phylogeny(ctx: dict) -> None:
    """代表配列から系統樹を作る (MAFFT 整列 → マスク → FastTree → midpoint root)。
    公式 moving-pictures ＝ マニュアル §6 ＝ アーカイブ。出力 4 つはすべて必須。
    入力は 03 の rep-seqs.qza (改名は table 側だけで、ASV の ID は変わらない)。"""
    log, out, cfg = ctx["log"], ctx["out"], ctx["cfg"]
    log.head("05. Phylogeny (系統樹の生成)")
    C.check(ctx, "05", 1, 1, "系統樹の生成", "aligned-rep-seqs.qza, masked-aligned-rep-seqs.qza, unrooted-tree.qza, rooted-tree.qza")   # 条件の照合（261014 E1a）
    # 261009 本人決定（D11b）: threads に関係なく 1 スレッド。2 以上だと QIIME 2 は非決定的な FastTreeMP に切り替わり、
    # 作り直すたびに木（UniFrac・Faith PD）が変わる（q2-phylogeny plugin_setup.py の説明・点検 F060/F178 の実測）
    threads = "1"
    d = out / "05_phylogeny"
    d.mkdir(parents=True, exist_ok=True)   # ステップ冒頭で明示（アーカイブ準拠）
    seqs = merged(ctx)["rep-seqs"]
    outs = [d / "aligned-rep-seqs.qza", d / "masked-aligned-rep-seqs.qza",
            d / "unrooted-tree.qza", d / "rooted-tree.qza"]
    Q.run(log, "05", 1, 1, "系統樹の生成", outs, [
        "phylogeny", "align-to-tree-mafft-fasttree",
        "--i-sequences", str(seqs),
        "--p-n-threads", threads,
        "--o-alignment", str(outs[0]),
        "--o-masked-alignment", str(outs[1]),
        "--o-tree", str(outs[2]),
        "--o-rooted-tree", str(outs[3]),
    ])


# ============================================================ 06 Rarefaction
def decide_sampling_depth(ctx: dict, d: Path) -> int:
    """B-4: 下限 floor 以上の検体のうち最小の non-chimeric 本数を depth に採用し、
    除外される検体を excluded_samples.tsv に明示する。決めた値は sampling_depth.txt に残し、
    07 diversity も同じ値を読む (再実行しても揺れない)。"""
    log, cfg, out = ctx["log"], ctx["cfg"], ctx["out"]
    floor = cfg["diversity"]["sampling_depth_floor"]   # 必須 (261021。以前は書かないと 10000)
    t = read_denoising_stats(merged(ctx)["stats"])
    md = pd.read_csv(out / "01_metadata" / "sample-metadata.txt", sep="\t", dtype=str)
    t = t.merge(md[["#SampleID", "newID"]], left_on="sample-id", right_on="#SampleID", how="left")
    nc = "non-chimeric"
    above = t[t[nc] >= floor]
    if above.empty:
        log(f"  [エラー] 下限 {floor:,} 以上の検体が 1 つもありません。config の sampling_depth_floor を見直してください。")
        raise SystemExit(1)
    depth = int(above[nc].min())
    at = above.loc[above[nc].idxmin()]
    low = t[t[nc] < floor].sort_values(nc)
    (d / "sampling_depth.txt").write_text(f"{depth}\n", encoding="utf-8")
    exc = low[["sample-id", "newID", nc]].rename(columns={nc: "non-chimeric"})
    exc["reason"] = f"non-chimeric < sampling_depth_floor ({floor})"
    exc.to_csv(d / "excluded_samples.tsv", sep="\t", index=False)
    log.step("06", 1, 3, False, "sampling-depth の動的決定", f"depth = {depth:,}")
    log("-" * 65)
    log(f"下限 {floor:,} 以上の検体のうち最小: {at['newID']}")
    if len(low):
        log(f"下限 {floor:,} 未満: 06_rarefaction/excluded_samples.tsv:")
        for _, r in low.iterrows():
            log(f"  {r['sample-id']:<16} {str(r['newID']):<26} non-chimeric {int(r[nc]):>8,}")
        mn = int(low[nc].min())            # 全検体を残すには、除外された検体の最小値以下にする必要がある (261021・F6。以前は最大値を案内していた)
        if mn > 0:
            log(f"  ※ 全検体を残したいなら、新しい実験フォルダで config の sampling_depth_floor を {mn:,} 以下にしてください")
        else:
            log("  ※ リード数 0 の検体があるため、全検体を残すことはできません")
        log("-" * 65)
    else:
        log(f"下限 {floor:,} 未満: なし ({len(t)} 検体すべてが depth 以上)")
        log("-" * 65)
    return depth


def step06_rarefaction(ctx: dict) -> None:
    """rarefaction curve。公式 moving-pictures 準拠 (既定 metrics ＋ 系統樹あり)。
    [注意] マニュアル §8・アーカイブにある chao1 版 (旧 3/5・5/5)は **公式に無いため採らない** (260914 本人判断)。
    その結果、段数は 5 → 3 になった。"""
    log, out = ctx["log"], ctx["out"]
    log.head("06. Rarefaction (sampling-depth の決定)")
    C.check(ctx, "06", 1, 3, "sampling-depth の動的決定", "sampling_depth.txt, excluded_samples.tsv")   # 条件の照合（261014 E1a）
    d = out / "06_rarefaction"
    d.mkdir(parents=True, exist_ok=True)   # ステップ冒頭で明示（アーカイブ準拠）
    depth = decide_sampling_depth(ctx, d)
    table_cn = out / "04_relabel" / "table_cn.qza"
    tree = out / "05_phylogeny" / "rooted-tree.qza"
    md_cn = out / "01_metadata" / "sample-metadata_cn.txt"

    qzv1 = d / "alpha-rarefaction.qzv"
    Q.run(log, "06", 2, 3, "rarefaction curve の可視化", [qzv1], [
        "diversity", "alpha-rarefaction",
        "--i-table", str(table_cn),
        "--i-phylogeny", str(tree),
        "--p-max-depth", str(depth),
        "--m-metadata-file", str(md_cn),
        "--o-visualization", str(qzv1),
    ])
    Q.export_qzv(log, "06", 3, 3, "rarefaction の HTML 化", qzv1, out / "HTML_06_alpha-rarefaction.qzv")


# ============================================================ 07 Diversity
def step07_diversity(ctx: dict) -> None:
    """α・β多様性 (core-metrics-phylogenetic)。公式 moving-pictures ＝ マニュアル §7 ＝ アーカイブ。
    depth は 06 が決めた sampling_depth.txt を読む (再実行しても揺れない)。
    `--output-dir` は既存ディレクトリを拒否するので、代表的な出力の有無でスキップを判定する。"""
    log, out, cfg = ctx["log"], ctx["out"], ctx["cfg"]
    log.head("07. Diversity (多様性の解析)")
    C.check(ctx, "07", 1, 2, "α・β多様性ファイルの生成", "rarefied_table.qza, shannon_vector.qza, …")   # 条件の照合（261014 E1a）
    d = out / "07_diversity"
    depth_file = out / "06_rarefaction" / "sampling_depth.txt"
    if not depth_file.exists():
        log("  [エラー] 06_rarefaction/sampling_depth.txt がありません。Step 06 を先に実行してください。")
        raise SystemExit(1)
    depth = int(depth_file.read_text().strip())
    threads = threads_arg(cfg, "auto")
    seed = rarefaction_seed(cfg)

    sentinel = [d / "rarefied_table.qza", d / "shannon_vector.qza",
                d / "weighted_unifrac_distance_matrix.qza", d / "bray_curtis_pcoa_results.qza"]
    # `--output-dir` は存在しないディレクトリを要求する。空なら消し、成果物がそろっていればスキップ判定に任せる。
    # 途中で失敗して成果物が欠けたフォルダが残っていたら止めて案内する（261015: --force を廃止したため。自動では消さない）
    if d.exists() and not any(d.iterdir()):
        d.rmdir()
    if d.exists() and not all(x.exists() for x in sentinel):
        log("  [エラー] 07_diversity が途中までしか作られていません (前回の実行が途中で止まった可能性)。")
        log("  ※ 07_diversity をゴミ箱へ移してから実行してください")
        raise SystemExit(1)
    # [注意] 07 だけはフォルダを作らない。`--output-dir` が自分で作る（アーカイブも mkdir していない）
    Q.run(log, "07", 1, 2, "α・β多様性ファイルの生成", sentinel, [
        "diversity", "core-metrics-phylogenetic",
        "--i-phylogeny", str(out / "05_phylogeny" / "rooted-tree.qza"),
        "--i-table", str(out / "04_relabel" / "table_cn.qza"),
        "--p-sampling-depth", str(depth),
        "--p-n-jobs-or-threads", threads,
        "--p-random-seed", str(seed),          # 261009 本人決定（D11a）: 希薄化の種を固定
        "--m-metadata-file", str(out / "01_metadata" / "sample-metadata_cn.txt"),
        "--output-dir", str(d),
    ])

    # 2/2 作図側へ渡す TSV（アーカイブ準拠の改名）
    todo = [("*_vector.qza", "alpha-diversity.tsv"),
            ("*_distance_matrix.qza", "distance-matrix.tsv"),
            ("*_pcoa_results.qza", "ordination.txt")]
    targets = [q.with_suffix(".tsv") for pat, _ in todo for q in sorted(d.glob(pat))]
    if targets and all(t.exists() for t in targets):
        log.step("07", 2, 2, True, "α・β多様性ファイルの .tsv 化", f"既存) {len(targets)} ファイル")
        return
    log.step("07", 2, 2, False, "α・β多様性ファイルの .tsv 化", f"{len(targets)} ファイル")
    import subprocess
    tmp = d / "_export"
    for pat, produced in todo:
        for q in sorted(d.glob(pat)):
            tmp.mkdir(exist_ok=True)
            Q.run_plain(log, ["qiime", "tools", "export", "--input-path", str(q), "--output-path", str(tmp)])
            (tmp / produced).replace(q.with_suffix(".tsv"))
    if tmp.exists():
        import shutil
        shutil.rmtree(tmp)
    for t in targets:                                  # 書き出したファイルを箇条書き（260925 本人指定）
        log(f"  - {t.name}")


# ============================================================ 09 Taxonomy
def read_classifier_facts(path: Path) -> dict[str, str]:
    """分類器に付属するファイルから素性を読み取る (config に書き写さない)。

    - `<name>.qza.doi`  … SILVA が配布する DOI (引用用)
    - `README`          … SILVA 同梱。`Creation Details` に版・領域・habitat がある
    どちらも無ければ空のまま返す (実行は止めない)。
    """
    facts: dict[str, str] = {}
    doi = path.with_suffix(path.suffix + ".doi")
    if doi.exists():
        facts["DOI"] = doi.read_text(encoding="utf-8").strip()
    readme = path.parent / "README"
    if readme.exists():
        txt = readme.read_text(encoding="utf-8", errors="replace")
        for key, label in [("Source Database", "Source"), ("SILVA Version", "SILVA"),
                           ("QIIME2 Version", "QIIME 2"), ("Clawback Version", "q2-clawback")]:
            m = re.search(rf"^{re.escape(key)}:\s*(.+)$", txt, re.M)
            if m:
                facts[label] = m.group(1).strip()
        for key in ("Region", "Taxonomic Ranks", "Propagation", "Dereplication Mode", "Habitat"):
            m = re.search(rf"^- {re.escape(key)}:\s*(.+)$", txt, re.M)
            if m:
                facts[key] = m.group(1).strip()
        hab, claw = facts.get("Habitat", "--"), facts.get("q2-clawback", "--")
        facts["Weighting"] = "uniform" if hab in ("--", "") and claw in ("--", "") else f"weighted / {hab}"
    return facts


def verify_classifier(ctx: dict) -> Path:
    """B-5: 分類器は自作しない。config が指すファイルの存在を確認し、素性をログに残す。

    [注意] md5 照合は行わない (260915 本人判断で廃止)。完全性はダウンロード時に人が 1 度確認する運用。
    [注意] SILVA 版・領域・habitat・DOI は **config に書かず、付属ファイルから読む** (260915 本人判断)。
       config に書き写すと実ファイルと食い違いうるため、実物を単一の情報源にする。
    """
    log, cfg = ctx["log"], ctx["cfg"]
    try:
        path = classifier_path(cfg)              # 01 の前にも同じ確認を通している (261021・C7)。途中で消えたときだけここで止まる
    except SystemExit as e:
        log(str(e))
        raise SystemExit(1)
    log("-" * 65)
    log(f"分類器: {path.name} ({path.stat().st_size / 1024 / 1024:,.0f} MB)")
    facts = read_classifier_facts(path)
    if facts:
        for k in ("Source", "SILVA", "QIIME 2", "Region", "Weighting",
                  "Taxonomic Ranks", "Propagation", "Dereplication Mode", "q2-clawback", "DOI"):
            if k in facts:
                log(f"  {k:<18} {facts[k]}")
    else:
        log(f"  [注意] 付属の README / .doi が見つからず、分類器の素性を記録できません ({path.parent} に SILVA 同梱の README と .doi を置いてください)")
    log("-" * 65)
    return path


def step09_taxonomy(ctx: dict) -> None:
    """細菌の同定。番号・出力名はアーカイブ準拠 (8 段)。公式 = マニュアル = アーカイブで形は一致。"""
    log, out, cfg = ctx["log"], ctx["out"], ctx["cfg"]
    log.head("09. Taxonomy (細菌の同定)")
    C.check(ctx, "09", 1, 8, "代表配列分類結果ファイルの生成", "rep-seqs_classified.qza")   # 条件の照合（261014 E1a）
    d = out / "09_taxonomy"
    d.mkdir(parents=True, exist_ok=True)   # ステップ冒頭で明示（アーカイブ準拠）
    tcfg = cfg["taxonomy"]                 # 必須 (261021。以前は書かないと order・family・genus と 0.7)
    ranks = [str(x).strip().lower() for x in tcfg["collapse_ranks"]]
    conf = str(tcfg["confidence"])
    n_jobs = str(cfg["classifier"]["n_jobs"])
    classifier = ctx.get("classifier") or verify_classifier(ctx)   # 08 で確認済み（260928: 分類器の記録は 08 に出す）

    seqs = merged(ctx)["rep-seqs"]
    table_cn = out / "04_relabel" / "table_cn.qza"
    md_cn = out / "01_metadata" / "sample-metadata_cn.txt"
    cls = d / "rep-seqs_classified.qza"

    # 1/8 classify-sklearn（ピークメモリを計測: n_jobs の見直しに使う）
    ran, peak = Q.run_measured(log, "09", 1, 8, "代表配列分類結果ファイルの生成", [cls], [
        "feature-classifier", "classify-sklearn",
        "--i-classifier", str(classifier),
        "--i-reads", str(seqs),
        "--p-confidence", conf,
        "--p-n-jobs", n_jobs,
        "--o-classification", str(cls),
    ])
    if ran:
        (d / "classify_peak_rss_mb.txt").write_text(f"n_jobs={n_jobs}\tpeak_rss_mb={peak:.0f}\n", encoding="utf-8")

    # 2/8 可視化
    cls_qzv = d / "rep-seqs_classified.qzv"
    Q.run(log, "09", 2, 8, "代表配列分類結果ファイルの可視化", [cls_qzv], [
        "metadata", "tabulate", "--m-input-file", str(cls), "--o-visualization", str(cls_qzv),
    ])

    # 3/8 .tsv 化（export は taxonomy.tsv を吐くので改名）
    cls_tsv = d / "rep-seqs_classified.tsv"
    if not cls_tsv.exists():
        log.step("09", 3, 8, False, "代表配列分類結果ファイルの .tsv 化", cls_tsv.name)   # 改名後の名前を出す (261019・F258)
        Q.run(log, "09", 3, 8, "代表配列分類結果ファイルの .tsv 化", [d / "taxonomy.tsv"], [
            "tools", "export", "--input-path", str(cls), "--output-path", str(d),
        ], force=True, quiet_head=True)
        (d / "taxonomy.tsv").replace(cls_tsv)
    else:
        log.step("09", 3, 8, True, "代表配列分類結果ファイルの .tsv 化", "既存) rep-seqs_classified.tsv")

    # 4/8 HTML 化
    Q.export_qzv(log, "09", 4, 8, "代表配列分類結果ファイルの HTML 化", cls_qzv,
                 out / "HTML_09_rep-seqs_classified.qzv")

    # 5/8 taxa barplot
    bar = d / "taxa-bar-plots.qzv"
    Q.run(log, "09", 5, 8, "細菌叢組成ファイルの可視化", [bar], [
        "taxa", "barplot",
        "--i-table", str(table_cn), "--i-taxonomy", str(cls),
        "--m-metadata-file", str(md_cn), "--o-visualization", str(bar),
    ])

    # 6/8 HTML 化
    Q.export_qzv(log, "09", 6, 8, "細菌叢組成ファイルの HTML 化", bar, out / "HTML_09_taxa-bar-plots.qzv")

    # 7/8 collapse（階級ごと）。--p-level の番号は分類結果から実測（SILVA 144 は k__ が入り 138 と番号がずれる）
    layout = rank_layout(cls_tsv)
    resolved = resolve_ranks(layout, ranks, log)
    for r, lv in resolved:
        nm = RANK_JA[RANK_PREFIX[r]]
        qza = d / f"collapsed_table_{r}.qza"
        Q.run(log, "09", 7, 8, f"{nm} (level {lv}) における組成ファイルの生成", [qza], [
            "taxa", "collapse",
            "--i-table", str(table_cn), "--i-taxonomy", str(cls),
            "--p-level", str(lv), "--o-collapsed-table", str(qza),
        ])

    # 8/8 .tsv 化（export → biom → tsv。アーカイブは .txt だが中身はタブ区切りなので .tsv とする）
    import subprocess, shutil
    for r, lv in resolved:
        nm = RANK_JA[RANK_PREFIX[r]]
        tsv = d / f"collapsed_table_{r}.tsv"
        if tsv.exists():
            log.step("09", 8, 8, True, f"{nm} における組成ファイルの .tsv 化", f"既存) {tsv.name}")
            continue
        log.step("09", 8, 8, False, f"{nm} における組成ファイルの .tsv 化", tsv.name)
        tmp = d / f"_export_{r}"
        Q.run_plain(log, ["qiime", "tools", "export", "--input-path", str(d / f"collapsed_table_{r}.qza"),
                          "--output-path", str(tmp)])
        Q.run_plain(log, ["biom", "convert", "-i", str(tmp / "feature-table.biom"), "-o", str(tsv), "--to-tsv"])
        shutil.rmtree(tmp)


def rank_layout(cls_tsv: Path) -> dict[int, str]:
    """分類結果の文字列から「位置 → 接頭辞」を実測する (例 {1:'d',2:'k',3:'p',...})。
    SILVA 144 は 138 と違い位置 2 に k__ (Kingdom) が入るため、階級番号は分類器に依存する。"""
    import collections
    t = pd.read_csv(cls_tsv, sep="\t")
    pos: dict[int, collections.Counter] = collections.defaultdict(collections.Counter)
    for s in t["Taxon"].astype(str):
        for i, r in enumerate([x.strip() for x in s.split(";")]):
            if len(r) >= 3 and r[1:3] == "__":
                pos[i + 1][r[0]] += 1
    return {i: c.most_common(1)[0][0] for i, c in sorted(pos.items())}


RANK_JA = {"d": "Domain", "k": "Kingdom", "p": "Phylum", "c": "Class", "o": "Order",
           "f": "Family", "g": "Genus", "s": "Species"}
RANK_PREFIX = {v.lower(): k for k, v in RANK_JA.items()}   # "genus" -> "g"


def resolve_ranks(layout: dict[int, str], ranks: list[str], log) -> list[tuple[str, int]]:
    """config の階級名 (order/family/genus …)を、実測レイアウトで --p-level の番号に変換する。
    分類結果に無い階級 (例: genus 分類器で species)を指定したら止める。"""
    inv = {pre: pos for pos, pre in layout.items()}
    out = []
    for r in ranks:
        key = str(r).strip().lower()
        pre = RANK_PREFIX.get(key)
        if pre is None:
            log(f"  [エラー] config の collapse_ranks に知らない階級名があります: '{r}' (使える名前: {', '.join(RANK_PREFIX)})")
            raise SystemExit(1)
        if pre not in inv:
            log(f"  [エラー] 階級 '{r}' ({pre}__) はこの分類結果に存在しません。"
                f" 実測の階級: {' '.join(f'{i}={layout[i]}__' for i in sorted(layout))}")
            raise SystemExit(1)
        out.append((key, inv[pre]))
    return out
