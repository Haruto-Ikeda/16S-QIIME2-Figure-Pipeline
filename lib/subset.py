#!/usr/bin/env python3
"""subset.py — 絞り込みの形 (イ。config の parent:・select:)。親の実験の結果から検体を選んだ写しを、子の Qiime2/ に作る。

10 (作図) はこの写しを、01〜09 の結果と同じように読む (figures.py は変えない。261026 決定・案 a)。
- 親から使うもの (7-2): 01 (検体の表)・03 の denoising-stats.tsv (#1)・06 (depth・除外・#2 のレアファクションカーブの CSV)・07 (α の値・距離行列)・09 (階級ごとの表)
- α の値と距離はそのまま使い (行・列を絞るだけ)、PCoA は絞った距離行列から計算し直す。09 の表はどの検体にも出ない分類群の行を除く
- 写しは実行するたびに作り直す (select を書き換えたら、そのまま反映する。261026 決定・案 a)
- 親が作り直されていれば止める (7-4): 初回に親の 01〜09 の条件の記録の要約を _parent.json に残し、次から照らし合わせる
"""
from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import pandas as pd

from . import conditions as C
from .steps import ID_RE, abs_path, is_none, merged, runs_of

PARENT_FILE = "_parent.json"


def _die(msg: str) -> None:
    raise SystemExit(f"  [エラー] {msg}")


def _values(v) -> str:
    return "all" if isinstance(v, str) else (", ".join(map(str, v)) or "なし")


def _check_parent(log, pexp: Path, out: Path) -> str:
    """親の 01〜09 の条件の記録の要約を、初回は残し、次からは照らし合わせる (7-4)。"""
    pout = pexp / "output/16S/Qiime2"
    cur = {}
    for step, d in C.STEP_DIRS.items():
        f = pout / d / C.FILE
        if not f.exists():
            _die(f"親の実験に {step} の条件の記録 ({d}/{C.FILE}) がありません: {pexp} (親を新しい形で解析し直してください)")
        cur[step] = hashlib.md5(f.read_bytes()).hexdigest()
    rec = {"parent": str(pexp), "conditions": cur}
    f = out / PARENT_FILE
    if not f.exists():
        f.write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
        return "01〜09 の条件の記録を _parent.json に残しました (初回)"
    prev = json.loads(f.read_text(encoding="utf-8"))
    if prev.get("parent") != rec["parent"]:
        _die(f"この実験フォルダは、別の親の実験から作られています: {prev.get('parent')} (今の config は {pexp})\n"
             "  ※ 親を変えるなら、新しい実験フォルダを作って実行してください")
    changed = [s for s in cur if prev.get("conditions", {}).get(s) != cur[s]]
    if changed:
        _die(f"親の実験が、この絞り込みを作ったあとに作り直されています ({', '.join(changed)})\n"
             "  ※ 新しい親の結果で描くなら、新しい実験フォルダを作って実行してください")
    return "前回と同じです (親は作り直されていません)"


def _select(md: pd.DataFrame, sel: dict) -> pd.DataFrame:
    """select: で検体を選ぶ。書いた値・ID がデータに無ければ止める (書き間違いに気づくため)。"""
    keep = pd.Series(True, index=md.index)
    for col in ("Group_name", "Subgroup_name"):
        v = sel[col]
        if isinstance(v, list):
            want = [str(x).strip() for x in v]
            miss = [x for x in want if x not in set(md[col])]
            if miss:
                _die(f"config の select.{col} に、親の実験に無い値があります: {', '.join(miss)} (ある値: {', '.join(sorted(set(md[col])))})")
            keep &= md[col].isin(want)
    exc = [] if is_none(sel["exclude"]) else [x.strip() for x in sel["exclude"]]
    miss = [x for x in exc if x not in set(md["#SampleID"])]
    if miss:
        _die(f"config の select.exclude に、親の実験に無い ID があります: {', '.join(miss)}")
    keep &= ~md["#SampleID"].isin(exc)
    if not keep.any():
        _die("select: で選んだ検体が 0 件です (select の書き方を見直してください)")
    return md[keep]


def _write_tsv_rows(src: Path, dst: Path, keep: set, comment_rows: int = 0) -> None:
    """1 列目が ID の tsv を、keep の行だけにして写す (先頭の comment_rows 行と列名の行はそのまま)。"""
    lines = src.read_text(encoding="utf-8").rstrip("\n").split("\n")
    head, body = lines[:comment_rows + 1], lines[comment_rows + 1:]
    dst.write_text("\n".join(head + [ln for ln in body if ln.split("\t")[0] in keep or ln.startswith("#")]) + "\n", encoding="utf-8")


def build(ctx: dict) -> dict:
    """子の Qiime2/ に絞った写しを作り、10 に渡す cfg (親のランなどを借りたもの) を返す。"""
    from skbio import DistanceMatrix
    from skbio.stats.ordination import pcoa
    log, cfg, out = ctx["log"], ctx["cfg"], ctx["out"]
    pexp = abs_path(cfg["parent"], "parent")
    pout = pexp / "output/16S/Qiime2"
    import yaml
    pcfg = yaml.safe_load((pexp / "input/16S/config.yaml").read_text(encoding="utf-8")) or {}
    # 10 に渡す cfg: 作図の設定 (grouping など)・parent・select は子のもの、ラン・多様性・分類器・分類は親のもの
    ecfg = {**cfg, **{k: pcfg.get(k) for k in ("threads", "runs", "sample_info", "diversity", "classifier", "taxonomy")}}

    log.head("絞り込み (親の実験の結果から検体を選ぶ)")
    log("-" * 65)
    log(f"親の実験     : {pexp}")
    sel = cfg["select"]
    log(f"select       : Group_name {_values(sel['Group_name'])} / Subgroup_name {_values(sel['Subgroup_name'])}"
        f" / exclude {'なし' if is_none(sel['exclude']) else _values(sel['exclude'])}")
    log(f"親の条件     : {_check_parent(log, pexp, out)}")

    md = pd.read_csv(pout / "01_metadata/sample-metadata.txt", sep="\t", dtype=str, keep_default_na=False)
    sub = _select(md, sel)
    ids, names = set(sub["#SampleID"]), set(sub["newID"])

    # 01: 検体の表 3 つ
    d = out / "01_metadata"
    d.mkdir(parents=True, exist_ok=True)
    sub.to_csv(d / "sample-metadata.txt", sep="\t", index=False)
    cn = pd.read_csv(pout / "01_metadata/sample-metadata_cn.txt", sep="\t", dtype=str, keep_default_na=False)
    cn[cn["#SampleID"].isin(names)].to_csv(d / "sample-metadata_cn.txt", sep="\t", index=False)
    si = pd.read_csv(pout / "01_metadata/sample_information.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)
    si[si["#SampleID"].isin(ids)].to_csv(d / "sample_information.csv", index=False, encoding="utf-8-sig")

    # 03: #1 と検体情報の表のリード数 (親のランの数に合わせて、子の同じ場所に置く)
    src, dst = merged({"cfg": pcfg, "out": pout})["stats"], merged({"cfg": ecfg, "out": out})["stats"]
    dst.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv_rows(src, dst, ids, comment_rows=1)          # 2 行目 (#q2:types) はそのまま

    # 06: depth はそのまま・除外は選んだ検体の行だけ
    d = out / "06_rarefaction"
    d.mkdir(parents=True, exist_ok=True)
    (d / "sampling_depth.txt").write_text((pout / "06_rarefaction/sampling_depth.txt").read_text(encoding="utf-8"), encoding="utf-8")
    _write_tsv_rows(pout / "06_rarefaction/excluded_samples.tsv", d / "excluded_samples.tsv", ids)
    n_excl = len(pd.read_csv(d / "excluded_samples.tsv", sep="\t", dtype=str))
    # 06: #2 のレアファクションカーブの値 (HTML_06_alpha-rarefaction.qzv/<指標>.csv。1 列目が検体名)。CSV だけを写す (HTML は写さない)
    d = out / "HTML_06_alpha-rarefaction.qzv"
    d.mkdir(parents=True, exist_ok=True)
    for f in sorted((pout / "HTML_06_alpha-rarefaction.qzv").glob("*.csv")):
        lines = f.read_text(encoding="utf-8").rstrip("\n").split("\n")
        (d / f.name).write_text("\n".join(lines[:1] + [ln for ln in lines[1:] if ln.split(",")[0] in names]) + "\n", encoding="utf-8")

    # 07: α の値 (行を絞る)・距離行列 (行と列を絞る)・PCoA (絞った距離行列から計算し直す。QIIME 2 の diversity pcoa と同じ scikit-bio の pcoa)
    d = out / "07_diversity"
    d.mkdir(parents=True, exist_ok=True)
    for f in sorted((pout / "07_diversity").glob("*_vector.tsv")):
        _write_tsv_rows(f, d / f.name, names)
    n_pcoa = 0
    for f in sorted((pout / "07_diversity").glob("*_distance_matrix.tsv")):
        D = DistanceMatrix.read(str(f))
        D = D.filter([x for x in D.ids if x in names])
        D.write(str(d / f.name))
        if len(D.ids) >= 3:
            with warnings.catch_warnings():            # 全次元の計算・負の固有値 (Bray-Curtis などでは普通) の知らせは出さない
                warnings.simplefilter("ignore", RuntimeWarning)
                res = pcoa(D)
            res.write(str(d / f.name.replace("_distance_matrix.tsv", "_pcoa_results.tsv")))
            n_pcoa += 1
        else:
            _die(f"希薄化のあとに残った検体が {len(D.ids)} 件しかなく、PCoA を計算できません ({f.name})")

    # 09: 階級ごとの表 (検体の列を絞り、どの検体にも出ない分類群の行を除く)
    d = out / "09_taxonomy"
    d.mkdir(parents=True, exist_ok=True)
    ranks = []
    for f in sorted((pout / "09_taxonomy").glob("collapsed_table_*.tsv")):
        first = f.read_text(encoding="utf-8").split("\n", 1)[0]
        t = pd.read_csv(f, sep="\t", skiprows=1, index_col=0)
        t = t[[c for c in t.columns if c in names]]
        t = t[(t != 0).any(axis=1)]
        with open(d / f.name, "w", encoding="utf-8") as w:
            w.write(first + "\n")
            t.to_csv(w, sep="\t")
        ranks.append(f.stem.replace("collapsed_table_", ""))

    log(f"検体         : 親の {len(md)} 件のうち {len(sub)} 件" + (f" (うち希薄化で外れた検体 {n_excl} 件)" if n_excl else ""))
    log(f"写し         : 01_metadata・{dst.parent.relative_to(out)}・06_rarefaction・HTML_06_alpha-rarefaction.qzv (#2 の CSV だけ)・07_diversity (PCoA {n_pcoa} つを計算し直し)"
        f"・09_taxonomy ({len(ranks)} 階級)")
    log("-" * 65)
    return ecfg
