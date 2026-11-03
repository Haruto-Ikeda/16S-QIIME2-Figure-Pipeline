"""各ステップの条件の記録と照合 (261014 本人決定・E1a 案 b)。

01〜09 は出力があれば飛ばす。そのままだと config や検体情報を変えても古い結果で進み、記録と中身が食い違う。
そこで各ステップのフォルダに、そのステップの条件を `_conditions.json` として残し、飛ばす前に今の条件と照らし合わせる。
- 違えば止める。解析の条件を変えるときは実験フォルダ自体を新しく作る (261015 本人方針)ので、その旨を案内する。自動では消さない。
- 記録が無い既存の出力 (このしくみを入れる前の出力)は、今の条件で記録して進む ([注記] に出す。照合はしていない)。
- 上流のステップの記録の要約 (digest)も一緒に残す。上流だけ作り直されて下流が古いまま、という食い違いも止められる。
threads は結果を変えないので照らし合わせない (05 の系統樹は常に 1 スレッド。D11)。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

FILE = "_conditions.json"
STEP_DIRS = {"01": "01_metadata", "02": "02_import", "03": "03_dada2", "04": "04_relabel",
             "05": "05_phylogeny", "06": "06_rarefaction", "07": "07_diversity", "09": "09_taxonomy"}
DEPS = {"01": [], "02": ["01"], "03": ["02"], "04": ["03", "01"], "05": ["04"], "06": ["05"], "07": ["06"],
        "09": ["04"]}                              # そのステップが入力に使う上流のステップ


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return "md5:" + h.hexdigest()


def own(ctx: dict, step: str) -> dict:
    """そのステップ自身の条件 (キー = ログに出す名前)。md5: で始まる値は中身の要約 (値そのものはログに出さない)。"""
    cfg, exp, base = ctx["cfg"], ctx["exp"], ctx["base"]
    if step == "01":                               # ラン・SampleSheet・検体情報 CSV・exclude・absent_ok (261022 に新しい設計へ)
        from .steps import runs_of, sample_info_of
        c = {}
        for r in runs_of(cfg):
            fq = sorted(p.name for p in r["fastq"].glob("*.fastq.gz"))
            c[f"{r['name']}: SampleSheet.txt"] = _md5(r["sheet"])
            c[f"{r['name']}: FASTQ のファイル名 ({len(fq)} ファイル)"] = "md5:" + hashlib.md5("\n".join(fq).encode()).hexdigest()
            c[f"{r['name']}: exclude"] = list(r["exclude"])
        for p, v in sorted(sample_info_of(cfg).items()):
            c[f"検体情報 CSV P{p:04d} ({v['path'].name})"] = _md5(v["path"])
            c[f"P{p:04d}: absent_ok"] = sorted(v["absent_ok"])
        return c
    if step == "03":                               # ランごとの quality_control (261026)
        from .steps import runs_of
        return {f"{r['name']}: quality_control.{k}": v for r in runs_of(cfg) for k, v in r["qc"].items()}
    if step == "06":
        return {"diversity.sampling_depth_floor": (cfg.get("diversity") or {}).get("sampling_depth_floor")}
    if step == "07":
        return {"diversity.random_seed": (cfg.get("diversity") or {}).get("random_seed")}
    if step == "09":
        from .steps import classifier_path
        p = classifier_path(cfg)
        t = cfg.get("taxonomy") or {}
        return {"classifier (ファイル名・大きさ)": f"{p.name} ({p.stat().st_size if p.exists() else 0:,} bytes)",
                "taxonomy.confidence": t.get("confidence"),
                "taxonomy.collapse_ranks": list(t.get("collapse_ranks") or [])}
    return {}


def _digest(out: Path, step: str) -> str | None:
    f = out / STEP_DIRS[step] / FILE
    return hashlib.md5(f.read_bytes()).hexdigest() if f.exists() else None


def current(ctx: dict, step: str) -> dict:
    return {"step": step, "own": own(ctx, step),
            "upstream": {dep: _digest(ctx["out"], dep) for dep in DEPS[step]}}


def downstream(step: str) -> list[str]:
    """そのステップと、その結果を (間接的にも)使う後ろのステップ。"""
    res, todo = [step], [step]
    while todo:
        s = todo.pop()
        for k, deps in DEPS.items():
            if s in deps and k not in res:
                res.append(k)
                todo.append(k)
    return sorted(res)


def _has_outputs(d: Path) -> bool:
    return d.is_dir() and any(p.name != FILE for p in d.iterdir())


def check(ctx: dict, step: str, n: int, total: int, title: str, outs: str) -> None:
    """飛ばす前の照合。出力が無いときは何もしない (実行後に record で記録する)。"""
    log, out = ctx["log"], ctx["out"]
    d = out / STEP_DIRS[step]
    if not _has_outputs(d):
        return
    cur = current(ctx, step)
    f = d / FILE
    if not f.exists():
        f.write_text(json.dumps(cur, ensure_ascii=False, indent=1), encoding="utf-8")
        log("  [注意] 条件の記録が無かったので、今の config で記録しました (前回の条件とは照らし合わせていません)")
        return
    prev = json.loads(f.read_text(encoding="utf-8"))
    lines = []
    keys = list(cur["own"]) + [k for k in prev.get("own", {}) if k not in cur["own"]]
    for k in keys:
        a, b = prev.get("own", {}).get(k, "なし"), cur["own"].get(k, "なし")
        if a != b:
            if str(a).startswith("md5:") or str(b).startswith("md5:"):
                lines.append(f"    {k}: 前回と中身が違います")
            else:
                lines.append(f"    {k}: {a} (前回) → {b} (今回)")
    for dep, dg in cur["upstream"].items():
        if prev.get("upstream", {}).get(dep) != dg:
            lines.append(f"    上流の {dep} ({STEP_DIRS[dep]}): 前回と違う条件で作り直されています")
    if not lines:
        return
    # 261015 本人方針: QIIME 2 の解析（01〜09）の条件を変えるときは、実験フォルダ自体を新しく作る（同じフォルダで一部を作り直さない）
    log(f"{step}  {n}/{total}|　中止　|{title}: 既存) {outs}")
    log("  [エラー] 既存の出力は、今の config と違う条件で作られています。")
    for ln in lines:
        log(ln)
    if any(not ln.lstrip().startswith("上流の") for ln in lines):
        log("  ※ 条件を変えて解析し直すなら、新しい実験フォルダを作って実行してください")
        log("  ※ 前回の条件のまま進めるなら、config を前回の値に戻してください")
    else:                                          # 上流の手順だけが作り直されている（フォルダを一部だけ消して実行したなど）
        log("  ※ 途中の手順だけが作り直されています。output/16S/Qiime2 をすべてゴミ箱へ移して最初から実行するか、新しい実験フォルダで実行してください")
    raise SystemExit(1)


def record(ctx: dict, step: str) -> None:
    """ステップを実行 (または照合を通過)したあとに、今の条件を記録する。"""
    d = ctx["out"] / STEP_DIRS[step]
    if d.is_dir():
        (d / FILE).write_text(json.dumps(current(ctx, step), ensure_ascii=False, indent=1), encoding="utf-8")
