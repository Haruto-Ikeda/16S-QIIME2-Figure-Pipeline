#!/usr/bin/env python3
"""qiime.py — QIIME 2 コマンドの実行・ログ・スキップ判定・エラー処理。

設計の正本は ../README.md。ここは「各ステップが共通で使う土台」だけを持つ。

- 実行環境は **利用者が先に `conda activate rachis-qiime2-2026.7` している前提** (B-2)。
  絶対パス起動は R_HOME が張られず dada2 が segfault するため、env の外にいたら止める。
- 出力が既にあればスキップする (アーカイブから受け継いだ再開設計)。飛ばす前に各ステップが条件を照合する (lib/conditions.py・261014 E1a)。
- ログの書式もアーカイブ準拠: `01  1/3|　実行　|demux.qza`
"""
from __future__ import annotations

import os
import subprocess
import re
import sys
from datetime import datetime
from pathlib import Path

ENV_NAME = "rachis-qiime2-2026.7"        # README の手順で作る env の名前（案内に使うだけ。判定は版で行う）
QIIME_VERSION = "2026.7"                 # この版の QIIME 2 でだけ動かす（261015 E4b）


class Log:
    """標準出力と実行ログへ同時に書き出す (長時間ジョブの全量をログに残す)。"""

    def __init__(self, path: Path | None = None):
        self.f = open(path, "a", encoding="utf-8") if path else None
        self.where = ""                                  # いまの手順（警告の見出しに使う）
        self.notes: dict[str, list[str]] = {}            # 手順 → 警告・[注記]（途中では出さず、最後にまとめて出す。260927 本人決定）

    def __call__(self, msg: str = "") -> None:
        if ("[注意]" in msg or "[注記]" in msg) and "[エラー]" not in msg:
            for item in _split_note(msg):
                lst = self.notes.setdefault(self.where, [])
                if item not in lst:                      # 同じ文は 1 回だけ
                    lst.append(item)
            return
        m = re.match(r"^(\d\d)\s+\d+/\d+\|[^|]*\|([^:]+)", msg)
        if m:
            self.where = f"{m.group(1)} {m.group(2).strip()}"
        elif re.match(r"^\d\d\. ", msg):
            self.where = msg
        print(msg, flush=True)
        if self.f:
            self.f.write(msg + "\n")
            self.f.flush()

    def step(self, step: str, n: int, total: int, done: bool, name: str, detail: str = "") -> None:
        """アーカイブ準拠の進捗行。done=True なら「スキップ」。"""
        mark = "スキップ" if done else "　実行　"
        self(f"{step}  {n}/{total}|{mark}|{name}" + (f": {detail}" if detail else ""))

    def head(self, title: str) -> None:
        self("")
        self(title)

    def summary(self) -> list[str]:
        """最後の「注記: N件」と一覧 (260927 本人指定)。
        続く手順で注記がまったく同じなら 1 つにまとめ、見出しを「10 #4 PCoA PCo1×PCo2, PCoA 交差形式 PCo1〜4」とする。
        N = まとめたあとの見出しの数。"""
        blocks: list[list] = []                          # [[見出し, …], 注記の並び]
        for where, items in self.notes.items():
            if blocks and blocks[-1][1] == items:
                blocks[-1][0].append(where)
            else:
                blocks.append([[where], items])
        lines = [f"注記: {len(blocks)}件"]                  # 「警告」でなく「注記」（260927 本人指示）
        for heads, items in blocks:
            head = heads[0] or "(手順の前)"
            if len(heads) > 1 and all(re.match(r"^\d\d\. ", h or "") for h in heads):   # 手順の見出しが並ぶときは番号だけ（261015）
                lines.append(", ".join(h[:2] for h in heads))
                lines += [f"  - {it}" for it in items]
                lines.append("")
                continue
            for h in heads[1:]:                          # 同じ「10 #4」などは 2 つ目から外す
                a, b = head.split(" "), h.split(" ")
                k = 0
                while k < min(len(a), len(b)) - 1 and a[k] == b[k] and (k == 0 or a[k].startswith("#")):
                    k += 1
                head += ", " + " ".join(b[k:])
            lines.append(head)
            lines += [f"  - {it}" for it in items]
            lines.append("")
        if lines[-1] == "":
            lines.pop()
        return lines

    def close(self) -> None:
        if self.f:
            self.f.close()


def _split_note(msg: str) -> list[str]:
    """警告の 1 行を一覧用に整える (印と字下げを外す。見切れの警告は「N群 (最長 …)」の 1 行)。"""
    s = msg.replace("[注意]", "").replace("[注記]", "").strip()
    m = re.match(r"^(.*?)で見切れている可能性があります \(群名を短くしてください\): (.*)$", s)
    if m:
        items = [(name.strip(), float(mm)) for name, mm in re.findall(r"(.+?) \(約 ([\d.]+) mm\)(?:, )?", m.group(2))]
        name, mm = max(items, key=lambda x: x[1])      # 数と一番長いものだけ（260928 本人決定・案 B）
        return [f"{m.group(1)}で見切れている群: {len(items)}群 (最長 {name} 約 {mm:.1f} mm)"]
    return [s]


def run_plain(log: Log, cmd: list[str], timeout: int = 3600) -> None:
    """Q.run を通さない補助のコマンド (qiime tools export・biom convert)。失敗したら Q.run と同じ形でログに残して止める
    (261019・点検 F131。以前は check=True で落ち、原因の stderr がログに残らなかった)。"""
    t0 = datetime.now()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        log(f"  [エラー] 時間切れで止めました ({timeout} 秒): {' '.join(cmd[:3])} …")
        sys.exit(1)
    if p.returncode != 0:
        dt = (datetime.now() - t0).total_seconds()
        log(f"  [エラー] 失敗しました (exit={p.returncode}, {dt:.1f} 秒): {' '.join(cmd[:3])} …")
        log("  --- stderr ---")
        for ln in (p.stderr or "").rstrip().split("\n"):
            log("    " + ln)
        sys.exit(1)


def check_env(log: Log) -> None:
    """QIIME 2 の版と R_HOME を確かめる。合わなければ原因と対処を示して止める (B-2)。
    261015 本人決定 (E4b): env の名前ではなく QIIME 2 の版 (QIIME_VERSION.x)で判定する (別名の env・コンテナでも動くように)。"""
    import importlib.metadata as md
    lines = []
    try:
        ver = md.version("qiime2")
    except md.PackageNotFoundError:
        ver = None
    if ver is None or not ver.startswith(QIIME_VERSION + "."):
        lines += [f"  [エラー] QIIME 2 {QIIME_VERSION} の環境に入っていません (今の Python: {sys.prefix}・QIIME 2: {ver or 'なし'})。",
                  f"  ※ conda activate {ENV_NAME} (または QIIME 2 {QIIME_VERSION} の環境)をしてから、もう一度起動してください"]
    else:
        r_home = os.environ.get("R_HOME", "")
        if not r_home.startswith(sys.prefix):
            lines += [f"  [エラー] R_HOME が env の外を指しています ({r_home or '未設定'})。",
                      "    dada2 は R で動くため、このままだと Segmentation fault で落ちます。",
                      "  ※ conda activate をやり直してください (activate フックが R_HOME を設定します)"]
        else:
            lines += _check_r()
    if lines:
        log("")
        for ln in lines:
            log(ln)
        sys.exit(1)


_R_CHECKED: list[str] | None = None


def _check_r() -> list[str]:
    """env の R と、10 の #11・#12 (ANCOM-BC2) が使う R のパッケージ ANCOMBC・phyloseq があるかを起動時に確かめる
    (無いと 01〜09 のあと 10 の #11 で初めて止まるため。261027)。結果は 1 回だけ調べて覚えておく。"""
    global _R_CHECKED
    if _R_CHECKED is not None:
        return _R_CHECKED
    import subprocess
    rscript = Path(sys.prefix) / "bin" / "Rscript"
    if not rscript.exists():
        _R_CHECKED = [f"  [エラー] env の中に R (Rscript) がありません: {rscript}",
                      f"  ※ QIIME 2 {QIIME_VERSION} の env を公式の手順で作り直してください (R は env に含まれます)"]
        return _R_CHECKED
    try:
        r = subprocess.run([str(rscript), "-e", "suppressMessages({library(ANCOMBC); library(phyloseq)})"],
                           capture_output=True, text=True, timeout=300)
        err = [ln.strip() for ln in (r.stderr or "").strip().split("\n") if ln.strip()]
        hit = [ln for ln in err if "package" in ln or "パッケージ" in ln]   # どのパッケージが無いかの行 (無ければ最後の行)
        ok, why = r.returncode == 0, (hit or err or ["理由不明"])[0]
    except (OSError, subprocess.TimeoutExpired) as e:
        ok, why = False, f"{type(e).__name__}"
    _R_CHECKED = [] if ok else [
        f"  [エラー] R のパッケージ ANCOMBC・phyloseq を読み込めません ({why})",
        f"  ※ QIIME 2 {QIIME_VERSION} の env には両方が入っているはずです。env を公式の手順で作り直してください"]
    return _R_CHECKED


def env_label() -> str:
    """開始の枠に出す env の名前と QIIME 2 の版。"""
    import importlib.metadata as md
    return f"{Path(sys.prefix).name} (QIIME 2 {md.version('qiime2')})"


def pipeline_label() -> str:
    """開始の枠に出すパイプラインの版と実行した場所 (261102 本人決定・案 1)。
    VERSION が無くても止めない。更新のたびに図の版が増えないよう、図の版の判定 (VERSION_KEYS) には入れない。"""
    root = Path(__file__).resolve().parent.parent
    try:
        ver = (root / "VERSION").read_text(encoding="utf-8").splitlines()[0].strip()
    except (OSError, UnicodeDecodeError, IndexError):
        ver = "不明 (VERSION がありません)"
    return f"{ver}  ({root})"


def lib_versions() -> str:
    """作図で使うライブラリの版 (開始の枠に 1 行。261021 本人決定・F202/F206 案 b)。"""
    import subprocess
    import matplotlib
    import numpy
    import pandas
    import scipy
    import skbio
    parts = [f"scipy {scipy.__version__}", f"scikit-bio {skbio.__version__}", f"numpy {numpy.__version__}",
             f"pandas {pandas.__version__}", f"matplotlib {matplotlib.__version__}"]
    try:
        r = subprocess.run([str(Path(sys.prefix) / "bin" / "Rscript"), "-e",
                            'cat(paste0("R ", R.version$major, ".", R.version$minor, " / R ANCOMBC ", packageVersion("ANCOMBC")))'],
                           capture_output=True, text=True, timeout=120)
        parts.append(r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else "R ANCOMBC (版を読めませんでした)")
    except (OSError, subprocess.TimeoutExpired):
        parts.append("R ANCOMBC (版を読めませんでした)")
    return " / ".join(parts)


def save_env_list(path: Path) -> str | None:
    """env の全パッケージの一覧 (conda list --explicit) をファイルに残す (261021 本人決定・F202/F206 案 b)。
    同じ種類の機械なら conda create --name <env> --file <このファイル> で同じ env を作り直せる。
    残せなかったときは理由を返す。"""
    import os
    import subprocess
    conda = os.environ.get("CONDA_EXE")
    if not conda:
        return "conda の場所 (CONDA_EXE) が分かりません"
    try:
        r = subprocess.run([conda, "list", "--explicit", "--prefix", sys.prefix], capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as e:
        return f"conda list が失敗しました ({type(e).__name__})"
    if r.returncode != 0 or "@EXPLICIT" not in r.stdout:
        return f"conda list が失敗しました (exit={r.returncode})"
    path.write_text(r.stdout, encoding="utf-8")
    return None


def _log_cmd(log: Log, cmd: list[str]) -> None:
    """コマンドを「$ qiime <サブコマンド>」＋ オプションごとに 1 行 (字下げ 2)で書き、最後に空行 (260925 本人指定)。"""
    head, rest = [], list(cmd)
    while rest and not rest[0].startswith("--"):
        head.append(rest.pop(0))
    log("$ " + " ".join(_q(a) for a in head))
    line = []
    for a in rest:
        if a.startswith("--") and line:
            log("  " + " ".join(line))
            line = []
        line.append(_q(a))
    if line:
        log("  " + " ".join(line))
    log("")


def run(log: Log, step: str, n: int, total: int, name: str,
        outputs: list[Path], args: list[str], force: bool = False, quiet_head: bool = False) -> bool:
    """qiime を 1 回実行する。outputs が全部あればスキップ。

    戻り値: 実際に実行したら True、スキップしたら False。
    [注意] 出力フォルダはここでは作らない。**各ステップの冒頭で明示的に作る** (アーカイブ準拠)。
       `--output-dir` を使うコマンド (07)は「存在しないディレクトリ」を要求するため、暗黙の mkdir は事故になる。
    """
    if not force and outputs and all(o.exists() for o in outputs):
        log.step(step, n, total, True, name, "既存) " + ", ".join(o.name for o in outputs))
        return False
    if not quiet_head:
        log.step(step, n, total, False, name, ", ".join(o.name for o in outputs))
    cmd = ["qiime"] + args
    _log_cmd(log, cmd)
    t0 = datetime.now()
    p = subprocess.run(cmd, capture_output=True, text=True)
    dt = (datetime.now() - t0).total_seconds()
    if p.returncode != 0:
        log("")
        log(f"  [エラー] 失敗しました (exit={p.returncode}, {dt:.1f} 秒)")
        log("  --- stderr ---")
        for ln in (p.stderr or "").rstrip().split("\n"):
            log("    " + ln)
        if p.stdout.strip():
            log("  --- stdout ---")
            for ln in p.stdout.rstrip().split("\n"):
                log("    " + ln)
        sys.exit(1)
    missing = [o for o in outputs if not o.exists()]
    if missing:
        log(f"  [エラー] 正常終了しましたが出力がありません: {[m.name for m in missing]}")
        sys.exit(1)
    return True


def _tree_rss_mb(root_pid: int) -> float:
    """root_pid とその子孫プロセスの RSS 合計 (MB)。ps で ppid を辿る (macOS/Linux 共通)。"""
    try:
        out = subprocess.run(["ps", "-eo", "pid=,ppid=,rss="], capture_output=True, text=True).stdout
    except Exception:
        return 0.0
    kids: dict[int, list[int]] = {}
    rss: dict[int, int] = {}
    for ln in out.split("\n"):
        parts = ln.split()
        if len(parts) != 3:
            continue
        pid, ppid, r = int(parts[0]), int(parts[1]), int(parts[2])
        kids.setdefault(ppid, []).append(pid)
        rss[pid] = r
    total, stack = 0, [root_pid]
    while stack:
        pid = stack.pop()
        total += rss.get(pid, 0)
        stack.extend(kids.get(pid, []))
    return total / 1024.0


def run_measured(log: Log, step: str, n: int, total: int, name: str,
                 outputs: list[Path], args: list[str], force: bool = False,
                 interval: float = 2.0) -> tuple[bool, float]:
    """run() と同じだが、実行中のプロセス木の RSS 合計を定期サンプリングし、ピーク (MB) を返す。
    classify-sklearn の n_jobs を決めるための実測に使う。"""
    if not force and outputs and all(o.exists() for o in outputs):
        log.step(step, n, total, True, name, "既存) " + ", ".join(o.name for o in outputs))
        return False, 0.0
    log.step(step, n, total, False, name, ", ".join(o.name for o in outputs))
    cmd = ["qiime"] + args
    _log_cmd(log, cmd)
    t0 = datetime.now()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    peak = 0.0
    import time
    while p.poll() is None:
        peak = max(peak, _tree_rss_mb(p.pid))
        time.sleep(interval)
    out_s, err_s = p.communicate()
    dt = (datetime.now() - t0).total_seconds()
    if p.returncode != 0:
        log("")
        log(f"  [エラー] 失敗しました (exit={p.returncode}, {dt:.1f} 秒, ピーク RSS ≈ {peak:,.0f} MB)")
        log("  --- stderr ---")
        for ln in (err_s or "").rstrip().split("\n"):
            log("    " + ln)
        sys.exit(1)
    missing = [o for o in outputs if not o.exists()]
    if missing:
        log(f"  [エラー] 正常終了しましたが出力がありません: {[m.name for m in missing]}")
        sys.exit(1)
    return True, peak


def export_qzv(log: Log, step: str, n: int, total: int, name: str,
               qzv: Path, out_dir: Path, force: bool = False) -> bool:
    """qzv を out_dir に展開して、ブラウザで index.html を直接開けるようにする。
    出力先の名前は呼び出し側が決める (規則は README「HTML フォルダの命名規則」)。"""
    if not force and out_dir.is_dir() and (out_dir / "index.html").exists():
        log.step(step, n, total, True, name, f"既存) {out_dir.name}")
        return False
    log.step(step, n, total, False, name, out_dir.name)      # 既存のときと同じくフォルダ名を出す（260925 本人指示）
    return run(log, step, n, total, name, [out_dir / "index.html"], [
        "tools", "export", "--input-path", str(qzv), "--output-path", str(out_dir),
    ], force=True, quiet_head=True)


def _q(s: str) -> str:
    return f"'{s}'" if (" " in s or "[" in s) else s
