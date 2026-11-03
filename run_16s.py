#!/usr/bin/env python3
"""run_16s.py — メタ16S 解析パイプライン (エントリポイント)

設計の単一正本は ./README.md。判断の経緯 (A・B-1〜B-9)もそこにある。

使い方
------
    conda activate rachis-qiime2-2026.7        # ※必須。絶対パス起動は不可 (B-2)
    python ~/Claude/Scripts/16S/run_16s.py     # 実験フォルダの相対パスを対話で尋ねる

    --stage 1   import で止める (demux.qzv を見たいとき。既定は通し)
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:                                           # env の外（QIIME 2 の無い Python）で起動したら、Python のエラーでなく案内を出して止める（261015 E4b）
    import yaml  # noqa: E402

    from lib import metadata as MD  # noqa: E402
    from lib import qiime as Q  # noqa: E402
    from lib import conditions as C  # noqa: E402
    from lib import steps as ST  # noqa: E402
    from lib import figures as FG  # noqa: E402
    from lib import plotcfg as PC  # noqa: E402
    from lib import subset as SB  # noqa: E402
except ModuleNotFoundError as e:
    print(f"  [エラー] QIIME 2 2026.7 の環境に入っていません (今の Python: {sys.prefix}・見つからない部品: {e.name})。")
    print("  ※ conda activate rachis-qiime2-2026.7 (または QIIME 2 2026.7 の環境)をしてから、もう一度起動してください")
    sys.exit(1)

REQUIRED = ["input/16S/config.yaml"]            # FASTQ・SampleSheet・検体情報 CSV は config の runs:・sample_info: で指す (261022)


def prompt_experiment_dir() -> Path:
    """起動時に実験フォルダのパスを尋ね、検証して Path を返す。

    絶対パス (~ も可) で受け取る。相対パスも受け付けるが、今いるフォルダが基準になる。必要ファイルが揃わなければ
    直さずに聞き直す (存在しないパスのまま走り出さない)。
    """
    print("実験フォルダのパスを入力してください (絶対パスか ~ で始まるパス。例: ~/experiments/X01_example)")
    while True:
        raw = re.sub(r"\\(.)", r"\1", input("実験フォルダ> ").strip().strip('"').strip("'"))   # ドラッグ入力の「\ 」「\(」などの逃がしを戻す (点検 F284)
        if not raw:
            print("  ※ 空です。もう一度入力してください")
            continue
        d = Path(raw).expanduser()
        if not d.is_dir():
            print(f"  [エラー] フォルダが見つかりません: {d}")
            continue
        missing = [r for r in REQUIRED if not (d / r).exists()]
        if missing:
            print(f"  [エラー] 必要なファイルがありません: {', '.join(missing)}")
            continue
        print(f"  実験フォルダ: {d.resolve()}")
        return d.resolve()                         # ログのコマンドに絶対パスを残す (点検 F195・F282)


def main() -> None:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("exp_dir", nargs="?", help="実験フォルダ (省略時は対話で尋ねる)")
    ap.add_argument("--stage", type=int, choices=[1], help="1 = import で止める")
    a = ap.parse_args()

    log = Q.Log()
    Q.check_env(log)

    if a.exp_dir:
        exp = Path(a.exp_dir).expanduser().resolve()   # 相対パスで渡されても絶対パスにそろえる (点検 F195・F282)
        missing = [r for r in REQUIRED if not (exp / r).exists()]
        if missing:
            log(f"  [エラー] 必要なファイルがありません: {', '.join(missing)}")
            sys.exit(1)
    else:
        exp = prompt_experiment_dir()

    base = exp / "input/16S"
    cfg = yaml.safe_load(open(base / "config.yaml", encoding="utf-8")) or {}
    root16 = (exp / "output/16S").resolve()
    out = root16 / "Qiime2"                        # QIIME 2 の出力（01〜09・HTML_*）。作図は figure/・記録は logs/（260928 本人決定）
    out.mkdir(parents=True, exist_ok=True)
    logs = root16 / "logs"
    logs.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now()
    log = Q.Log(logs / f"{stamp:%y%m%d_%H%M%S}_16S.log")   # 秒まで (同じ分に 2 回起動しても混ざらない。261021 本人決定・F288)
    Q.check_env(log)
    log("=" * 72)
    log(f"開始日時     : {stamp:%Y-%m-%d %H:%M:%S}")
    log(f"パイプライン : {Q.pipeline_label()}")
    log(f"実験フォルダ : {exp.resolve()}")
    log(f"config       : {(base / 'config.yaml').resolve()}")
    log(f"出力先       : {root16}")
    log(f"QIIME 2 env  : {Q.env_label()}")
    log(f"作図の版     : {Q.lib_versions()}")
    env_file = logs / f"{stamp:%y%m%d_%H%M%S}_env.txt"
    why = Q.save_env_list(env_file)
    if why is None:
        log(f"env の一覧   : logs/{env_file.name}")
    else:
        log(f"  [注意] env の全パッケージの一覧を残せませんでした ({why})")
    log("=" * 72)

    try:
        mode = ST.check_config(cfg)             # 01〜09 の config (すべて必須・乱数の種 D11a・分類器のパス E4a を含む。261021)
        PC.check_config(ST.with_parent(cfg))    # 作図の config のうち、読むだけで分かるもの (261021 本人決定・案 B)
    except SystemExit as e:
        log(str(e))
        sys.exit(1)
    except PC.Stop as e:
        log(f"  [エラー] {e}")
        sys.exit(1)

    # 01〜10 を包み、止まったときも終了の枠と注記を必ず書く (261019 本人決定・点検 F006・F141・F224・F241)
    status = ""                                    # 空 = 最後まで進んだ
    try:
        if mode == "parent":                       # 絞り込み (イ): 01〜09 は行わず、親の結果から選んだ写しを作って 10 だけ (261026・7-1)
            ctx = dict(exp=exp, base=base, out=out, cfg=cfg, log=log)
            try:
                ctx["cfg"] = SB.build(ctx)
            except SystemExit as e:                # 止めた理由 ([エラー] の文) をログに残す
                if isinstance(e.code, str):
                    log(e.code)
                raise SystemExit(1)
            FG.step10_figures(ctx)
            return

        # ---------------- 01 Metadata ----------------
        log.head("01. Metadata (メタデータの生成)")
        md_dir = out / "01_metadata"
        md_dir.mkdir(parents=True, exist_ok=True)   # ステップ冒頭で明示（アーカイブ準拠）
        p1, p2, p3 = md_dir / "sample-metadata.txt", md_dir / "sample-metadata_cn.txt", md_dir / "sample_information.csv"
        ctx = dict(exp=exp, base=base, out=out, cfg=cfg, log=log)
        listed = {r["path"].resolve() for r in ST.runs_of(cfg)}   # 自分の runs/ にあって config の runs: に無いランは使わないことを知らせる (261026 本人決定・設計 1-1)
        extra = sorted(d.name for d in (base / "runs").glob("*") if d.is_dir() and d.resolve() not in listed) if (base / "runs").is_dir() else []
        if extra:
            log(f"  [注記] input/16S/runs/ に、config の runs: に書いていないランがあります (使いません): {', '.join(extra)}")
        C.check(ctx, "01", 1, 1, "サンプルメタデータの生成", "sample-metadata.txt, sample-metadata_cn.txt, sample_information.csv")   # 条件の照合（261014 E1a）
        if p1.exists() and p2.exists() and p3.exists():
            log.step("01", 1, 1, True, "サンプルメタデータの生成",
                     "既存) sample-metadata.txt, sample-metadata_cn.txt, sample_information.csv")
        else:
            log.step("01", 1, 1, False, "サンプルメタデータの生成",
                     "sample-metadata.txt, sample-metadata_cn.txt, sample_information.csv")
            try:
                before, after, sinfo = MD.build(ST.runs_of(cfg), ST.sample_info_of(cfg), log)
            except MD.Stop as e:
                log("")
                log(f"  [エラー] {e}")
                sys.exit(1)
            MD.write(before, after, sinfo, md_dir)
        C.record(ctx, "01")

        # ---------------- 02 Import ----------------
        ST.step02_import(ctx)
        C.record(ctx, "02")
        if a.stage == 1:
            log("")
            log("━━ --stage 1 のためここで停止します ━━")
            log("  ※ Qiime2/HTML_02_demux_<M 番号>.qzv/index.html (ランごと) を開き、Interactive Quality Plot を確認してください")
            log("  ※ trunc_len を変えるなら、input/16S/config.yaml の runs: の各ランの quality_control を編集して、再実行してください")
            status = "--stage 1 のため import で止めました"
            return

        # ---------------- 03 DADA2 ----------------
        ST.step03_dada2(ctx)
        C.record(ctx, "03")

        # ---------------- 04 Relabel ----------------
        ST.step04_relabel(ctx)
        C.record(ctx, "04")

        # ---------------- 05 Phylogeny ----------------
        ST.step05_phylogeny(ctx)
        C.record(ctx, "05")

        # ---------------- 06 Rarefaction ----------------
        ST.step06_rarefaction(ctx)
        C.record(ctx, "06")

        # ---------------- 07 Diversity ----------------
        ST.step07_diversity(ctx)
        C.record(ctx, "07")

        # ---------------- 08 Classifier ----------------
        # 分類器は作らない（B-5: SILVA 公式の学習済み分類器を使う）。番号の抜けが分かるようにスキップと記録する（260928 本人指示）
        log.head("08. Classifier (分類器の作成)")
        log.step("08", 1, 1, True, "分類器の作成", "既存) (config の classifier.path)")
        ctx["classifier"] = ST.verify_classifier(ctx)   # 使う分類器の確認と記録（260928: 09 から移した）

        # ---------------- 09 Taxonomy ----------------
        ST.step09_taxonomy(ctx)
        C.record(ctx, "09")

        # ---------------- 10 Figures ----------------
        FG.step10_figures(ctx)
    except SystemExit as e:
        if e.code not in (None, 0):
            status = "途中で止まりました (上の [エラー] を参照してください)"
        raise
    except Exception as e:                         # 想定外のエラーは詳細 (traceback) もログに残す
        import traceback
        for ln in traceback.format_exc().rstrip().split("\n"):
            log("    " + ln)
        status = f"想定外のエラーで止まりました: {type(e).__name__}: {e}"
        raise
    finally:
        fin = datetime.now()
        sec = int((fin - stamp).total_seconds())
        end = ["=" * 72,                              # 冒頭の「開始日時」のブロックと対にする（260927 本人指定）
               f"終了日時     : {fin:%Y-%m-%d %H:%M:%S}",
               f"所要時間     : {sec // 60} 分 {sec % 60} 秒"]
        if status:
            end.append(f"状態         : {status}")
        end += ["-" * 65] + log.summary() + ["=" * 72]
        for line in end:
            log(line)


if __name__ == "__main__":
    main()
