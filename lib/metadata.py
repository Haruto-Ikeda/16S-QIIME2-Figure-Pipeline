#!/usr/bin/env python3
"""metadata.py — ラン・SampleSheet・検体情報 CSV を突き合わせて sample-metadata を作る (01。261022 に新しい設計へ)。

入力 (すべて人が置く。config の runs:・sample_info: で指す。設計の正本 README「FASTQ の置き方」)
  runs/<YYMMDD>_M<4 桁>/fastq/          FASTQ (名前は M<4 桁>P<4 桁>N<4 桁>_S…_L001_R1_001.fastq.gz)
  runs/<YYMMDD>_M<4 桁>/SampleSheet.txt そのランの検体表 (タブ区切り・1 列目 #SampleID)
  <YYMMDD>_<P 番号 4 桁>_….csv          検体情報 (2 行目 Project No.・7 行目が列名・No. / Sample_name / Group_name / Subgroup_name)

出力 (output/16S/Qiime2/01_metadata/)
  sample-metadata.txt     1 列目 #SampleID = FASTQ の ID。Run・Project・SampleSheet の列・群・newID (右端)
  sample-metadata_cn.txt  1 列目 #SampleID = newID (検体名)。ほかは上と同じ並び (Sample_name・newID を除く)
  sample_information.csv  解析した検体の検体情報 CSV の全列 (左は上と同じ列・No.)。10 がリード数を足して 00_info/ に写す
                          (CSV を読むのは 01 の 1 回だけにして、図と表の食い違いを防ぐ。261024 決定・案 B)
  行は ID の順 (M → P → N)

CSV の読み取り仕様はコードに固定する (config に書かない)。形が変わったら黙って追随せず、止めて気づけるようにするため。
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

# ---- CSV の固定仕様（config に出さない。変わったら止める）----
CSV_HEADER_ROW = 7          # 1 始まり（260925 に 6 → 7: 修正した CSV が正本。6 行目は空行）
CSV_PROJECT_ROW = 2         # 2 行目 = 「Project No.,<番号>」(261021 決定・3-4。ファイル名の番号と照らし合わせる)
COL_NO = "No."
COL_NAME = "Sample_name"
COL_GROUP = "Group_name"
COL_SUBGROUP = "Subgroup_name"

ID_RE = re.compile(r"^M(\d{4})P(\d{4})N(\d{4})$")   # 検体の ID (261021 決定・②)

# ---- QIIME 2 の識別子ルール ----
ALLOWED = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.-_:|")
MAX_LEN = 36
RESERVED_CI = {"id", "sampleid", "sample id", "sample-id",
               "featureid", "feature id", "feature-id"}
RESERVED_CS = {"#SampleID", "#Sample ID", "#OTUID", "#OTU ID", "sample_name"}
# QIIME 2 が数として読む名前 (qiime2.Metadata.load で 36 通りを照らした。float() は nan・inf・1_000 でずれるので使わない。261102)
NUMERIC_RE = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")
UTF8_SHOW_MAX = 5           # UTF-8 でない行を、エラーに最大何行出すか (261102 決定・案 2)
UTF8_SHOW_WIDTH = 30        # [ここ] の前後に何文字出すか


class Stop(Exception):
    """検証に失敗して止めるときに投げる。呼び出し側がメッセージを出す。"""


def _die(msg: str) -> None:
    raise Stop(msg)


def _at(p: Path) -> str:
    """場所の示し方: ラン名/ファイル (例 250101_M0001/SampleSheet.txt)。無いものを指すときだけ絶対パスを使う (261023 決定・案 a)。"""
    return f"{p.parent.name}/{p.name}"


def _list(items) -> str:
    """値を「a, b, c」と並べる (リストの [ ] や ' ' を出さない。261023 決定)。"""
    return ", ".join(str(x) for x in items)


def _lines(items) -> str:
    return "".join(f"       - {x}\n" for x in items)


def _not_utf8_where(path: Path) -> str:
    """UTF-8 でない文字のある行を、行番号と前後の文字で示す (最大 UTF8_SHOW_MAX 行と合計の行数。261102 決定・案 2)。
    行番号はファイルの行で数える (セルの中の改行があると Excel の行番号とずれる)。"""
    found: list[str] = []
    for no, raw in enumerate(re.split(rb"\r\n|\r|\n", path.read_bytes().removeprefix(b"\xef\xbb\xbf")), 1):
        try:
            raw.decode("utf-8")
        except UnicodeDecodeError as e:
            head = raw[:e.start].decode("utf-8", "replace")[-UTF8_SHOW_WIDTH:]
            tail = raw[e.end:].decode("utf-8", "replace")[:UTF8_SHOW_WIDTH]
            found.append(f"     {no} 行目: …{head}[ここ]{tail}…\n")
    return "".join(found[:UTF8_SHOW_MAX]) + f"     合計 {len(found)} 行\n"


# ---------------------------------------------------------------- FASTQ
def fastq_ids(fastq_dir: Path) -> dict[str, dict[str, Path]]:
    """fastq/ の FASTQ を {ID: {"R1": ファイル, "R2": ファイル}} にする。ID は QIIME 2 の Casava 規則と同じく右から '_' で区切る (点検 F191)。
    FASTQ 以外・Mac の隠しファイル (._*) があれば止める (点検 F245。読み込みに混ざると分かりにくく止まるため)。"""
    if not fastq_dir.is_dir():
        _die(f"FASTQ のフォルダがありません: {fastq_dir}")
    stray = sorted(p.name for p in fastq_dir.iterdir() if not p.name.endswith(".fastq.gz") or p.name.startswith("."))
    if stray:
        _die(f"FASTQ のフォルダに FASTQ 以外のファイルがあります ({_at(fastq_dir)}):\n" + _lines(stray)
             + "     ._ で始まる Mac の隠しファイルも含めて、ゴミ箱へ移してから、もう一度実行してください (自動では消しません)。")
    ids: dict[str, dict[str, Path]] = {}
    for f in sorted(fastq_dir.glob("*.fastq.gz")):
        parts = f.name.removesuffix(".fastq.gz").rsplit("_", 4)
        if len(parts) != 5 or not re.fullmatch(r"L\d{3}", parts[2]) or parts[3] not in ("R1", "R2") or parts[4] != "001":
            _die(f"FASTQ の名前が <ID>_S…_L001_R1_001.fastq.gz の形ではありません: {f.name} ({_at(fastq_dir)})")
        if parts[3] in ids.setdefault(parts[0], {}):
            _die(f"同じ検体の {parts[3]} の FASTQ が 2 つあります (レーンが複数なら今は扱えません): {ids[parts[0]][parts[3]].name} と {f.name} ({_at(fastq_dir)})")
        ids[parts[0]][parts[3]] = f
    if not ids:
        _die(f"FASTQ が 1 つもありません ({_at(fastq_dir)})")
    lacking = sorted(i for i, d in ids.items() if set(d) != {"R1", "R2"})
    if lacking:
        _die(f"R1 と R2 がそろっていない FASTQ があります ({_at(fastq_dir)}):\n" + _lines(lacking))
    return ids


# ---------------------------------------------------------------- SampleSheet
def read_sample_sheet(path: Path) -> pd.DataFrame:
    """ランの検体表を文字のまま読む。同じ #SampleID が 2 行あれば、使わない行 (相乗りしたほかの方の検体) でも止める
    (261021 決定・④-2・点検 F209。261023 に「使う行だけ」も検討したが、シートの壊れに気づけるよう全行のまま)。"""
    try:
        ss = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8-sig")
    except UnicodeDecodeError:                                    # 261102: 場所を出す (A2)
        _die(f"SampleSheet.txt を読めません: UTF-8 でない文字があります ({_at(path)})\n" + _not_utf8_where(path)
             + "     この文字を消すか打ち直して、UTF-8 で保存し直してください (自動で文字コードを切り替えません)")
    except pd.errors.ParserError as e:                            # 261027: CSV と同じく案内して止める
        _die(f"SampleSheet.txt を読めません: {type(e).__name__} ({_at(path)})\n"
             "     UTF-8 で保存し直してください (自動で文字コードを切り替えません)")
    if ss.columns[0] != "#SampleID":
        _die(f"SampleSheet.txt の 1 列目が '#SampleID' ではありません: '{ss.columns[0]}' ({_at(path)})")
    ss["#SampleID"] = ss["#SampleID"].str.strip()
    dup = sorted(set(ss["#SampleID"][ss["#SampleID"].duplicated()]))
    if dup:
        _die(f"SampleSheet.txt に同じ #SampleID が 2 行以上あります ({_at(path)}):\n" + _lines(dup))
    return ss


# ---------------------------------------------------------------- 検体情報 CSV
def _read_csv(path: Path, **kw) -> pd.DataFrame:
    try:
        return pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False, **kw)
    except UnicodeDecodeError:                                    # 261102: 場所を出す (A2)
        _die(f"検体情報 CSV を読めません: UTF-8 でない文字があります ({path.name})\n" + _not_utf8_where(path)
             + "     この文字を消すか打ち直して、UTF-8 (Excel の「CSV UTF-8」) で保存し直してください (自動で文字コードを切り替えません)")
    except pd.errors.ParserError as e:                            # 点検 F190・F192・F246
        _die(f"検体情報 CSV を読めません: {type(e).__name__} ({path.name})\n"
             "     CSV を UTF-8 (Excel の「CSV UTF-8」) で保存し直してください (自動で文字コードを切り替えません)")


def read_sample_info(path: Path, pnum: int, log) -> pd.DataFrame:
    """検体情報 CSV を固定仕様で読む。2 行目の Project No. がファイル名の番号と違えば止める (261021 決定・3-4)。
    すべて文字のまま読み (点検 F008・F211)、No. の重複・整数でない No.・空欄の名前は止める (点検 F210・F223・F211)。"""
    head = _read_csv(path, header=None, nrows=CSV_PROJECT_ROW)
    row = [str(x).strip() for x in head.iloc[CSV_PROJECT_ROW - 1].tolist()]
    if row[0] != "Project No.":
        _die(f"検体情報 CSV の {CSV_PROJECT_ROW} 行目が「Project No.」ではありません: {_list(row[:2])} ({path.name})")
    if not row[1].isdigit() or int(row[1]) != pnum:
        _die(f"検体情報 CSV のファイル名の P 番号 ({pnum:04d}) と、{CSV_PROJECT_ROW} 行目の Project No. ({row[1]}) が違います ({path.name})")
    df = _read_csv(path, header=CSV_HEADER_ROW - 1)
    df.columns = [str(c).strip() for c in df.columns]
    need = [COL_NO, COL_NAME, COL_GROUP, COL_SUBGROUP]
    lack = [c for c in need if c not in df.columns]
    if lack:
        _die(f"CSV に必要な列がありません: {_list(lack)} ({path.name})\n"
             f"     {CSV_HEADER_ROW} 行目を列名として読みました。実際の列: {_list(list(df.columns)[:8])} …\n"
             f"     列名が変わった場合は、勝手に追随せず仕様を見直してください。")
    raw = df.copy()
    df = df[df[COL_NO].str.strip() != ""].copy()                    # No. が空の行 (末尾の空行など) は読まない
    bad = sorted(df.loc[~df[COL_NO].str.strip().str.fullmatch(r"\d+"), COL_NO])
    if bad:
        _die(f"検体情報 CSV の No. に整数でない値があります: {_list(bad)} ({path.name})")
    df[COL_NO] = df[COL_NO].str.strip().astype(int)
    dup = sorted(set(df[COL_NO][df[COL_NO].duplicated()]))
    if dup:
        _die(f"検体情報 CSV に同じ No. が 2 行以上あります: No. {_list(dup)} ({path.name})")
    warn_whitespace(raw[raw[COL_NO].str.strip() != ""], [COL_NAME, COL_GROUP, COL_SUBGROUP], log)
    for c in (COL_NAME, COL_GROUP, COL_SUBGROUP):
        df[c] = df[c].str.strip()                                     # QIIME 2 は前後の空白を無視する。合わせる
    empty = df.loc[(df[[COL_NAME, COL_GROUP, COL_SUBGROUP]] == "").any(axis=1), COL_NO].tolist()
    if empty:
        _die(f"検体情報 CSV に、検体名・群・サブグループが空欄の行があります: No. {_list(empty)} ({path.name})")
    return df.assign(Project=f"P{pnum:04d}")                       # 全列 (残りの列は sample_information.csv に使う。261024)


# ---------------------------------------------------------------- 検証
def validate(df: pd.DataFrame, columns: list[str], log, id_col: str = "#SampleID") -> None:
    """識別子になりうる値を検証する。違反があれば止める。"""
    bad: list[str] = []
    for col in columns:
        vals = sorted(df[col].astype(str).unique())
        for v in vals:
            ng = sorted({c for c in v if c not in ALLOWED})
            if ng:
                bad.append(f"{col}「{v}」: 使えない文字 {_list(ng)}")
            if len(v) > MAX_LEN:
                bad.append(f"{col}「{v}」: {len(v)} 文字 ({MAX_LEN} 文字以内)")
            if v.startswith("#"):
                bad.append(f"{col}「{v}」: # で始まっています")
    # 検体名は 04 で ID になるので、予約語も使えない (点検 F179)
    hits = sorted(v for v in df[COL_NAME] if v.lower() in RESERVED_CI or v in RESERVED_CS)
    if hits:
        bad.append(f"{COL_NAME}: QIIME 2 の予約語です: {_list(hits)}")
    # 全部が数だと 04 の feature-table group が受け付けない。混在は 04 を通るが、数の名前だけが残ると作図で名前が合わなくなる
    # (pandas が 727 を数に読む)。1 つでも数なら止める (261102 決定・案 2)
    hits = sorted(v for v in df[COL_NAME] if NUMERIC_RE.fullmatch(v))
    if hits:
        bad.append(f"{COL_NAME}: 数として読める検体名があります: {_list(hits)} (数字以外の文字を 1 文字以上入れてください。例 _727・M727)")
    # 予約語（列名）。ID 列だけは予約語を「使わなければならない」ので対象外にする。
    hits = [c for c in df.columns if c != id_col and (c.lower() in RESERVED_CI or c in RESERVED_CS)]
    if hits:
        bad.append(f"ID 列以外に QIIME 2 の予約語が使われています: {_list(hits)}")
    if bad:
        _die("QIIME 2 が識別子に使えない値があります。CSV を修正してください。\n" + _lines(bad)
             + f"     使える文字: {''.join(sorted(ALLOWED))}\n"
             + f"     {MAX_LEN} 文字以内・# で始めない・重複なし")


def warn_whitespace(raw: pd.DataFrame, columns: list[str], log) -> None:
    """前後の空白は注意のみ (QIIME 2 は無視するが pandas は無視しないため)。"""
    for col in columns:
        v = raw[col].astype(str)
        ws = sorted(set(v[v != v.str.strip()]))
        if ws:
            log(f"  [注意] {col} の前後に空白がある値がありました: {_list(ws)} (読み込み時に除去しました)")


# ---------------------------------------------------------------- 生成
def build(runs: list[dict], infos: dict[int, dict], log) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ラン・SampleSheet・検体情報 CSV を突き合わせて (変更前, 変更後) の metadata を返す (261021〜22 決定・①〜⑥)。
    runs = steps.runs_of(cfg)、infos = steps.sample_info_of(cfg) ({P 番号: {path, absent_ok}})。"""
    log("-" * 65)
    # 検体情報 CSV を先に読む (exclude・absent_ok の検体名をログに出すため。261022 決定)
    info = pd.concat([read_sample_info(v["path"], p, log) for p, v in sorted(infos.items())], ignore_index=True)
    info["_P"] = info["Project"].str[1:].astype(int)

    def name_of(p: int, n: int) -> str:
        hit = info.loc[(info["_P"] == p) & (info[COL_NO] == n), COL_NAME]
        return hit.iloc[0] if len(hit) else "(検体情報 CSV に無い検体)"

    runs = sorted(runs, key=lambda r: r["m"])                   # ランは M の順 (ログも表の ID の順とそろえる。261023 決定)
    rows, seen_fastq = [], []                                    # seen_fastq = FASTQ がある (P, N) (exclude で外したものも含む)
    for r in runs:
        ids = fastq_ids(r["fastq"])
        bad = sorted(i for i in ids if not ID_RE.fullmatch(i))
        if bad:
            _die(f"FASTQ の ID が M<4 桁>P<4 桁>N<4 桁> の形ではありません ({_at(r['fastq'])}):\n" + _lines(bad))   # ④-3
        other = sorted(i for i in ids if ID_RE.fullmatch(i).group(1) != r["m"][1:])
        if other:
            _die(f"ほかのランの FASTQ ({r['m']} 以外) があります ({_at(r['fastq'])}):\n" + _lines(other))   # ④-4
        miss_ex = sorted(set(r["exclude"]) - set(ids))
        if miss_ex:
            _die(f"config の exclude に、このランに無い ID があります ({r['name']}):\n" + _lines(miss_ex))
        ss = read_sample_sheet(r["sheet"])
        not_in = sorted(set(ids) - set(ss["#SampleID"]))
        if not_in:
            _die(f"SampleSheet.txt に載っていない FASTQ があります ({_at(r['sheet'])}):\n" + _lines(not_in))   # ④-1
        seen_fastq += [tuple(int(x) for x in ID_RE.fullmatch(i).groups()[1:]) for i in ids]
        use = sorted(set(ids) - set(r["exclude"]))
        sel = ss[ss["#SampleID"].isin(use)].copy()
        sel.insert(1, "Run", r["m"])
        rows.append(sel)
        log(f"{r['name']}: FASTQ {sum(len(v) for v in ids.values())} ファイル (検体 {len(ids)} 件)"
            f"・SampleSheet {len(ss)} 行のうち {len(sel)} 行を使う")
        if r["exclude"]:                                         # ID と検体名を 1 行ずつ (261022 決定)
            log(f"  exclude {len(r['exclude'])} 件:")
            for x in sorted(r["exclude"]):
                p, n = (int(v) for v in ID_RE.fullmatch(x).groups()[1:])
                log(f"    {x}  {name_of(p, n)}")
        else:
            log("  exclude: なし")
    m = pd.concat(rows, ignore_index=True)
    pn = m["#SampleID"].str.extract(ID_RE)
    m["_P"], m["_N"] = pn[1].astype(int), pn[2].astype(int)

    # 同じ検体 (P・N) が 2 つのランに残っていないか (読み足し。6-7)
    dup = m[m.duplicated(["_P", "_N"], keep=False)].sort_values(["_P", "_N"])
    if len(dup):
        _die("同じ検体 (P・N) が 2 つのランにあります (読み足しなら、どちらかのランの exclude に書いてください):\n"
             + _lines(f"P{p:04d}N{n:04d}: " + ", ".join(sorted(g["#SampleID"])) for (p, n), g in dup.groupby(["_P", "_N"])))

    # 検体情報 CSV を P 番号で選び、N で No. を引く (③)
    noinfo = sorted(set(m["_P"]) - set(infos))
    if noinfo:
        _die("FASTQ の P 番号に当たる検体情報 CSV が config の sample_info: にありません: " + ", ".join(f"P{p:04d}" for p in noinfo))
    # ログ: CSV ごとの absent_ok を、No. と検体名で 1 行ずつ (261022 決定)
    for p, v in sorted(infos.items()):
        log(f"検体情報 CSV (P{p:04d}): {v['path'].name}・{int((info['_P'] == p).sum())} 行")
        if v["absent_ok"]:
            log(f"  absent_ok {len(v['absent_ok'])} 件 (FASTQ が無くてよい検体):")
            for n in sorted(v["absent_ok"]):
                log(f"    N{n:04d}  {name_of(p, n)}")
        else:
            log("  absent_ok: なし")
    names = info[COL_NAME]
    dupn = sorted(set(names[names.duplicated()]))
    if dupn:
        _die(f"検体情報 CSV に同じ検体名 ({COL_NAME}) が 2 行以上あります (CSV をまたぐものを含む。04 で別の検体が合算されるため):\n"
             + _lines(dupn))                                         # 3-5
    m = m.merge(info[["_P", COL_NO, COL_NAME, COL_GROUP, COL_SUBGROUP]], left_on=["_P", "_N"], right_on=["_P", COL_NO], how="left")
    lost = m[m[COL_NAME].isna()]
    if len(lost):
        _die("検体情報 CSV に対応する No. が無い FASTQ があります:\n"
             + _lines(f"{r['#SampleID']} (P{r['_P']:04d} の No. {r['_N']})" for _, r in lost.iterrows()))

    # CSV に載っているのに FASTQ が無い検体は、absent_ok に書いたものだけ許す (6-8)
    have = set(seen_fastq)
    for p, v in sorted(infos.items()):
        nos = set(info.loc[info["_P"] == p, COL_NO])
        unknown = sorted(v["absent_ok"] - nos)
        if unknown:
            _die("absent_ok に、検体情報 CSV に無い No. があります: " + ", ".join(f"N{n:04d}" for n in unknown) + f" ({v['path'].name})")
        both = sorted(n for n in v["absent_ok"] if (p, n) in have)
        if both:
            _die("absent_ok に書いた検体に FASTQ があります: " + ", ".join(f"N{n:04d}" for n in both) + f" ({v['path'].name})\n"
                 "     使わないなら absent_ok から消し、ランの exclude に書いてください")
        absent = sorted(n for n in nos if (p, n) not in have and n not in v["absent_ok"])
        if absent:
            _die(f"検体情報 CSV に載っているのに、どのランにも FASTQ が無い検体があります ({v['path'].name}):\n"
                 + _lines(f"N{n:04d} ({info.loc[(info['_P'] == p) & (info[COL_NO] == n), COL_NAME].iloc[0]})" for n in absent)
                 + "     FASTQ の入れ忘れでなければ、sample_info: のその CSV の absent_ok に書いてください")

    m = m.sort_values("#SampleID").reset_index(drop=True)          # ID の順 (M → P → N。6-9)
    m.insert(2, "Project", m["_P"].map(lambda p: f"P{p:04d}"))
    validate(m, [COL_NAME, COL_GROUP, COL_SUBGROUP], log)
    log(f"検体: {len(m)} 件 (ラン {len(runs)}・プロジェクト {m['Project'].nunique()})・文字種・文字数・先頭文字・重複・予約語の検証: クリア")
    log("-" * 65)

    # 変更前: #SampleID・Run・Project・SampleSheet の列・群・newID (右端)
    before = m.drop(columns=["_P", "_N", COL_NO]).copy()
    before["newID"] = m[COL_NAME]
    # 変更後: 1 列目が newID。ほかは変更前と同じ並び (Sample_name・newID は 1 列目と同じなので外す。261022 決定・案 a)
    after = before.drop(columns=[COL_NAME, "newID"]).assign(**{"#SampleID": before["newID"]})
    # 検体情報: 変更前の列 (newID を除く・Sample_name の前に No.) ＋ CSV の残りの列。解析した検体だけ・ID の順 (261024 決定)
    left = [c for c in before.columns if c != "newID"]
    left.insert(left.index(COL_NAME), COL_NO)
    rest = [c for c in info.columns if c not in left and c not in ("_P", "Project")]
    sinfo = m[left + ["_P"]].merge(info[["_P", COL_NO] + rest], on=["_P", COL_NO], how="left").drop(columns="_P")
    return before, after, sinfo.fillna("")


def write(before: pd.DataFrame, after: pd.DataFrame, sinfo: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    before.to_csv(out_dir / "sample-metadata.txt", sep="\t", index=False, na_rep="")
    after.to_csv(out_dir / "sample-metadata_cn.txt", sep="\t", index=False, na_rep="")
    sinfo.to_csv(out_dir / "sample_information.csv", index=False, encoding="utf-8-sig")   # Excel でそのまま開けるように BOM 付き
