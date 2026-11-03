#!/usr/bin/env python3
"""groupstats.py — 群間の検定 (stats.method に従う)。#1・#3・#10 で使う (#8 の星も)。

メタボローム (Scripts/figure/metabolome/lib/stats.py)と同じ規則:
  tukey   = 全ペア。Tukey HSD (scipy.stats.tukey_hsd)
  dunnett = 各群 vs 対照群。群が 3 つ以上なら Dunnett (scipy.stats.dunnett)、2 つなら Welch の t 検定
  welch   = 群が 2 つのとき。Welch の t 検定 (等分散を仮定しない。261009 本人決定)。
            config で群が 2 つのときは welch しか選べない (plotcfg で止める)。dunnett の 2 群は、n < 2 の群を外して
            結果的に 2 群になったときだけ起きる
どちらも手法の中で多重比較を補正するので、BH などは重ねない。
Dunnett は多変量 t 分布を準モンテカルロで積分するため、種を固定して決定論的にする。
n < min_n (stats.min_n・必須)の群は検定から外し、外したことを返り値と rows (p = NaN・test = not_tested)で知らせる。
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy import stats

DUNNETT_SEED = 0


def compare(values: dict[str, np.ndarray], order: list[str], method: str, control: str | None = None,
            min_n: int = 2) -> dict:
    """返り値: {"rows": [{groupA, groupB, n_A, n_B, p, test, reason}], "skipped": [群], "test": 検定名, "min_n": min_n}
    n < min_n の群は検定から外す (stats.min_n。必須。261009 本人決定)。外した群が関わる組も rows に残し、
    p は NaN・test は not_tested・reason は "n < {min_n}" にする (図・CSV で「検定していない」と分かるように)。
    tukey で外した結果 2 群しか残らないときは、2 群の Tukey (= Student の t)ではなく Welch の t に切り替える (D1 とそろえる)。"""
    usable = [g for g in order if len(values.get(g, [])) >= min_n]
    skipped = [g for g in order if g not in usable]
    rows = []
    reason = f"n < {min_n}"

    def n_of(g):
        return len(values.get(g, []))

    def untested(pairs):
        for a, b in pairs:
            rows.append(dict(groupA=a, groupB=b, n_A=n_of(a), n_B=n_of(b), p=np.nan, test="not_tested", reason=reason))

    def welch(a, b):
        r = stats.ttest_ind(np.asarray(values[a], float), np.asarray(values[b], float), equal_var=False)
        rows.append(dict(groupA=a, groupB=b, n_A=n_of(a), n_B=n_of(b), p=float(r.pvalue), test="welch_t", reason=""))

    out = {"rows": rows, "skipped": skipped, "min_n": min_n}
    if method == "tukey":
        pairs = list(itertools.combinations(order, 2))
        untested([(a, b) for a, b in pairs if a in skipped or b in skipped])
        if len(usable) == 2:                       # 外した結果 2 群だけ残った → Welch
            welch(*usable)
            return {**out, "test": "welch_t", "fallback": True}
        if len(usable) > 2:
            res = stats.tukey_hsd(*[np.asarray(values[g], float) for g in usable])
            for i, j in itertools.combinations(range(len(usable)), 2):
                rows.append(dict(groupA=usable[i], groupB=usable[j], n_A=n_of(usable[i]), n_B=n_of(usable[j]),
                                 p=float(res.pvalue[i, j]), test="tukey_hsd", reason=""))
        return {**out, "test": "tukey_hsd"}

    if method == "welch":
        a, b = order[0], order[1]
        if skipped:
            untested([(a, b)])
        else:
            welch(a, b)
        return {**out, "test": "welch_t"}

    if method == "dunnett":
        others = [g for g in order if g != control]
        if control in skipped:
            untested([(control, g) for g in others])
            return {**out, "test": "dunnett"}
        untested([(control, g) for g in others if g in skipped])
        ok = [g for g in others if g in usable]
        cvals = np.asarray(values[control], float)
        if len(ok) == 1:
            welch(control, ok[0])
            test = "welch_t"
        elif ok:
            r = stats.dunnett(*[np.asarray(values[g], float) for g in ok], control=cvals, rng=DUNNETT_SEED)
            for g, p in zip(ok, r.pvalue):
                rows.append(dict(groupA=control, groupB=g, n_A=len(cvals), n_B=n_of(g), p=float(p), test="dunnett", reason=""))
            test = "dunnett"
        else:
            test = "dunnett"
        return {**out, "test": test}

    raise ValueError(f"未知の検定法: {method}")


def letters(groups: list[str], rows: list[dict], alpha: float = 0.05) -> dict[str, str]:
    """全ペアの検定結果から compact letter display (CLD)を作る (Piepho 2004 の insert-and-absorb)。
    同じ文字を 1 つでも共有する群どうしは有意差なし (p >= alpha)。p が NaN の組 (検定していない)は使わない。
    文字は group の並び順で、最初に現れる群から a, b, c … の順に付ける。返り値: {群: 文字列}。"""
    cols = [set(groups)]
    sig = [(r["groupA"], r["groupB"]) for r in rows
           if r["groupA"] in groups and r["groupB"] in groups and np.isfinite(r["p"]) and r["p"] < alpha]
    for a, b in sig:
        new = []
        for c in cols:
            if a in c and b in c:
                new += [c - {a}, c - {b}]
            else:
                new.append(c)
        cols = [c for c in new if c and not any(c < o for o in new)]          # 他の列に含まれる列を吸収
        cols = [c for i, c in enumerate(cols) if c not in cols[:i]]          # 重複を除く
    pos = {g: i for i, g in enumerate(groups)}
    cols.sort(key=lambda c: min(pos[g] for g in c))
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    out = {g: "" for g in groups}
    for i, c in enumerate(cols):
        for g in c:
            out[g] += alphabet[i] if i < len(alphabet) else "?"
    return out
