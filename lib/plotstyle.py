#!/usr/bin/env python3
"""plotstyle.py — 作図の体裁 (Skills/figure-plotting の規約に揃える)。

設計の正本は ../README.md「作図パイプラインの設計」節。
**コードは Scripts/figure/ と共有しない** (260917 決定・規約の例外)が、
数値と描き方は Scripts/figure/plot_summary.py に合わせてある:
  固定 mm レイアウト (tight_layout / bbox_inches="tight" は使わない)・軸高 48 mm・前景 #222222・
  フォント 基準 9 / 軸ラベル 10 / 星 9・箱の幅 0.55・群色の透明度 0.55・点は前景色。
"""
from __future__ import annotations

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from pathlib import Path  # noqa: E402

MM = 1 / 25.4
GREY = "#222222"
PLOT_H_MM = 48                 # 軸高（全図共通。260923 決定: 46 → 48 mm。上位 12 タクソン × 4 mm = 48 mm でヒートマップとそろう）
# 文字の大きさ（260921 決定・規約 Skills/figure-plotting §3）: 軸ラベル 8・目盛 8・星 9・群名 8・タイトル 9・凡例 8・カラーバー 8・注記 7・マスの中 8
FS_TICK, FS_LABEL, FS_STAR, FS_NOTE = 8, 8, 9, 7
FS_LEGEND, FS_CBAR, FS_CELL = 8, 8, 8     # 注記 7・カラーバー 8（260923 に 7 → 8: 軸ラベル・目盛と同じ）・行列のマスの中 8
# 余白は全図で固定（260921 決定・規約全体）: 左 15・右 4・上 13・下 20 mm（上は 260923 に 5 → 13: #10 の長いタクソン名を 3 行まで折り返して収めるため。下は 260921 に 15 → 20）。凡例は図に入れず別ファイルに出す
MARGIN_L, MARGIN_R, MARGIN_T, MARGIN_B = 15, 4, 13, 20
# 横軸は「目盛 1 つあたりの幅」を固定（260921 決定・確定。規約 Skills/figure-plotting §3）
MM_PER_GROUP = 8              # 群の図: 1 群 = 8 mm。群名は 45° に傾ける
READS_PER_UNIT, MM_PER_UNIT = 5000, 16  # #2 の横軸: 5,000 リード = 16 mm・目盛は 5,000 刻み（260929 決定で 2,000 リード 8 mm から変更。目盛を少なくするため）
XLABEL_ROT = 45
XLABEL_FS = 8
STAR_PAD_PT = -4.0             # 星と横棒の隙間（plot_summary と同じ）
JITTER_SEED = 0                # 点の横ずらしの乱数（決定論的にする）

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
    "font.size": FS_TICK, "axes.labelsize": FS_LABEL, "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "axes.linewidth": 0.8,
    "text.color": GREY, "axes.edgecolor": GREY, "axes.labelcolor": GREY,
    "xtick.color": GREY, "ytick.color": GREY,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.default": "regular",   # 下付きなどの数式部分も本文と同じ Helvetica で描く
})


def p_stars(p: float) -> str:
    """*** p<0.001 / ** p<0.01 / * p<0.05 / それ以外は空 (n.s. は書かない)。"""
    if p is None or not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def text_mm(s: str, fontsize: float) -> float:
    """文字列の幅を実測して mm で返す (ラベルの重なりを実測で避けるため)。"""
    fig = plt.figure()
    t = fig.text(0, 0, s, fontsize=fontsize)
    fig.canvas.draw()
    w = t.get_window_extent().width / fig.dpi * 25.4
    plt.close(fig)
    return w


def nice_step(span: float) -> float:
    """span を 4〜7 目盛で割り切れる丸い刻み (1/2/2.5/5 ×10^k)。"""
    if span <= 0:
        return 1.0
    mag = 10.0 ** np.floor(np.log10(span / 5.0))
    best = None
    for m in (1, 2, 2.5, 5, 10):
        st = m * mag
        if 4 <= span / st <= 7:
            return st
        if best is None or abs(span / st - 5) < abs(span / best - 5):
            best = st
    return best


def choose_ylim(dlo: float, dhi: float, n_levels: int) -> tuple[float, float, float]:
    """y 軸の (下端, 上端, 刻み)。Scripts/figure/plot_summary.py の choose_ylim と同じ規則:
    正値のみなら 0 基底 / 有意差ブラケットの段数ぶん上端を軸高の一定割合で確保 / 目盛は丸い刻み。"""
    span = max(dhi - dlo, 1e-9)
    lo = dlo - span * 0.08 if dlo < 0 else 0.0
    frac = 0.0 if n_levels == 0 else min(0.24, 0.035 * n_levels + 0.03)
    hi0 = dhi + span * 0.08
    hi = hi0 if frac == 0 else (hi0 - frac * lo) / (1.0 - frac)
    step = nice_step(hi - lo)
    if lo != 0.0:
        lo = np.floor(lo / step) * step
        hi = hi0 if frac == 0 else (hi0 - frac * lo) / (1.0 - frac)
        step = nice_step(hi - lo)
    return lo, hi, step


def draw_brackets(ax, xpos: dict, sig_pairs: list, data_top: float, lo: float, hi: float) -> float:
    """有意ペアのブラケットを data_top より上に積む。収まらなければ上端を伸ばして返す。
    (Scripts/figure/plot_summary.py の draw_brackets と同じ描き方)
    要素は (a, b, 文字) か (a, b, 文字, 線を引くか)。線を引かない要素は「n < 3」など検定しなかった印で、
    文字だけをブラケットと同じ段の高さに書く (a == b なら群の真上。261009 決定)。"""
    if not sig_pairs:
        return hi
    n = len(sig_pairs)
    ax.figure.canvas.draw()
    ax_pt = ax.get_window_extent().height * 72.0 / ax.figure.dpi
    for _ in range(3):
        per_pt = (hi - lo) / max(ax_pt, 1e-9)
        h = min((hi - data_top) / (n + 0.9) * 0.18, 2.2 * per_pt)
        star_pt = FS_STAR * 1.05 + STAR_PAD_PT
        step = max((hi - data_top) / (n + 0.9), h + star_pt * per_pt)
        need = data_top + step * (n - 1 + 0.55) + h + (star_pt + 1.0) * per_pt
        if need <= hi:
            break
        hi = need
    ax.set_ylim(lo, hi)
    per_pt = (hi - lo) / max(ax_pt, 1e-9)
    ax._brackets = dict(xpos=xpos, pairs=sig_pairs, data_top=data_top, ax_pt=ax_pt,
                        step_pt=step / per_pt, h_pt=h / per_pt, artists=[])   # 軸を伸ばしたあと同じ見た目で描き直すため (261020・F122)
    _place_brackets(ax, step, h)
    return hi


def _place_brackets(ax, step: float, h: float) -> None:
    b = ax._brackets
    for art in b["artists"]:
        art.remove()
    b["artists"] = []
    xpos = b["xpos"]
    for level, item in enumerate(sorted(b["pairs"], key=lambda s: abs(xpos[s[0]] - xpos[s[1]]))):
        a, bb, star = item[:3]
        line = item[3] if len(item) > 3 else True
        x1, x2 = sorted([xpos[a], xpos[bb]])
        y = b["data_top"] + step * (level + 0.55)
        if line:
            b["artists"] += ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.8, color=GREY, clip_on=False)
        b["artists"].append(ax.annotate(star, xy=((x1 + x2) / 2, y + h), textcoords="offset points", xytext=(0, STAR_PAD_PT),
                                        ha="center", va="bottom", fontsize=FS_STAR if line else FS_TICK, color=GREY,
                                        annotation_clip=False))


def keep_brackets(ax) -> None:
    """目盛をそろえるために縦軸を伸ばしたあと、ブラケットを描いたときと同じ見た目 (pt) で描き直す。
    ブラケットの間隔の決め方は変えず、伸ばした分は一番上のブラケットの上の余白になる (261020 決定・F122)。"""
    b = getattr(ax, "_brackets", None)
    if not b:
        return
    lo, hi = ax.get_ylim()
    per_pt = (hi - lo) / max(b["ax_pt"], 1e-9)
    _place_brackets(ax, b["step_pt"] * per_pt, b["h_pt"] * per_pt)


def draw_control_stars(ax, xpos: dict, stars: dict, tops: dict, lo: float, hi: float) -> float:
    """Dunnett / Welch (各群 vs 対照群)の星を、各群のデータの上に置く。上端が足りなければ伸ばす。"""
    if not any(stars.values()):
        return hi
    ax.figure.canvas.draw()
    ax_pt = ax.get_window_extent().height * 72.0 / ax.figure.dpi
    per_pt = (hi - lo) / max(ax_pt, 1e-9)
    need = max(tops[g] for g, s in stars.items() if s) + (FS_STAR * 1.4) * per_pt
    hi = max(hi, need)
    ax.set_ylim(lo, hi)
    for g, s in stars.items():
        if s:
            ax.annotate(s, xy=(xpos[g], tops[g]), textcoords="offset points", xytext=(0, 1),
                        ha="center", va="bottom", fontsize=FS_STAR, color=GREY, annotation_clip=False)
    return hi


def fixed_axes(ax_w: float, ax_h: float = PLOT_H_MM):
    """固定の余白 (左 15・右 4・上 13・下 20 mm)で軸を 1 つ置いた図を作る。"""
    fig_w = MARGIN_L + ax_w + MARGIN_R
    fig_h = MARGIN_B + ax_h + MARGIN_T
    fig = plt.figure(figsize=(fig_w * MM, fig_h * MM))
    ax = fig.add_axes([MARGIN_L / fig_w, MARGIN_B / fig_h, ax_w / fig_w, ax_h / fig_h])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return fig, ax


def rotated_label_mm(s: str, fs: float = XLABEL_FS, rot: float = XLABEL_ROT) -> float:
    """45° に傾けた群名が、軸から下に何 mm はみ出すかを実測する。"""
    fig = plt.figure()
    t = fig.text(0.5, 0.5, s, fontsize=fs, rotation=rot, ha="right", va="top", rotation_mode="anchor")
    fig.canvas.draw()
    h = t.get_window_extent().height / fig.dpi * 25.4
    plt.close(fig)
    return h


def labels_overflow(labels: list[str]) -> list[tuple[str, float]]:
    """下余白 (20 mm)に収まらない可能性のある群名を返す (目盛の長さ＋目盛と文字の隙間を含めて判定)。"""
    gap_mm = (plt.rcParams["xtick.major.size"] + plt.rcParams["xtick.major.pad"]) / 72 * 25.4
    return [(s, h) for s in labels if (h := rotated_label_mm(s) + gap_mm) > MARGIN_B]


def box_panel(labels: list[str], ylabel: str):
    """群ごとの箱ひげ＋点の枠を作る。横幅 = 1 群 8 mm × 群の数 (260921 決定)。
    返り値: (fig, ax, xpos)"""
    fig, ax = fixed_axes(MM_PER_GROUP * len(labels))
    ax.set_ylabel(ylabel)
    ax._group_labels = labels            # 目盛ラベルは dots_mean_sd が点を描くときに貼る
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return fig, ax, {i: i for i in range(len(labels))}


DOT_S = 10                     # 群の図の点の大きさ（面積）。14 と 10 を実物で比べて 10 に決定（260921）
DOT_JITTER = 0.2               # 点の横の散らし幅（±。1 群の幅 = 1）。乱数の種は JITTER_SEED で固定
MEAN_HALF = 0.25               # 平均の横線の半分の長さ（1 群の幅の約半分 = ±0.25）
MEAN_LW = 0.8                  # 平均の横線の太さ (pt)。エラーバー・軸と同じ 0.8 pt（260921 決定。1.2 と比較）


def dots_mean_sd(ax, data: list, colors: list[str], hollow: list | None = None) -> list[tuple[float, float]]:
    """群ごとの点 (群色・ランダムに横へ散らす・種は固定)＋平均の横線＋SD のエラーバー (前景色・キャップあり)。
    16S パイプラインの箱ひげ図をこの形に置き換えた (260921 指示。n = 3 の実験が多いため)。
    hollow[i] が True の点は白抜き (中を白・外枠は群色)。返り値: 群ごとの (平均, SD)。"""
    rng = np.random.default_rng(JITTER_SEED)
    labels = getattr(ax, "_group_labels", None)
    if labels is not None:
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=XLABEL_FS, rotation=XLABEL_ROT, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=plt.rcParams["ytick.major.size"], width=plt.rcParams["axes.linewidth"])
    ax.set_xlim(-0.5, len(data) - 0.5)
    stats = []
    for i, (v, c) in enumerate(zip(data, colors)):
        v = np.asarray(v, float)
        x = i + rng.uniform(-DOT_JITTER, DOT_JITTER, len(v))
        ho = np.asarray(hollow[i], bool) if hollow is not None else np.zeros(len(v), bool)
        ax.scatter(x[~ho], v[~ho], s=DOT_S, color=c, lw=0.3, edgecolors="white", zorder=2)
        if ho.any():
            ax.scatter(x[ho], v[ho], s=DOT_S, facecolors=OPEN_FACE, edgecolors=c, linewidths=OPEN_LW, zorder=2)
        if len(v) == 0:
            stats.append((np.nan, np.nan))
            continue
        m = float(v.mean())
        sd = float(v.std(ddof=1)) if len(v) > 1 else 0.0
        ax.plot([i - MEAN_HALF, i + MEAN_HALF], [m, m], color=GREY, lw=MEAN_LW, zorder=4, solid_capstyle="butt")
        if len(v) > 1:
            ax.errorbar(i, m, yerr=sd, fmt="none", ecolor=GREY, elinewidth=0.8, capsize=2.5, capthick=0.8, zorder=4)
        stats.append((m, sd))
    return stats


# ---------------- サブグループの記号（メタボローム lib/plots.py と同じ規則） ----------------
SUBGROUP_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">", "h"]
OPEN_LW = 0.6                  # 白抜きの記号の外枠の太さ (pt)。260921 決定（0.9/0.6/0.4/0.3 と白塗りの有無を実物で比較）
OPEN_FACE = "white"            # 白抜きの中身を白で塗る（重なっても下の線が透けず、白抜きと分かる。260921 決定）


def subgroup_markers(subgroups: dict) -> dict:
    """サブグループ → 記号。subgroups.<名前>.marker で指定し、auto のものは定義順に自動で割り当てる (重複は避ける)。"""
    override = {k: v["marker"] for k, v in subgroups.items() if str(v.get("marker", "auto")).strip().lower() != "auto"}
    used = set(override.values())
    pool = [m for m in SUBGROUP_MARKERS if m not in used] or SUBGROUP_MARKERS
    out, i = {}, 0
    for x in subgroups:
        if x in override:
            out[x] = override[x]
        else:
            out[x] = pool[i % len(pool)]
            i += 1
    return out


def scatter_subgroup(ax, x, y, color, spec: dict, marker: str, s: float = 14, **kw):
    """サブグループの記号で散布する。open: true なら白抜き (輪郭だけ群の色)。scale は面積の倍率。"""
    s = float(s) * float(spec.get("scale", 1.0))
    if spec.get("open", False):
        return ax.scatter(x, y, marker=marker, s=s, facecolors=OPEN_FACE, edgecolors=color, linewidths=OPEN_LW, **kw)
    return ax.scatter(x, y, marker=marker, s=s, color=color, lw=0.3, edgecolors="white", **kw)


def save_legend(path, G, markers: dict | None = None, present_subs: list | None = None, log=None) -> None:
    """凡例を別ファイルに出す (260921 決定。規約「同一凡例は 1 つに集約」)。
    群 (色)だけの版と、群＋サブグループ (記号)の版を作れる。"""
    from matplotlib.lines import Line2D
    if G.color_by_subgroup:                         # サブグループごとに色を付けるとき（use_subgroup: true・show_group_subgroup: false）
        h_g = [Line2D([], [], ls="", marker="o", ms=5, color=G.color[x], label=x) for x in G.order]
        rows = [("Subgroup", h_g)]
    else:
        h_g = [Line2D([], [], ls="", marker="o", ms=5, color=G.group_color[x], label=x) for x in G.group_order]
        rows = [("Group", h_g)]
    if markers:
        h_s = []
        for x in present_subs:
            spec = G.subgroups.get(x, {})
            style = ({"markerfacecolor": OPEN_FACE, "markeredgecolor": "#666666", "markeredgewidth": OPEN_LW}
                     if spec.get("open") else {})
            h_s.append(Line2D([], [], ls="", marker=markers[x], color="#666666", label=x,
                              ms=5 * float(spec.get("scale", 1.0)) ** 0.5, **style))
        rows.append(("Subgroup" + (f"\n({G.open_note})" if G.open_note else ""), h_s))
    # 高さは他の図と同じ（下 20 + 軸 48 + 上 13 = 81 mm。Keynote でグリッドにそろえるため。260921 指示）。
    # 横幅は中身を実測して決める。上端は軸の上端（上余白 5 mm）にそろえ、Group と Subgroup のあいだは 1 行ぶん空ける。
    fig_h = MARGIN_B + PLOT_H_MM + MARGIN_T
    line_mm = FS_LEGEND / 72 * 25.4 * 1.2

    def _place(fig_w):
        """Group と Subgroup を上から順に置き、(右端 mm, 下端 mm) を返す。"""
        y_mm, right_mm = MARGIN_T, 0.0
        for lg in legends:
            lg.set_bbox_to_anchor((2 / fig_w, 1 - y_mm / fig_h), transform=fig.transFigure)
            fig.canvas.draw()
            bb = lg.get_window_extent()
            y_mm += bb.height / fig.dpi * 25.4 + line_mm                  # 次の凡例は 1 行ぶん空けて下に置く
            right_mm = max(right_mm, 2 + bb.width / fig.dpi * 25.4)
        return right_mm, y_mm - line_mm

    # 高さに収まらないときは列を増やす（高さ 81 mm は変えず、横に広げる。260925 決定・案 B）
    for ncol in range(1, max(len(h) for _, h in rows) + 1):   # 列数に上限は設けない（群が多いほど横に伸びる。260925 指示）
        fig = plt.figure(figsize=(400 * MM, fig_h * MM))
        legends = []
        for title, handles in rows:
            lg = fig.legend(handles=handles, fontsize=FS_LEGEND, frameon=False, loc="upper left", borderaxespad=0,
                            borderpad=0, handletextpad=0.4, labelspacing=0.4, title=title, title_fontsize=FS_LEGEND,
                            alignment="left", ncol=ncol, columnspacing=1.0)
            fig.add_artist(lg)
            legends.append(lg)
        right_mm, bottom_mm = _place(400)
        if bottom_mm <= fig_h - 2:
            break
        plt.close(fig)
    fig_w = right_mm + 2                                                  # 横幅は中身を実測（右に 2 mm）
    fig.set_size_inches(fig_w * MM, fig_h * MM)
    _place(fig_w)
    if bottom_mm > fig_h - 2:
        msg = (f"  [注意] 凡例が高さ {fig_h:.0f} mm に収まっていない可能性があります: {Path(path).name}"
               f" (必要な高さ 約 {bottom_mm + 2:.0f} mm)")
        (log or print)(msg)                       # ログと注記の一覧に残す (261021 決定・F3。以前は print だけ。見切れの恐れなので [注意])
    fig.savefig(path)
    plt.close(fig)


# ---------------- 散布図（PCA・PCoA）の規約（260921 確定・Skills/figure-plotting §3） ----------------
SCATTER_AX_MM = PLOT_H_MM      # 軸は 48 × 48 mm（軸高と同じ。交差形式の 1 マスも同じ）
SCATTER_S = 14                 # 点の大きさ（面積）
ORIGIN_LINE = dict(color="#DDDDDD", lw=0.5, zorder=0)
ELLIPSE_LEVEL = 0.95           # 群ごとの 95% データ楕円（χ² 近似・半径一定。261009 D8 決定 b）
TITLE_FS = 9                   # 上余白に書く「何の指標か」のタイトル（260921 指定: 9 pt）


def data_ellipse(x, y, level: float = ELLIPSE_LEVEL):
    """群の 95% データ楕円 (χ² 近似。二変量正規ならその群の検体の level が入る範囲)。n < 3 なら None。
    半径 r = sqrt(χ²(level; 2)) ≈ 2.45 SD で n によらず一定 (vegan ordiellipse(kind="sd", conf=0.95) と同じ)。
    261009 D8 で Hotelling T² (n が小さいと SD の数倍〜28 倍に膨らむ)から変更。
    返り値: dict(cx, cy, width, height, angle, n)"""
    from scipy import stats as _st
    pts = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    n = len(pts)
    if n < 3:
        return None
    vals, vecs = np.linalg.eigh(np.cov(pts.T))
    o = np.argsort(-vals)
    vals, vecs = vals[o], vecs[:, o]
    if np.any(vals <= 0):
        return None
    r = np.sqrt(_st.chi2.ppf(level, 2))
    return dict(cx=pts[:, 0].mean(), cy=pts[:, 1].mean(), width=2 * r * np.sqrt(vals[0]),
                height=2 * r * np.sqrt(vals[1]), angle=float(np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))), n=n)


def draw_ellipse(ax, e: dict, color: str) -> None:
    from matplotlib.patches import Ellipse
    ax.add_patch(Ellipse((e["cx"], e["cy"]), e["width"], e["height"], angle=e["angle"], facecolor=color,
                         alpha=0.10, edgecolor=color, lw=0.6, zorder=1))
    ax.add_patch(Ellipse((e["cx"], e["cy"]), e["width"], e["height"], angle=e["angle"], facecolor="none",
                         edgecolor=color, lw=0.6, zorder=1))


def symmetric_ticks(lim: float, inner: float = 0.8) -> list[float]:
    """原点中心の軸に、±lim の inner 以内の丸い値で目盛を振る (隣のマスの数字と重ならないように)。"""
    step = nice_step(2 * lim)
    k = int(np.floor(lim * inner / step))
    return [i * step for i in range(-k, k + 1)]


def wrap_title(name: str, width_mm: float, fs: float = TITLE_FS) -> str:
    """タイトルを軸の幅に収まるよう折り返す (260923 決定)。
    切れ目は `;` の直後と、1 つだけの `_` の直後 (`__` = 空の階級の印は割らない。261021・F4)。収まらない長い語は文字の途中でも切る。"""
    import re
    lines, cur = [], ""
    for part in re.split(r"(?<=;)|(?<=[^_]_)(?=[^_])", name):
        while text_mm(part, fs) > width_mm:            # 切れ目の無い長い語は文字の途中で切る
            k = len(part)
            while k > 1 and text_mm(part[:k], fs) > width_mm:
                k -= 1
            if cur:
                lines.append(cur)
                cur = ""
            lines.append(part[:k])
            part = part[k:]
        if cur and text_mm(cur + part, fs) > width_mm:
            lines.append(cur)
            cur = part
        else:
            cur += part
    if cur:
        lines.append(cur)
    return "\n".join(lines)


# 値の軸は目盛 5 個（4 等分）にそろえ、刻みは下の候補から「収まる一番小さいもの」を選ぶ（260924 決定・全図共通）。
# 候補を 1/2/2.5/5 系に限ると上端の無駄が平均 39%（16S #10 の上位 60 タクソンで実測）。広げると 18% に下がる
TICK_N = 4
TICK_CAND = (1, 1.5, 2, 2.5, 3, 4, 5, 7.5)


def tick_step(hi: float, lo: float = 0.0, n: int = TICK_N) -> float:
    """lo から n 等分して hi 以上になる、一番小さい丸い刻みを返す。"""
    span = hi - lo
    best = None
    for p in range(-8, 9):
        for c in TICK_CAND:
            st = c * 10.0 ** p
            if st * n >= span and (best is None or st < best):
                best = st
    return best


def fixed_ticks(ax, lo: float, hi: float, axis: str = "y") -> tuple[float, float]:
    """値の軸に目盛を 5 個置く (上端は刻み × 4)。返り値: (下端, 新しい上端)。"""
    st = tick_step(hi, lo)
    top = lo + st * TICK_N
    ticks = [lo + st * k for k in range(TICK_N + 1)]
    if axis == "y":
        ax.set_ylim(lo, top)
        ax.set_yticks(ticks)
    else:
        ax.set_xlim(lo, top)
        ax.set_xticks(ticks)
    return lo, top


VIRIDIS_TOP = 0.95     # 連続の色は viridis の 0〜95% を使う（一番明るい黄色は使わない。全図共通・261009 決定）


def viridis():
    """全図共通の連続の色: viridis の 0〜VIRIDIS_TOP。マスの中の文字は白 (261009 決定)。"""
    import matplotlib.colors as mc
    import matplotlib.pyplot as plt
    return mc.LinearSegmentedColormap.from_list("viridis_095", plt.get_cmap("viridis")(np.linspace(0, VIRIDIS_TOP, 256)))
