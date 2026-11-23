import os

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURES = os.path.join(ROOT, "figures")
BUILDS = ("canonical", "mixed")
LABEL = {"canonical": "Canonical homodimer", "mixed": "Mixed dimer"}
COL = {"canonical": "#4d4d4d", "mixed": "#b2182b"}
BURN_COL = "#efefef"


def stop(msg):
    raise SystemExit("STOP: " + msg)


def picture(ax, path, note):
    import matplotlib.image
    import numpy as np
    from render_md_struct import strip_credentials
    ax.axis("off")
    if os.path.isfile(path):
        gone = strip_credentials(path)
        if gone:
            print("stripped %s from %s" % (", ".join(gone), path))
    else:
        ax.text(0.5, 0.5, note, ha="center", va="center", fontsize=6.4,
                color="#8a8a8a", wrap=True, transform=ax.transAxes)
        return False
    try:
        from PIL import Image
    except ImportError:
        ax.imshow(matplotlib.image.imread(path))
        return True
    im = Image.open(path).convert("RGBA")
    box = im.getbbox()
    if box:
        im = im.crop(box)
    ax.imshow(np.asarray(im))
    return True


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42, "figure.dpi": 100,
    "svg.hashsalt": "kynu-truncated-isoform",
})

_LETTERS = []


def letter(ax, s):
    _LETTERS.append((ax, s))


def _groups(vals, tol=0.03):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out, cur = [], [order[0]]
    for i in order[1:]:
        if vals[i] - vals[cur[-1]] <= tol:
            cur.append(i)
        else:
            out.append(cur)
            cur = [i]
    out.append(cur)
    return out


def _place_letters(fig):
    items = [(ax, s) for ax, s in _LETTERS if ax.figure is fig]
    if not items:
        return
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    bbs = [ax.get_tightbbox(rend, bbox_extra_artists=[]).transformed(inv)
           for ax, _ in items]
    pos = [ax.get_position() for ax, _ in items]
    ys, xs = [b.y1 for b in bbs], [b.x0 for b in bbs]
    for grp in _groups([p.y1 for p in pos]):
        top = max(bbs[i].y1 for i in grp)
        for i in grp:
            ys[i] = top
    for grp in _groups([b.x0 for b in bbs], tol=0.02):
        left = min(bbs[i].x0 for i in grp)
        for i in grp:
            xs[i] = left
    for (ax, s), x, y in zip(items, xs, ys):
        fig.text(x - 0.012, y + 0.004, s, fontsize=9, fontweight="bold",
                 va="top", ha="left")
    _LETTERS[:] = [it for it in _LETTERS if it[0].figure is not fig]


def save(fig, name):
    _place_letters(fig)
    os.makedirs(FIGURES, exist_ok=True)
    fig.savefig(os.path.join(FIGURES, name + ".png"), dpi=400, bbox_inches="tight",
                metadata={"Software": None})
    fig.savefig(os.path.join(FIGURES, name + ".svg"), bbox_inches="tight",
                metadata={"Date": None, "Creator": None})
    plt.close(fig)
    print("wrote", os.path.join(FIGURES, name + ".png"), "and .svg")
