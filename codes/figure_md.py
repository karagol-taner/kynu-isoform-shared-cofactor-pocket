import json
import os
import sys

import matplotlib
import matplotlib.image
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from render_md_struct import strip_credentials

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
MD = os.path.join(RESULTS, "md")
SERIES = os.path.join(MD, "series")
from figure_style import (FIGURES, BUILDS, LABEL, COL, BURN_COL, letter,
                          save, stop, picture as _picture)

def xvg(build, name, column=1, scale=1.0):
    p = os.path.join(SERIES, build, name)
    if not os.path.isfile(p):
        stop("no %s. The series each run wrote live under %s, one folder per build."
             % (p, SERIES))
    D = np.loadtxt(p, comments=["#", "@"])
    if D.ndim == 1:
        D = D.reshape(1, -1)
    if column >= D.shape[1]:
        stop("%s has %d columns and column %d was asked for." % (p, D.shape[1], column))
    return D[:, 0] / 1000.0, D[:, column] * scale


def picture(ax, name, note):
    return _picture(ax, os.path.join(MD, name), note)


REP = {b: json.load(open(os.path.join(MD, b + "_report.json"))) for b in BUILDS}
DRIFT = json.load(open(os.path.join(MD, "cofactor_drift.json")))
CRIT = json.load(open(os.path.join(MD, "kill_criterion.json")))

for k in ("ns", "frames", "burn"):
    if REP["canonical"][k] != REP["mixed"][k]:
        stop("the two runs disagree on %s, %s, and the panels put them on one axis."
             % (k, {b: REP[b][k] for b in BUILDS}))
NS = REP["canonical"]["ns"]
FRAMES = REP["canonical"]["frames"]
BURN = REP["canonical"]["burn"]
BURN_NS = BURN * NS / (FRAMES - 1)
HOLD = CRIT["thresholds"]["hold"]

S = {}
for b in BUILDS:
    t, ab = xvg(b, "sasa_dimer.xvg")
    _, s1 = xvg(b, "sasa_sub1.xvg")
    _, s2 = xvg(b, "sasa_sub2.xvg")
    n = min(len(ab), len(s1), len(s2))
    S[b] = {
        "t": t[:n],
        "buried": (s1[:n] + s2[:n] - ab[:n]) * 100.0 / 2.0,
        "cofactor": xvg(b, "sasa_cofactor.xvg", column=2, scale=100.0)[1],
        "contacts": xvg(b, "numcont_dimer.xvg")[1],
        "mindist": xvg(b, "mindist_dimer.xvg", scale=10.0)[1],
        "drift": xvg(b, "cofactor_drift.xvg", scale=10.0)[1],
    }

for b in BUILDS:
    for key, want in (("buried", REP[b]["buried_per_subunit_A2"]["mean"]),
                      ("cofactor", REP[b]["cofactor_sasa_A2"]["mean"]),
                      ("contacts", REP[b]["contacts_0_5nm"]["mean"]),
                      ("mindist", REP[b]["mindist_subunits_A"]["mean"]),
                      ("drift", DRIFT[b]["mean"])):
        v = np.asarray(S[b][key], dtype=float)
        if len(v) != FRAMES:
            stop("the %s %s series has %d frames and the report counts %d."
                 % (b, key, len(v), FRAMES))
        got = float(v[BURN:].mean())
        if abs(got - float(want)) > max(1e-6, abs(float(want)) * 1e-6):
            stop("the %s %s series averages %.6f over the scored window and the report "
                 "says %.6f, so they are not the same run." % (b, key, got, want))
print("every series reproduces the statistic its report carries")

fig = plt.figure(figsize=(7.2, 4.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1.85, 1.0, 1.0], hspace=0.20, wspace=0.52,
                      left=0.055, right=0.985, top=0.93, bottom=0.10)

ax_a = fig.add_subplot(gs[0, 0])
picture(ax_a, "render_systems.png",
        "structural panel needs PyMOL\nrun codes/render_md_struct.py")

ax_b = fig.add_subplot(gs[0, 1])
ax_b.axvspan(0, BURN_NS, color=BURN_COL, lw=0, zorder=0)
for b in BUILDS:
    ax_b.plot(S[b]["t"], S[b]["buried"], lw=0.45, color=COL[b], alpha=0.9, zorder=2)
thr = float(np.asarray(S["canonical"]["buried"])[BURN:].mean()) * HOLD
ax_b.axhline(thr, color="black", lw=0.7, ls="--", zorder=3)
ax_b.annotate("%.2f of canonical" % HOLD, xy=(0.98, thr),
              xycoords=("axes fraction", "data"), xytext=(0, 3),
              textcoords="offset points", fontsize=6, ha="right", va="bottom")
ax_b.set_ylabel("Buried area per\nsubunit, $\\mathrm{\\AA}^2$")

ax_c = fig.add_subplot(gs[0, 2])
ax_c.axvspan(0, BURN_NS, color=BURN_COL, lw=0, zorder=0)
for b in BUILDS:
    ax_c.plot(S[b]["t"], S[b]["contacts"], lw=0.45, color=COL[b], alpha=0.9, zorder=2)
ax_c.set_ylabel("Inter-subunit contacts")

ax_d = fig.add_subplot(gs[1, 0])
picture(ax_d, "render_cofactor.png",
        "cofactor panel needs the end frames\nrun kynu_MD_3c_frames.ipynb, then\n"
        "codes/render_md_struct.py")

ax_e = fig.add_subplot(gs[1, 1])
ax_e.axvspan(0, BURN_NS, color=BURN_COL, lw=0, zorder=0)
for b in BUILDS:
    ax_e.plot(S[b]["t"], S[b]["drift"], lw=0.45, color=COL[b], alpha=0.9, zorder=2)
ax_e.set_ylabel("Cofactor displacement, $\\mathrm{\\AA}$")

for ax in (ax_b, ax_c, ax_e):
    ax.set_xlim(0, NS)
    ax.set_xlabel("Time, ns")
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))

ax_f = fig.add_subplot(gs[1, 2])
allv = np.concatenate([np.asarray(S[b]["cofactor"])[BURN:] for b in BUILDS])
grid = np.linspace(max(0.0, allv.min() - 10), allv.max() + 10, 400)
for b in BUILDS:
    v = np.asarray(S[b]["cofactor"])[BURN:]
    k = stats.gaussian_kde(v)
    ax_f.plot(grid, k(grid), lw=1.0, color=COL[b])
    ax_f.fill_between(grid, k(grid), color=COL[b], alpha=0.10, lw=0)
    ax_f.axvline(v.mean(), color=COL[b], lw=0.7, ls=":")
ax_f.set_xlabel("Cofactor exposure over\nthe scored window, $\\mathrm{\\AA}^2$")
ax_f.set_ylabel("Density")
ax_f.set_xlim(grid[0], grid[-1])

for ax in (ax_b, ax_c, ax_e, ax_f):
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4))
    box = ax.get_position()
    short = box.height * 0.82
    ax.set_position([box.x0, box.y0 + (box.height - short) * 0.5,
                     box.width, short])

for ax, s in zip((ax_a, ax_b, ax_c, ax_d, ax_e, ax_f), "abcdef"):
    letter(ax, s)

handles = [plt.Line2D([], [], color=COL[b], lw=1.4, label=LABEL[b]) for b in BUILDS]
handles.append(matplotlib.patches.Patch(facecolor=BURN_COL, edgecolor="none",
                                        label="Discarded before scoring"))
fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, 1.02))
save(fig, "Figure5")
