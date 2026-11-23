import json
import os
import sys

import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import COL, letter, save, stop

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DCS = os.path.join(RESULTS, "dcs")
SUMMARY = os.path.join(RESULTS, "dcs_summary.json")
COPY_COL = "#a8a8a8"

if not os.path.isfile(SUMMARY):
    stop("no %s. Run codes/dcs_analysis.py first." % SUMMARY)
S = json.load(open(SUMMARY))
SPAN_LO, SPAN_HI = S["method"]["scored_span"]


def table(build, chain):
    p = os.path.join(DCS, "%s_%s.tsv" % (build, chain))
    if not os.path.isfile(p):
        stop("no %s. Run codes/dcs_analysis.py first." % p)
    rows = []
    with open(p) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for L in fh:
            r = dict(zip(head, L.rstrip("\n").split("\t")))
            rows.append({"position": int(r["position"]), "residue": r["residue"],
                         "block": r["block"], "rmsf": float(r["rmsf_A"]),
                         "nd": float(r["norm_dynamics"]),
                         "cons": float(r["conservation"]) if r["conservation"] else None,
                         "dcs": float(r["dcs_0_10"]) if r["dcs_0_10"] else None,
                         "rcs": float(r["rcs_0_10"]) if r["rcs_0_10"] else None})
    return rows


def scored(rows):
    return [r for r in rows if r["dcs"] is not None]


def series(rows, lo=None, hi=None):
    sel = [r for r in rows if (lo is None or r["position"] >= lo)
           and (hi is None or r["position"] <= hi)]
    return (np.array([r["position"] for r in sel], float),
            np.array([r["rmsf"] for r in sel], float))


T = {(b, c): table(b, c) for b in ("canonical", "mixed") for c in ("sub1", "sub2")}

iso = S["isoform_against_full_partner"]
lo, hi = iso["first_position"], iso["last_position"]
a_iso = {r["position"]: r["rmsf"] for r in scored(T[("mixed", "sub2")])}
a_par = {r["position"]: r["rmsf"] for r in scored(T[("canonical", "sub2")])}
shared = [p for p in sorted(set(a_iso) & set(a_par)) if lo <= p <= hi]
checks = [("the isoform retained block", np.mean([a_iso[p] for p in shared]),
           iso["isoform_rmsf_A_mean"]),
          ("a full length partner", np.mean([a_par[p] for p in shared]),
           iso["partner_rmsf_A_mean"]),
          ("positions compared", float(len(shared)), float(iso["positions"]))]
for name, got, want in checks:
    if abs(got - want) > 5e-3:
        stop("%s comes to %.4f from the score tables and the summary reports %.4f, "
             "so the two disagree." % (name, got, want))
print("every series reproduces the statistic the score summary carries")

fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1.0], hspace=0.42, wspace=0.28,
                      left=0.075, right=0.985, top=0.90, bottom=0.095)
ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
ax_c, ax_d = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

x, y = series(scored(T[("mixed", "sub2")]), lo, hi)
ax_a.plot(x, y, lw=0.8, color=COL["mixed"], label="Isoform chain")
x2, y2 = series(scored(T[("canonical", "sub2")]), lo, hi)
ax_a.plot(x2, y2, lw=0.8, color=COL["canonical"], label="Full length partner")
ax_a.set_yscale("log")
ax_a.set_xlabel("Residue, retained block")
ax_a.set_ylabel("C$_\\alpha$ fluctuation, $\\mathrm{\\AA}$")
ax_a.set_xlim(lo, hi)
ax_a.set_ylim(min(y.min(), y2.min()) * 0.85, max(y.max(), y2.max()) * 1.30)
ax_a.legend(frameon=False, loc="upper left", handlelength=1.4)
ax_a.annotate("%.2f$\\times$ at %d of %d positions"
              % (iso["ratio"], iso["positions_isoform_higher"], iso["positions"]),
              xy=(0.985, 0.97), xycoords="axes fraction", ha="right", va="top",
              fontsize=6.4, color=COL["mixed"])

for key, col, lab, z in ((("canonical", "sub2"), COPY_COL, "Second copy, canonical", 1),
                         (("canonical", "sub1"), COL["canonical"], "Canonical dimer", 2),
                         (("mixed", "sub1"), COL["mixed"], "Mixed dimer", 3)):
    x, y = series(scored(T[key]), SPAN_LO, SPAN_HI)
    ax_b.plot(x, y, lw=0.6, color=col, label=lab, zorder=z)
ax_b.set_yscale("log")
ax_b.set_xlabel("Residue, canonical subunit")
ax_b.set_ylabel("C$_\\alpha$ fluctuation, $\\mathrm{\\AA}$")
ax_b.set_xlim(SPAN_LO, SPAN_HI)
_bv = np.concatenate([series(scored(T[k]), SPAN_LO, SPAN_HI)[1]
                      for k in (("canonical", "sub1"), ("canonical", "sub2"),
                                ("mixed", "sub1"))])
ax_b.set_ylim(_bv.min() * 0.85, _bv.max() * 1.55)
ax_b.legend(frameon=False, loc="upper left", handlelength=1.4)

rows = scored(T[("canonical", "sub1")])
nd = np.array([r["rmsf"] for r in rows])
cons = np.array([r["cons"] for r in rows])
ax_c.scatter(nd, cons, s=3.5, color="#cfcfcf", lw=0, zorder=1)
ax_c.set_xscale("log")
top_d = {r["position"] for r in S["chains"]["canonical_sub1"]["dynamic_conserved"]}
top_r = {r["position"] for r in S["chains"]["canonical_sub1"]["rigid_conserved"]}
for sel, col, lab in ((top_r, COL["canonical"], "Rigid conserved"),
                      (top_d, COL["mixed"], "Dynamic conserved")):
    pts = [(r["rmsf"], r["cons"]) for r in rows if r["position"] in sel]
    ax_c.scatter([p[0] for p in pts], [p[1] for p in pts], s=6, color=col, lw=0,
                 label="%s, top %d" % (lab, len(pts)), zorder=2)
lys = [r for r in rows if r["position"] == int(S["method"]["schiff_base_residue"])]
if lys:
    r = lys[0]
    ax_c.annotate("%s%d, the catalytic lysine" % (r["residue"], r["position"]),
                  xy=(r["rmsf"], r["cons"]), xytext=(18, -26),
                  textcoords="offset points", fontsize=6.2,
                  color=COL["canonical"], ha="left",
                  bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none",
                            alpha=0.85),
                  arrowprops=dict(arrowstyle="-", lw=0.5, color=COL["canonical"]))
ax_c.set_xlabel("C$_\\alpha$ fluctuation, $\\mathrm{\\AA}$")
ax_c.set_ylabel("Orthologue identity")
ax_c.set_ylim(-0.03, 1.06)
ax_c.legend(frameon=False, loc="lower right", handlelength=1.0, borderpad=0.2)

c1 = {r["position"]: r["rmsf"] for r in scored(T[("canonical", "sub1")])}
c2 = {r["position"]: r["rmsf"] for r in scored(T[("canonical", "sub2")])}
m1 = {r["position"]: r["rmsf"] for r in scored(T[("mixed", "sub1")])}
both = sorted(set(c1) & set(m1) & set(c2))
d_build = np.array([m1[p] - c1[p] for p in both])
d_copy = np.array([c2[p] - c1[p] for p in both])
edges = np.linspace(-0.6, 0.6, 41)
for vals, col, lab in ((d_copy, COPY_COL, "Second copy, no truncation"),
                       (d_build, COL["mixed"], "Mixed against canonical")):
    ax_d.hist(vals, bins=edges, histtype="step", lw=1.0, color=col, label=lab)
ax_d.axvline(0.0, color="#000000", lw=0.5, ls=":")
ax_d.set_xlabel("Change in fluctuation per residue, $\\mathrm{\\AA}$")
ax_d.set_ylabel("Residues")
ax_d.set_xlim(edges[0], edges[-1])
ax_d.legend(frameon=False, loc="upper left", handlelength=1.4, borderpad=0.2)

for ax in (ax_c, ax_d):
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
for ax in (ax_a, ax_b):
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.5, 1, 2, 5]))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
ax_c.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.5, 1, 2, 4]))
ax_c.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax_c.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
for ax, s in zip((ax_a, ax_b, ax_c, ax_d), "abcd"):
    letter(ax, s)
save(fig, "Figure6")
