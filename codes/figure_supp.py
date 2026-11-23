import csv
import glob
import json
import os
import sys

import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_common as cc
from figure_style import BUILDS, BURN_COL, COL, LABEL, letter, save, stop

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
MD = os.path.join(RESULTS, "md")
SERIES = os.path.join(MD, "series")
AF3 = json.load(open(os.path.join(RESULTS, "af3_scores.json")))
DCS = json.load(open(os.path.join(RESULTS, "dcs_summary.json")))
REP = {b: json.load(open(os.path.join(MD, b + "_report.json"))) for b in BUILDS}
BLOCK_COL = {cc.BLOCKS[0][0]: "#bfbfbf", cc.BLOCKS[1][0]: "#b2182b",
             cc.BLOCKS[2][0]: "#4d4d4d"}

JOB_LABEL = {"canonical__canonical": "canonical", "canonical__isoform2": "isoform",
             "canonical__trunc_cterm": "tail removed",
             "canonical__trunc_nterm": "N face removed",
             "canonical__tail_scramble_01": "permuted 1",
             "canonical__tail_scramble_02": "permuted 2",
             "canonical__tail_scramble_03": "permuted 3",
             "isoform2__isoform2": "isoform pair"}
JOB_ORDER = ("canonical__canonical", "canonical__isoform2", "canonical__trunc_cterm",
             "canonical__trunc_nterm", "canonical__tail_scramble_01",
             "canonical__tail_scramble_02", "canonical__tail_scramble_03",
             "isoform2__isoform2")
if set(AF3["detail"]) != set(JOB_ORDER):
    stop("af3_scores.json carries the jobs %s and this figure is laid out for %s."
         % (sorted(AF3["detail"]), sorted(JOB_ORDER)))
for k in ("ns", "frames", "burn"):
    if REP["canonical"][k] != REP["mixed"][k]:
        stop("the two runs disagree on %s, %s, and the panels put them on one axis."
             % (k, {b: REP[b][k] for b in BUILDS}))
NS, FRAMES, BURN = (REP["canonical"][k] for k in ("ns", "frames", "burn"))
BURN_NS = NS * BURN / float(FRAMES - 1)


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


fig = plt.figure(figsize=(7.2, 4.4))
gs = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.30, left=0.065, right=0.90,
                      top=0.92, bottom=0.10)
img = None
for i, job in enumerate(JOB_ORDER):
    ax = fig.add_subplot(gs[i // 4, i % 4])
    row = AF3["detail"][job]
    name = row["model_file"].replace("_model_", "_full_data_")[:-4] + ".json"
    hits = glob.glob(os.path.join(ROOT, "data", "af3", "**", name), recursive=True)
    if len(hits) != 1:
        stop("%d copies of %s under data/af3, and one is needed." % (len(hits), name))
    pae = np.asarray(json.load(open(hits[0]))["pae"], dtype=float)
    cut = row["interface"]["chain_residues"][0]
    if pae.shape[0] != sum(row["interface"]["chain_residues"]):
        stop("the error matrix of %s is %d wide and its chains hold %d residues."
             % (job, pae.shape[0], sum(row["interface"]["chain_residues"])))
    img = ax.imshow(pae, cmap="Greys", vmin=0, vmax=31.75, interpolation="nearest")
    for at in (cut,):
        ax.axhline(at - 0.5, color=COL["mixed"], lw=0.5)
        ax.axvline(at - 0.5, color=COL["mixed"], lw=0.5)
    ax.set_title("%s\nipTM %.2f" % (JOB_LABEL[job], float(row["summary"]["iptm"])),
                 fontsize=6.2, pad=2)
    ax.set_xticks([0, cut, pae.shape[0]])
    ax.set_yticks([0, cut, pae.shape[0]])
    ax.tick_params(labelsize=5.4)
    if i % 4 == 0:
        ax.set_ylabel("Residue")
    if i // 4 == 1:
        ax.set_xlabel("Residue")
    letter(ax, "abcdefgh"[i])
cb = fig.colorbar(img, cax=fig.add_axes([0.915, 0.12, 0.014, 0.74]))
cb.set_label("Predicted aligned error, $\\mathrm{\\AA}$", fontsize=6.4)
cb.ax.tick_params(labelsize=5.8)
save(fig, "FigureS1")


fig = plt.figure(figsize=(7.2, 4.2))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.28, left=0.085, right=0.985,
                      top=0.95, bottom=0.11)
PANELS = (("rmsd.xvg", 1, 10.0, "Backbone deviation, $\\mathrm{\\AA}$"),
          ("mindist_dimer.xvg", 1, 10.0, "Subunit closest approach, $\\mathrm{\\AA}$"),
          ("sasa_dimer.xvg", 1, 1.0, "Dimer surface, nm$^2$"),
          ("mindist_cofactor_partner.xvg", 1, 10.0, "Cofactor to partner, "
                                                    "$\\mathrm{\\AA}$"))
for i, (name, col, scale, ylab) in enumerate(PANELS):
    ax = fig.add_subplot(gs[i // 2, i % 2])
    ax.axvspan(0, BURN_NS, color=BURN_COL, lw=0)
    for build in BUILDS:
        t, v = xvg(build, name, col, scale)
        if abs(t[-1] - NS) > 0.5:
            stop("%s of the %s run ends at %.1f ns and the report says %.1f."
                 % (name, build, t[-1], NS))
        ax.plot(t, v, lw=0.5, color=COL[build], label=LABEL[build])
    ax.set_xlabel("Time, ns")
    ax.set_ylabel(ylab)
    ax.set_xlim(0, NS)
    if i == 0:
        ax.legend(frameon=False, fontsize=6.0, handlelength=1.2, loc="lower right")
    letter(ax, "abcd"[i])
save(fig, "FigureS2")


def chain_table(key):
    p = os.path.join(RESULTS, "dcs", key + ".tsv")
    if not os.path.isfile(p):
        stop("no %s. Run codes/dcs_analysis.py first." % p)
    return list(csv.DictReader(open(p), delimiter="\t"))


fig = plt.figure(figsize=(7.2, 5.0))
gs = fig.add_gridspec(3, 1, hspace=0.52, left=0.085, right=0.985, top=0.92,
                      bottom=0.085)
ax_a, ax_b, ax_c = (fig.add_subplot(gs[i, 0]) for i in range(3))

STYLE = {("canonical", "sub1"): (COL["canonical"], "-"),
         ("canonical", "sub2"): (COL["canonical"], ":"),
         ("mixed", "sub1"): (COL["mixed"], ":"),
         ("mixed", "sub2"): (COL["mixed"], "-")}
for key, rec in DCS["chains"].items():
    rows = chain_table(key)
    if len(rows) != rec["residues"]:
        stop("%s holds %d rows and dcs_summary.json records %d residues."
             % (key, len(rows), rec["residues"]))
    colour, dash = STYLE[(rec["build"], rec["chain"])]
    ax_a.plot([int(r["position"]) for r in rows],
              [float(r["rmsf_A"]) for r in rows], lw=0.6, color=colour, ls=dash,
              label="%s, segment %s" % (LABEL[rec["build"]], rec["segment"]))
ax_a.set_ylabel("Fluctuation, $\\mathrm{\\AA}$")
ax_a.set_yscale("log")
ax_a.legend(frameon=False, fontsize=5.8, ncol=4, handlelength=1.6, loc="lower center",
            bbox_to_anchor=(0.5, 1.02), columnspacing=1.4)

ref = chain_table("canonical_sub1")
pos = [int(r["position"]) for r in ref]
ax_b.bar(pos, [100.0 * float(r["conservation"]) for r in ref], width=1.0,
         color=[BLOCK_COL[r["block"]] for r in ref], lw=0)
ax_b.set_ylabel("Orthologue identity, per cent")
ax_b.set_ylim(0, 105)

scored = [r for r in ref if r["dcs_0_10"]]
if not scored:
    stop("no position of canonical_sub1 carries a score, so the profiles are empty.")
sx = [int(r["position"]) for r in scored]
ax_c.fill_between(sx, [float(r["rcs_0_10"]) for r in scored], color=COL["canonical"],
                  alpha=0.35, lw=0, label="rigid conserved score")
ax_c.plot(sx, [float(r["dcs_0_10"]) for r in scored], lw=0.7, color=COL["mixed"],
          label="dynamic conserved score")
ax_c.set_ylabel("Score, 0 to 10")
ax_c.legend(frameon=False, fontsize=5.8, ncol=2, handlelength=1.6, loc="upper left")

for ax in (ax_a, ax_b, ax_c):
    ax.set_xlabel("Residue")
    ax.set_xlim(0, cc.CHAIN_LENGTH + 1)
for ax, s in zip((ax_a, ax_b, ax_c), "abc"):
    letter(ax, s)
save(fig, "FigureS3")
