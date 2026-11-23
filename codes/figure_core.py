import csv
import json
import os
import sys

import matplotlib.patches
import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_common as cc
from figure_style import COL, letter, picture, save, stop

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
RENDERS = os.path.join(RESULTS, "renders")
ENTRIES = ("2HZP", "3E9K")
SCLY_ENTRIES = ("3GZC", "3GZD")
(RET, RLO, RHI), (SUB, SLO, SHI), (DEL, DLO, DHI) = cc.BLOCKS
BLOCK_COL = {RET: "#bfbfbf", SUB: "#b2182b", DEL: "#4d4d4d"}
SPAN_ALPHA = {RET: 0.0, SUB: 0.16, DEL: 0.09}


def J(name):
    p = os.path.join(RESULTS, name)
    if not os.path.isfile(p):
        stop("no %s. The pipeline in the README writes it." % p)
    return json.load(open(p))


def T(name):
    p = os.path.join(RESULTS, name)
    if not os.path.isfile(p):
        stop("no %s. The pipeline in the README writes it." % p)
    return list(csv.DictReader(open(p), delimiter="\t"))


def blocks_behind(ax):
    for name, lo, hi in cc.BLOCKS:
        if SPAN_ALPHA[name] > 0:
            ax.axvspan(lo - 0.5, hi + 0.5, color=BLOCK_COL[name],
                       alpha=SPAN_ALPHA[name], lw=0)


def stack(ax, x, parts, width=0.62, label_above=0.0):
    bottom = np.zeros(len(x))
    for name in (RET, SUB, DEL):
        vals = np.asarray(parts[name], dtype=float)
        ax.bar(x, vals, bottom=bottom, width=width, color=BLOCK_COL[name], lw=0,
               label=name)
        if label_above:
            for i, v in enumerate(vals):
                if v >= label_above:
                    ax.text(x[i], bottom[i] + v / 2.0, "%.1f" % v, ha="center",
                            va="center", fontsize=6.0,
                            color="white" if name == DEL else "#222222")
        bottom = bottom + vals
    return bottom


PART = J("interface_partition.json")["entries"]["KYNU"]
POCKET = J("cofactor_pocket.json")["entries"]
AF3 = J("af3_scores.json")
BSA = T("interface_residue_bsa.tsv")
SWEEP = T("interface_parameter_sweep.tsv")
CONS = T("per_position_conservation.tsv")
HB = J("phosphate_hbonds.json")
EVID = J("existence_evidence.json")
SEQ = J("orthologue_sequences.json")[cc.ACCESSION]

if len(SEQ) != cc.CHAIN_LENGTH:
    stop("the recorded human sequence runs to %d residues and conservation_common "
         "sets the chain at %d." % (len(SEQ), cc.CHAIN_LENGTH))
for e in ENTRIES:
    got = sum(float(r["bsa_angstrom2"]) for r in BSA
              if r["protein"] == "KYNU" and r["entry"] == e)
    want = PART[e]["total_bsa"]
    if abs(got - want) > 1.0:
        stop("the per residue areas of %s add to %.1f and the partition reports %.1f."
             % (e, got, want))
    share = sum(PART[e]["fraction_" + n] for n in (RET, SUB, DEL))
    if abs(share - 1.0) > 0.002:
        stop("the three block shares of %s add to %.4f." % (e, share))
SWEEP_KEY = ("entry", "implementation", "algorithm", "probe_radius", "bsa_cutoff")
if len(SWEEP) != len({tuple(r[k] for k in SWEEP_KEY) for r in SWEEP}):
    stop("the parameter sweep carries repeated combinations.")
if len(SWEEP) != len({r["entry"] for r in SWEEP}) * len(
        {(r["implementation"], r["algorithm"]) for r in SWEEP}) * len(
        {r["probe_radius"] for r in SWEEP}) * len({r["bsa_cutoff"] for r in SWEEP}):
    stop("the parameter sweep is not the full grid of entry, method, probe and cutoff.")
for key, rec in POCKET.items():
    if rec["cofactor_sasa_dimer"] > rec["cofactor_sasa_monomer"] + 0.01:
        stop("%s buries less cofactor as a subunit than as a dimer." % key)
    parts = rec["buried_by_retained_block"] + rec["buried_by_lost_block"]
    if abs(parts - rec["buried_by_whole_partner"]) > 0.05:
        stop("the two blocks of %s bury %.2f and the whole partner buries %.2f."
             % (key, parts, rec["buried_by_whole_partner"]))
SOURCES = ("same_subunit", "partner_retained", "partner_lost")
for key, rec in HB["summary"].items():
    if sum(rec[s] for s in SOURCES) != rec["total"]:
        stop("the hydrogen bonds of %s by source do not add to its total." % key)
    seen = sum(1 for b in HB["bonds"] if "%s %s" % (b["protein"], b["entry"]) == key)
    if seen != rec["total"]:
        stop("%s lists %d hydrogen bonds and its summary counts %d."
             % (key, seen, rec["total"]))


fig = plt.figure(figsize=(7.2, 5.5))
gs = fig.add_gridspec(3, 2, height_ratios=[0.72, 1.22, 0.92], hspace=0.56,
                      wspace=0.22, left=0.085, right=0.985, top=0.95, bottom=0.08)
ax_a = fig.add_subplot(gs[0, :])
ax_b, ax_c = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
ax_d = fig.add_subplot(gs[2, :])

rec = EVID["summary"]["unique_peptide_record"]
width = SHI - SLO + 1
if int(rec["end"]) != SHI or EVID["summary"]["isoform_length"] != SHI:
    stop("the recorded isoform ends at %s and its peptide at %s, and the substituted "
         "block of conservation_common ends at %d."
         % (EVID["summary"]["isoform_length"], rec["end"], SHI))
hexa, tail = SEQ[SLO - 1:SHI], rec["peptide"][-width:]
if len(hexa) != width or len(tail) != width:
    stop("the substituted block spans %d residues and the two hexapeptides read %d "
         "and %d." % (width, len(hexa), len(tail)))

TRACKS = ((0.56, "canonical", cc.BLOCKS),
          (0.12, "isoform", ((RET, RLO, RHI), (SUB, SLO, SHI))))
for y0, name, blocks in TRACKS:
    for block, lo, hi in blocks:
        ax_a.add_patch(matplotlib.patches.Rectangle(
            (lo - 0.5, y0), hi - lo + 1, 0.30, facecolor=BLOCK_COL[block],
            edgecolor="white", lw=0.8))
        if hi - lo > 20:
            ax_a.text((lo + hi) / 2.0, y0 + 0.15, "%d to %d" % (lo, hi), ha="center",
                      va="center", fontsize=6.4,
                      color="white" if block == DEL else "#222222")
    ax_a.text(-6, y0 + 0.15, name, ha="right", va="center", fontsize=6.8,
              color="#333333")
ax_a.annotate("%d to %d  %s" % (SLO, SHI, hexa), xy=((SLO + SHI) / 2.0, 0.86),
              xytext=(SHI + 26, 0.99), ha="left", va="center", fontsize=6.4,
              color=BLOCK_COL[SUB],
              arrowprops=dict(arrowstyle="-", lw=0.6, color=BLOCK_COL[SUB]))
ax_a.annotate("%d to %d  %s" % (SLO, SHI, tail), xy=((SLO + SHI) / 2.0, 0.12),
              xytext=(SHI + 26, -0.01), ha="left", va="center", fontsize=6.4,
              color=BLOCK_COL[SUB],
              arrowprops=dict(arrowstyle="-", lw=0.6, color=BLOCK_COL[SUB]))
ax_a.set_xlim(-40, cc.CHAIN_LENGTH + 1)
ax_a.set_ylim(-0.1, 1.08)
ax_a.set_yticks([])
ax_a.set_xlabel("Residue")
ax_a.set_xticks([1, 100, 200, 300, 400, cc.CHAIN_LENGTH])
ax_a.spines["bottom"].set_bounds(1, cc.CHAIN_LENGTH)
for side in ("left", "top", "right"):
    ax_a.spines[side].set_visible(False)

picture(ax_b, os.path.join(RENDERS, "dimer_blocks.png"),
        "codes/render_core_struct.py draws the dimer coloured by block")

ax_c.axis("off")
WRAP_AT, LINE, PAD, GAP = 34, 1.0, 0.45, 0.30
cards = []
for lay in EVID["layers"]:
    text = lay["identifier"]
    if len(text) > WRAP_AT and ", " in text:
        head, _, rest = text.rpartition(", ")
        text = head + ",\n" + rest
    cards.append((lay["layer"], text, text.count("\n") + 1))
units = [2 * PAD + LINE * (1 + n) for _, _, n in cards]
scale = 1.0 / (sum(units) + GAP * (len(cards) - 1))
top = 1.0
for (title, text, _), height in zip(cards, units):
    box = height * scale
    ax_c.add_patch(matplotlib.patches.FancyBboxPatch(
        (0.02, top - box), 0.96, box, boxstyle="round,pad=0.004",
        facecolor="#f4f4f4", edgecolor="#d5d5d5", lw=0.6,
        transform=ax_c.transAxes))
    ax_c.text(0.055, top - (PAD + LINE * 0.5) * scale, title, fontsize=6.6,
              va="center", transform=ax_c.transAxes, color="#222222")
    ax_c.text(0.055, top - (PAD + LINE) * scale, text, fontsize=5.9, va="top",
              linespacing=1.25, transform=ax_c.transAxes, color=COL["mixed"])
    top -= box + GAP * scale

win = [r for r in CONS if SLO - 12 <= int(r["position"]) <= SHI + 13]
pos = [int(r["position"]) for r in win]
blocks_behind(ax_d)
ax_d.bar(pos, [float(r["percent_identical"]) for r in win], width=0.82,
         color=[BLOCK_COL[r["block"]] for r in win], lw=0)
ax_d.set_xlabel("Residue")
ax_d.set_ylabel("Orthologue identity, per cent")
ax_d.set_xlim(min(pos) - 0.6, max(pos) + 0.6)
ax_d.set_ylim(0, 105)
ax_d.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))

for ax, s in zip((ax_a, ax_b, ax_c, ax_d), "abcd"):
    letter(ax, s)
save(fig, "Figure1")


fig = plt.figure(figsize=(7.2, 4.5))
gs = fig.add_gridspec(2, 3, height_ratios=[0.90, 1.0], hspace=0.52, wspace=0.40,
                      left=0.085, right=0.985, top=0.94, bottom=0.14)
ax_a = fig.add_subplot(gs[0, :])
ax_b, ax_c, ax_d = (fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]),
                    fig.add_subplot(gs[1, 2]))

blocks_behind(ax_a)
for e, sign, col in zip(ENTRIES, (1.0, -1.0), (COL["canonical"], COL["mixed"])):
    rows = [r for r in BSA if r["protein"] == "KYNU" and r["entry"] == e]
    ax_a.bar([int(r["residue"]) for r in rows],
             [sign * float(r["bsa_angstrom2"]) for r in rows], width=1.7, color=col,
             lw=0)
    ax_a.text(0.005, 0.95 if sign > 0 else 0.05, e, transform=ax_a.transAxes,
              ha="left", va="top" if sign > 0 else "bottom", fontsize=6.6, color=col)
ax_a.axhline(0, color="#333333", lw=0.6)
ax_a.set_xlabel("Residue")
ax_a.set_ylabel("Buried area, $\\mathrm{\\AA}^2$")
ax_a.set_xlim(0, cc.CHAIN_LENGTH + 1)
ax_a.yaxis.set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % abs(v)))
ax_a.spines["bottom"].set_visible(False)
ax_a.tick_params(axis="x", length=0)

x = np.arange(len(ENTRIES))
stack(ax_b, x, {n: [100.0 * PART[e]["fraction_" + n] for e in ENTRIES]
                for n in (RET, SUB, DEL)}, width=0.56, label_above=6.0)
ax_b.set_xticks(x)
ax_b.set_xticklabels(ENTRIES)
ax_b.set_xlim(-0.62, len(ENTRIES) - 0.38)
ax_b.set_ylabel("Share of buried area, per cent")
ax_b.set_ylim(0, 100)
ax_b.legend(frameon=False, fontsize=6.0, loc="upper center", ncol=3,
            bbox_to_anchor=(0.5, 1.20), handlelength=1.0, columnspacing=0.9)

picture(ax_c, os.path.join(RENDERS, "interface_patch.png"),
        "codes/render_core_struct.py draws the interface patch")

methods = sorted({(r["implementation"], r["algorithm"]) for r in SWEEP})
for i, meth in enumerate(methods):
    for e, col, off in zip(ENTRIES, (COL["canonical"], COL["mixed"]), (-0.16, 0.16)):
        v = [100.0 * float(r["f_removed"]) for r in SWEEP
             if (r["implementation"], r["algorithm"]) == meth and r["entry"] == e]
        ax_d.scatter(np.full(len(v), i + off) + np.linspace(-0.09, 0.09, len(v)), v,
                     s=5, color=col, lw=0, alpha=0.8, label=e if i == 0 else None)
allv = [100.0 * float(r["f_removed"]) for r in SWEEP]
for edge in (min(allv), max(allv)):
    ax_d.axhline(edge, color="#999999", lw=0.6, ls=":")
ax_d.set_xticks(range(len(methods)))
ax_d.set_xticklabels(["%s\n%s" % m for m in methods], fontsize=5.6)
ax_d.set_ylabel("Removed share, per cent")
ax_d.set_xlim(-0.5, len(methods) - 0.5)
ax_d.set_ylim(min(allv) - 0.5, max(allv) + 1.1)
ax_d.legend(frameon=False, fontsize=6.0, handlelength=0.8, loc="upper left",
            ncol=2, columnspacing=0.9)
ax_d.annotate("%d settings, %.1f to %.1f per cent"
              % (len(SWEEP), min(allv), max(allv)), xy=(0.5, 0.03),
              xycoords="axes fraction", ha="center", fontsize=6.0, color="#555555")


for ax, s in zip((ax_a, ax_b, ax_c, ax_d), "abcd"):
    letter(ax, s)
save(fig, "Figure2")


fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 3, width_ratios=[1.22, 1.0, 1.0], hspace=0.46, wspace=0.42,
                      left=0.055, right=0.985, top=0.94, bottom=0.13)
ax_a = fig.add_subplot(gs[:, 0])
ax_b, ax_c = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
ax_d, ax_e = fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])

picture(ax_a, os.path.join(RENDERS, "pocket_environment.png"),
        "codes/render_core_struct.py draws the cofactor environment")

STATES = (("cofactor_sasa_monomer", "Subunit\nalone"),
          ("cofactor_sasa_dimer", "Full\ndimer"),
          ("cofactor_sasa_isoform_mimic", "Isoform\nmimic"))
ALL_KEYS = ["KYNU " + e for e in ENTRIES] + ["SCLY " + e for e in SCLY_ENTRIES]


def states_panel(ax, keys, title):
    pos = np.arange(len(STATES))
    for i, (key, col) in enumerate(zip(keys, (COL["canonical"], COL["mixed"]))):
        ax.bar(pos + (i - 0.5) * 0.36, [POCKET[key][k] for k, _ in STATES],
               width=0.34, color=col, lw=0, label=key.split()[1])
    ax.set_xticks(pos)
    ax.set_xticklabels([n for _, n in STATES], fontsize=6.0)
    ax.set_ylabel("Cofactor exposure, $\\mathrm{\\AA}^2$")
    ax.set_title(title, fontsize=7.0, pad=3)
    ax.set_ylim(0, max(POCKET[k][s] for k in ALL_KEYS for s, _ in STATES) * 1.22)
    ax.legend(frameon=False, fontsize=6.0, handlelength=1.0, loc="upper right")


states_panel(ax_b, ["KYNU " + e for e in ENTRIES], "Kynureninase")
states_panel(ax_c, ["SCLY " + e for e in SCLY_ENTRIES], "Selenocysteine lyase")

pos = np.arange(len(ALL_KEYS))
ax_d.bar(pos, [POCKET[k]["buried_by_retained_block"] for k in ALL_KEYS], width=0.6,
         color=BLOCK_COL[RET], lw=0, label="retained block")
ax_d.bar(pos, [POCKET[k]["buried_by_lost_block"] for k in ALL_KEYS], width=0.6,
         bottom=[POCKET[k]["buried_by_retained_block"] for k in ALL_KEYS],
         color=BLOCK_COL[DEL], lw=0, label="removed block")
ax_d.set_xticks(pos)
ax_d.set_xticklabels([k.split()[1] for k in ALL_KEYS], fontsize=5.8, rotation=28,
                     ha="right")
ax_d.set_ylabel("Buried by the partner, $\\mathrm{\\AA}^2$")
ax_d.legend(frameon=False, fontsize=5.8, handlelength=1.0, loc="upper right")
ax_d.set_ylim(0, max(POCKET[k]["buried_by_whole_partner"] for k in ALL_KEYS) * 1.30)

SOURCE_LOOK = (("same_subunit", "#bfbfbf", "same subunit"),
               ("partner_retained", "#7f7f7f", "partner, retained block"),
               ("partner_lost", BLOCK_COL[DEL], "partner, removed block"))
bottom = np.zeros(len(ALL_KEYS))
for key, col, lab in SOURCE_LOOK:
    vals = np.array([HB["summary"][k][key] for k in ALL_KEYS], dtype=float)
    ax_e.bar(pos, vals, bottom=bottom, width=0.6, color=col, lw=0, label=lab)
    bottom += vals
ax_e.set_xticks(pos)
ax_e.set_xticklabels([k.split()[1] for k in ALL_KEYS], fontsize=5.8, rotation=28,
                     ha="right")
ax_e.set_ylabel("Bonds to the cofactor phosphate")
ax_e.set_ylim(0, bottom.max() * 1.62)
ax_e.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True, nbins=5))
ax_e.legend(frameon=False, fontsize=5.6, handlelength=1.0, loc="upper right",
            labelspacing=0.3)

for ax, s in zip((ax_a, ax_b, ax_c, ax_d, ax_e), "abcde"):
    letter(ax, s)
save(fig, "Figure3")


JOB_ORDER = ("canonical__canonical", "canonical__isoform2", "canonical__trunc_cterm",
             "canonical__trunc_nterm", "canonical__tail_scramble_01",
             "canonical__tail_scramble_02", "canonical__tail_scramble_03",
             "isoform2__isoform2")
JOB_LABEL = {"canonical__canonical": "canonical", "canonical__isoform2": "isoform",
             "canonical__trunc_cterm": "tail removed",
             "canonical__trunc_nterm": "N face removed",
             "canonical__tail_scramble_01": "permuted 1",
             "canonical__tail_scramble_02": "permuted 2",
             "canonical__tail_scramble_03": "permuted 3",
             "isoform2__isoform2": "isoform pair"}
CLAIM = "canonical__isoform2"
rows = {r["job"]: r for r in AF3["scores"]}
if set(rows) != set(JOB_ORDER):
    stop("af3_scores.json carries the jobs %s and this figure is laid out for %s."
         % (sorted(rows), sorted(JOB_ORDER)))

fig = plt.figure(figsize=(7.2, 5.1))
gs = fig.add_gridspec(2, 3, height_ratios=[1.22, 1.0], width_ratios=[1.32, 1.0, 1.0],
                      hspace=0.68, wspace=0.40, left=0.065, right=0.985, top=0.94,
                      bottom=0.14)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1:])
ax_c, ax_d, ax_e = (fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]),
                    fig.add_subplot(gs[1, 2]))

picture(ax_a, os.path.join(RENDERS, "predicted_heterodimer.png"),
        "codes/render_core_struct.py draws the predicted heterodimer")

x = np.arange(len(JOB_ORDER))
claim_at = JOB_ORDER.index(CLAIM)
labels = [JOB_LABEL[j] for j in JOB_ORDER]
for ax in (ax_b, ax_c):
    ax.axvspan(claim_at - 0.5, claim_at + 0.5, color=BLOCK_COL[SUB], alpha=0.10, lw=0)
for i, (key, col) in enumerate(zip(("iptm", "pDockQ", "pDockQ2"),
                                   ("#bfbfbf", COL["canonical"], COL["mixed"]))):
    ax_b.bar(x + (i - 1) * 0.27, [float(rows[j][key]) for j in JOB_ORDER], width=0.26,
             color=col, lw=0, label=key)
ax_b.set_xticks(x)
ax_b.set_xticklabels(labels, fontsize=5.8, rotation=26, ha="right")
ax_b.set_ylabel("Interface score")
ax_b.set_ylim(0, 1.0)
ax_b.legend(frameon=False, fontsize=6.0, ncol=3, handlelength=1.0, loc="upper center",
            bbox_to_anchor=(0.5, 1.16), columnspacing=0.9)

stack(ax_c, x, {n: [float(rows[j]["partner_iface_" + n]) for j in JOB_ORDER]
                for n in (RET, SUB, DEL)}, width=0.58)
ax_c.set_xticks(x)
ax_c.set_xticklabels(labels, fontsize=5.2, rotation=30, ha="right")
ax_c.set_ylabel("Partner interface residues")
ax_c.legend(frameon=False, fontsize=5.6, handlelength=1.0, loc="upper right")

pred = AF3["detail"]["canonical__canonical"]["interface"]["partner_interface_by_block"]
total = float(sum(pred.values()))
if not total:
    stop("the predicted homodimer records no partner interface residue by block.")
names = ["predicted"] + list(ENTRIES)
parts = {n: [100.0 * pred[n] / total] + [100.0 * PART[e]["fraction_" + n]
                                         for e in ENTRIES] for n in (RET, SUB, DEL)}
stack(ax_d, np.arange(len(names)), parts, width=0.58)
ax_d.set_xticks(np.arange(len(names)))
ax_d.set_xticklabels(names, fontsize=6.0)
ax_d.set_ylabel("Share of the interface, per cent")
ax_d.set_ylim(0, 100)

for job, col, lab in ((CLAIM, COL["mixed"], JOB_LABEL[CLAIM]),
                      ("canonical__trunc_nterm", COL["canonical"],
                       JOB_LABEL["canonical__trunc_nterm"])):
    v = AF3["detail"][job]["interface"].get("interface_pae_pairs")
    if not v:
        stop("af3_scores.json carries no per pair error for %s. Run "
             "codes/score_af3.py again." % job)
    ax_e.hist(v, bins=np.linspace(0, 16, 33), histtype="step", lw=1.0, color=col,
              label=lab)
ax_e.set_xlabel("Interface pair error, $\\mathrm{\\AA}$")
ax_e.set_ylabel("Residue pairs")
ax_e.set_xlim(0, 16)
ax_e.legend(frameon=False, fontsize=5.8, handlelength=1.0)

for ax, s in zip((ax_a, ax_b, ax_c, ax_d, ax_e), "abcde"):
    letter(ax, s)
save(fig, "Figure4")
