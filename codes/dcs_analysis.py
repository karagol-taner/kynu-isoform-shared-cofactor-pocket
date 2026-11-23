import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_common as cc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
MD = os.path.join(RESULTS, "md")
RMSF_DIR = os.path.join(MD, "rmsf")
OUT_DIR = os.path.join(RESULTS, "dcs")
BUILDS = ("canonical", "mixed")
CHAINS = ("sub1", "sub2")
TOP_FRACTION = 0.10

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "HSD": "H", "HSE": "H", "HSP": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
SIBLING_ENTRIES = ("3GZC", "3GZD")
HINGE = ("G", "P")
HYDROPHOBIC = ("A", "V", "L", "I", "M", "F", "W", "C")


def need(path, what):
    if not os.path.isfile(path):
        raise SystemExit("STOP: no %s, so %s. " % (path, what))
    return path


def tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, L.rstrip("\n").split("\t"))) for L in fh if L.strip()]


def deposition_span():
    import glob
    paths = sorted(glob.glob(os.path.join(RESULTS, "assemblies", "*_A.pdb")))
    seen = {}
    for path in paths:
        entry = os.path.basename(path).split("_")[0]
        if entry in SIBLING_ENTRIES:
            continue
        res = set()
        for L in open(path, errors="replace"):
            if L.startswith(("ATOM", "HETATM")):
                try:
                    res.add(int(L[22:26]))
                except ValueError:
                    continue
        if res:
            seen[entry] = (min(res), max(res))
    if not seen:
        raise SystemExit("STOP: no deposited chain of this protein under %s, so the "
                         "span to score over cannot be taken from the structures. "
                         "The files it looked for are named <entry>_A.pdb and the "
                         "entries of the sibling protein, %s, are skipped."
                         % (os.path.join(RESULTS, "assemblies"),
                            ", ".join(sorted(SIBLING_ENTRIES))))
    lo = max(v[0] for v in seen.values())
    hi = min(v[1] for v in seen.values())
    if not (1 <= lo < hi <= cc.CHAIN_LENGTH):
        raise SystemExit("STOP: the deposited chains resolve %d to %d, which is not "
                         "a span inside a chain of %d." % (lo, hi, cc.CHAIN_LENGTH))
    return lo, hi, dict(sorted(seen.items()))


def minmax(vals):
    lo, hi = min(vals), max(vals)
    span = hi - lo
    return [0.0 for _ in vals] if span == 0 else [(v - lo) / span for v in vals]


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return 0.0 if sxx == 0 or syy == 0 else sxy / (sxx * syy) ** 0.5


cons_rows = tsv(need(os.path.join(RESULTS, "per_position_conservation.tsv"),
                     "there is no conservation to score against. Run "
                     "codes/per_position_conservation.py first"))
CONS = {int(r["position"]): r for r in cons_rows}

seqs = json.load(open(need(os.path.join(RESULTS, "orthologue_sequences.json"),
                           "the human sequence the residue numbering is checked "
                           "against is not on disk. Run "
                           "codes/per_position_conservation.py first")))
if cc.ACCESSION not in seqs:
    raise SystemExit("STOP: %s carries no sequence for %s."
                     % (os.path.join(RESULTS, "orthologue_sequences.json"),
                        cc.ACCESSION))
HUMAN = seqs[cc.ACCESSION]
if len(HUMAN) != cc.CHAIN_LENGTH:
    raise SystemExit("STOP: the cached %s sequence is %d residues and this project "
                     "works with %d." % (cc.ACCESSION, len(HUMAN), cc.CHAIN_LENGTH))

SPAN_LO, SPAN_HI, SPAN_ENTRIES = deposition_span()

gate = json.load(open(need(os.path.join(MD, "kill_criterion.json"),
                           "whether the cofactor site may be scored is not known. "
                           "Run the part 3 notebook first")))
if "pocket_evaluable" not in gate:
    raise SystemExit("STOP: %s carries no pocket_evaluable flag, so it was written "
                     "by an earlier part 3. Run the part 3 notebook again."
                     % os.path.join(MD, "kill_criterion.json"))
POCKET_EVALUABLE = bool(gate["pocket_evaluable"])

reports, rmsf = {}, {}
for build in BUILDS:
    reports[build] = json.load(open(need(
        os.path.join(MD, "%s_report.json" % build),
        "which chain of %s is the canonical one is not known" % build)))
    for chain in CHAINS:
        rmsf[(build, chain)] = tsv(need(
            os.path.join(RMSF_DIR, build, "rmsf_%s.tsv" % chain),
            "there is no per residue fluctuation for %s of %s. Run "
            "codes/colab/kynu_MD_4_rmsf.ipynb" % (chain, build)))

summary = {"method": {
    "definition": "DCS = conservation * norm_dynamics, RCS = conservation * "
                  "(1 - norm_dynamics), each rescaled 0 to 10 over the chain",
    "dynamics": "C-alpha RMSF per residue, each chain superposed on itself",
    "replicas": 1,
    "replicas_in_published_method": 3,
    "conservation": "orthologue identity fraction at the aligned position",
    "orthologues": len(json.load(open(
        os.path.join(RESULTS, "variants_and_conservation.json")))["orthologues"]),
    "top_fraction": TOP_FRACTION,
    "schiff_base_residue": int(json.load(open(
        os.path.join(MD, "canonical_report.json")))["roles"]["scored_copy"][1]),
    "scored_span": [SPAN_LO, SPAN_HI],
    "scored_span_source": "resolved in every deposited chain of this protein",
    "scored_span_entries": SPAN_ENTRIES,
    "pocket_evaluable": POCKET_EVALUABLE}, "chains": {}}

os.makedirs(OUT_DIR, exist_ok=True)
scored = {}
for build in BUILDS:
    roles = reports[build]["roles"]
    built_residue = int(roles["scored_copy"][1])
    for chain in CHAINS:
        seg = roles["canonical" if chain == "sub1" else "partner"]
        rows = rmsf[(build, chain)]
        want = int(roles["lengths"][seg])
        if len(rows) != want:
            raise SystemExit("STOP: the fluctuation file for %s of %s carries %d "
                             "residues and the run report calls that chain %d long."
                             % (chain, build, len(rows), want))

        retained_end = cc.BLOCKS[0][2]
        is_isoform = want < cc.CHAIN_LENGTH
        novel = sorted(int(r["residue"]) for r in rows
                       if is_isoform and int(r["residue"]) > retained_end)

        mismatch = []
        for r in rows:
            pos = int(r["residue"])
            one = THREE_TO_ONE.get(r["resname"])
            if pos == built_residue and one is None:
                continue
            if pos in novel:
                if one is None:
                    mismatch.append("%d is %s in the isoform tail, which is not a "
                                    "residue this reads" % (pos, r["resname"]))
                continue
            if pos > cc.CHAIN_LENGTH:
                mismatch.append("%d is past the end of the human chain" % pos)
            elif one is None:
                mismatch.append("%d is %s, which is not a residue this reads"
                                % (pos, r["resname"]))
            elif one != HUMAN[pos - 1]:
                mismatch.append("%d is %s and the human chain has %s there"
                                % (pos, one, HUMAN[pos - 1]))
        if mismatch:
            raise SystemExit("STOP: %s of %s does not line up with the human "
                             "sequence at %d positions, so the fluctuation cannot "
                             "be paired with conservation. The first few are:\n  %s"
                             % (chain, build, len(mismatch),
                                "\n  ".join(mismatch[:5])))
        if is_isoform:
            tail = "".join(THREE_TO_ONE.get(r["resname"], "?") for r in rows
                           if int(r["residue"]) in novel)
            canon = HUMAN[retained_end:want]
            if len(tail) != want - retained_end:
                raise SystemExit("STOP: the isoform chain of %s runs to %d and its "
                                 "retained block ends at %d, which leaves %d "
                                 "positions, and %d of them were read."
                                 % (build, want, retained_end,
                                    want - retained_end, len(tail)))
            if tail == canon:
                raise SystemExit("STOP: the tail of the isoform chain of %s reads %s, "
                                 "which is the canonical sequence at those positions. "
                                 "That chain is then not the isoform the construct "
                                 "describes." % (build, tail))

        out = []
        for r in rows:
            pos, val = int(r["residue"]), float(r["rmsf_A"])
            c = CONS.get(pos)
            depth = int(c["aligned"]) if c else 0
            letter = THREE_TO_ONE.get(r["resname"])
            if letter is None and pos == built_residue:
                letter = HUMAN[pos - 1]
            out.append({"position": pos, "resname": r["resname"],
                        "residue": letter or r["resname"],
                        "rmsf_A": val,
                        "block": c["block"] if c else "novel",
                        "aligned": depth,
                        "conservation": (float(c["percent_identical"]) / 100.0
                                         if c and depth else None),
                        "novel_sequence": pos in novel})

        nd = minmax([r["rmsf_A"] for r in out])
        for r, v in zip(out, nd):
            r["norm_dynamics"] = v
        usable = [r for r in out if r["conservation"] is not None
                  and not r["novel_sequence"]
                  and SPAN_LO <= r["position"] <= SPAN_HI]
        if len(usable) < 2:
            raise SystemExit("STOP: %s of %s has %d positions carrying both a "
                             "fluctuation and a conservation value, which is not "
                             "enough to score." % (chain, build, len(usable)))
        raw_d = [r["conservation"] * r["norm_dynamics"] for r in usable]
        raw_r = [r["conservation"] * (1.0 - r["norm_dynamics"]) for r in usable]
        for r, d, s in zip(usable, minmax(raw_d), minmax(raw_r)):
            r["dcs"] = d * 10.0
            r["rcs"] = s * 10.0

        cd = pearson([r["norm_dynamics"] for r in usable], [r["dcs"] for r in usable])
        cr = pearson([r["norm_dynamics"] for r in usable], [r["rcs"] for r in usable])
        if not (cd > 0 and cr < 0):
            raise SystemExit("STOP: the dynamic score correlates with the dynamics "
                            "at %.3f and the rigid score at %.3f in %s of %s, and "
                            "the definitions require the first positive and the "
                            "second negative." % (cd, cr, chain, build))

        path = os.path.join(OUT_DIR, "%s_%s.tsv" % (build, chain))
        with open(path, "w") as fh:
            fh.write("position\tresidue\tresname\tblock\taligned\tconservation\t"
                     "rmsf_A\tnorm_dynamics\tdcs_0_10\trcs_0_10\n")
            for r in out:
                fh.write("%d\t%s\t%s\t%s\t%d\t%s\t%.4f\t%.6f\t%s\t%s\n"
                         % (r["position"], r["residue"], r["resname"], r["block"],
                            r["aligned"],
                            "" if r["conservation"] is None
                            else "%.4f" % r["conservation"],
                            r["rmsf_A"], r["norm_dynamics"],
                            "%.3f" % r["dcs"] if "dcs" in r else "",
                            "%.3f" % r["rcs"] if "rcs" in r else ""))

        vals = [r["rmsf_A"] for r in out]
        k = max(1, int(round(TOP_FRACTION * len(usable))))
        top_d = sorted(usable, key=lambda r: -r["dcs"])[:k]
        top_r = sorted(usable, key=lambda r: -r["rcs"])[:k]
        summary["chains"]["%s_%s" % (build, chain)] = {
            "build": build, "chain": chain, "segment": seg,
            "residues": len(out), "scored": len(usable),
            "novel_sequence_positions": sorted(novel),
            "no_conservation": sorted(r["position"] for r in out
                                      if r["conservation"] is None),
            "outside_scored_span": sorted(r["position"] for r in out
                                          if not SPAN_LO <= r["position"] <= SPAN_HI),
            "rmsf_A": {"mean": sum(vals) / len(vals), "min": min(vals),
                       "max": max(vals)},
            "top_n": k,
            "dynamic_conserved": [
                {"position": r["position"], "residue": r["residue"],
                 "block": r["block"], "dcs": r["dcs"], "rmsf_A": r["rmsf_A"],
                 "conservation": r["conservation"]} for r in top_d],
            "rigid_conserved": [
                {"position": r["position"], "residue": r["residue"],
                 "block": r["block"], "rcs": r["rcs"], "rmsf_A": r["rmsf_A"],
                 "conservation": r["conservation"]} for r in top_r],
            "hinge_fraction_in_dynamic": sum(1 for r in top_d
                                             if r["residue"] in HINGE) / float(k),
            "hinge_fraction_in_chain": sum(1 for r in usable
                                           if r["residue"] in HINGE) / float(len(usable)),
            "hydrophobic_fraction_in_rigid": sum(1 for r in top_r
                                                 if r["residue"] in HYDROPHOBIC) / float(k),
            "hydrophobic_fraction_in_chain": sum(1 for r in usable
                                                 if r["residue"] in HYDROPHOBIC) / float(len(usable))}
        scored[(build, chain)] = {r["position"]: r for r in usable}

A, B = scored[("canonical", "sub1")], scored[("mixed", "sub1")]
both = sorted(set(A) & set(B))
if len(both) < 2:
    raise SystemExit("STOP: the canonical chain is scored at %d positions in the "
                     "canonical build and %d in the mixed, with %d in common, which "
                     "is not enough to compare." % (len(A), len(B), len(both)))

iface = {}
for r in tsv(os.path.join(RESULTS, "interface_residue_bsa.tsv")):
    if r["protein"] != "KYNU":
        continue
    iface.setdefault(int(r["residue"]), set()).add(r["entry"])
entries = sorted({e for v in iface.values() for e in v})
consensus = sorted(p for p, es in iface.items() if len(es) == len(entries))

def ranks(d, positions):
    order = sorted(positions, key=lambda p: d[p]["dcs"])
    return {p: i + 1 for i, p in enumerate(order)}

ra, rb = ranks(A, both), ranks(B, both)
d_rmsf = [B[p]["rmsf_A"] - A[p]["rmsf_A"] for p in both]


def over(positions, label):
    ps = [p for p in positions if p in A and p in B]
    if not ps:
        return {"label": label, "positions": 0}
    dr = [B[p]["rmsf_A"] - A[p]["rmsf_A"] for p in ps]
    return {"label": label, "positions": len(ps),
            "rmsf_A_canonical": sum(A[p]["rmsf_A"] for p in ps) / len(ps),
            "rmsf_A_mixed": sum(B[p]["rmsf_A"] for p in ps) / len(ps),
            "delta_rmsf_A_mean": sum(dr) / len(dr),
            "delta_rmsf_A_mean_abs": sum(abs(v) for v in dr) / len(dr),
            "delta_rmsf_A_max": max(dr, key=abs),
            "rank_shift_mean": sum(abs(rb[p] - ra[p]) for p in ps) / float(len(ps))}

comparison = {
    "chain": "sub1, the canonical subunit in both builds",
    "positions": len(both),
    "dcs_rank_correlation": pearson([ra[p] for p in both], [rb[p] for p in both]),
    "rmsf_correlation": pearson([A[p]["rmsf_A"] for p in both],
                                [B[p]["rmsf_A"] for p in both]),
    "delta_rmsf_A_mean": sum(d_rmsf) / len(d_rmsf),
    "delta_rmsf_A_mean_abs": sum(abs(v) for v in d_rmsf) / len(d_rmsf),
    "blocks": [over([p for p in both if A[p]["block"] == name], name)
               for name, lo, hi in cc.BLOCKS],
    "interface_entries": entries,
    "interface_consensus": over(consensus, "interface in all %d structures"
                                % len(entries)),
    "interface_consensus_positions": consensus,
    "moved_most": [{"position": p, "residue": A[p]["residue"],
                    "block": A[p]["block"],
                    "rmsf_A_canonical": A[p]["rmsf_A"],
                    "rmsf_A_mixed": B[p]["rmsf_A"],
                    "delta_rmsf_A": B[p]["rmsf_A"] - A[p]["rmsf_A"]}
                   for p in sorted(both, key=lambda q: -abs(B[q]["rmsf_A"]
                                                            - A[q]["rmsf_A"]))[:10]]}

contacts = [p for p in cc.CROSS_SUBUNIT_COFACTOR if p in A and p in B]
if POCKET_EVALUABLE:
    comparison["cofactor_contacts"] = over(contacts, "cross subunit cofactor contact")
else:
    comparison["cofactor_contacts_withheld"] = {
        "positions": contacts,
        "reason": "the cofactor left its site in both runs, which %s records as "
                  "leaving the pocket arm unscored, so the fluctuation of the "
                  "residues that coordinate it is measured against a cofactor in an "
                  "undetermined position" % os.path.basename(
                      os.path.join(MD, "kill_criterion.json"))}
summary["canonical_against_mixed"] = comparison

C1, C2 = scored[("canonical", "sub1")], scored[("canonical", "sub2")]
pair = sorted(set(C1) & set(C2))
if pair:
    dp = [C2[p]["rmsf_A"] - C1[p]["rmsf_A"] for p in pair]
    rp1, rp2 = ranks(C1, pair), ranks(C2, pair)

    def over_pair(positions, label):
        ps = [p for p in positions if p in C1 and p in C2]
        if not ps:
            return {"label": label, "positions": 0}
        d = [C2[p]["rmsf_A"] - C1[p]["rmsf_A"] for p in ps]
        return {"label": label, "positions": len(ps),
                "delta_rmsf_A_mean": sum(d) / len(d),
                "delta_rmsf_A_mean_abs": sum(abs(v) for v in d) / len(d)}

    summary["two_copies_in_one_build"] = {
        "build": "canonical",
        "what": "the two copies of the same chain in one box, as the size of "
                "difference a second copy gives with no truncation",
        "positions": len(pair),
        "dcs_rank_correlation": pearson([rp1[p] for p in pair],
                                        [rp2[p] for p in pair]),
        "rmsf_correlation": pearson([C1[p]["rmsf_A"] for p in pair],
                                    [C2[p]["rmsf_A"] for p in pair]),
        "delta_rmsf_A_mean": sum(dp) / len(dp),
        "delta_rmsf_A_mean_abs": sum(abs(v) for v in dp) / len(dp),
        "blocks": [over_pair([p for p in pair if C1[p]["block"] == name], name)
                   for name, lo, hi in cc.BLOCKS],
        "interface_consensus": over_pair(consensus, "interface in all %d structures"
                                         % len(entries))}

ISO, PARTNER = scored[("mixed", "sub2")], scored[("canonical", "sub2")]
shared = sorted(p for p in set(ISO) & set(PARTNER)
                if ISO[p]["block"] == cc.BLOCKS[0][0])
if len(shared) < 2:
    raise SystemExit("STOP: the isoform chain and a full length partner share %d "
                     "scored positions in the retained block, which is not enough "
                     "to compare." % len(shared))
iso_v = [ISO[p]["rmsf_A"] for p in shared]
par_v = [PARTNER[p]["rmsf_A"] for p in shared]
own_v = [scored[("canonical", "sub1")][p]["rmsf_A"] for p in shared]
higher = sum(1 for a, b in zip(iso_v, par_v) if a > b)
baseline_higher = sum(1 for a, b in zip(own_v, par_v) if a > b)
summary["isoform_against_full_partner"] = {
    "what": "the retained block of the isoform chain against the same residues of "
            "a full length partner, each the partner copy of its own dimer",
    "positions": len(shared),
    "first_position": shared[0], "last_position": shared[-1],
    "isoform_rmsf_A_mean": sum(iso_v) / len(iso_v),
    "partner_rmsf_A_mean": sum(par_v) / len(par_v),
    "ratio": (sum(iso_v) / len(iso_v)) / (sum(par_v) / len(par_v)),
    "positions_isoform_higher": higher,
    "baseline": {
        "what": "the other copy of the canonical build over the same residues",
        "rmsf_A_mean": sum(own_v) / len(own_v),
        "ratio": (sum(own_v) / len(own_v)) / (sum(par_v) / len(par_v)),
        "positions_higher": baseline_higher}}

json.dump(summary, open(os.path.join(RESULTS, "dcs_summary.json"), "w"), indent=2)
print("wrote %s and %d per chain tables in %s"
      % (os.path.join(RESULTS, "dcs_summary.json"), len(scored), OUT_DIR))
for key, s in sorted(summary["chains"].items()):
    print("  %-18s %3d residues, %3d scored, RMSF %.2f A mean, hinge in top DCS "
          "%.2f against %.2f in the chain"
          % (key, s["residues"], s["scored"], s["rmsf_A"]["mean"],
             s["hinge_fraction_in_dynamic"], s["hinge_fraction_in_chain"]))
c = summary["canonical_against_mixed"]
print()
print("the canonical subunit in both builds, over %d shared positions:" % c["positions"])
print("  DCS ranking correlates at %.3f, the fluctuation itself at %.3f"
      % (c["dcs_rank_correlation"], c["rmsf_correlation"]))
print("  mean change in fluctuation %+.3f A, mean absolute change %.3f A"
      % (c["delta_rmsf_A_mean"], c["delta_rmsf_A_mean_abs"]))
for b in c["blocks"]:
    if b["positions"]:
        print("  %-12s %3d positions, %.2f to %.2f A, mean change %+.3f A"
              % (b["label"], b["positions"], b["rmsf_A_canonical"],
                 b["rmsf_A_mixed"], b["delta_rmsf_A_mean"]))
i = c["interface_consensus"]
if i["positions"]:
    print("  %-12s %3d positions, %.2f to %.2f A, mean change %+.3f A"
          % ("interface", i["positions"], i["rmsf_A_canonical"],
             i["rmsf_A_mixed"], i["delta_rmsf_A_mean"]))
t = summary.get("two_copies_in_one_build")
if t:
    print()
    print("the two copies of the same chain in the canonical build, as the size a "
          "second copy gives:")
    print("  DCS ranking correlates at %.3f, the fluctuation itself at %.3f"
          % (t["dcs_rank_correlation"], t["rmsf_correlation"]))
    print("  mean change in fluctuation %+.3f A, mean absolute change %.3f A"
          % (t["delta_rmsf_A_mean"], t["delta_rmsf_A_mean_abs"]))
    for b in t["blocks"]:
        if b["positions"]:
            print("  %-12s %3d positions, mean change %+.3f A, mean absolute %.3f A"
                  % (b["label"], b["positions"], b["delta_rmsf_A_mean"],
                     b["delta_rmsf_A_mean_abs"]))
    i = t["interface_consensus"]
    if i["positions"]:
        print("  %-12s %3d positions, mean change %+.3f A, mean absolute %.3f A"
              % ("interface", i["positions"], i["delta_rmsf_A_mean"],
                 i["delta_rmsf_A_mean_abs"]))
    print()
i = summary["isoform_against_full_partner"]
print("the retained block of the isoform chain against the same residues of a full "
      "length partner,")
print("  %d positions, %d to %d: %.3f A against %.3f A, a ratio of %.2f"
      % (i["positions"], i["first_position"], i["last_position"],
         i["isoform_rmsf_A_mean"], i["partner_rmsf_A_mean"], i["ratio"]))
print("  the isoform copy is the more mobile at %d of %d positions"
      % (i["positions_isoform_higher"], i["positions"]))
print("  the other copy of the canonical build gives a ratio of %.2f at %d of %d"
      % (i["baseline"]["ratio"], i["baseline"]["positions_higher"], i["positions"]))
print()
if "cofactor_contacts_withheld" in c:
    print("  the cofactor contact positions %s are not reported, since the cofactor "
          "left its site" % c["cofactor_contacts_withheld"]["positions"])
