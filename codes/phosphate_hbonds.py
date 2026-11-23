import ast, csv, json, os

import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(ROOT, "structures")
RESULTS = os.path.join(ROOT, "results")

_conf = open(os.path.join(ROOT, "codes", "interface_partition.py")).read()
SPECS = next(ast.literal_eval(n.value) for n in ast.parse(_conf).body
             if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", None) == "SPECS")
COFACTOR_NAMES = next(ast.literal_eval(n.value) for n in ast.parse(_conf).body
                      if isinstance(n, ast.Assign)
                      and getattr(n.targets[0], "id", None) == "COFACTOR_NAMES")
HBOND_CUTOFF = 3.5
PHOSPHATE_OXYGENS = ("O1P", "O2P", "O3P", "OP1", "OP2", "OP3")
DONOR_ELEMENTS = ("N", "O")

rows, summary = [], {}
for protein, spec in SPECS.items():
    retained_end = spec["retained_end"]
    lost_from = (spec["substituted"][0] if spec["substituted"]
                 else spec["deleted"][0])
    lost_to = spec["deleted"][1]
    for entry in spec["primary"]:
        path = os.path.join(STRUCT, entry + ".cif")
        if not os.path.isfile(path):
            raise SystemExit("STOP: no %s. Run the fetch and check script first." % path)
        st = gemmi.read_structure(path)
        st.setup_entities()
        st.transform_to_assembly("1", gemmi.HowToNameCopiedChain.AddNumber)
        st.remove_waters()
        atoms = [(ch.name, r.name, r.seqid.num, a.name, a.pos, a.element.name)
                 for ch in st[0] for r in ch for a in r]
        cofactors = sorted({(c, n, i) for c, n, i, _an, _p, _e in atoms
                            if n in COFACTOR_NAMES})
        if not cofactors:
            raise SystemExit("STOP: %s holds no cofactor of %s" % (entry, COFACTOR_NAMES))
        host_chain, cof_name, cof_num = cofactors[0]
        phosphate = [x for x in atoms if x[0] == host_chain and x[1] == cof_name
                     and x[3] in PHOSPHATE_OXYGENS]
        if not phosphate:
            raise SystemExit("STOP: %s %s%d carries none of the phosphate oxygens %s"
                             % (entry, cof_name, cof_num, PHOSPHATE_OXYGENS))
        best = {}
        for _c, _n, _i, oxygen, pos, _e in phosphate:
            for ch2, rn2, ri2, an2, pos2, el2 in atoms:
                if rn2 in COFACTOR_NAMES or el2 not in DONOR_ELEMENTS:
                    continue
                d = pos.dist(pos2)
                if d > HBOND_CUTOFF:
                    continue
                key = (ch2, rn2, ri2, an2, oxygen)
                if key not in best or d < best[key]:
                    best[key] = d
        counts = {"same_subunit": 0, "partner_retained": 0, "partner_lost": 0}
        for (ch2, rn2, ri2, an2, oxygen), d in sorted(best.items(), key=lambda x: x[1]):
            if ch2 == host_chain:
                where = "same_subunit"
            elif ri2 <= retained_end:
                where = "partner_retained"
            elif lost_from <= ri2 <= lost_to:
                where = "partner_lost"
            else:
                where = "partner_retained"
            counts[where] += 1
            rows.append(dict(protein=protein, entry=entry, cofactor=cof_name,
                             phosphate_oxygen=oxygen, donor_chain=ch2,
                             donor_residue=rn2, donor_number=ri2, donor_atom=an2,
                             distance_A=round(d, 2), source=where))
        summary["%s %s" % (protein, entry)] = dict(
            cofactor=cof_name, host_chain=host_chain, total=len(best), **counts)
        print("%-5s %-5s %s%-4d  %d hydrogen bonds to the phosphate: %d from the same "
              "subunit, %d from the partner's retained block, %d from the span the "
              "isoform lacks"
              % (protein, entry, cof_name, cof_num, len(best), counts["same_subunit"],
                 counts["partner_retained"], counts["partner_lost"]))
        for r in rows[-len(best):]:
            print("      %-4s %s%-4d %-5s to %-4s  %.2f A   %s"
                  % (r["donor_chain"], r["donor_residue"], r["donor_number"],
                     r["donor_atom"], r["phosphate_oxygen"], r["distance_A"],
                     r["source"].replace("_", " ")))

with open(os.path.join(RESULTS, "phosphate_hbonds.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
json.dump(dict(cutoff_A=HBOND_CUTOFF, donor_elements=list(DONOR_ELEMENTS),
               phosphate_oxygens=list(PHOSPHATE_OXYGENS), summary=summary, bonds=rows),
          open(os.path.join(RESULTS, "phosphate_hbonds.json"), "w"), indent=1)
print("\nwrote results/phosphate_hbonds.tsv and .json")
