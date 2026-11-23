import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "codes", "cofactor_pocket.py")
MARK = "rows, summary = [], {}"
text = open(SRC).read()
if MARK not in text:
    raise SystemExit("STOP: %s no longer begins its run at %r, so the definitions "
                     "this script borrows cannot be separated from the run" % (SRC, MARK))
G = {"__name__": "cofactor_pocket_defs", "__file__": SRC}
exec(compile(text[:text.index(MARK)], SRC, "exec"), G)
SPECS, COFACTOR_NAMES = G["SPECS"], G["COFACTOR_NAMES"]
WORK, RESULTS = G["WORK"], G["RESULTS"]
assembly, write, cofactor_sasa, partner_of = (G["assembly"], G["write"],
                                              G["cofactor_sasa"], G["partner_of"])

rows, summary = [], {}
for protein, spec in SPECS.items():
    end = spec["retained_end"]
    for entry in spec["entries"]:
        st = assembly(entry, end)
        partner = partner_of(st)
        states = {
            "own_whole_alone": dict(keep={"A"}, limit=None, who=partner),
            "own_truncated_alone": dict(keep={"A"}, limit=end, who="A"),
            "dimer_full": dict(keep={"A", partner}, limit=None, who=partner),
            "canonical_site_in_mixed": dict(keep={"A", partner}, limit=end, who=partner),
            "isoform_site_in_mixed": dict(keep={"A", partner}, limit=end, who="A"),
        }
        sasa = {}
        for label, s in states.items():
            path = os.path.join(WORK, "%s_own_%s.pdb" % (entry, label))
            write(st, s["keep"], s["limit"], path, keep_partner_cofactor=True,
                  partner=s["who"])
            got = cofactor_sasa(path)
            on_a = [k for k in got if k[0] == "A"]
            if len(on_a) != 1:
                raise SystemExit("STOP: %s %s puts %d cofactors on chain A"
                                 % (entry, label, len(on_a)))
            sasa[label] = got[on_a[0]][0]
        own_lost = sasa["own_truncated_alone"] - sasa["own_whole_alone"]
        if own_lost < -0.01:
            raise SystemExit("STOP: %s buries %.2f A2 more of its cofactor with part "
                             "of its own chain removed, which cannot be"
                             % (entry, -own_lost))
        rec = dict(protein=protein, entry=entry, partner_chain=partner,
                   retained_end=end,
                   own_whole_alone=round(sasa["own_whole_alone"], 2),
                   own_truncated_alone=round(sasa["own_truncated_alone"], 2),
                   buried_by_own_lost_block=round(own_lost, 2),
                   dimer_full=round(sasa["dimer_full"], 2),
                   canonical_site_in_mixed=round(sasa["canonical_site_in_mixed"], 2),
                   isoform_site_in_mixed=round(sasa["isoform_site_in_mixed"], 2),
                   canonical_site_excess=round(sasa["canonical_site_in_mixed"]
                                               - sasa["dimer_full"], 2),
                   isoform_site_excess=round(sasa["isoform_site_in_mixed"]
                                             - sasa["dimer_full"], 2))
        rows.append(rec)
        summary["%s %s" % (protein, entry)] = rec
        print("%-5s %-5s  own chain beyond %d buries %6.2f A2 of its own cofactor"
              % (protein, entry, end, own_lost))
        print("        canonical site %6.2f A2 in the full dimer, %6.2f in a mixed "
              "dimer  (+%.2f)"
              % (rec["dimer_full"], rec["canonical_site_in_mixed"],
                 rec["canonical_site_excess"]))
        print("        isoform's own site %6.2f A2 in a mixed dimer              "
              "  (+%.2f)" % (rec["isoform_site_in_mixed"], rec["isoform_site_excess"]))

with open(os.path.join(RESULTS, "own_subunit_pocket.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
json.dump(dict(note="The isoform's own chain is modelled as the canonical chain cut "
                    "at the retained end. Its six novel residues are not in any "
                    "structure, so what it supplies to its own site is a lower bound.",
               entries=summary),
          open(os.path.join(RESULTS, "own_subunit_pocket.json"), "w"), indent=1)
print("\nwrote results/own_subunit_pocket.tsv and .json")
