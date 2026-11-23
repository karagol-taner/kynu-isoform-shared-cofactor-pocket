import os, json, urllib.request
import sys

import gemmi, freesasa

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structures_common import assembly

SPECS = {
    "KYNU": {
        "accession": "Q16719",
        "length": 465,
        "entries": ("2HZP", "3E9K", "7S3V"),
        "primary": ("2HZP", "3E9K"),
        "retained_end": 301,
        "substituted": (302, 307),
        "deleted": (308, 465),
        "cofactor_lysine": 276,
        "published_cross_subunit_hbond": (303, 305, 318, 333),
        "published_reported_only": (102, 110, 306, 310, 313, 314, 316),
    },
    "SCLY": {
        "accession": "Q96I15",
        "length": 445,
        "entries": ("3GZC", "3GZD"),
        "primary": ("3GZC", "3GZD"),
        "retained_end": 307,
        "substituted": None,
        "deleted": (308, 445),
        "cofactor_lysine": 259,
        "published_cross_subunit_hbond": (),
        "published_reported_only": (),
    },
}
COFACTOR_NAMES = ("PLP", "LLP", "PLR")
BSA_CUTOFF = 1.0
PATCH_LINKAGE = 8.0
RETAINED_AREA_THRESHOLD = 800.0
COFACTOR_CONTACT_CUTOFF = 4.0

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(ROOT, "structures")
RESULTS = os.path.join(ROOT, "results")
WORK = os.path.join(RESULTS, "assemblies")
for d in (STRUCT, RESULTS, WORK):
    os.makedirs(d, exist_ok=True)
freesasa.setVerbosity(freesasa.silent)


def fetch(entry):
    path = os.path.join(STRUCT, entry + ".cif")
    if not os.path.exists(path):
        with urllib.request.urlopen("https://files.rcsb.org/download/%s.cif" % entry, timeout=60) as r:
            open(path, "wb").write(r.read())
    if gemmi.read_structure(path).name.upper() != entry.upper():
        raise SystemExit("entry mismatch in %s" % path)
    return path


def residue_areas(path):
    result = freesasa.calc(freesasa.Structure(path))
    out = {}
    for chain, residues in result.residueAreas().items():
        for number, area in residues.items():
            if number.strip().isdigit():
                out[(chain, int(number))] = area.total
    return out


def write_subset(source_pdb, keep, target):
    st = gemmi.read_structure(source_pdb)
    st.setup_entities()
    for name in [c.name for c in st[0]]:
        if name not in keep:
            st[0].remove_chain(name)
    st.write_pdb(target)
    return sum(1 for line in open(target) if line.startswith("ATOM"))


def buried_against(entry, ref, partner, full_pdb):
    pair = os.path.join(WORK, "%s_%s%s.pdb" % (entry, ref, partner))
    mono = os.path.join(WORK, "%s_%s.pdb" % (entry, ref))
    n_pair = write_subset(full_pdb, {ref, partner}, pair)
    n_mono = write_subset(full_pdb, {ref}, mono)
    if n_pair <= n_mono:
        raise SystemExit("%s pair %s%s has no partner atoms" % (entry, ref, partner))
    bound = residue_areas(pair)
    free = residue_areas(mono)
    bsa = {}
    for (chain, number), area in free.items():
        if chain != ref or (chain, number) not in bound:
            continue
        value = area - bound[(chain, number)]
        if value > BSA_CUTOFF:
            bsa[number] = value
    return bsa, mono


def patches(monomer, numbers):
    st = gemmi.read_structure(monomer)
    st.setup_entities()
    ca = {}
    for residue in st[0][0]:
        atom = residue.find_atom("CA", "*")
        if atom:
            ca[residue.seqid.num] = atom.pos
    members = sorted(n for n in numbers if n in ca)
    parent = {n: n for n in members}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(members):
        for b in members[i + 1:]:
            if ca[a].dist(ca[b]) < PATCH_LINKAGE:
                ra, rb = root(a), root(b)
                if ra != rb:
                    parent[ra] = rb
    groups = {}
    for n in members:
        groups.setdefault(root(n), []).append(n)
    return sorted(groups.values(), key=len, reverse=True)


def cofactor_contacts(entry, ref, partner, spec):
    st = assembly(entry, keep_ligands=True)
    chains = {c.name: c for c in st[0]}
    a, b = chains[ref], chains[partner]
    cofactor = [atom.pos for res in b if res.name in COFACTOR_NAMES for atom in res]
    names = sorted({res.name for res in b if res.name in COFACTOR_NAMES})
    lysine = spec["cofactor_lysine"]
    if not cofactor:
        res = next((r for r in b if r.seqid.num == lysine), None)
        raise SystemExit("%s partner %s has no cofactor; residue %d is %s"
                         % (entry, partner, lysine, res.name if res else "absent"))
    out = {"cofactor_names": names, "contacts": {}}
    for res in a:
        if res.het_flag != "A" and res.name not in COFACTOR_NAMES:
            continue
        d = min((at.pos.dist(p) for at in res for p in cofactor), default=None)
        if d is not None and d < COFACTOR_CONTACT_CUTOFF:
            out["contacts"][res.seqid.num] = {"residue": res.name, "distance": d}
    return out


def coverage(entry, length):
    st = gemmi.read_structure(os.path.join(STRUCT, entry + ".cif"))
    st.setup_entities()
    out = {}
    for chain in st[0]:
        nums = sorted({res.seqid.num for res in chain.get_polymer()})
        if not nums:
            continue
        segs, start, prev = [], nums[0], nums[0]
        for n in nums[1:]:
            if n == prev + 1:
                prev = n
                continue
            segs.append((start, prev))
            start = prev = n
        segs.append((start, prev))
        out[chain.name] = {"observed": len(nums), "fraction": len(nums) / length,
                           "segments": ";".join("%d-%d" % s for s in segs)}
    return out


def block_of(number, spec):
    if number <= spec["retained_end"]:
        return "retained"
    if spec["substituted"] and spec["substituted"][0] <= number <= spec["substituted"][1]:
        return "substituted"
    return "deleted"


report, rows, lines = {}, [], []
for protein, spec in SPECS.items():
    report[protein] = {}
    for entry in spec["entries"]:
        fetch(entry)
        st = assembly(entry, keep_ligands=False)
        full = os.path.join(WORK, entry + "_assembly.pdb")
        st.write_pdb(full)
        names = [c.name for c in st[0]]
        ref = names[0]
        per_partner = {}
        for partner in names[1:]:
            bsa, mono = buried_against(entry, ref, partner, full)
            per_partner[partner] = sum(bsa.values())
        best = max(per_partner, key=per_partner.get)
        bsa, mono = buried_against(entry, ref, best, full)
        total = sum(bsa.values())
        blocks = {"retained": 0.0, "substituted": 0.0, "deleted": 0.0}
        for n, v in bsa.items():
            blocks[block_of(n, spec)] += v
        groups = patches(mono, [n for n in bsa if block_of(n, spec) == "retained"])
        cof = cofactor_contacts(entry, ref, best, spec)
        lost_cof = {n: v for n, v in cof["contacts"].items() if block_of(n, spec) != "retained"}
        kept_cof = {n: v for n, v in cof["contacts"].items() if block_of(n, spec) == "retained"}
        r = {
            "chains_in_assembly": names, "reference": ref, "partner": best,
            "bsa_per_partner": {k: round(v, 1) for k, v in per_partner.items()},
            "coverage": coverage(entry, spec["length"]),
            "interface_residues": len(bsa), "total_bsa": total,
            "bsa_retained": blocks["retained"], "bsa_substituted": blocks["substituted"],
            "bsa_deleted": blocks["deleted"],
            "fraction_retained": blocks["retained"] / total,
            "fraction_substituted": blocks["substituted"] / total,
            "fraction_deleted": blocks["deleted"] / total,
            "f_removed": (blocks["substituted"] + blocks["deleted"]) / total,
            "retained_patches": len(groups),
            "largest_patch_residues": len(groups[0]) if groups else 0,
            "largest_patch_bsa": sum(bsa[n] for n in groups[0]) if groups else 0.0,
            "largest_patch_span": [min(groups[0]), max(groups[0])] if groups else None,
            "cofactor_names": cof["cofactor_names"],
            "cross_subunit_cofactor_contacts_lost": {k: (v["residue"], round(v["distance"], 2)) for k, v in sorted(lost_cof.items())},
            "cross_subunit_cofactor_contacts_kept": {k: (v["residue"], round(v["distance"], 2)) for k, v in sorted(kept_cof.items())},
            "published_hbond_bsa": {n: round(bsa.get(n, 0.0), 1) for n in spec["published_cross_subunit_hbond"]},
            "published_reported_only_bsa": {n: round(bsa.get(n, 0.0), 1) for n in spec["published_reported_only"]},
        }
        report[protein][entry] = r
        for n in sorted(bsa):
            rows.append((protein, entry, n, round(bsa[n], 2), block_of(n, spec)))
        lines.append("%s %s  chains %s  reference %s  partner %s  BSA per partner %s"
                     % (protein, entry, names, ref, best, r["bsa_per_partner"]))
        for chain, cov in sorted(r["coverage"].items()):
            lines.append("    coverage chain %s: %d of %d (%.3f), segments %s"
                         % (chain, cov["observed"], spec["length"], cov["fraction"], cov["segments"]))
        lines.append("    interface residues %d  total BSA %.0f A2" % (r["interface_residues"], total))
        lines.append("    retained 1-%d %.0f A2 %.3f | substituted %s %.0f A2 %.3f | replaced-or-deleted %d-%d %.0f A2 %.3f"
                     % (spec["retained_end"], r["bsa_retained"], r["fraction_retained"],
                        spec["substituted"], r["bsa_substituted"], r["fraction_substituted"],
                        spec["deleted"][0], spec["deleted"][1], r["bsa_deleted"], r["fraction_deleted"]))
        lines.append("    retained patches %d  largest %d residues %.0f A2 span %s"
                     % (r["retained_patches"], r["largest_patch_residues"], r["largest_patch_bsa"], r["largest_patch_span"]))
        lines.append("    partner cofactor %s  cross-subunit contacts under %.1f A lost: %s  kept: %s"
                     % (r["cofactor_names"], COFACTOR_CONTACT_CUTOFF,
                        r["cross_subunit_cofactor_contacts_lost"] or "none",
                        r["cross_subunit_cofactor_contacts_kept"] or "none"))
        if spec["published_cross_subunit_hbond"]:
            lines.append("    published cross-subunit hydrogen-bond residues, BSA: %s" % r["published_hbond_bsa"])
            lines.append("    residues reported without an explicit contact, BSA: %s" % r["published_reported_only_bsa"])

gates = {}
for protein, spec in SPECS.items():
    p = [report[protein][e] for e in spec["primary"]]
    spread = max(x["f_removed"] for x in p) - min(x["f_removed"] for x in p)
    gates[protein] = {
        "A_published_hbond_recovered": (all(x["published_hbond_bsa"][n] > BSA_CUTOFF
                                            for x in p for n in spec["published_cross_subunit_hbond"])
                                        if spec["published_cross_subunit_hbond"] else None),
        "B_informative": all(0.10 < x["f_removed"] < 0.90 for x in p),
        "C_retained_area": all(x["bsa_retained"] >= RETAINED_AREA_THRESHOLD and x["retained_patches"] == 1 for x in p),
        "C_retained_bsa": [round(x["bsa_retained"]) for x in p],
        "D_spread_points": spread * 100,
        "D_pass": spread * 100 <= 5.0,
        "cofactor_contacts_lost": [len(x["cross_subunit_cofactor_contacts_lost"]) for x in p],
    }

with open(os.path.join(RESULTS, "interface_residue_bsa.tsv"), "w") as fh:
    fh.write("protein\tentry\tresidue\tbsa_angstrom2\tblock\n")
    for row in rows:
        fh.write("\t".join(str(v) for v in row) + "\n")
json.dump({"entries": report, "gates": gates}, open(os.path.join(RESULTS, "interface_partition.json"), "w"),
          indent=1, default=str)
lines.append("")
for protein, g in gates.items():
    lines.append("%s gates: %s" % (protein, json.dumps(g)))
open(os.path.join(RESULTS, "gate_report.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
