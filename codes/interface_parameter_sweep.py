import os, csv, itertools
import freesasa
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

ENTRIES = ("2HZP", "3E9K")
RETAINED_END = 301
SUBSTITUTED = (302, 307)
DELETED = (308, 465)
PROBE_RADII = (1.2, 1.4, 1.6)
ALGORITHMS = {"LeeRichards": freesasa.LeeRichards, "ShrakeRupley": freesasa.ShrakeRupley}
BSA_CUTOFFS = (0.0, 1.0, 2.5, 5.0, 10.0)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "assemblies")
RESULTS = os.path.join(ROOT, "results")
freesasa.setVerbosity(freesasa.silent)


def areas(path, algorithm, probe):
    params = freesasa.Parameters({"algorithm": algorithm, "probe-radius": probe})
    result = freesasa.calc(freesasa.Structure(path), params)
    out = {}
    for chain, residues in result.residueAreas().items():
        for number, area in residues.items():
            if number.strip().isdigit():
                out[(chain, int(number))] = area.total
    return out


def biopython_areas(path, probe):
    structure = PDBParser(QUIET=True).get_structure("x", path)
    ShrakeRupley(probe_radius=probe).compute(structure[0], level="R")
    out = {}
    for chain in structure[0]:
        for residue in chain:
            out[(chain.id, residue.id[1])] = residue.sasa
    return out


def partition(bsa):
    total = sum(bsa.values())
    retained = sum(v for n, v in bsa.items() if n <= RETAINED_END)
    substituted = sum(v for n, v in bsa.items() if SUBSTITUTED[0] <= n <= SUBSTITUTED[1])
    deleted = sum(v for n, v in bsa.items() if DELETED[0] <= n <= DELETED[1])
    return total, retained, substituted, deleted


rows = []
for entry in ENTRIES:
    pair = os.path.join(WORK, entry + "_AB.pdb")
    mono = os.path.join(WORK, entry + "_A.pdb")
    for name, algorithm in ALGORITHMS.items():
        for probe in PROBE_RADII:
            bound = areas(pair, algorithm, probe)
            free = areas(mono, algorithm, probe)
            raw = {}
            for (chain, number), area in free.items():
                if chain != "A" or (chain, number) not in bound:
                    continue
                raw[number] = area - bound[(chain, number)]
            for cutoff in BSA_CUTOFFS:
                bsa = {n: v for n, v in raw.items() if v > cutoff}
                total, retained, substituted, deleted = partition(bsa)
                rows.append({
                    "entry": entry, "implementation": "freesasa", "algorithm": name,
                    "probe_radius": probe, "bsa_cutoff": cutoff,
                    "interface_residues": len(bsa),
                    "total_bsa": round(total, 1), "retained_bsa": round(retained, 1),
                    "f_retained": round(retained / total, 4),
                    "f_substituted": round(substituted / total, 4),
                    "f_removed": round((substituted + deleted) / total, 4),
                })
    for probe in PROBE_RADII:
        bound = biopython_areas(pair, probe)
        free = biopython_areas(mono, probe)
        raw = {}
        for (chain, number), area in free.items():
            if chain != "A" or (chain, number) not in bound:
                continue
            raw[number] = area - bound[(chain, number)]
        for cutoff in BSA_CUTOFFS:
            bsa = {n: v for n, v in raw.items() if v > cutoff}
            total, retained, substituted, deleted = partition(bsa)
            rows.append({
                "entry": entry, "implementation": "biopython", "algorithm": "ShrakeRupley",
                "probe_radius": probe, "bsa_cutoff": cutoff,
                "interface_residues": len(bsa),
                "total_bsa": round(total, 1), "retained_bsa": round(retained, 1),
                "f_retained": round(retained / total, 4),
                "f_substituted": round(substituted / total, 4),
                "f_removed": round((substituted + deleted) / total, 4),
            })

target = os.path.join(RESULTS, "interface_parameter_sweep.tsv")
with open(target, "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

f = [r["f_removed"] for r in rows]
print("combinations: %d" % len(rows))
print("f_removed over all settings: min %.4f  max %.4f  spread %.1f points" % (min(f), max(f), (max(f) - min(f)) * 100))
for impl in ("freesasa", "biopython"):
    sub = [r["f_removed"] for r in rows if r["implementation"] == impl]
    print("  %-10s min %.4f max %.4f  n=%d" % (impl, min(sub), max(sub), len(sub)))
for entry in ENTRIES:
    sub = [r["f_removed"] for r in rows if r["entry"] == entry]
    print("  %-6s min %.4f max %.4f" % (entry, min(sub), max(sub)))
ret = [r["retained_bsa"] for r in rows]
print("retained BSA over all settings: min %.0f  max %.0f  (threshold 800)" % (min(ret), max(ret)))
print("wrote %s" % target)
