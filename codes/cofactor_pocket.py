import csv, json, os
import gemmi
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

SPECS = {
    "KYNU": {"entries": ("2HZP", "3E9K"), "retained_end": 301, "deleted": (302, 465), "cofactor_lysine": 276},
    "SCLY": {"entries": ("3GZC", "3GZD"), "retained_end": 307, "deleted": (308, 445), "cofactor_lysine": 259},
}
COFACTOR_NAMES = ("PLP", "LLP", "PLR")
PHOSPHATE_ATOMS = ("P", "O1P", "O2P", "O3P", "OP1", "OP2", "OP3")
PROBE_RADIUS = 1.40
N_POINTS = 960
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(ROOT, "structures")
RESULTS = os.path.join(ROOT, "results")
WORK = os.path.join(RESULTS, "pocket_work")
os.makedirs(WORK, exist_ok=True)


def assembly(entry, retained_end):
    st = gemmi.read_structure(os.path.join(STRUCT, entry + ".cif"))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.transform_to_assembly("1", gemmi.HowToNameCopiedChain.AddNumber)
    st.remove_waters()
    model = st[0]
    if len(model) < 2:
        raise SystemExit("%s assembly 1 has %d chains" % (entry, len(model)))
    for i, chain in enumerate(model):
        chain.name = "ABCDEFGH"[i]
    st.setup_entities()
    return st


def write(st, keep_chains, partner_limit, path, keep_partner_cofactor=False, partner="B"):
    s = gemmi.Structure()
    s.add_model(gemmi.Model("1"))
    for chain in st[0]:
        if chain.name not in keep_chains:
            continue
        out = gemmi.Chain(chain.name)
        for res in chain:
            if chain.name == partner and partner_limit is not None:
                if res.name in COFACTOR_NAMES:
                    if not keep_partner_cofactor:
                        continue
                elif res.seqid.num > partner_limit:
                    continue
                elif res.het_flag == "H":
                    continue
            out.add_residue(res)
        s[0].add_chain(out)
    s.setup_entities()
    s.write_pdb(path)
    return sum(1 for line in open(path) if line.startswith(("ATOM", "HETATM")))


def cofactor_sasa(path):
    structure = PDBParser(QUIET=True).get_structure("x", path)
    ShrakeRupley(probe_radius=PROBE_RADIUS, n_points=N_POINTS).compute(structure[0], level="A")
    total = {}
    for chain in structure[0]:
        for res in chain:
            if res.get_resname().strip() in COFACTOR_NAMES:
                per_atom = {a.get_id(): round(a.sasa, 2) for a in res}
                total[(chain.id, res.id[1], res.get_resname().strip())] = (sum(a.sasa for a in res), per_atom)
    return total


def partner_of(st, reference="A"):
    best, best_n = None, -1
    ref = [a.pos for ch in st[0] if ch.name == reference for r in ch for a in r]
    for ch in st[0]:
        if ch.name == reference:
            continue
        n = sum(1 for r in ch for a in r if any(a.pos.dist(p) < 5.0 for p in ref[::40]))
        if n > best_n:
            best, best_n = ch.name, n
    if best is None:
        raise SystemExit("no partner chain found")
    return best


rows, summary = [], {}
for protein, spec in SPECS.items():
  RETAINED_END = spec["retained_end"]
  for entry in spec["entries"]:
    st = assembly(entry, RETAINED_END)
    partner = partner_of(st)
    present = sorted({res.name for chain in st[0] for res in chain if res.name in COFACTOR_NAMES})
    if not present:
        raise SystemExit("%s carries no cofactor" % entry)
    paths = {
        "monomer_A": write(st, {"A"}, None, os.path.join(WORK, "%s_A_cofactor.pdb" % entry), partner=partner),
        "dimer_full": write(st, {"A", partner}, None, os.path.join(WORK, "%s_AB_cofactor.pdb" % entry), partner=partner),
        "dimer_isoform_mimic": write(st, {"A", partner}, RETAINED_END,
                                     os.path.join(WORK, "%s_AB_iso_mimic.pdb" % entry), partner=partner),
        "dimer_isoform_mimic_with_partner_cofactor": write(
            st, {"A", partner}, RETAINED_END,
            os.path.join(WORK, "%s_AB_iso_mimic_cof.pdb" % entry), keep_partner_cofactor=True, partner=partner),
    }
    sasa = {}
    for label in paths:
        name = {"monomer_A": "%s_A_cofactor.pdb", "dimer_full": "%s_AB_cofactor.pdb",
                "dimer_isoform_mimic": "%s_AB_iso_mimic.pdb",
                "dimer_isoform_mimic_with_partner_cofactor": "%s_AB_iso_mimic_cof.pdb"}[label] % entry
        sasa[label] = cofactor_sasa(os.path.join(WORK, name))
    key_a = [k for k in sasa["monomer_A"] if k[0] == "A"]
    if len(key_a) != 1:
        raise SystemExit("%s: expected one cofactor on chain A, found %d" % (entry, len(key_a)))
    k = key_a[0]
    alone = sasa["monomer_A"][k][0]
    full = sasa["dimer_full"][k][0]
    mimic = sasa["dimer_isoform_mimic"][k][0]
    buried_by_partner = alone - full
    buried_by_retained = alone - mimic
    buried_by_lost = mimic - full
    control = sasa["dimer_isoform_mimic_with_partner_cofactor"][k][0]
    if abs(control - mimic) > 0.01:
        raise SystemExit("%s: the partner cofactor itself buries %.2f A2, so burial cannot be attributed to "
                         "protein residues alone" % (entry, mimic - control))
    atoms_f = sasa["dimer_full"][k][1]
    atoms_m = sasa["dimer_isoform_mimic"][k][1]
    phosphate = sum(atoms_m.get(a, 0) - atoms_f.get(a, 0) for a in atoms_f if a in PHOSPHATE_ATOMS)
    if buried_by_partner <= 0:
        raise SystemExit("%s: the partner buries no cofactor surface, so the comparison is void" % entry)
    summary["%s %s" % (protein, entry)] = {
        "protein": protein, "partner_chain": partner, "retained_end": RETAINED_END,
        "retained_block": "1-%d" % RETAINED_END,
        "lost_block": "%d-%d" % (RETAINED_END + 1, spec["deleted"][1]),
        "isoform_mimic": "partner chain truncated to residues 1 to %d, its cofactor removed" % RETAINED_END,
        "cofactor": k[2], "cofactor_residue": k[1], "atoms_in_models": paths,
        "cofactor_sasa_monomer": round(alone, 2),
        "cofactor_sasa_dimer": round(full, 2),
        "cofactor_sasa_isoform_mimic": round(mimic, 2),
        "buried_by_whole_partner": round(buried_by_partner, 2),
        "buried_by_retained_block": round(buried_by_retained, 2),
        "buried_by_lost_block": round(buried_by_lost, 2),
        "fraction_of_partner_burial_lost": round(buried_by_lost / buried_by_partner, 4),
        "partner_cofactor_control_sasa": round(control, 2),
        "phosphate_group_burial_lost": round(phosphate, 2),
        "phosphate_share_of_lost_burial": round(phosphate / buried_by_lost, 4) if buried_by_lost else None,
    }
    atoms_full = sasa["dimer_full"][k][1]
    atoms_mimic = sasa["dimer_isoform_mimic"][k][1]
    for atom in sorted(atoms_full):
        rows.append({"protein": protein, "entry": entry, "cofactor": k[2], "atom": atom,
                     "sasa_monomer": sasa["monomer_A"][k][1].get(atom, ""),
                     "sasa_dimer": atoms_full[atom],
                     "sasa_isoform_mimic": atoms_mimic.get(atom, ""),
                     "exposed_when_partner_truncated": round(atoms_mimic.get(atom, 0) - atoms_full[atom], 2)})

with open(os.path.join(RESULTS, "cofactor_pocket_atoms.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
json.dump({"probe_radius": PROBE_RADIUS, "n_points": N_POINTS,
           "isoform_mimic": "per entry, see each entry's own isoform_mimic field",
           "entries": summary}, open(os.path.join(RESULTS, "cofactor_pocket.json"), "w"), indent=1)

for entry, s in summary.items():
    print("%s  cofactor %s%d  partner chain %s  retained block 1-%d"
          % (entry, s["cofactor"], s["cofactor_residue"], s["partner_chain"], s["retained_end"]))
    print("   cofactor solvent accessibility, chain A alone          %8.2f A2" % s["cofactor_sasa_monomer"])
    print("   with the whole partner present                         %8.2f A2" % s["cofactor_sasa_dimer"])
    print("   with the partner truncated to the isoform              %8.2f A2" % s["cofactor_sasa_isoform_mimic"])
    print("   buried by the whole partner                            %8.2f A2" % s["buried_by_whole_partner"])
    print("   buried by the partner's retained block                 %8.2f A2" % s["buried_by_retained_block"])
    print("   buried by the span the isoform lacks                   %8.2f A2  (%.1f per cent of the partner's contribution)"
          % (s["buried_by_lost_block"], 100 * s["fraction_of_partner_burial_lost"]))
    if s["phosphate_share_of_lost_burial"] is None:
        print("   no burial is lost, so no phosphate share applies")
    else:
        print("   of that, the cofactor phosphate group accounts for      %8.2f A2  (%.1f per cent)"
              % (s["phosphate_group_burial_lost"], 100 * s["phosphate_share_of_lost_burial"]))
    print("   control, partner cofactor retained in the mimic         %8.2f A2, unchanged" % s["partner_cofactor_control_sasa"])
    print()
