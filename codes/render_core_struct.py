import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_common as cc
from render_md_struct import start_pymol, strip_credentials, WIDTH, HEIGHT, say
from structures_common import assembly

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
OUT_DIR = os.path.join(RESULTS, "renders")

RETAINED = "0xBFBFBF"
SUBSTITUTED = "0xB2182B"
DELETED = "0x4D4D4D"
COFACTOR = "0xE0A33C"
PARTNER_PALE = "0xD8D8D8"
SURFACE_BASE = "0xF7F7F7"

(RET_NAME, RET_LO, RET_HI), (SUB_NAME, SUB_LO, SUB_HI), (DEL_NAME, DEL_LO, DEL_HI) = cc.BLOCKS


def pocket_record(entry):
    path = os.path.join(RESULTS, "cofactor_pocket.json")
    if not os.path.isfile(path):
        raise SystemExit("STOP: no %s, so which residue the cofactor is cannot be "
                         "read. Run codes/cofactor_pocket.py first." % path)
    rec = json.load(open(path))["entries"]
    hit = [v for k, v in rec.items() if k.endswith(entry)]
    if len(hit) != 1:
        raise SystemExit("STOP: %d records for %s in cofactor_pocket.json."
                         % (len(hit), entry))
    return hit[0]


def entry_path(entry):
    out = os.path.join(RESULTS, "assemblies", "%s_AB_ligands.pdb" % entry)
    if not os.path.isfile(out):
        st = assembly(entry, keep_ligands=True)
        st.remove_waters()
        st.write_pdb(out)
    return out


def colour_blocks(c, selection):
    c.color(RETAINED, "%s and resi %d-%d" % (selection, RET_LO, RET_HI))
    c.color(SUBSTITUTED, "%s and resi %d-%d" % (selection, SUB_LO, SUB_HI))
    c.color(DELETED, "%s and resi %d-%d" % (selection, DEL_LO, DEL_HI))


def dimer_blocks(p, out, entry="2HZP"):
    c = p.cmd
    c.delete("all")
    c.load(entry_path(entry), "dim")
    c.hide("everything")
    c.show("cartoon", "polymer")
    colour_blocks(c, "polymer")
    cof = "resn %s" % pocket_record(entry)["cofactor"]
    if c.count_atoms(cof) == 0:
        raise SystemExit("STOP: the assembly of %s holds no %s, so the cofactor "
                         "cannot be drawn." % (entry, cof))
    c.show("spheres", "%s and not elem H" % cof)
    c.color(COFACTOR, cof)
    c.color("atomic", "%s and not elem C" % cof)
    c.orient("polymer")
    c.turn("x", -15)
    c.zoom("polymer", 2)
    c.ray(WIDTH, HEIGHT)
    c.png(out, dpi=400)
    return os.path.isfile(out)


def face_towards(c, target, whole, orient_on=None):
    c.orient(orient_on or target)
    m = c.get_view()[:9]
    a, b = c.centerofmass(target), c.centerofmass(whole)
    d = [a[i] - b[i] for i in range(3)]
    z = m[2] * d[0] + m[5] * d[1] + m[8] * d[2]
    if z < 0:
        c.turn("y", 180)


def interface_patch(p, out, entry="2HZP"):
    rows = [r for r in csv.DictReader(
        open(os.path.join(RESULTS, "interface_residue_bsa.tsv")), delimiter="\t")
        if r["protein"] == "KYNU" and r["entry"] == entry]
    if not rows:
        raise SystemExit("STOP: no interface residues recorded for %s in "
                         "interface_residue_bsa.tsv." % entry)
    by_block = {}
    for r in rows:
        by_block.setdefault(r["block"], []).append(int(r["residue"]))
    unknown = sorted(set(by_block) - {RET_NAME, SUB_NAME, DEL_NAME})
    if unknown:
        raise SystemExit("STOP: interface_residue_bsa.tsv names the blocks %s, which "
                         "conservation_common does not define." % unknown)
    c = p.cmd
    c.delete("all")
    c.load(entry_path(entry), "dim")
    c.remove("not chain A")
    c.remove("not polymer")
    c.hide("everything")
    c.show("surface", "polymer")
    c.color(SURFACE_BASE, "polymer")
    parts = []
    for block, colour in ((RET_NAME, RETAINED), (SUB_NAME, SUBSTITUTED),
                          (DEL_NAME, DELETED)):
        resi = sorted(by_block.get(block, []))
        if not resi:
            continue
        sel = "resi %s" % "+".join(str(q) for q in resi)
        c.color(colour, sel)
        parts.append(sel)
    face_towards(c, "(%s)" % ") or (".join(parts), "polymer")
    c.zoom("polymer", 2)
    c.ray(WIDTH, HEIGHT)
    c.png(out, dpi=400)
    return os.path.isfile(out)


def pocket_environment(p, out, entry="2HZP"):
    c = p.cmd
    c.delete("all")
    c.load(entry_path(entry), "dim")
    c.hide("everything")
    name = pocket_record(entry)["cofactor"]
    cof = "resn %s and chain A" % name
    if c.count_atoms(cof) == 0:
        raise SystemExit("STOP: chain A of the %s assembly holds no %s."
                         % (entry, name))
    near = "byres (polymer within 5 of (%s))" % cof
    host, partner = "(%s) and chain A" % near, "(%s) and chain B" % near
    if c.count_atoms(partner) == 0:
        raise SystemExit("STOP: no residue of the partner reaches the %s of %s, so "
                         "the pocket panel has nothing to colour." % (name, entry))
    c.show("cartoon", "byres (polymer within 12 of (%s))" % cof)
    c.set("cartoon_transparency", 0.45)
    c.color("0xE4E4E4", "chain A")
    c.color("0xA8A8A8", "chain B")
    c.show("lines", "(%s) and not elem H and (sidechain or name CA)" % host)
    c.color("0xB4B4B4", host)
    c.show("sticks", "(%s) and not elem H and (sidechain or name CA)" % partner)
    colour_blocks(c, partner)
    c.show("sticks", "%s and not elem H" % cof)
    c.color(COFACTOR, cof)
    c.color("atomic", "(%s or (%s)) and not elem C and not elem H" % (cof, partner))
    c.set("stick_radius", 0.19)
    c.set("line_width", 1.0)
    c.orient(cof)
    c.zoom("%s or (%s)" % (cof, partner), 4)
    c.ray(WIDTH, HEIGHT)
    c.png(out, dpi=400)
    return os.path.isfile(out)


def predicted_heterodimer(p, out):
    scores = json.load(open(os.path.join(RESULTS, "af3_scores.json")))
    job = next((k for k in scores["detail"] if "isoform2" in k and "canonical" in k), None)
    if job is None:
        raise SystemExit("STOP: af3_scores.json carries no canonical with isoform job, "
                         "so the predicted heterodimer cannot be drawn.")
    row = scores["detail"][job]
    name = row["model_file"]
    hits = glob.glob(os.path.join(ROOT, "data", "af3", "**", name), recursive=True)
    if len(hits) != 1:
        raise SystemExit("STOP: %d copies of %s under data/af3, and one is needed."
                         % (len(hits), name))
    c = p.cmd
    c.delete("all")
    c.load(hits[0], "pred")
    c.hide("everything")
    want = dict(zip(row["interface"]["chains"], row["interface"]["chain_residues"]))
    for chain, n in sorted(want.items()):
        got = c.count_atoms("chain %s and polymer and name CA" % chain)
        if got != n:
            raise SystemExit("STOP: chain %s of %s holds %d residues and af3_scores.json "
                             "records %d, so the chains are not the ones scored."
                             % (chain, name, got, n))
    tail = "chain B and resi %d-%d" % (SUB_LO, SUB_HI)
    c.show("cartoon", "polymer")
    c.set("cartoon_loop_radius", 0.35)
    c.color("0xC4C4C4", "chain A")
    c.color(DELETED, "chain B")
    c.color(SUBSTITUTED, tail)
    c.show("sticks", "%s and not elem H" % tail)
    c.set("stick_radius", 0.28)
    c.show("spheres", "not polymer and not solvent and not elem H")
    c.color(COFACTOR, "not polymer and not solvent")
    face_towards(c, tail, "polymer", orient_on="polymer")
    c.turn("x", -15)
    c.zoom("polymer", 1)
    c.ray(WIDTH, HEIGHT)
    c.png(out, dpi=400)
    return os.path.isfile(out)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    p = start_pymol()
    if p is None:
        raise SystemExit("STOP: no pymol2 in this interpreter, so the structural "
                         "panels cannot be drawn.")
    jobs = [("dimer_blocks.png", dimer_blocks),
            ("interface_patch.png", interface_patch),
            ("pocket_environment.png", pocket_environment),
            ("predicted_heterodimer.png", predicted_heterodimer)]
    for name, fn in jobs:
        out = os.path.join(OUT_DIR, name)
        if fn(p, out):
            gone = strip_credentials(out)
            say("wrote %s%s" % (out, " (stripped %s)" % ", ".join(gone) if gone else ""))
        else:
            raise SystemExit("STOP: %s reported success but %s is not there."
                             % (fn.__name__, out))
    p.stop()


if __name__ == "__main__":
    main()
