import csv, glob, json, math, os, re, sys

import gemmi
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
CONTACT_A = 8.0
PAE_D0 = 10.0

_conf = open(os.path.join(ROOT, "codes", "interface_partition.py")).read()
import ast as _ast
_SPECS = next(_ast.literal_eval(n.value) for n in _ast.parse(_conf).body
              if isinstance(n, _ast.Assign) and getattr(n.targets[0], "id", None) == "SPECS")
RETAINED_END = _SPECS["KYNU"]["retained_end"]
SUBSTITUTED_END = _SPECS["KYNU"]["substituted"][1]
CANONICAL_END = _SPECS["KYNU"]["deleted"][1]

JOB_FOLDER = {
    "canonical__canonical": "canonical_canonical",
    "isoform2__isoform2": "isoform2_isoform2",
    "canonical__isoform2": "canonical_isoform2",
    "canonical__trunc_cterm": "canonical_trunc_cterm",
    "canonical__trunc_nterm": "canonical_trunc_nterm",
    "canonical__tail_scramble_01": "canonical_tail_scramble_01",
    "canonical__tail_scramble_02": "canonical_tail_scramble_02",
    "canonical__tail_scramble_03": "canonical_tail_scramble_03",
}

def fatal(msg):
    raise SystemExit("STOP: " + msg)

def find_models_root(given):
    if given:
        if not os.path.isdir(given):
            fatal("%s is not a folder" % given)
        return given
    here = os.path.join(DATA, "af3")
    hits = set()
    for dirpath, dirnames, filenames in os.walk(here):
        if any(f.endswith("_model_0.cif") for f in filenames):
            hits.add(os.path.dirname(dirpath))
    hits = sorted(hits)
    if not hits:
        fatal("no folder under %s holds a file ending in _model_0.cif. Unpack the "
              "AlphaFold3 results there and run this again." % here)
    if len(hits) > 1:
        fatal("%d folders under %s hold model files, %s. Which run to score is not "
              "guessed; pass one as an argument." % (len(hits), here, hits))
    return hits[0]

def top_model(folder):
    best = None
    for name in sorted(os.listdir(folder)):
        m = re.match(r".*_summary_confidences_(\d+)\.json$", name)
        if not m:
            continue
        d = json.load(open(os.path.join(folder, name)))
        if "ranking_score" not in d:
            fatal("%s carries no ranking_score" % os.path.join(folder, name))
        if best is None or d["ranking_score"] > best[1]["ranking_score"]:
            best = (int(m.group(1)), d)
    if best is None:
        fatal("%s holds no summary_confidences file" % folder)
    k, summary = best
    cif = [f for f in os.listdir(folder) if f.endswith("_model_%d.cif" % k)]
    pae = [f for f in os.listdir(folder) if f.endswith("_full_data_%d.json" % k)]
    if len(cif) != 1 or len(pae) != 1:
        fatal("%s holds %d model and %d data files for sample %d, and one of each is "
              "needed" % (folder, len(cif), len(pae), k))
    return k, summary, os.path.join(folder, cif[0]), os.path.join(folder, pae[0])

def residue_frame(cif_path, data):
    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    atoms = [(ch.name, r.seqid.num, r.name, a.name, a.pos)
             for ch in st[0] for r in ch for a in r]
    plddt = data["atom_plddts"]
    if len(atoms) != len(plddt):
        fatal("%s holds %d atoms and its data file carries %d pLDDT values"
              % (os.path.basename(cif_path), len(atoms), len(plddt)))
    rows = {}
    for (chain, num, resname, aname, pos), p in zip(atoms, plddt):
        want = "CA" if resname == "GLY" else "CB"
        if aname == want:
            rows[(chain, num)] = dict(chain=chain, num=num, resname=resname,
                                      pos=pos, plddt=float(p))
    if not rows:
        fatal("%s yielded no representative atoms" % cif_path)
    return rows

def pae_index(data):
    chains, nums = data["token_chain_ids"], data["token_res_ids"]
    if len(chains) != len(nums):
        fatal("the data file carries %d token chains and %d token residues"
              % (len(chains), len(nums)))
    index = {}
    for i, (c, n) in enumerate(zip(chains, nums)):
        index.setdefault((c, int(n)), i)
    return index

def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line and name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

CONSTRUCTS = read_fasta(os.path.join(DATA, "constructs.fasta"))
CANONICAL = CONSTRUCTS["canonical"]
def canonical_offset(seq, name):
    probe = seq[:50]
    at = CANONICAL.find(probe)
    if at < 0:
        fatal("the first 50 residues of %s are not in the canonical sequence, so its "
              "numbering cannot be placed" % name)
    if CANONICAL.find(probe, at + 1) >= 0:
        fatal("the first 50 residues of %s occur more than once in the canonical "
              "sequence, so its offset is ambiguous" % name)
    return at

def chain_sequence(st, chain_name):
    out = []
    for ch in st[0]:
        if ch.name != chain_name:
            continue
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info is not None and info.is_amino_acid():
                out.append(gemmi.find_tabulated_residue(r.name).one_letter_code.upper())
    return "".join(out)

def score(cif_path, pae_path):
    data = json.load(open(pae_path))
    res = residue_frame(cif_path, data)
    idx = pae_index(data)
    pae = np.asarray(data["pae"], dtype=float)
    keys = sorted(res)
    chains = sorted({k[0] for k in keys})
    if len(chains) != 2:
        fatal("%s holds %d chains and this scores one pair"
              % (os.path.basename(cif_path), len(chains)))
    a = [k for k in keys if k[0] == chains[0]]
    b = [k for k in keys if k[0] == chains[1]]
    pa = np.array([[res[k]["pos"].x, res[k]["pos"].y, res[k]["pos"].z] for k in a])
    pb = np.array([[res[k]["pos"].x, res[k]["pos"].y, res[k]["pos"].z] for k in b])
    d = np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=2)
    close = d <= CONTACT_A
    npairs = int(close.sum())
    if npairs == 0:
        return dict(interface_pairs=0, interface_residues=0, interface_plddt=0.0,
                    pDockQ=0.0, pDockQ2=0.0, chain_residues=[len(a), len(b)],
                    chains=chains)
    iface = ([a[i] for i in np.where(close.any(axis=1))[0]]
             + [b[j] for j in np.where(close.any(axis=0))[0]])
    mean_plddt = float(np.mean([res[k]["plddt"] for k in iface]))
    x1 = mean_plddt * math.log10(npairs)
    pdockq = 0.724 / (1 + math.exp(-0.052 * (x1 - 152.611))) + 0.018
    missing = [k for k in a + b if k not in idx]
    if missing:
        fatal("%d residues of %s are absent from the predicted aligned error index, "
              "first %s" % (len(missing), os.path.basename(cif_path), missing[:3]))
    ia = np.array([idx[k] for k in a])
    ib = np.array([idx[k] for k in b])
    sub = pae[np.ix_(ia, ib)]
    ptm_like = 1.0 / (1.0 + (sub / PAE_D0) ** 2)
    mean_ptm = float(ptm_like[close].mean())
    x2 = mean_plddt * mean_ptm
    pdockq2 = 1.31 / (1 + math.exp(-0.075 * (x2 - 84.733))) + 0.005
    iface_a = [a[i][1] for i in np.where(close.any(axis=1))[0]]
    iface_b = [b[j][1] for j in np.where(close.any(axis=0))[0]]
    return dict(interface_pairs=npairs, interface_residues=len(iface),
                interface_plddt=round(mean_plddt, 2),
                interface_pae=round(float(sub[close].mean()), 2),
                pDockQ=round(pdockq, 4), pDockQ2=round(pdockq2, 4),
                chain_residues=[len(a), len(b)], chains=chains,
                interface_numbers={chains[0]: sorted(iface_a),
                                   chains[1]: sorted(iface_b)},
                interface_pae_pairs=[round(float(v), 2) for v in sub[close]])

MODELS = find_models_root(sys.argv[1] if len(sys.argv) > 1 else None)
print("models root %s" % MODELS)
print()
ROLE = {}
with open(os.path.join(DATA, "pairings.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        ROLE[r["job"]] = r["role"]

rows, detail = [], {}
for job, folder_name in JOB_FOLDER.items():
    folder = os.path.join(MODELS, folder_name)
    if not os.path.isdir(folder):
        fatal("no folder %s under %s for job %s" % (folder_name, MODELS, job))
    k, summary, cif, pae = top_model(folder)
    s = score(cif, pae)
    st = gemmi.read_structure(cif)
    st.setup_entities()
    partner = s["chains"][1]
    offset = canonical_offset(chain_sequence(st, partner), job)
    blocks = dict(retained=0, substituted=0, deleted=0, outside=0)
    for n in s["interface_numbers"][partner]:
        c = n + offset
        if c <= RETAINED_END:
            blocks["retained"] += 1
        elif c <= SUBSTITUTED_END:
            blocks["substituted"] += 1
        elif c <= CANONICAL_END:
            blocks["deleted"] += 1
        else:
            blocks["outside"] += 1
    s["partner_offset_in_canonical"] = offset
    s["partner_interface_by_block"] = blocks
    row = dict(job=job, role=ROLE.get(job, ""), model=k,
               chain_a_residues=s["chain_residues"][0],
               chain_b_residues=s["chain_residues"][1],
               iptm=summary["iptm"], ptm=summary["ptm"],
               ranking_score=summary["ranking_score"],
               has_clash=summary["has_clash"],
               fraction_disordered=summary["fraction_disordered"],
               interface_residues=s["interface_residues"],
               interface_pairs=s["interface_pairs"],
               interface_plddt=s["interface_plddt"],
               interface_pae=s.get("interface_pae"),
               pDockQ=s["pDockQ"], pDockQ2=s["pDockQ2"],
               partner_offset=s["partner_offset_in_canonical"],
               partner_iface_retained=s["partner_interface_by_block"]["retained"],
               partner_iface_substituted=s["partner_interface_by_block"]["substituted"],
               partner_iface_deleted=s["partner_interface_by_block"]["deleted"])
    rows.append(row)
    detail[job] = dict(summary=summary, model_file=os.path.basename(cif),
                       pae_file=os.path.basename(pae), interface=s)
    print("%-28s %4d+%-4d ipTM %.2f pDockQ %.3f pDockQ2 %.3f iface %3d plDDT %.1f  "
          "partner iface by block %3d/%d/%-3d"
          % (job, row["chain_a_residues"], row["chain_b_residues"], row["iptm"],
             row["pDockQ"], row["pDockQ2"], row["interface_residues"],
             row["interface_plddt"], row["partner_iface_retained"],
             row["partner_iface_substituted"], row["partner_iface_deleted"]))

tsv = os.path.join(RESULTS, "af3_scores.tsv")
with open(tsv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
settings, seeds, templates = {}, set(), set()
for req in sorted(glob.glob(os.path.join(MODELS, "*", "fold_*_job_request.json"))):
    for job in json.load(open(req)):
        for chain in job["sequences"]:
            for body in chain.values():
                if "useStructureTemplate" in body:
                    templates.add(bool(body["useStructureTemplate"]))
        for sd in job["modelSeeds"]:
            seeds.add(str(sd))
        settings[job["name"]] = dict(
            dialect=job.get("dialect"), version=job.get("version"),
            seeds=[str(x) for x in job["modelSeeds"]])
if not settings:
    raise SystemExit("STOP: no job request beside the predictions under %s, so what "
                     "the server was asked for is not recorded." % MODELS)
if len(templates) != 1 or len(seeds) != 1:
    raise SystemExit("STOP: the jobs do not agree on their settings. Structure "
                     "templates %s and seeds %s appear across %d jobs, so one "
                     "sentence cannot describe them."
                     % (sorted(templates), sorted(seeds), len(settings)))

json.dump(dict(models_root=os.path.relpath(MODELS, ROOT), contact_cutoff_A=CONTACT_A, pae_d0=PAE_D0,
               server=dict(jobs=len(settings), use_structure_template=templates.pop(),
                           seed=seeds.pop(), per_job=settings),
               scores=rows, detail=detail),
          open(os.path.join(RESULTS, "af3_scores.json"), "w"), indent=1)
print("\nwrote %s and af3_scores.json" % tsv)
