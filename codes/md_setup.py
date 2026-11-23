import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHARMM = os.path.join(ROOT, "md", "charmmgui")
BUILDS = ("canonical", "mixed")
ION_NAMES = ("POT", "CLA", "SOD", "CAL", "MG", "ZN")
WATER_NAMES = ("TIP3", "TIP4", "SOL", "HOH", "SPC")
NOTEBOOK = os.path.join(ROOT, "codes", "colab", "kynu_MD_2_canonical.ipynb")
SPEC_NAME = "build_spec.json"


def mdp(path):
    out = {}
    for line in open(path, errors="replace"):
        line = line.split(";")[0]
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip().replace("_", "-").lower()] = v.strip()
    return out


def restraints(define):
    out = {}
    for m in re.finditer(r"-DPOSRES_FC_([A-Z]+)=([0-9.eE+-]+)", define or ""):
        out[{"BB": "backbone", "SC": "sidechain"}.get(m.group(1), m.group(1).lower())] \
            = float(m.group(2))
    return out


def builder_record(run):
    up = os.path.dirname(run)
    cfg = os.path.join(up, "input.config.dat")
    if not os.path.isfile(cfg):
        raise SystemExit("STOP: no %s, so what the builder was asked for cannot be "
                         "read back." % cfg)
    rec = json.loads(open(cfg, errors="replace").read())
    ver = None
    readme = os.path.join(run, "README")
    if os.path.isfile(readme):
        m = re.search(r"CHARMM-GUI[^\n]*?\bv([0-9][0-9.]*)",
                      open(readme, errors="replace").read())
        ver = m.group(1) if m else None
    return dict(force_field=(rec.get("forcefield") or {}).get("type"),
                hmr_requested=bool(rec.get("hmr")),
                box_requested=rec.get("dimensions"),
                charmm_gui_version=ver)


def find_build(build):
    hits = [d for d, _sub, files in os.walk(CHARMM)
            if "step3_input.gro" in files and "step5_production.mdp" in files
            and build in os.path.relpath(d, CHARMM)]
    other = [b for b in BUILDS if b != build][0]
    hits = [d for d in hits if other not in os.path.relpath(d, CHARMM)]
    if len(hits) != 1:
        raise SystemExit("STOP: %d folders under %s carry the word %s and a finished "
                         "GROMACS build, %s. One is needed."
                         % (len(hits), CHARMM, build, sorted(hits)))
    return hits[0]


def gro_residue_at(path, resid):
    names = set()
    with open(path, errors="replace") as fh:
        fh.readline()
        n = int(fh.readline())
        for _ in range(n):
            L = fh.readline()
            if L[:5].strip() == str(resid):
                names.add(L[5:10].strip())
    return sorted(names)


def itp_masses(path):
    rows, inside = [], False
    for line in open(path, errors="replace"):
        line = line.split(";")[0].rstrip()
        if line.strip().startswith("["):
            inside = line.strip().replace(" ", "") == "[atoms]"
            continue
        if not inside or not line.strip():
            continue
        bits = line.split()
        if len(bits) >= 8:
            rows.append((bits[1], float(bits[7])))
    return rows


def hydrogen_mass(path):
    seen = {}
    for atype, mass in itp_masses(path):
        if atype.upper().startswith("H"):
            seen[round(mass, 4)] = seen.get(round(mass, 4), 0) + 1
    if not seen:
        raise SystemExit("STOP: no hydrogens found in %s, so whether the build "
                         "repartitioned their mass cannot be read back." % path)
    return max(seen.items(), key=lambda kv: kv[1])[0]


spec = None
for d, _sub, files in os.walk(os.path.join(ROOT, "md")):
    if SPEC_NAME in files:
        if spec is not None:
            raise SystemExit("STOP: more than one %s under md/, so which one part 1 "
                             "wrote is not guessed here." % SPEC_NAME)
        spec = json.load(open(os.path.join(d, SPEC_NAME)))
if spec is None or "schiff_base_residue" not in spec:
    raise SystemExit("STOP: no %s under md/ carrying schiff_base_residue. Part 1 "
                     "writes it and the residue the cofactor sits on is read from it "
                     "rather than typed here." % SPEC_NAME)
SCHIFF = spec["schiff_base_residue"]

out = {}
for build in BUILDS:
    run = find_build(build)
    prod = mdp(os.path.join(run, "step5_production.mdp"))
    eq = mdp(os.path.join(run, "step4.1_equilibration.mdp"))
    mini = mdp(os.path.join(run, "step4.0_minimization.mdp"))
    with open(os.path.join(run, "step3_input.gro"), errors="replace") as fh:
        fh.readline()
        atoms = int(fh.readline())
        box = None
        for line in fh:
            box = line
        box_nm = [float(x) for x in box.split()[:3]]
    top = open(os.path.join(run, "topol.top"), errors="replace").read()
    mols = {}
    block = top.split("[ molecules ]")[-1]
    for line in block.splitlines():
        line = line.split(";")[0].strip()
        if not line or line.startswith("["):
            continue
        bits = line.split()
        if len(bits) == 2 and bits[1].isdigit():
            mols[bits[0]] = int(bits[1])
    ff = sorted({m.group(1) for m in re.finditer(r'#include\s+"([^"]+)"', top)})
    water = [f for f in ff if "water" in f.lower() or "tip" in f.lower()]
    protein = [f for f in ff if os.path.basename(f).split(".")[0].upper()
               in {k.upper() for k, v in mols.items()
                   if k.upper() not in ION_NAMES + WATER_NAMES}]
    if not protein or not water:
        raise SystemExit("STOP: %s lists %s, and the protein and water includes have "
                         "to be told apart to read the hydrogen masses back."
                         % (os.path.join(run, "topol.top"), ff))
    solute_h = hydrogen_mass(os.path.join(run, protein[0]))
    water_h = hydrogen_mass(os.path.join(run, water[0]))
    out[build] = dict(
        folder=os.path.relpath(run, ROOT),
        atoms=atoms, box_nm=[round(v, 5) for v in box_nm],
        box_A=[round(v * 10.0, 2) for v in box_nm],
        cubic=max(box_nm) - min(box_nm) < 1e-6,
        dt_ps=float(prod["dt"]), integrator=prod.get("integrator"),
        constraints=prod.get("constraints"),
        constraint_algorithm=prod.get("constraint-algorithm"),
        cutoff_scheme=prod.get("cutoff-scheme"), nstlist=prod.get("nstlist"),
        vdwtype=prod.get("vdwtype"), vdw_modifier=prod.get("vdw-modifier"),
        rvdw_switch_nm=prod.get("rvdw-switch"),
        ref_t_K=sorted({float(v) for v in prod["ref-t"].split()}),
        tcoupl=prod.get("tcoupl"), pcoupl=prod.get("pcoupl"),
        pcoupltype=prod.get("pcoupltype"),
        tau_t_ps=sorted({float(v) for v in prod["tau-t"].split()}),
        tau_p_ps=float(prod["tau-p"]),
        compressibility_bar=prod.get("compressibility"),
        ref_p_bar=prod.get("ref-p"), coulombtype=prod.get("coulombtype"),
        rvdw_nm=prod.get("rvdw"), rcoulomb_nm=prod.get("rcoulomb"),
        minimisation=dict(integrator=mini.get("integrator"),
                          steps=int(mini["nsteps"]),
                          emtol_kJ_mol_nm=float(mini["emtol"]),
                          restraint_kJ_mol_nm2=restraints(mini.get("define"))),
        equilibration_ns=round(float(eq["dt"]) * int(eq["nsteps"]) / 1000.0, 4),
        equilibration_dt_ps=float(eq["dt"]),
        equilibration_ensemble="NVT" if eq.get("pcoupl", "no").lower() == "no" else "NPT",
        equilibration_restraint_kJ_mol_nm2=restraints(eq.get("define")),
        velocities_generated_at_K=(float(eq["gen-temp"])
                                   if eq.get("gen-vel", "no").lower() == "yes" else None),
        molecules=mols,
        waters=sum(v for k, v in mols.items() if k.upper() in WATER_NAMES),
        ions={k: v for k, v in mols.items() if k.upper() in ION_NAMES},
        toppar=[f for f in ff if f.endswith(".itp") or "toppar" in f],
        water_model=water,
        water_model_name=os.path.basename(water[0]).split(".")[0],
        solute_hydrogen_mass_amu=solute_h,
        water_hydrogen_mass_amu=water_h,
        substituted_residue=gro_residue_at(os.path.join(run, "step3_input.gro"),
                                           SCHIFF),
        **builder_record(run))

nb = open(NOTEBOOK, errors="replace").read()
m = re.search(r'GMX_VER\s*=\s*\\"([0-9][0-9.]*)\\"', nb)
if not m:
    raise SystemExit("STOP: no GMX_VER in %s, so the version these ran under cannot "
                     "be read back." % NOTEBOOK)

chunk = re.search(r"NS_PER_CHUNK\s*=\s*([0-9][0-9.]*)", nb)
if not chunk:
    raise SystemExit("STOP: no NS_PER_CHUNK in %s. Production was carried in chunks "
                     "and the length of one is read from the notebook that ran them."
                     % NOTEBOOK)

shared = {"gromacs_version": m.group(1), "schiff_base_residue": SCHIFF,
          "chunk_ns": float(chunk.group(1))}
for key in ("dt_ps", "ref_t_K", "tcoupl", "pcoupl", "pcoupltype", "tau_t_ps",
            "tau_p_ps", "compressibility_bar", "ref_p_bar", "coulombtype",
            "rvdw_nm", "rcoulomb_nm", "rvdw_switch_nm", "vdwtype", "vdw_modifier",
            "cutoff_scheme", "nstlist", "constraints", "constraint_algorithm",
            "integrator", "box_A", "minimisation",
            "equilibration_ns", "equilibration_dt_ps", "equilibration_ensemble",
            "equilibration_restraint_kJ_mol_nm2", "velocities_generated_at_K",
            "water_model_name", "force_field", "hmr_requested",
            "charmm_gui_version", "solute_hydrogen_mass_amu",
            "water_hydrogen_mass_amu", "substituted_residue"):
    vals = {json.dumps(out[b][key]) for b in BUILDS}
    if len(vals) != 1:
        raise SystemExit("STOP: the builds disagree on %s, %s. The two have to share "
                         "every condition that is not the difference between them."
                         % (key, {b: out[b][key] for b in BUILDS}))
    shared[key] = out[BUILDS[0]][key]
excess = {b: out[b]["ions"].get("POT", 0) - out[b]["ions"].get("CLA", 0)
          for b in BUILDS}
shared["neutralising_cation_excess"] = excess
shared["solute_charge"] = {b: -v for b, v in excess.items()}
shared["hydrogen_mass_repartitioning"] = (shared["solute_hydrogen_mass_amu"]
                                          > shared["water_hydrogen_mass_amu"])
AVOGADRO = 6.02214076e23
for b in BUILDS:
    litres = (out[b]["box_nm"][0] * out[b]["box_nm"][1] * out[b]["box_nm"][2]) * 1e-24
    out[b]["salt_M"] = round(out[b]["ions"].get("CLA", 0) / (litres * AVOGADRO), 4)
shared["salt_M"] = {b: out[b]["salt_M"] for b in BUILDS}
if shared["hmr_requested"] != shared["hydrogen_mass_repartitioning"]:
    raise SystemExit("STOP: the builder record says repartitioning was %s and the "
                     "topology carries solute hydrogens of %g amu against water at "
                     "%g. The two have to agree."
                     % ("asked for" if shared["hmr_requested"] else "not asked for",
                        shared["solute_hydrogen_mass_amu"],
                        shared["water_hydrogen_mass_amu"]))
if shared["dt_ps"] >= 0.004 and not shared["hydrogen_mass_repartitioning"]:
    raise SystemExit("STOP: the timestep is %g fs and the solute hydrogens weigh "
                     "%g amu, so the build was not repartitioned and the step is not "
                     "safe." % (shared["dt_ps"] * 1000.0,
                                shared["solute_hydrogen_mass_amu"]))
shared["timestep_fs"] = shared["dt_ps"] * 1000.0

doc = dict(shared=shared, builds=out)
path = os.path.join(ROOT, "results", "md_setup.json")
json.dump(doc, open(path, "w"), indent=2)
for b in BUILDS:
    print("%-10s %7d atoms, box %s A, %d waters, ions %s, solute charge %+d"
          % (b, out[b]["atoms"], out[b]["box_A"][0], out[b]["waters"],
             out[b]["ions"], -excess[b]))
print("shared: %s from CHARMM-GUI v%s, GROMACS %s, %g fs timestep, %s by %s, "
      "%s K, %s, %s, residue %s at %d"
      % (shared["force_field"], shared["charmm_gui_version"],
         shared["gromacs_version"], shared["timestep_fs"], shared["constraints"],
         shared["constraint_algorithm"], shared["ref_t_K"], shared["tcoupl"],
         shared["pcoupl"], shared["substituted_residue"], SCHIFF))
print("        KCl at %s M in the box"
      % ", ".join("%.3f" % out[b]["salt_M"] for b in BUILDS))
print("        solute H %g amu against water H %g amu, %g ns equilibration in %s, "
      "production in %g ns chunks"
      % (shared["solute_hydrogen_mass_amu"], shared["water_hydrogen_mass_amu"],
         shared["equilibration_ns"], shared["equilibration_ensemble"],
         shared["chunk_ns"]))
print("wrote", path)
