import os

import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(ROOT, "structures")


def assembly(entry, keep_ligands):
    st = gemmi.read_structure(os.path.join(STRUCT, entry + ".cif"))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if "1" not in [a.name for a in st.assemblies]:
        raise SystemExit("%s has no assembly 1" % entry)
    st.transform_to_assembly("1", gemmi.HowToNameCopiedChain.AddNumber)
    if not keep_ligands:
        st.remove_ligands_and_waters()
    model = st[0]
    if len(model) < 2:
        raise SystemExit("%s assembly 1 has %d chains" % (entry, len(model)))
    if len(model) > 8:
        raise SystemExit("%s assembly 1 has %d chains, beyond this analysis" % (entry, len(model)))
    for i, chain in enumerate(model):
        chain.name = "ABCDEFGH"[i]
    st.setup_entities()
    return st
