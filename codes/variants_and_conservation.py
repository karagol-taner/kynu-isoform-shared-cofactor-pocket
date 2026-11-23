import json, os, sys, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from conservation_common import (ACCESSION, FAMILY, BLOCKS, CHAIN_LENGTH,
                                 INTERFACE_IN_LOST_SPAN, CROSS_SUBUNIT_COFACTOR,
                                 get, block_of, aligner)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)


payload = json.loads(get("https://www.ebi.ac.uk/proteins/api/variation/%s" % ACCESSION,
                         accept="application/json"))
record = payload[0] if isinstance(payload, list) else payload
reference = record["sequence"]
if len(reference) != CHAIN_LENGTH:
    raise SystemExit("variation endpoint serves a %d residue sequence" % len(reference))

usable, malformed = [], []
for f in record.get("features", []):
    begin = f.get("begin")
    if not begin or not str(begin).isdigit():
        continue
    position = int(begin)
    if position < 1 or position > len(reference):
        continue
    if f.get("wildType") != reference[position - 1]:
        malformed.append({"position": position, "stated_wild_type": f.get("wildType"),
                          "reference_residue": reference[position - 1],
                          "alternative": f.get("alternativeSequence"),
                          "sources": [x.get("name") for x in f.get("xrefs") or []]})
        continue
    usable.append((position, f))

by_block = collections.Counter(block_of(p) for p, _ in usable)
sig_by_block = collections.defaultdict(collections.Counter)
for p, f in usable:
    for s in f.get("clinicalSignificances") or []:
        sig_by_block[block_of(p)][s.get("type")] += 1

with open(os.path.join(RESULTS, "variants_by_block.tsv"), "w") as fh:
    fh.write("block\tstart\tend\tresidues\tvariants\tvariants_per_residue\tclinical_significance_tally\n")
    for name, lo, hi in BLOCKS:
        n = hi - lo + 1
        fh.write("%s\t%d\t%d\t%d\t%d\t%.3f\t%s\n" % (
            name, lo, hi, n, by_block[name], by_block[name] / n,
            dict(sig_by_block[name]) or "none"))

at_contact = []
for p, f in usable:
    if p in INTERFACE_IN_LOST_SPAN:
        at_contact.append({"position": p, "wild_type": f["wildType"],
                           "alternative": f.get("alternativeSequence"),
                           "consequence": f.get("consequenceType"),
                           "clinical": [s.get("type") for s in f.get("clinicalSignificances") or []] or None,
                           "sources": [x.get("name") for x in f.get("xrefs") or []],
                           "cross_subunit_cofactor_contact": p in CROSS_SUBUNIT_COFACTOR})

family = json.loads(get("https://rest.uniprot.org/uniprotkb/search?query=xref:interpro-%s+AND+reviewed:true"
                        "&fields=accession,organism_name,sequence&format=json&size=200" % FAMILY))
human = None
orthologues = {}
for item in family.get("results", []):
    acc = item["primaryAccession"]
    org = item.get("organism", {}).get("scientificName", "")
    seq = item["sequence"]["value"]
    if acc == ACCESSION:
        human = seq
        continue
    if org not in orthologues or len(seq) > len(orthologues[org][1]):
        orthologues[org] = (acc, seq)
if human is None:
    raise SystemExit("human %s absent from the %s family query" % (ACCESSION, FAMILY))
if human != reference:
    raise SystemExit("family query sequence disagrees with the variation endpoint sequence")

aligner = aligner()

positions = sorted(set(INTERFACE_IN_LOST_SPAN))
table = {p: {} for p in positions}
for org, (acc, seq) in sorted(orthologues.items()):
    alignment = aligner.align(human, seq)[0]
    mapping = {}
    for (hs, he), (os_, oe) in zip(alignment.aligned[0], alignment.aligned[1]):
        for offset in range(he - hs):
            mapping[hs + offset + 1] = seq[os_ + offset]
    for p in positions:
        table[p][org] = mapping.get(p, "-")

with open(os.path.join(RESULTS, "contact_conservation.tsv"), "w") as fh:
    organisms = sorted(orthologues)
    fh.write("position\thuman\tcross_subunit_cofactor_contact\tidentical\taligned\tpercent_identical\t%s\n"
             % "\t".join(organisms))
    for p in positions:
        row = [table[p][o] for o in organisms]
        aligned = sum(1 for r in row if r != "-")
        identical = sum(1 for r in row if r == human[p - 1])
        fh.write("%d\t%s\t%s\t%d\t%d\t%.1f\t%s\n" % (
            p, human[p - 1], p in CROSS_SUBUNIT_COFACTOR, identical, aligned,
            100.0 * identical / aligned if aligned else 0.0, "\t".join(row)))

json.dump({"variants_total": len(record.get("features", [])), "variants_usable": len(usable),
           "variants_malformed": malformed, "by_block": dict(by_block),
           "significance_by_block": {k: dict(v) for k, v in sig_by_block.items()},
           "variants_at_interface_residues": at_contact,
           "orthologues": {o: a for o, (a, _) in orthologues.items()}},
          open(os.path.join(RESULTS, "variants_and_conservation.json"), "w"), indent=1)

print("variants served %d | usable %d | rejected for a wild type disagreeing with the reference %d"
      % (len(record.get("features", [])), len(usable), len(malformed)))
for m in malformed:
    print("   rejected: position %d stated %s against reference %s, sources %s"
          % (m["position"], m["stated_wild_type"], m["reference_residue"], m["sources"]))
print()
for name, lo, hi in BLOCKS:
    n = hi - lo + 1
    print("%-12s %3d-%3d  %3d residues  %4d variants  %.3f per residue  %s"
          % (name, lo, hi, n, by_block[name], by_block[name] / n, dict(sig_by_block[name]) or "none"))
print("\nvariants at the interface residues inside the lost span:")
for v in at_contact:
    print("   %d %s>%s  %s  cofactor contact: %s  sources %s  clinical %s"
          % (v["position"], v["wild_type"], v["alternative"], v["consequence"],
             v["cross_subunit_cofactor_contact"], v["sources"], v["clinical"]))
print("\northologues aligned: %d" % len(orthologues))
print(open(os.path.join(RESULTS, "contact_conservation.tsv")).read())
