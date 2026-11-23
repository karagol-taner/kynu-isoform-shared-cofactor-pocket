import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_common as cc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
RECORD = os.path.join(RESULTS, "variants_and_conservation.json")
CONTACTS = os.path.join(RESULTS, "contact_conservation.tsv")
CACHE = os.path.join(RESULTS, "orthologue_sequences.json")
OUT = os.path.join(RESULTS, "per_position_conservation.tsv")

for p in (RECORD, CONTACTS):
    if not os.path.isfile(p):
        raise SystemExit("STOP: no %s, so the orthologue set this has to use is not "
                         "known. Run codes/variants_and_conservation.py first." % p)

record = json.load(open(RECORD))
if "orthologues" not in record:
    raise SystemExit("STOP: %s carries no orthologues key, so the set cannot be "
                     "pinned to the one the paper reports." % RECORD)
recorded = dict(record["orthologues"])
if not recorded:
    raise SystemExit("STOP: %s records an empty orthologue set." % RECORD)

cache = json.load(open(CACHE)) if os.path.isfile(CACHE) else {}
if not isinstance(cache, dict):
    raise SystemExit("STOP: %s is not an object of accession to sequence." % CACHE)


def fasta_sequence(text, accession):
    lines = [L.strip() for L in text.strip().split("\n")]
    if not lines or not lines[0].startswith(">"):
        raise SystemExit("STOP: the reply for %s does not begin with a FASTA header, "
                         "so it is not the sequence that was asked for. It begins "
                         "%r." % (accession, text[:80]))
    if accession not in lines[0]:
        raise SystemExit("STOP: the FASTA header returned for %s names something "
                         "else: %s" % (accession, lines[0][:120]))
    seq = "".join(lines[1:])
    if not seq or not seq.isalpha():
        raise SystemExit("STOP: the sequence returned for %s is empty or carries "
                         "characters that are not residues." % accession)
    return seq


def sequence(accession):
    if accession not in cache:
        cache[accession] = fasta_sequence(
            cc.get("https://rest.uniprot.org/uniprotkb/%s.fasta" % accession),
            accession)
    return cache[accession]


human = sequence(cc.ACCESSION)
lo, hi = cc.BLOCKS[0][1], cc.BLOCKS[-1][2]
if len(human) != hi:
    raise SystemExit("STOP: %s is %d residues long and the block table this project "
                     "works with ends at %d, so the two disagree about the chain."
                     % (cc.ACCESSION, len(human), hi))

missing = sorted(org for org, acc in recorded.items() if not acc)
if missing:
    raise SystemExit("STOP: %d of the recorded orthologues carry no accession, "
                     "beginning with %s." % (len(missing), missing[0]))

aligner = cc.aligner()
table = {}
for n, (org, acc) in enumerate(sorted(recorded.items()), 1):
    mapping = cc.human_to_orthologue(human, sequence(acc), aligner)
    table[org] = mapping
    if n % 10 == 0 or n == len(recorded):
        print("  aligned %d of %d" % (n, len(recorded)), flush=True)
json.dump(cache, open(CACHE, "w"), indent=0, sort_keys=True)

organisms = sorted(recorded)
rows = []
for p in range(lo, hi + 1):
    res = [table[o].get(p, "-") for o in organisms]
    aligned = sum(1 for r in res if r != "-")
    identical = sum(1 for r in res if r == human[p - 1])
    rows.append({"position": p, "human": human[p - 1], "aligned": aligned,
                 "identical": identical,
                 "percent_identical": 100.0 * identical / aligned if aligned else 0.0,
                 "block": cc.block_of(p),
                 "cross_subunit_cofactor_contact": p in cc.CROSS_SUBUNIT_COFACTOR})

published = {}
with open(CONTACTS) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    for L in fh:
        f = L.rstrip("\n").split("\t")
        r = dict(zip(head, f))
        published[int(r["position"])] = (r["human"], int(r["identical"]),
                                         int(r["aligned"]), float(r["percent_identical"]))
if not published:
    raise SystemExit("STOP: %s carries no rows to check against." % CONTACTS)
by_pos = {r["position"]: r for r in rows}
bad = []
for p, (hres, ident, align, pct) in sorted(published.items()):
    r = by_pos.get(p)
    if r is None:
        bad.append("%d is absent from this table" % p)
        continue
    if (r["human"], r["identical"], r["aligned"]) != (hres, ident, align) \
            or abs(r["percent_identical"] - pct) > 0.05:
        bad.append("%d: this gives %s %d of %d (%.1f) and %s gives %s %d of %d (%.1f)"
                   % (p, r["human"], r["identical"], r["aligned"],
                      r["percent_identical"], os.path.basename(CONTACTS),
                      hres, ident, align, pct))
if bad:
    raise SystemExit("STOP: the per position identity disagrees with the published "
                     "contact table at %d of %d positions, so one of the two "
                     "alignments is not the one the paper reports:\n  %s"
                     % (len(bad), len(published), "\n  ".join(bad)))

with open(OUT, "w") as fh:
    fh.write("position\thuman\tblock\tcross_subunit_cofactor_contact\tidentical\t"
             "aligned\tpercent_identical\n")
    for r in rows:
        fh.write("%d\t%s\t%s\t%s\t%d\t%d\t%.1f\n"
                 % (r["position"], r["human"], r["block"],
                    r["cross_subunit_cofactor_contact"], r["identical"],
                    r["aligned"], r["percent_identical"]))

depth = [r["aligned"] for r in rows]
print()
print("orthologues pinned to the recorded set: %d" % len(recorded))
print("positions written: %d, %d to %d" % (len(rows), rows[0]["position"],
                                           rows[-1]["position"]))
print("alignment depth per position: %d to %d, median %d"
      % (min(depth), max(depth), sorted(depth)[len(depth) // 2]))
print("agrees with %s at all %d of its positions"
      % (os.path.basename(CONTACTS), len(published)))
print("wrote %s" % OUT)
print("sequences cached in %s" % CACHE)
