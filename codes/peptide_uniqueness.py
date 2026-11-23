import gzip, json, os, time, urllib.error, urllib.request

ACCESSION = "Q16719"
ISOFORM = "Q16719-2"
PROTEOME = "UP000005640"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
for d in (DATA, RESULTS):
    os.makedirs(d, exist_ok=True)
FASTA = os.path.join(DATA, "human_reviewed_with_isoforms.fasta.gz")
MISSED_CLEAVAGES = 2
MIN_LEN, MAX_LEN = 7, 40
RETRIES = 4
BACKOFF = 5


def fetch(url, target=None):
    last = None
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=600) as r:
                payload = r.read()
            if target:
                open(target, "wb").write(payload)
                return target
            return payload.decode()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            code = getattr(exc, "code", None)
            if code is not None and code not in (429, 500, 502, 503, 504):
                raise
            last = exc
            time.sleep(BACKOFF * (attempt + 1))
    raise SystemExit("%s unreachable after %d attempts: %s" % (url, RETRIES, last))


def tryptic(seq):
    sites = [0] + [i + 1 for i in range(len(seq) - 1) if seq[i] in "KR" and seq[i + 1] != "P"] + [len(seq)]
    sites = sorted(set(sites))
    out = set()
    for i in range(len(sites) - 1):
        for m in range(MISSED_CLEAVAGES + 1):
            j = i + 1 + m
            if j < len(sites):
                p = seq[sites[i]:sites[j]]
                if MIN_LEN <= len(p) <= MAX_LEN:
                    out.add(p)
    return out


canonical = json.loads(fetch("https://rest.uniprot.org/uniprotkb/%s.json" % ACCESSION))["sequence"]["value"]
isoform = "".join(fetch("https://rest.uniprot.org/uniprotkb/%s.fasta" % ISOFORM).split("\n")[1:]).strip()
diagnostic = sorted(p for p in tryptic(isoform) if p not in canonical)
if not diagnostic:
    raise SystemExit("no isoform-diagnostic tryptic peptide was derived")

if not os.path.exists(FASTA):
    fetch("https://rest.uniprot.org/uniprotkb/stream?query=%28proteome%3A" + PROTEOME +
          "%29+AND+%28reviewed%3Atrue%29&includeIsoform=true&format=fasta&compressed=true", FASTA)

entries, seqs = 0, {}
name, buf = None, []
with gzip.open(FASTA, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
            entries += 1
        else:
            buf.append(line.strip())
if name:
    seqs[name] = "".join(buf)
if entries < 20000:
    raise SystemExit("the proteome file holds only %d entries" % entries)
if not any(k.endswith("|%s|" % ACCESSION) or ("|%s|" % ACCESSION) in k for k in seqs):
    matching = [k for k in seqs if ACCESSION in k]
    if not matching:
        raise SystemExit("%s is absent from the downloaded proteome" % ACCESSION)

rows = []
for peptide in diagnostic:
    carriers = sorted(k for k, s in seqs.items() if peptide in s)
    accessions = sorted({k.split("|")[1] for k in carriers if "|" in k})
    genes = sorted({k.split("|")[1].split("-")[0] for k in carriers if "|" in k})
    rows.append({"peptide": peptide, "length": len(peptide),
                 "protein_entries_carrying_it": len(carriers),
                 "distinct_accessions": ";".join(accessions) or "none",
                 "distinct_base_accessions": len(genes),
                 "observed_in_peptideatlas": peptide == "HAHTIKPARSEFFN"})

with open(os.path.join(RESULTS, "peptide_uniqueness.tsv"), "w") as fh:
    fh.write("peptide\tlength\tprotein_entries_carrying_it\tdistinct_base_accessions\tdistinct_accessions\tobserved_in_peptideatlas\n")
    for r in rows:
        fh.write("%s\t%d\t%d\t%d\t%s\t%s\n" % (r["peptide"], r["length"], r["protein_entries_carrying_it"],
                                               r["distinct_base_accessions"], r["distinct_accessions"],
                                               r["observed_in_peptideatlas"]))
json.dump({"proteome": PROTEOME, "entries_searched": entries,
           "missed_cleavages": MISSED_CLEAVAGES, "length_window": [MIN_LEN, MAX_LEN],
           "diagnostic_peptides": rows}, open(os.path.join(RESULTS, "peptide_uniqueness.json"), "w"), indent=1)

controls = {}
for peptide in sorted(tryptic(canonical)):
    if peptide in isoform:
        controls.setdefault("shared", []).append(peptide)
    elif canonical.find(peptide) + 1 >= 302:
        controls.setdefault("canonical_specific", []).append(peptide)
panel = os.path.join(DATA, "peptide_panel.tsv")
with open(panel, "w") as fh:
    fh.write("peptide\tlength\trole\tstart_in_its_own_sequence\tproteome_entries_carrying_it\n")
    for r in rows:
        fh.write("%s\t%d\tisoform_diagnostic\t%d\t%d\n"
                 % (r["peptide"], r["length"], isoform.find(r["peptide"]) + 1, r["protein_entries_carrying_it"]))
    for peptide in controls.get("canonical_specific", [])[:8]:
        fh.write("%s\t%d\tcanonical_specific_control\t%d\t%d\n"
                 % (peptide, len(peptide), canonical.find(peptide) + 1,
                    sum(1 for s2 in seqs.values() if peptide in s2)))
    for peptide in controls.get("shared", [])[:8]:
        fh.write("%s\t%d\tshared_detection_control\t%d\t%d\n"
                 % (peptide, len(peptide), canonical.find(peptide) + 1,
                    sum(1 for s2 in seqs.values() if peptide in s2)))
print("panel written: %s" % panel)
print("searched %d reviewed human entries including isoforms" % entries)
print("isoform-diagnostic tryptic peptides, up to %d missed cleavages:" % MISSED_CLEAVAGES)
for r in rows:
    print("  %-28s %2d aa  carried by %d entry/entries  %s  observed: %s"
          % (r["peptide"], r["length"], r["protein_entries_carrying_it"],
             r["distinct_accessions"], r["observed_in_peptideatlas"]))
