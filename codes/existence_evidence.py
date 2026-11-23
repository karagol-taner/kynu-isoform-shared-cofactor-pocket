import json, os, time, urllib.error, urllib.request

ACCESSION = "Q16719"
ISOFORM = "Q16719-2"
REFSEQ_PROTEIN = "NP_001028170.1"
ENSEMBL_TRANSCRIPT = "ENST00000375773"
UNIQUE_PEPTIDE = "HAHTIKPARSEFFN"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
for d in (DATA, RESULTS):
    os.makedirs(d, exist_ok=True)


RETRIES = 4
BACKOFF = 5


def get(url, accept=None):
    req = urllib.request.Request(url, headers={"Accept": accept} if accept else {})
    last = None
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            code = getattr(exc, "code", None)
            if code is not None and code not in (429, 500, 502, 503, 504):
                raise
            last = exc
            time.sleep(BACKOFF * (attempt + 1))
    raise SystemExit("%s unreachable after %d attempts: %s" % (url, RETRIES, last))


layers = []
entry = json.loads(get("https://rest.uniprot.org/uniprotkb/%s.json" % ACCESSION))
canonical = entry["sequence"]["value"]
isoform = "".join(get("https://rest.uniprot.org/uniprotkb/%s.fasta" % ISOFORM).split("\n")[1:]).strip()
if len(isoform) != 307:
    raise SystemExit("isoform length is %d" % len(isoform))

xrefs = {}
for x in entry.get("uniProtKBCrossReferences", []):
    xrefs.setdefault(x["database"], []).append((x["id"], x.get("isoformId"),
                                                {p["key"]: p["value"] for p in x.get("properties", [])}))

ccds = [i for i, iso, _ in xrefs.get("CCDS", []) if iso == ISOFORM]
if not ccds:
    raise SystemExit("no CCDS cross-reference on %s" % ISOFORM)
layers.append({"layer": "Consensus CDS", "identifier": ";".join(ccds),
               "what_it_establishes": "NCBI, Ensembl and UCSC agree on a protein-coding CDS for this isoform",
               "verified_how": "cross-reference on the UniProt entry, isoform-resolved"})

refseq_curated = sorted({i for i, iso, _ in xrefs.get("RefSeq", []) if iso == ISOFORM and i.startswith("NP_")})
refseq_model = sorted({i for i, iso, _ in xrefs.get("RefSeq", []) if iso == ISOFORM and i.startswith("XP_")})
fasta = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            "?db=protein&id=%s&rettype=fasta&retmode=text" % REFSEQ_PROTEIN)
refseq_seq = "".join(fasta.split("\n")[1:]).strip()
if refseq_seq != isoform:
    raise SystemExit("%s does not match the %s sequence" % (REFSEQ_PROTEIN, ISOFORM))
layers.append({"layer": "Curated RefSeq", "identifier": "%s from %s" % (
    REFSEQ_PROTEIN, [p.get("NucleotideSequenceId") for i, iso, p in xrefs["RefSeq"] if i == REFSEQ_PROTEIN][0]),
    "what_it_establishes": "a manually curated RefSeq protein, not a model prediction, whose sequence is identical to the isoform",
    "verified_how": "sequence fetched from NCBI and compared residue by residue against the UniProt isoform"})

ens = json.loads(get("https://rest.ensembl.org/lookup/id/%s?expand=1" % ENSEMBL_TRANSCRIPT,
                     accept="application/json"))
if ens.get("biotype") != "protein_coding":
    raise SystemExit("%s biotype is %s" % (ENSEMBL_TRANSCRIPT, ens.get("biotype")))
translation = ens.get("Translation", {}).get("length")
if translation != len(isoform):
    raise SystemExit("%s translates to %s residues against %d" % (ENSEMBL_TRANSCRIPT, translation, len(isoform)))
layers.append({"layer": "Ensembl transcript", "identifier": "%s.%s" % (ens["id"], ens.get("version", "")),
               "what_it_establishes": "biotype %s, translating to %d residues, transcript support level %s" % (
                   ens["biotype"], translation, ens.get("transcript_support_level", "not given")),
               "verified_how": "Ensembl REST lookup, translation length compared against the isoform"})

prot = json.loads(get("https://www.ebi.ac.uk/proteins/api/proteomics?accession=%s" % ISOFORM,
                      accept="application/json"))
feats = prot[0]["features"]
unique = [f for f in feats if f.get("unique") and f.get("peptide") == UNIQUE_PEPTIDE]
if len(unique) != 1:
    raise SystemExit("expected one unique peptide record, found %d" % len(unique))
rec = unique[0]
if isoform[int(rec["begin"]) - 1:int(rec["end"])] != UNIQUE_PEPTIDE:
    raise SystemExit("the unique peptide does not sit at %s-%s of the isoform" % (rec["begin"], rec["end"]))
if UNIQUE_PEPTIDE in canonical:
    raise SystemExit("the unique peptide occurs in the canonical sequence")
sources = sorted({ev["source"]["name"] for f in feats for ev in (f.get("evidences") or [])
                  if ev.get("source", {}).get("name")})
unique_sources = sorted({ev["source"]["name"] for ev in rec["evidences"]})
props = rec["evidences"][0]["source"]["properties"]
layers.append({"layer": "Observed peptide", "identifier": "%s at isoform %s-%s, %s" % (
    UNIQUE_PEPTIDE, rec["begin"], rec["end"], rec["evidences"][0]["source"]["id"]),
    "what_it_establishes": "best probability %s, %s observations in %s sample, %s protein mapping; the only isoform-unique peptide of %d mapped" % (
        props.get("Best Probability"), props.get("Times of observations"), props.get("Number of sample"),
        props.get("Number of protein mappings"), len(feats)),
    "verified_how": "peptide located in the isoform and confirmed absent from the canonical sequence here"})

summary = {
    "accession": ACCESSION, "isoform": ISOFORM,
    "canonical_length": len(canonical), "isoform_length": len(isoform),
    "mane_select_isoform": [p.get("ProteinId") for i, iso, p in xrefs.get("MANE-Select", [])],
    "ccds": ccds, "refseq_curated": refseq_curated, "refseq_model_predicted": refseq_model,
    "ensembl_transcript": ens["id"], "ensembl_biotype": ens["biotype"],
    "peptides_mapped": len(feats), "unique_peptides": sum(1 for f in feats if f.get("unique")),
    "peptide_sources_all": sources, "peptide_sources_for_the_unique_peptide": unique_sources,
    "unique_peptide_record": rec,
}
json.dump({"layers": layers, "summary": summary}, open(os.path.join(RESULTS, "existence_evidence.json"), "w"), indent=1)
with open(os.path.join(DATA, "existence_evidence.tsv"), "w") as fh:
    fh.write("layer\tidentifier\twhat_it_establishes\tverified_how\n")
    for l in layers:
        fh.write("%s\t%s\t%s\t%s\n" % (l["layer"], l["identifier"], l["what_it_establishes"], l["verified_how"]))

for l in layers:
    print("%-20s %s" % (l["layer"], l["identifier"]))
    print("%-20s %s" % ("", l["what_it_establishes"]))
print("\nMANE-Select is on isoform 1: %s" % summary["mane_select_isoform"])
print("model-predicted RefSeqs for the isoform: %d (%s)" % (len(refseq_model), ",".join(refseq_model)))
print("peptide sources across all %d peptides: %s" % (len(feats), sources))
print("sources for the unique peptide: %s" % unique_sources)
