import json, os, time, urllib.error, urllib.request, csv

GENE = "ENSG00000115919"
CANONICAL_LENGTH = 465
ISOFORM_LENGTH = 307
ISOFORM_TRANSCRIPT = "ENST00000375773"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)
RETRIES = 4
BACKOFF = 5


def get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    last = None
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode())
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            code = getattr(exc, "code", None)
            if code is not None and code not in (429, 500, 502, 503, 504):
                raise
            last = exc
            time.sleep(BACKOFF * (attempt + 1))
    raise SystemExit("%s unreachable after %d attempts: %s" % (url, RETRIES, last))


gene = get("https://rest.ensembl.org/lookup/id/%s?expand=1" % GENE)
transcripts = gene.get("Transcript", [])
if not transcripts:
    raise SystemExit("%s returned no transcripts" % GENE)

rows = []
for t in transcripts:
    aa = t.get("Translation", {}).get("length") if t.get("Translation") else None
    rows.append({
        "transcript": "%s.%s" % (t["id"], t.get("version", "")),
        "biotype": t.get("biotype", ""),
        "translation_aa": aa if aa is not None else "",
        "transcript_support_level": (t.get("transcript_support_level") or "not given").split()[0],
        "is_canonical": bool(t.get("is_canonical")),
        "encodes": ("canonical isoform 1" if aa == CANONICAL_LENGTH
                    else "isoform 2" if aa == ISOFORM_LENGTH else ""),
    })
rows.sort(key=lambda r: (r["biotype"], str(r["translation_aa"])))

iso = [r for r in rows if r["transcript"].startswith(ISOFORM_TRANSCRIPT)]
if len(iso) != 1:
    raise SystemExit("expected one %s transcript, found %d" % (ISOFORM_TRANSCRIPT, len(iso)))
if iso[0]["biotype"] != "protein_coding" or iso[0]["translation_aa"] != ISOFORM_LENGTH:
    raise SystemExit("%s is %s translating to %s" % (ISOFORM_TRANSCRIPT, iso[0]["biotype"], iso[0]["translation_aa"]))

nmd = [r for r in rows if "nmd" in r["biotype"].lower()]
coding = [r for r in rows if r["biotype"] == "protein_coding"]
with open(os.path.join(RESULTS, "transcript_support.tsv"), "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("gene %s, transcripts %d, protein-coding %d, NMD-biotype %d"
      % (GENE, len(rows), len(coding), len(nmd)))
print("the isoform 2 transcript is %s, biotype %s, translation %s aa, support level %s"
      % (iso[0]["transcript"], iso[0]["biotype"], iso[0]["translation_aa"], iso[0]["transcript_support_level"]))
print("transcripts encoding a 307 aa protein: %s"
      % [r["transcript"] for r in rows if r["translation_aa"] == ISOFORM_LENGTH])
print("transcripts with an NMD biotype: %s" % ([r["transcript"] for r in nmd] or "none"))
print()
print("%-22s %-28s %8s %5s %s" % ("transcript", "biotype", "aa", "tsl", "encodes"))
for r in rows:
    print("%-22s %-28s %8s %5s %s" % (r["transcript"], r["biotype"], r["translation_aa"],
                                      r["transcript_support_level"], r["encodes"]))
