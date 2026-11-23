import csv, json, os, time, urllib.error, urllib.request

TISSUES = ("liver", "kidney")
INSTRUMENTS = ("Astral", "Exploris", "HF-X", "Lumos", "timsTOF", "Eclipse", "Q Exactive")
RESULT_PATTERNS = ("peptide", "report", "msms", "evidence", "mztab", ".mzid", "spectronaut",
                   "diann", "dia-nn", "proteingroups", "psm", "result", "search")
TEXT_SUFFIXES = (".txt", ".tsv", ".csv", ".mztab", ".mzid", ".xml", ".gz", ".zip", ".xlsx")
SIZE_CAP_MB = 400
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
RETRIES = 3
BACKOFF = 4


def get(url):
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                body = r.read().decode()
            return json.loads(body) if body.strip() else []
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            code = getattr(exc, "code", None)
            if code is not None and code not in (429, 500, 502, 503, 504):
                return None
            last = exc
            time.sleep(BACKOFF * (attempt + 1))
    return None


candidates = []
with open(os.path.join(RESULTS, "public_dataset_candidates.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["data_independent_or_targeted"] != "True":
            continue
        if not any(t in r["tissues"] for t in TISSUES):
            continue
        if not any(i in r["instruments"] for i in INSTRUMENTS):
            continue
        candidates.append(r)
if not candidates:
    raise SystemExit("no candidates survived the filters")

rows = []
for i, c in enumerate(candidates, 1):
    files = get("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files?pageSize=100" % c["accession"])
    if files is None:
        rows.append({"accession": c["accession"], "tissues": c["tissues"], "files_listed": -1,
                     "raw": 0, "searchable_files": 0, "searchable_mb": 0.0,
                     "note": "file listing unavailable", "examples": ""})
        continue
    raw = sum(1 for f in files if (f.get("fileCategory") or {}).get("value") == "RAW")
    searchable = []
    for f in files:
        name = f["fileName"].lower()
        size_mb = (f.get("fileSizeBytes") or 0) / 1e6
        if (f.get("fileCategory") or {}).get("value") == "RAW":
            continue
        if not name.endswith(TEXT_SUFFIXES):
            continue
        if not any(p in name for p in RESULT_PATTERNS):
            continue
        if size_mb > SIZE_CAP_MB:
            continue
        urls = [l.get("value") for l in (f.get("publicFileLocations") or [])
                if isinstance(l, dict) and str(l.get("value", "")).startswith("http")]
        if not urls:
            continue
        searchable.append({"name": f["fileName"], "mb": round(size_mb, 2), "url": urls[0]})
    rows.append({"accession": c["accession"], "tissues": c["tissues"], "files_listed": len(files),
                 "raw": raw, "searchable_files": len(searchable),
                 "searchable_mb": round(sum(s["mb"] for s in searchable), 2),
                 "note": "", "examples": ";".join(s["name"] for s in searchable[:3])})
    if searchable:
        json.dump(searchable, open(os.path.join(RESULTS, "searchable_%s.json" % c["accession"]), "w"), indent=1)

with open(os.path.join(RESULTS, "dataset_searchability.tsv"), "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

usable = [r for r in rows if r["searchable_files"] > 0]
print("shortlist surveyed: %d projects" % len(rows))
print("file listing unavailable for: %d" % sum(1 for r in rows if r["files_listed"] == -1))
print("projects depositing a searchable peptide-level text file: %d" % len(usable))
print("total searchable payload: %.1f MB" % sum(r["searchable_mb"] for r in usable))
for r in sorted(usable, key=lambda x: -x["searchable_files"])[:20]:
    print("  %-11s %-16s files %3d  searchable %2d  %7.1f MB  %s"
          % (r["accession"], r["tissues"], r["files_listed"], r["searchable_files"],
             r["searchable_mb"], r["examples"][:70]))
