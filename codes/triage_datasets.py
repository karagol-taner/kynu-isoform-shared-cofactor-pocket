import csv, json, os, time, urllib.error, urllib.request, collections

TISSUES = ("liver", "kidney")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
RETRIES = 3
BACKOFF = 4


def get(url):
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
            time.sleep(BACKOFF * (attempt + 1))
    return None


def names(field):
    out = []
    for x in field or []:
        if isinstance(x, dict):
            v = x.get("name") or x.get("value")
            if v:
                out.append(str(v))
        elif x:
            out.append(str(x))
    return out


candidates = []
with open(os.path.join(RESULTS, "public_dataset_candidates.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["data_independent_or_targeted"] == "True" and any(t in r["tissues"] for t in TISSUES):
            candidates.append(r["accession"])
print("DIA or targeted candidates to triage: %d" % len(candidates), flush=True)

rows = []
for accession in candidates:
    meta = get("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s" % accession)
    if meta is None:
        rows.append({"accession": accession, "curated_tissue": "", "tissue_source": "metadata unavailable",
                     "organisms": "", "experiment_types": "", "files": 0, "total_gb": 0.0,
                     "mzml_files": 0, "raw_files": 0, "d_files": 0, "smallest_mzml_mb": "", "title": ""})
        continue
    parts = names(meta.get("organismsPart"))
    curated = sorted({t for t in TISSUES if any(t in p.lower() for p in parts)})
    freetext = sorted({t for t in TISSUES
                       if any(t in s.lower() for s in names(meta.get("keywords")) + [meta.get("title") or ""])})
    source = "curated organismsPart" if curated else ("free text only" if freetext else "neither")
    tissue = ";".join(curated or freetext)
    files = []
    if curated or freetext:
        for page in range(12):
            batch = get("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files?pageSize=100&page=%d"
                        % (accession, page))
            if not batch:
                break
            files.extend(batch)
            if len(batch) < 100:
                break
    ext = collections.Counter(f["fileName"].split(".")[-1].lower() for f in files)
    mzml = [f for f in files if f["fileName"].lower().endswith((".mzml", ".mzml.gz"))]
    rows.append({"accession": accession, "curated_tissue": tissue, "tissue_source": source,
                 "organisms": ";".join(names(meta.get("organisms"))),
                 "experiment_types": ";".join(names(meta.get("experimentTypes"))),
                 "files": len(files),
                 "total_gb": round(sum((f.get("fileSizeBytes") or 0) for f in files) / 1e9, 2),
                 "mzml_files": len(mzml), "raw_files": ext.get("raw", 0), "d_files": ext.get("d", 0),
                 "smallest_mzml_mb": (round(min((f.get("fileSizeBytes") or 0) for f in mzml) / 1e6, 1)
                                      if mzml else ""),
                 "title": (meta.get("title") or "").replace("\t", " ")[:120]})

rows.sort(key=lambda r: (r["tissue_source"] != "curated organismsPart", -r["mzml_files"], r["total_gb"]))
target = os.path.join(RESULTS, "dataset_triage.tsv")
with open(target, "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

curated = [r for r in rows if r["tissue_source"] == "curated organismsPart"]
print("curated liver or kidney: %d" % len(curated))
print("free text only: %d" % sum(1 for r in rows if r["tissue_source"] == "free text only"))
print("neither, so the earlier label was spurious: %d" % sum(1 for r in rows if r["tissue_source"] == "neither"))
print("curated projects depositing mzML: %d" % sum(1 for r in curated if r["mzml_files"]))
print()
print("%-11s %-14s %6s %9s %6s %6s %s" % ("project", "tissue", "files", "GB", "mzML", "raw", "title"))
for r in curated[:20]:
    print("%-11s %-14s %6d %9.1f %6d %6d %s"
          % (r["accession"], r["curated_tissue"], r["files"], r["total_gb"],
             r["mzml_files"], r["raw_files"], r["title"][:58]))
print("\nwrote %s" % target)
