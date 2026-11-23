import json, os, time, urllib.error, urllib.parse, urllib.request, csv

QUERIES = ("liver", "kidney", "lung", "placenta", "hepatocyte")
PAGES = 4
TISSUES = ("liver", "kidney", "lung", "placenta")
PAGE_SIZE = 100
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


projects, per_query = {}, {}
for query in QUERIES:
    total = 0
    for page in range(PAGES):
        url = ("https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects?keyword=%s&pageSize=%d&page=%d"
               % (urllib.parse.quote(query), PAGE_SIZE, page))
        found = get(url)
        if not found:
            break
        total += len(found)
        for p in found:
            projects.setdefault(p["accession"], p)
    per_query[query] = total
if not projects:
    raise SystemExit("the PRIDE search returned no projects for any query")


def human(p):
    return any("homo sapiens" in (o or "").lower() for o in p.get("organisms") or [])


def tissue_of(p):
    text = " ".join(str(x).lower() for x in
                    (p.get("organismsPart") or []) + (p.get("keywords") or []) + [p.get("title", "")])
    return sorted({t for t in TISSUES if t in text})


def is_dia(p):
    text = " ".join(str(x).lower() for x in
                    (p.get("experimentTypes") or []) + (p.get("keywords") or []) +
                    [p.get("title", ""), p.get("dataProcessingProtocol", "") or ""])
    return any(k in text for k in ("dia", "data-independent", "data independent", "swath", "prm", "parallel reaction"))


rows = []
for accession, p in projects.items():
    if not human(p):
        continue
    tissues = tissue_of(p)
    if not tissues:
        continue
    rows.append({
        "accession": accession,
        "tissues": ";".join(tissues),
        "data_independent_or_targeted": is_dia(p),
        "instruments": ";".join(p.get("instruments") or []),
        "experiment_types": ";".join(p.get("experimentTypes") or []),
        "submission_date": p.get("submissionDate", ""),
        "raw_files": sum(1 for f in p.get("projectFileNames") or []
                         if str(f).lower().endswith((".raw", ".d", ".wiff", ".mzml", ".mzxml"))),
        "title": (p.get("title") or "").replace("\t", " ")[:160],
    })
rows.sort(key=lambda r: (not r["data_independent_or_targeted"], -r["raw_files"], r["submission_date"]), reverse=False)

target = os.path.join(RESULTS, "public_dataset_candidates.tsv")
with open(target, "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("queries run, each against the PRIDE v3 project search:")
for query, n in per_query.items():
    print("   %-38s %d projects" % (query, n))
print("distinct projects returned: %d" % len(projects))
print("human projects naming one of %s: %d" % (", ".join(TISSUES), len(rows)))
print("of those, data-independent or targeted: %d" % sum(1 for r in rows if r["data_independent_or_targeted"]))
print("\ntop candidates:")
print("%-12s %-22s %-6s %-30s %s" % ("accession", "tissues", "DIA", "instruments", "title"))
for r in rows[:15]:
    print("%-12s %-22s %-6s %-30s %s" % (r["accession"], r["tissues"], r["data_independent_or_targeted"],
                                         r["instruments"][:30], r["title"][:70]))
print("\nwrote %s" % target)
