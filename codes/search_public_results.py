import csv, gzip, io, json, os, shutil, time, urllib.error, urllib.parse, urllib.request, zipfile

TISSUES = ("liver", "kidney")
PAGES = 4
PAGE_SIZE = 100
RESULT_PATTERNS = ("peptide", "report", "msms", "evidence", "mztab", "mzid", "spectronaut",
                   "diann", "dia-nn", "proteingroup", "psm", "result", "search", "combined")
TEXT_SUFFIXES = (".txt", ".tsv", ".csv", ".mztab", ".mzid", ".xml", ".gz", ".zip")
FILE_CAP_MB = 300
TOTAL_CAP_MB = 6000
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
SCRATCH = os.path.join(RESULTS, "download_scratch")
os.makedirs(SCRATCH, exist_ok=True)
RETRIES = 3
BACKOFF = 4

panel = {"isoform_diagnostic": [], "canonical_specific_control": [], "shared_detection_control": []}
with open(os.path.join(DATA, "peptide_panel.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        panel[r["role"]].append(r["peptide"])
if not panel["isoform_diagnostic"]:
    raise SystemExit("the peptide panel holds no diagnostic peptide")
ACCESSION_PROBES = ("Q16719-2", "Q16719")


def get_json(url):
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


def https_of(locations):
    for l in locations or []:
        if isinstance(l, dict):
            v = str(l.get("value", ""))
            if v.startswith("ftp://ftp.pride.ebi.ac.uk/"):
                return v.replace("ftp://ftp.pride.ebi.ac.uk/", "https://ftp.pride.ebi.ac.uk/")
            if v.startswith("https://"):
                return v
    return None


def texts_in(path):
    name = path.lower()
    try:
        if name.endswith(".gz"):
            with gzip.open(path, "rt", errors="ignore") as fh:
                yield fh.read()
        elif name.endswith(".zip"):
            with zipfile.ZipFile(path) as z:
                for inner in z.namelist():
                    if inner.lower().endswith((".txt", ".tsv", ".csv", ".mztab", ".mzid", ".xml")):
                        with z.open(inner) as fh:
                            yield io.TextIOWrapper(fh, errors="ignore").read()
        else:
            with open(path, "r", errors="ignore") as fh:
                yield fh.read()
    except Exception as exc:
        raise RuntimeError(str(exc))


projects = {}
for tissue in TISSUES:
    for page in range(PAGES):
        url = ("https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects"
               "?keyword=%s&pageSize=%d&page=%d&filter=submissionType%%3D%%3DCOMPLETE"
               % (urllib.parse.quote(tissue), PAGE_SIZE, page))
        found = get_json(url)
        if not found:
            break
        for p in found:
            if any("homo sapiens" in (o or "").lower() for o in p.get("organisms") or []):
                projects.setdefault(p["accession"], p)
print("COMPLETE human projects for %s: %d" % (", ".join(TISSUES), len(projects)), flush=True)

rows = []
downloaded_mb = 0.0
target = os.path.join(RESULTS, "public_result_file_search.tsv")
fields = ["project", "file", "size_mb", "diagnostic_hits", "canonical_control_hits",
          "shared_control_hits", "mentions_isoform_accession", "mentions_base_accession", "note"]
with open(target, "w") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    out.flush()
    for accession in sorted(projects):
        files = []
        for page in range(12):
            batch = get_json("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files?pageSize=%d&page=%d"
                             % (accession, PAGE_SIZE, page))
            if not batch:
                break
            files.extend(batch)
            if len(batch) < PAGE_SIZE:
                break
        picks = []
        for f in files:
            if (f.get("fileCategory") or {}).get("value") == "RAW":
                continue
            name = f["fileName"]
            low = name.lower()
            if not low.endswith(TEXT_SUFFIXES) or not any(p in low for p in RESULT_PATTERNS):
                continue
            size_mb = (f.get("fileSizeBytes") or 0) / 1e6
            if size_mb > FILE_CAP_MB:
                continue
            url = https_of(f.get("publicFileLocations"))
            if url:
                picks.append((name, size_mb, url))
        for name, size_mb, url in picks:
            if downloaded_mb + size_mb > TOTAL_CAP_MB:
                writer.writerow({"project": accession, "file": name, "size_mb": round(size_mb, 2),
                                 "diagnostic_hits": "", "canonical_control_hits": "",
                                 "shared_control_hits": "", "mentions_isoform_accession": "",
                                 "mentions_base_accession": "", "note": "skipped, total download cap reached"})
                out.flush()
                continue
            local = os.path.join(SCRATCH, os.path.basename(name))
            note = ""
            try:
                with urllib.request.urlopen(url, timeout=600) as r, open(local, "wb") as fh:
                    shutil.copyfileobj(r, fh)
                downloaded_mb += size_mb
                diag, canon, shared = set(), set(), set()
                iso_acc = base_acc = False
                for text in texts_in(local):
                    for p in panel["isoform_diagnostic"]:
                        if p in text:
                            diag.add(p)
                    for p in panel["canonical_specific_control"]:
                        if p in text:
                            canon.add(p)
                    for p in panel["shared_detection_control"]:
                        if p in text:
                            shared.add(p)
                    if ACCESSION_PROBES[0] in text:
                        iso_acc = True
                    if ACCESSION_PROBES[1] in text:
                        base_acc = True
                writer.writerow({"project": accession, "file": name, "size_mb": round(size_mb, 2),
                                 "diagnostic_hits": ";".join(sorted(diag)),
                                 "canonical_control_hits": len(canon), "shared_control_hits": len(shared),
                                 "mentions_isoform_accession": iso_acc, "mentions_base_accession": base_acc,
                                 "note": note})
                if diag:
                    print("HIT %s %s -> %s" % (accession, name, sorted(diag)), flush=True)
            except Exception as exc:
                writer.writerow({"project": accession, "file": name, "size_mb": round(size_mb, 2),
                                 "diagnostic_hits": "", "canonical_control_hits": "",
                                 "shared_control_hits": "", "mentions_isoform_accession": "",
                                 "mentions_base_accession": "", "note": "failed: %s" % str(exc)[:80]})
            finally:
                if os.path.exists(local):
                    os.remove(local)
                out.flush()

with open(target) as fh:
    done = list(csv.DictReader(fh, delimiter="\t"))
searched = [r for r in done if not r["note"]]
print("\nprojects with a searchable result file: %d" % len({r["project"] for r in done}))
print("files searched: %d, downloaded %.0f MB" % (len(searched), downloaded_mb))
print("files naming KYNU base accession: %d" % sum(1 for r in searched if r["mentions_base_accession"] == "True"))
print("files naming the isoform accession Q16719-2: %d" % sum(1 for r in searched if r["mentions_isoform_accession"] == "True"))
print("files reporting a shared KYNU control peptide: %d" % sum(1 for r in searched if r["shared_control_hits"] not in ("", "0")))
print("files reporting a canonical-specific control peptide: %d" % sum(1 for r in searched if r["canonical_control_hits"] not in ("", "0")))
print("FILES REPORTING AN ISOFORM-DIAGNOSTIC PEPTIDE: %d" % sum(1 for r in searched if r["diagnostic_hits"]))
for r in searched:
    if r["diagnostic_hits"]:
        print("   %s %s -> %s" % (r["project"], r["file"], r["diagnostic_hits"]))
