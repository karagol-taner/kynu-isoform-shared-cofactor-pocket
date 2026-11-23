import csv, gzip, io, json, os, re, shutil, time, urllib.error, urllib.parse, urllib.request, zipfile, collections

ACCESSION = "Q16719"
ISOFORM = "Q16719-2"
TISSUES = ("liver", "kidney")
PAGES = 4
PAGE_SIZE = 100
PREFERRED = ("peptides", "proteingroups", "report", "mztab", "evidence", "psm", "combined", "msms", "result")
TEXT_SUFFIXES = (".txt", ".tsv", ".csv", ".mztab", ".mzid", ".gz", ".zip")
FILES_PER_PROJECT = 2
FILE_CAP_MB = 250
TOTAL_CAP_MB = 9000
PAUSE = 1.0
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
SCRATCH = os.path.join(RESULTS, "sweep_scratch")
os.makedirs(SCRATCH, exist_ok=True)
RETRIES = 3
BACKOFF = 4


def get_json(url):
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                body = r.read().decode()
            return json.loads(body) if body.strip() else []
        except (OSError, json.JSONDecodeError) as exc:
            code = getattr(exc, "code", None)
            if code is not None and code not in (429, 500, 502, 503, 504):
                return None
            time.sleep(BACKOFF * (attempt + 1))
    return None


def get_text(url):
    with urllib.request.urlopen(url, timeout=300) as r:
        return r.read().decode()


entry = json.loads(get_text("https://rest.uniprot.org/uniprotkb/%s.json" % ACCESSION))
ensp = {}
for x in entry.get("uniProtKBCrossReferences", []):
    if x["database"] == "Ensembl":
        for p in x.get("properties", []):
            if p["key"] == "ProteinId" and p["value"].startswith("ENSP"):
                ensp[p["value"].split(".")[0]] = x.get("isoformId") or ACCESSION
iso_ensp = sorted(p for p, v in ensp.items() if v == ISOFORM)
can_ensp = sorted(p for p, v in ensp.items() if v != ISOFORM)
if not iso_ensp:
    raise SystemExit("no Ensembl protein identifier maps to %s" % ISOFORM)

panel = {}
with open(os.path.join(DATA, "peptide_panel.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        panel[r["peptide"]] = r["role"]
diagnostic = [p for p, r in panel.items() if r == "isoform_diagnostic"]
decoys = []
with open(os.path.join(DATA, "search_additions.fasta")) as fh:
    seq, name = [], None
    for line in fh:
        if line.startswith(">"):
            if name and "DECOY_TAIL" in name:
                decoys.append("".join(seq)[-14:])
            name, seq = line, []
        else:
            seq.append(line.strip())
    if name and "DECOY_TAIL" in name:
        decoys.append("".join(seq)[-14:])

print("isoform Ensembl ids: %s" % iso_ensp)
print("canonical Ensembl ids: %d" % len(can_ensp))
print("diagnostic peptides: %d | decoy probes: %d" % (len(diagnostic), len(decoys)), flush=True)

projects = {}
for tissue in TISSUES:
    for page in range(PAGES):
        found = get_json("https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects"
                         "?keyword=%s&pageSize=%d&page=%d&filter=submissionType%%3D%%3DCOMPLETE"
                         % (urllib.parse.quote(tissue), PAGE_SIZE, page))
        if not found:
            break
        for p in found:
            if any("homo sapiens" in str(o).lower() for o in p.get("organisms") or []):
                projects.setdefault(p["accession"], p)
print("COMPLETE human projects: %d" % len(projects), flush=True)


def https_of(locations):
    for l in locations or []:
        if isinstance(l, dict):
            v = str(l.get("value", ""))
            if v.startswith("ftp://ftp.pride.ebi.ac.uk/"):
                return v.replace("ftp://ftp.pride.ebi.ac.uk/", "https://ftp.pride.ebi.ac.uk/")
    return None


def texts_in(path):
    low = path.lower()
    if low.endswith(".gz"):
        with gzip.open(path, "rt", errors="ignore") as fh:
            yield os.path.basename(path), fh.read()
    elif low.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            for inner in z.namelist():
                if inner.lower().endswith((".txt", ".tsv", ".csv", ".mztab", ".mzid")):
                    with z.open(inner) as fh:
                        yield inner, io.TextIOWrapper(fh, errors="ignore").read()
    else:
        with open(path, "r", errors="ignore") as fh:
            yield os.path.basename(path), fh.read()


def sequence_column_hits(text, probes):
    lines = text.split("\n", 1)
    if len(lines) < 2 or "\t" not in lines[0]:
        return None
    header = lines[0].rstrip("\r").split("\t")
    idx = next((i for i, c in enumerate(header) if c.strip().lower() in
                ("sequence", "peptide", "peptide sequence", "stripped.sequence", "modified sequence")), None)
    if idx is None:
        return None
    found = set()
    for line in lines[1].split("\n"):
        f = line.split("\t")
        if len(f) > idx:
            v = re.sub(r"[^A-Z]", "", f[idx].upper())
            if v in probes:
                found.add(v)
    return found


target = os.path.join(RESULTS, "search_space_sweep.tsv")
fields = ["project", "file", "size_mb", "column_aware", "diagnostic_in_sequence_column",
          "diagnostic_substring", "decoy_substring", "isoform_ensp", "canonical_ensp",
          "uniprot_isoform", "uniprot_base", "gene_symbol", "note"]
downloaded = 0.0
done_projects = set()
if os.path.exists(target):
    with open(target) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            done_projects.add(r["project"])
print("resuming, %d projects already recorded" % len(done_projects), flush=True)
fresh = not os.path.exists(target)
with open(target, "a") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    if fresh:
        writer.writeheader()
    out.flush()
    for accession in sorted(projects):
        if accession in done_projects:
            continue
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
            low = f["fileName"].lower()
            if not low.endswith(TEXT_SUFFIXES):
                continue
            rank = next((i for i, k in enumerate(PREFERRED) if k in low), None)
            if rank is None:
                continue
            size_mb = (f.get("fileSizeBytes") or 0) / 1e6
            if size_mb > FILE_CAP_MB:
                continue
            url = https_of(f.get("publicFileLocations"))
            if url:
                picks.append((rank, size_mb, f["fileName"], url))
        picks.sort(key=lambda x: (x[0], x[1]))
        if not picks:
            writer.writerow({"project": accession, "file": "", "size_mb": 0, "note": "no peptide-level text file"})
            out.flush()
            continue
        for _, size_mb, name, url in picks[:FILES_PER_PROJECT]:
            if downloaded + size_mb > TOTAL_CAP_MB:
                writer.writerow({"project": accession, "file": name, "size_mb": round(size_mb, 2),
                                 "note": "skipped, download cap"})
                out.flush()
                continue
            local = os.path.join(SCRATCH, re.sub(r"[^A-Za-z0-9._-]", "_", name))
            row = {"project": accession, "file": name, "size_mb": round(size_mb, 2), "note": ""}
            try:
                time.sleep(PAUSE)
                with urllib.request.urlopen(url, timeout=600) as r, open(local, "wb") as fh:
                    shutil.copyfileobj(r, fh)
                downloaded += size_mb
                col_aware = False
                diag_col, diag_sub, dec_sub = set(), set(), set()
                iso_hit = can_hit = uni_iso = uni_base = gene = False
                for inner, text in texts_in(local):
                    hits = sequence_column_hits(text, set(diagnostic))
                    if hits is not None:
                        col_aware = True
                        diag_col |= hits
                    for p in diagnostic:
                        if p in text:
                            diag_sub.add(p)
                    for p in decoys:
                        if p in text:
                            dec_sub.add(p)
                    if any(e in text for e in iso_ensp):
                        iso_hit = True
                    if any(e in text for e in can_ensp):
                        can_hit = True
                    if ISOFORM in text:
                        uni_iso = True
                    if ACCESSION in text:
                        uni_base = True
                    if re.search(r"\bKYNU\b", text):
                        gene = True
                row.update({"column_aware": col_aware,
                            "diagnostic_in_sequence_column": ";".join(sorted(diag_col)),
                            "diagnostic_substring": ";".join(sorted(diag_sub)),
                            "decoy_substring": ";".join(sorted(dec_sub)),
                            "isoform_ensp": iso_hit, "canonical_ensp": can_hit,
                            "uniprot_isoform": uni_iso, "uniprot_base": uni_base, "gene_symbol": gene})
                if diag_col:
                    print("SEQUENCE-COLUMN HIT %s %s -> %s" % (accession, name, sorted(diag_col)), flush=True)
            except Exception as exc:
                row["note"] = "failed: %s" % str(exc)[:70]
            finally:
                if os.path.exists(local):
                    os.remove(local)
            writer.writerow(row)
            out.flush()

with open(target) as fh:
    done = list(csv.DictReader(fh, delimiter="\t"))
ok = [r for r in done if not r["note"]]
print("\nprojects touched: %d | files searched: %d | %.0f MB" % (len({r["project"] for r in done}), len(ok), downloaded))
print("files parsed column-aware: %d" % sum(1 for r in ok if r["column_aware"] == "True"))
print("DIAGNOSTIC IN A SEQUENCE COLUMN: %d" % sum(1 for r in ok if r["diagnostic_in_sequence_column"]))
print("diagnostic as a bare substring: %d" % sum(1 for r in ok if r["diagnostic_substring"]))
print("decoy as a bare substring: %d" % sum(1 for r in ok if r["decoy_substring"]))
print("files naming an isoform-2 Ensembl id: %d" % sum(1 for r in ok if r["isoform_ensp"] == "True"))
print("files naming a canonical Ensembl id: %d" % sum(1 for r in ok if r["canonical_ensp"] == "True"))
print("files naming Q16719-2: %d | Q16719: %d | KYNU: %d"
      % (sum(1 for r in ok if r["uniprot_isoform"] == "True"),
         sum(1 for r in ok if r["uniprot_base"] == "True"),
         sum(1 for r in ok if r["gene_symbol"] == "True")))
print("projects whose search space contained the isoform: %d"
      % len({r["project"] for r in ok if r["isoform_ensp"] == "True" or r["uniprot_isoform"] == "True"}))
