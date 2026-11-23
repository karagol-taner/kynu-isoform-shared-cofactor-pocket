import csv, io, json, os, re, sys, collections, urllib.request, zipfile

ACCESSION = "Q16719"
ISOFORM = "Q16719-2"
TABLES = ("txt/msms.txt", "txt/evidence.txt", "txt/peptides.txt")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")


def get(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode()


entry = json.loads(get("https://rest.uniprot.org/uniprotkb/%s.json" % ACCESSION))
canonical = entry["sequence"]["value"]
isoform = "".join(get("https://rest.uniprot.org/uniprotkb/%s.fasta" % ISOFORM).split("\n")[1:]).strip()
ensp = {}
for x in entry.get("uniProtKBCrossReferences", []):
    if x["database"] != "Ensembl":
        continue
    for p in x.get("properties", []):
        if p["key"] == "ProteinId" and p["value"].startswith("ENSP"):
            ensp[p["value"].split(".")[0]] = x.get("isoformId") or ACCESSION
if not ensp:
    raise SystemExit("no Ensembl protein identifiers cross-referenced on %s" % ACCESSION)

panel = {}
with open(os.path.join(DATA, "peptide_panel.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        panel[r["peptide"]] = r["role"]

print("KYNU Ensembl protein identifiers resolved from UniProt:")
for p, iso in sorted(ensp.items()):
    print("   %-18s %s" % (p, "isoform 2" if iso == ISOFORM else "isoform 1"))
print()

rows, summary = [], {}
for archive in sys.argv[1:]:
    label = os.path.basename(archive).replace(".zip", "")
    z = zipfile.ZipFile(archive)
    database = ""
    for name in ("txt/mqpar.xml", "txt/parameters.txt"):
        if name in z.namelist():
            text = z.read(name).decode("utf8", "ignore")
            found = re.findall(r"<fastaFilePath>(.*?)</fastaFilePath>", text)
            if not found:
                found = [l.split("\t")[-1] for l in text.split("\n") if l.lower().startswith("fasta file")]
            if found:
                database = os.path.basename(found[0].replace("\\", "/"))
                break
    if not database:
        raise SystemExit("%s does not record its search database" % archive)
    present = [t for t in TABLES if t in z.namelist()]
    if not present:
        raise SystemExit("%s holds none of %s" % (archive, TABLES))
    for table in present:
        with z.open(table) as fh:
            t = io.TextIOWrapper(fh, errors="ignore")
            header = t.readline().rstrip("\n").split("\t")
            si = header.index("Sequence")
            pi = [i for i, c in enumerate(header) if c in ("Proteins", "Leading razor protein")]
            n = 0
            panel_hits = collections.Counter()
            per_entry = collections.defaultdict(set)
            for line in t:
                f = line.rstrip("\n").split("\t")
                if len(f) <= si:
                    continue
                n += 1
                s = f[si]
                if s in panel:
                    panel_hits[s] += 1
                prot = " ".join(f[i] for i in pi if i < len(f))
                for p in ensp:
                    if p in prot:
                        per_entry[p].add(s)
            key = "%s::%s" % (label, table.replace("txt/", ""))
            summary[key] = {
                "search_database": database,
                "isoform_entry_in_results": any(ensp[p] == ISOFORM for p in per_entry),
                "data_rows": n,
                "diagnostic_peptides_found": sorted(p for p in panel_hits if panel[p] == "isoform_diagnostic"),
                "decoy_peptides_found": sorted(p for p in panel_hits if panel[p] == "scrambled_tail_decoy"),
                "canonical_specific_found": sorted(p for p in panel_hits if panel[p] == "canonical_specific_control"),
                "shared_found": sorted(p for p in panel_hits if panel[p] == "shared_detection_control"),
                "distinct_peptides_per_entry": {p: len(v) for p, v in sorted(per_entry.items())},
            }
            if table == "txt/peptides.txt":
                union = set().union(*per_entry.values()) if per_entry else set()
                for s in sorted(union):
                    pos = canonical.find(s) + 1
                    rows.append({"tissue": label, "search_database": database, "peptide": s, "length": len(s),
                                 "canonical_start": pos if pos else "",
                                 "within_isoform": s in isoform,
                                 "after_residue_307": bool(pos) and pos > 307,
                                 "entries": ";".join(sorted(p for p, v in per_entry.items() if s in v)),
                                 "panel_role": panel.get(s, "")})

with open(os.path.join(RESULTS, "wang2019_kynu_peptides.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
json.dump({"dataset": "PXD010154", "ensembl_protein_identifiers": ensp,
           "diagnostic_peptides_searched": sorted(p for p, r in panel.items() if r == "isoform_diagnostic"),
           "tables": summary}, open(os.path.join(RESULTS, "wang2019_kynu_summary.json"), "w"), indent=1)

for key, s in summary.items():
    print("%-34s db %-40s isoform entry searched %-5s rows %7d | diagnostic %d | decoy %d | after-307 controls %d"
          % (key, s["search_database"], s["isoform_entry_in_results"], s["data_rows"],
             len(s["diagnostic_peptides_found"]), len(s["decoy_peptides_found"]),
             len(s["canonical_specific_found"])))
print()
for tissue in sorted({r["tissue"] for r in rows}):
    sel = [r for r in rows if r["tissue"] == tissue]
    print("%s: %d KYNU peptides | within 1-307: %d | after 307: %d"
          % (tissue, len(sel), sum(1 for r in sel if not r["after_residue_307"]),
             sum(1 for r in sel if r["after_residue_307"])))
