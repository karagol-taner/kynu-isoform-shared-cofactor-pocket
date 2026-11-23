import csv, json, os, sys, collections
from pyteomics import mzml

FRAGMENT_PPM_HIGH_RES = 20.0
FRAGMENT_DA_LOW_RES = 0.5
LOW_RES_ANALYSERS = ("ion trap", "iontrap", "linear ion trap", "quadrupole ion trap")
MIN_FRAGMENTS = 3
DIA_WINDOW_MIN_WIDTH = 2.0
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)


def load_transitions():
    by_precursor = collections.defaultdict(list)
    meta = {}
    with open(os.path.join(DATA, "transitions.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["peptide"], int(r["precursor_charge"]))
            by_precursor[key].append((r["fragment"], float(r["fragment_mz"])))
            meta[key] = {"role": r["role"], "precursor_mz": float(r["precursor_mz"]),
                         "length": int(r["length"]), "note": r["note"]}
    if not by_precursor:
        raise SystemExit("transitions.tsv holds no transition")
    return by_precursor, meta


def analyser_of(spectrum):
    scan = ((spectrum.get("scanList") or {}).get("scan") or [{}])[0]
    for source in (scan, spectrum):
        for key in source:
            low = str(key).lower()
            if any(a in low for a in LOW_RES_ANALYSERS):
                return "low"
            if "orbitrap" in low or "fourier transform" in low or "time-of-flight" in low:
                return "high"
    filt = str(scan.get("filter string", "")).lower()
    if filt.startswith("itms") or " itms" in filt:
        return "low"
    if filt.startswith("ftms") or " ftms" in filt:
        return "high"
    return None


def tolerance_for(analyser, fragment_mz):
    if analyser == "low":
        return FRAGMENT_DA_LOW_RES
    if analyser == "high":
        return fragment_mz * FRAGMENT_PPM_HIGH_RES / 1e6
    return None


def window_of(spectrum):
    pre = (spectrum.get("precursorList") or {}).get("precursor") or [{}]
    w = pre[0].get("isolationWindow") or {}
    target = w.get("isolation window target m/z")
    lo = w.get("isolation window lower offset")
    hi = w.get("isolation window upper offset")
    if target is None:
        return None
    target = float(target)
    lo = float(lo) if lo is not None else 0.0
    hi = float(hi) if hi is not None else 0.0
    return target - lo, target + hi, lo + hi


def run(path):
    by_precursor, meta = load_transitions()
    widths = []
    hits = {k: {"scans_selected": 0, "best_fragments": 0, "best_intensity": 0.0,
                "best_rt": None, "best_scan": None, "matched": ""} for k in by_precursor}
    n_ms1 = n_ms2 = 0
    analysers = collections.Counter()
    unresolved = 0
    with mzml.read(path) as reader:
        for s in reader:
            level = s.get("ms level")
            if level == 1:
                n_ms1 += 1
                continue
            if level != 2:
                continue
            n_ms2 += 1
            analyser = analyser_of(s)
            analysers[analyser or "unresolved"] += 1
            if analyser is None:
                unresolved += 1
            win = window_of(s)
            if win is None:
                continue
            low, high, width = win
            widths.append(width)
            mzs = s.get("m/z array")
            ints = s.get("intensity array")
            if mzs is None or len(mzs) == 0:
                continue
            rt = ((s.get("scanList") or {}).get("scan") or [{}])[0].get("scan start time")
            rt = float(rt) if rt is not None else None
            for key, frags in by_precursor.items():
                pmz = meta[key]["precursor_mz"]
                if not (low <= pmz <= high):
                    continue
                hits[key]["scans_selected"] += 1
                matched, intensity = [], 0.0
                for name, fmz in frags:
                    tol = tolerance_for(analyser, fmz)
                    if tol is None:
                        continue
                    import bisect
                    i = bisect.bisect_left(mzs, fmz - tol)
                    best = 0.0
                    while i < len(mzs) and mzs[i] <= fmz + tol:
                        best = max(best, float(ints[i]))
                        i += 1
                    if best > 0:
                        matched.append(name)
                        intensity += best
                if len(matched) > hits[key]["best_fragments"] or (
                        len(matched) == hits[key]["best_fragments"] and intensity > hits[key]["best_intensity"]):
                    hits[key].update({"best_fragments": len(matched), "best_intensity": intensity,
                                      "best_rt": rt, "best_scan": s.get("id"), "matched": ";".join(matched)})
    median_width = sorted(widths)[len(widths) // 2] if widths else 0.0
    acquisition = "DIA" if median_width >= DIA_WINDOW_MIN_WIDTH else "DDA"
    if unresolved == n_ms2 and n_ms2:
        raise SystemExit("%s: the MS2 mass analyser could not be determined for any scan, so no fragment "
                         "tolerance can be chosen" % os.path.basename(path))
    return hits, meta, by_precursor, acquisition, median_width, n_ms1, n_ms2, analysers


if len(sys.argv) < 2:
    raise SystemExit("usage: extract_peptides.py <mzML> [<mzML> ...]")

rows = []
for path in sys.argv[1:]:
    if not os.path.exists(path):
        raise SystemExit("missing file: %s" % path)
    hits, meta, by_precursor, acquisition, median_width, n_ms1, n_ms2, analysers = run(path)
    print("%s | %s, median isolation width %.2f Da | MS1 %d MS2 %d | MS2 analysers %s"
          % (os.path.basename(path), acquisition, median_width, n_ms1, n_ms2, dict(analysers)), flush=True)
    for key in sorted(hits, key=lambda k: (meta[k]["role"], k[0], k[1])):
        peptide, charge = key
        h = hits[key]
        rows.append({"file": os.path.basename(path), "acquisition": acquisition,
                     "peptide": peptide, "role": meta[key]["role"], "charge": charge,
                     "precursor_mz": meta[key]["precursor_mz"],
                     "transitions_listed": len(by_precursor[key]),
                     "ms2_scans_covering_precursor": h["scans_selected"],
                     "best_fragments_matched": h["best_fragments"],
                     "best_summed_intensity": round(h["best_intensity"], 1),
                     "best_rt": h["best_rt"], "best_scan": h["best_scan"],
                     "ms2_analysers": ";".join("%s=%d" % kv for kv in sorted(analysers.items())),
                     "matched_fragments": h["matched"],
                     "called_detected": h["best_fragments"] >= MIN_FRAGMENTS})

target = os.path.join(RESULTS, "peptide_extraction.tsv")
with open(target, "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print()
by_role = collections.defaultdict(lambda: {"peptides": set(), "covered": set(), "detected": set()})
for r in rows:
    b = by_role[r["role"]]
    b["peptides"].add(r["peptide"])
    if r["ms2_scans_covering_precursor"]:
        b["covered"].add(r["peptide"])
    if r["called_detected"]:
        b["detected"].add(r["peptide"])
print("%-28s %9s %9s %9s" % ("role", "peptides", "precursor", "detected"))
for role in ("isoform_diagnostic", "scrambled_tail_decoy", "canonical_specific_control", "shared_detection_control"):
    b = by_role[role]
    print("%-28s %9d %9d %9d" % (role, len(b["peptides"]), len(b["covered"]), len(b["detected"])))
print("\ndetection rule: at least %d of the listed fragments in one MS2 scan whose isolation window covers the "
      "precursor, within %.0f ppm on a high-resolution analyser or %.2f Da on an ion trap"
      % (MIN_FRAGMENTS, FRAGMENT_PPM_HIGH_RES, FRAGMENT_DA_LOW_RES))
for r in rows:
    if r["called_detected"]:
        print("  DETECTED %-32s z=%d %s fragments %s rt %s"
              % (r["peptide"], r["charge"], r["role"], r["matched_fragments"], r["best_rt"]))
print("\nwrote %s" % target)
