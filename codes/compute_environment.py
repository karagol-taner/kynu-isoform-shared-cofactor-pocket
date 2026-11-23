import glob, io, json, os, platform, re, statistics, sys

import importlib.metadata as metadata

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECORDS = os.path.join(ROOT, "md", "run_records")
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
SETUP = os.path.join(RESULTS, "md_setup.json")
LIBS = ("numpy", "scipy", "matplotlib", "biopython", "gemmi", "freesasa",
        "openpyxl", "python-docx")
COLUMNS = ("segment", "gromacs", "simd", "gpus", "ranks", "openmp_threads", "ns_per_day")
SEGMENT = re.compile(r"^prod_\d+$")
PATTERNS = {
    "deffnm": re.compile(r"gmx mdrun\b[^\n]*-deffnm\s+(\S+)"),
    "gromacs": re.compile(r"GROMACS - gmx mdrun, (\S+) \(-:"),
    "simd": re.compile(r"Compiled SIMD is (\w+)"),
    "gpus": re.compile(r"^(\d+) GPU selected for this run"),
    "ranks": re.compile(r"^Using (\d+) MPI thread"),
    "openmp_threads": re.compile(r"^Using (\d+) OpenMP threads"),
    "ns_per_day": re.compile(r"^Performance:\s+([\d.]+)"),
    "card": re.compile(r"^(NVIDIA [^,]+), (\d+) MiB$"),
}


def outputs(path):
    nb = json.load(io.open(path, encoding="utf8"))
    text = []
    for cell in nb.get("cells", []):
        for out in cell.get("outputs", []) or []:
            got = "".join(out.get("text", []) or [])
            if not got and isinstance(out.get("data"), dict):
                got = "".join(out["data"].get("text/plain", []) or [])
            if got:
                text.append(got)
    return "\n".join(text).splitlines()


def extract(lines):
    card, rows, cur = None, [], None
    for line in lines:
        bare = line.strip()
        got = PATTERNS["card"].match(bare)
        if got:
            card = (got.group(1), int(got.group(2)))
            continue
        got = PATTERNS["deffnm"].search(bare)
        if got:
            if cur is None or cur.get("segment") != got.group(1):
                cur = {"segment": got.group(1)}
            continue
        if cur is None:
            continue
        for key in ("gromacs", "simd", "gpus", "ranks", "openmp_threads", "ns_per_day"):
            got = PATTERNS[key].search(bare)
            if got:
                cur[key] = got.group(1)
                if key == "ns_per_day" and SEGMENT.match(cur["segment"]):
                    rows.append(cur)
                    cur = None
                break
    return card, rows


def write_log(build, card, rows):
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, "md_run_log_%s.tsv" % build)
    out = ["# card\t%s\t%d MiB" % card, "\t".join(COLUMNS)]
    for r in rows:
        missing = [c for c in COLUMNS if c not in r]
        if missing:
            raise SystemExit("STOP: segment %s of the %s run printed no %s, so the "
                             "log extract would carry a hole."
                             % (r.get("segment"), build, ", ".join(missing)))
        out.append("\t".join(str(r[c]) for c in COLUMNS))
    io.open(path, "w", encoding="utf8").write("\n".join(out) + "\n")
    return path


def read_log(build):
    path = os.path.join(DATA, "md_run_log_%s.tsv" % build)
    if not os.path.isfile(path):
        raise SystemExit("STOP: %s is absent and no executed notebook for the %s run "
                         "is under %s, so the hardware cannot be read."
                         % (os.path.relpath(path, ROOT), build,
                            os.path.relpath(RECORDS, ROOT)))
    head, rows = None, []
    for line in io.open(path, encoding="utf8"):
        line = line.rstrip("\n")
        if line.startswith("# card"):
            _, name, mem = line.split("\t")
            head = (name, int(mem.split()[0]))
        elif line.startswith(COLUMNS[0] + "\t"):
            continue
        elif line.strip():
            rows.append(dict(zip(COLUMNS, line.split("\t"))))
    if head is None or not rows:
        raise SystemExit("STOP: %s carries no card line or no segment."
                         % os.path.relpath(path, ROOT))
    return head, rows


def one(rows, key, build):
    seen = sorted({r[key] for r in rows})
    if len(seen) != 1:
        raise SystemExit("STOP: the segments of the %s run disagree on %s: %s. One "
                         "figure cannot stand for the run." % (build, key, seen))
    return seen[0]


def environment():
    out = {"python": platform.python_version(),
           "platform": "%s %s" % (platform.system(), platform.machine()),
           "libraries": {}}
    for name in LIBS:
        try:
            out["libraries"][name] = metadata.version(name)
        except Exception:
            raise SystemExit("STOP: %s reports no version in this environment, and "
                             "the environment record names it." % name)
    return out


def build_record():
    setup = json.load(io.open(SETUP, encoding="utf8"))
    builds = list(setup["builds"])
    chunk_ns = float(setup["shared"]["chunk_ns"])
    doc = {"simulation": {}, "throughput": {}, "analysis": environment()}
    cards, shared = {}, {}
    for build in builds:
        found = glob.glob(os.path.join(RECORDS, "*%s*.ipynb" % build))
        if len(found) > 1:
            raise SystemExit("STOP: %s hold the %s run and only one can be its record."
                             % ([os.path.basename(f) for f in found], build))
        if found:
            card, rows = extract(outputs(found[0]))
            if card is None:
                raise SystemExit("STOP: %s reports no card, so the %s run has no "
                                 "hardware record." % (os.path.basename(found[0]), build))
            write_log(build, card, rows)
        card, rows = read_log(build)
        cards[build] = card
        for key in ("gromacs", "simd", "gpus", "ranks", "openmp_threads"):
            shared.setdefault(key, {})[build] = one(rows, key, build)
        if len({r["segment"] for r in rows}) != len(rows):
            raise SystemExit("STOP: the log of the %s run names a segment twice." % build)
        ns = [float(r["ns_per_day"]) for r in rows]
        median = statistics.median(ns)
        report = os.path.join(RESULTS, "md", "%s_report.json" % build)
        if not os.path.isfile(report):
            raise SystemExit("STOP: %s is absent, and the length of the %s run is read "
                             "from it." % (os.path.relpath(report, ROOT), build))
        asked = float(json.load(io.open(report, encoding="utf8"))["ns_requested"])
        doc["throughput"][build] = {
            "ns_per_day_median": round(median, 1),
            "ns_per_day_min": round(min(ns), 1),
            "ns_per_day_max": round(max(ns), 1),
            "segments_recorded": len(ns),
            "segments_in_run": int(round(asked / chunk_ns)),
            "hours_at_median": round(asked / median * 24.0, 1),
        }
    if len(set(cards.values())) != 1:
        raise SystemExit("STOP: the runs report different cards, %s, so one card cannot "
                         "stand for both." % cards)
    for key, per_build in shared.items():
        if len(set(per_build.values())) != 1:
            raise SystemExit("STOP: the runs disagree on %s: %s." % (key, per_build))
    name, mib = list(cards.values())[0]
    doc["simulation"] = {
        "gpu": name,
        "gpu_memory_MiB": mib,
        "gpu_memory_GB": int(round(mib * (2 ** 20) / 1e9)),
        "gpus_per_run": int(shared["gpus"][builds[0]]),
        "ranks": int(shared["ranks"][builds[0]]),
        "openmp_threads": int(shared["openmp_threads"][builds[0]]),
        "simd": shared["simd"][builds[0]],
        "gromacs": shared["gromacs"][builds[0]],
    }
    if doc["simulation"]["gromacs"] != setup["shared"]["gromacs_version"]:
        raise SystemExit("STOP: the runs printed GROMACS %s and results/md_setup.json "
                         "reports %s." % (doc["simulation"]["gromacs"],
                                          setup["shared"]["gromacs_version"]))
    return doc


def self_test():
    import tempfile
    lines = [
        "NVIDIA L4, 23034 MiB",
        "$ gmx mdrun -v -deffnm step4.1_equilibration -nsteps 125000",
        "                      :-) GROMACS - gmx mdrun, 2026.3 (-:",
        "Compiled SIMD is AVX2_256, but CPU also supports AVX_512 (see log).",
        "1 GPU selected for this run.",
        "Using 1 MPI thread",
        "Using 12 OpenMP threads",
        "Performance:       27.862        0.861",
        "$ gmx mdrun -v -deffnm prod_05 -nsteps 1250000",
        "                      :-) GROMACS - gmx mdrun, 2026.3 (-:",
        "Command line:",
        "  gmx mdrun -v -deffnm prod_05 -nsteps 1250000",
        "Compiled SIMD is AVX2_256, but CPU also supports AVX_512 (see log).",
        "1 GPU selected for this run.",
        "Using 1 MPI thread",
        "Using 12 OpenMP threads",
        "Performance:      113.154        0.212",
        "$ gmx mdrun -v -deffnm prod_06 -nsteps 1250000",
        "                      :-) GROMACS - gmx mdrun, 2026.3 (-:",
        "Compiled SIMD is AVX2_256, but CPU also supports AVX_512 (see log).",
        "1 GPU selected for this run.",
        "Using 1 MPI thread",
        "Using 12 OpenMP threads",
        "Performance:      103.425        0.232",
    ]
    card, rows = extract(lines)
    assert card == ("NVIDIA L4", 23034), card
    assert [r["segment"] for r in rows] == ["prod_05", "prod_06"], rows
    assert [r["ns_per_day"] for r in rows] == ["113.154", "103.425"], rows
    assert all(r["openmp_threads"] == "12" and r["ranks"] == "1" for r in rows)
    assert statistics.median([float(r["ns_per_day"]) for r in rows]) == 108.2895
    assert extract([l for l in lines if "NVIDIA" not in l])[0] is None
    for bad, why in (([{"segment": "prod_05", "ranks": "1"},
                       {"segment": "prod_06", "ranks": "4"}], "ranks"),):
        try:
            one(bad, why, "probe")
        except SystemExit as exc:
            assert why in str(exc), exc
        else:
            raise AssertionError("two values for %s passed" % why)
    probe = os.path.join(tempfile.mkdtemp(), "nb.ipynb")
    json.dump({"cells": [{"outputs": [{"text": [l + "\n" for l in lines]}]}]},
              io.open(probe, "w", encoding="utf8"))
    assert extract(outputs(probe))[0] == ("NVIDIA L4", 23034)
    print("self-test passed: %d segments, equilibration excluded, card read, "
          "disagreement refused" % len(rows))


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        self_test()
        raise SystemExit(0)
    doc = build_record()
    path = os.path.join(RESULTS, "compute_environment.json")
    json.dump(doc, io.open(path, "w", encoding="utf8"), indent=2)
    print("written", os.path.relpath(path, ROOT))
    print("card    %s, %d MiB, %d GB, %d per run" % (
        doc["simulation"]["gpu"], doc["simulation"]["gpu_memory_MiB"],
        doc["simulation"]["gpu_memory_GB"], doc["simulation"]["gpus_per_run"]))
    print("gromacs %s, SIMD %s, %d rank, %d OpenMP threads" % (
        doc["simulation"]["gromacs"], doc["simulation"]["simd"],
        doc["simulation"]["ranks"], doc["simulation"]["openmp_threads"]))
    for build, t in sorted(doc["throughput"].items()):
        print("%-10s median %.1f ns/day, %.1f to %.1f over %d of %d segments, "
              "%.1f hours at the median" % (
                  build, t["ns_per_day_median"], t["ns_per_day_min"],
                  t["ns_per_day_max"], t["segments_recorded"], t["segments_in_run"],
                  t["hours_at_median"]))
    print("analysis %s on %s" % (doc["analysis"]["python"], doc["analysis"]["platform"]))
