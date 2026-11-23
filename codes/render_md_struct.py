import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "results", "md")
FRAMES = os.path.join(MD, "frames")
BUILDS = ("canonical", "mixed")

GREY_A = "0xBFBFBF"
GREY_B = "0x8C8C8C"
RED = "0xB2182B"
DARK = "0x4D4D4D"
PALE = "0x8A8A8A"
POCKET = "0xD8D8D8"
COL = {"canonical": DARK, "mixed": RED}
WIDTH, HEIGHT = 1700, 1400


def say(msg):
    sys.stdout.write(msg + "\n")


_KEEP = {b"IHDR", b"PLTE", b"IDAT", b"IEND", b"tRNS", b"gAMA",
         b"cHRM", b"sRGB", b"iCCP", b"pHYs", b"bKGD", b"sBIT"}


def strip_credentials(path):
    import struct
    d = open(path, "rb").read()
    if d[:8] != b"\x89PNG\r\n\x1a\n":
        raise SystemExit("STOP: %s is not a PNG, so its metadata cannot be stripped"
                         % path)
    out, i, dropped = [d[:8]], 8, []
    while i < len(d):
        if i + 8 > len(d):
            raise SystemExit("STOP: %s is truncated at byte %d" % (path, i))
        ln = struct.unpack(">I", d[i:i + 4])[0]
        typ = d[i + 4:i + 8]
        end = i + 12 + ln
        if end > len(d):
            raise SystemExit("STOP: %s chunk %s overruns the file"
                             % (path, typ.decode("latin1")))
        if typ in _KEEP:
            out.append(d[i:end])
        else:
            dropped.append(typ.decode("latin1"))
        i = end
        if typ == b"IEND":
            break
    if not any(c[4:8] == b"IEND" for c in out[1:]):
        raise SystemExit("STOP: %s has no IEND chunk" % path)
    if dropped:
        open(path, "wb").write(b"".join(out))
    return dropped


AMINO = set("""ALA ARG ASN ASP CYS GLN GLU GLY HIS HSD HSE HSP ILE LEU LYS MET PHE
PRO SER THR TRP TYR VAL""".split())


def frame(build, which):
    return os.path.join(FRAMES, build, "frame_%s.pdb" % which)


def cofactor_resn(path):
    seen = set()
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            name = line[17:21].strip()
            if name and name not in AMINO:
                seen.add(name)
    if len(seen) != 1:
        raise SystemExit("STOP: %s holds %d residue names that are not amino acids, "
                         "%s. One is expected, the cofactor."
                         % (path, len(seen), sorted(seen)))
    return sorted(seen)[0]


def _refuse_marked_canvas(c):
    try:
        from PIL import Image
    except ImportError:
        raise SystemExit("STOP: Pillow is needed to check that this PyMOL leaves an empty "
                         "canvas empty. Install it into the interpreter running this.")
    import tempfile
    probe = os.path.join(tempfile.mkdtemp(), "probe.png")
    c.delete("all")
    c.ray(WIDTH, HEIGHT)
    c.png(probe, dpi=400)
    marked = Image.open(probe).convert("RGBA").getbbox() is not None
    os.remove(probe)
    os.rmdir(os.path.dirname(probe))
    if marked:
        raise SystemExit("STOP: this PyMOL puts its own text on an empty canvas, so every "
                         "render it makes carries that text. Use a licensed build or the "
                         "open source one.")


def start_pymol():
    try:
        import pymol2
    except ImportError:
        return None
    p = pymol2.PyMOL()
    p.start()
    c = p.cmd
    c.set("ray_opaque_background", 0)
    c.set("orthoscopic", 1)
    c.set("antialias", 2)
    c.set("ray_trace_mode", 0)
    c.set("cartoon_transparency", 0.0)
    c.set("ambient", 0.22)
    c.set("specular", 0.15)
    c.set("direct", 0.55)
    c.set("reflect", 0.35)
    c.set("depth_cue", 0)
    c.set("ray_shadows", 0)
    c.set("cartoon_sampling", 12)
    c.set("surface_quality", 1)
    c.set("stick_radius", 0.20)
    c.set("sphere_scale", 0.38)
    _refuse_marked_canvas(c)
    return p


def render_systems(p, out):
    c = p.cmd
    c.delete("all")
    f = frame("mixed", "first")
    if not os.path.isfile(f):
        say("skipped %s: no %s" % (os.path.basename(out), f))
        return False
    cof = cofactor_resn(f)
    c.load(f, "mixed")
    c.hide("everything")
    c.show("cartoon", "polymer")
    c.color(GREY_A, "chain A")
    c.color(RED, "chain B")
    c.set("cartoon_transparency", 0.18)
    c.show("spheres", "resn %s and not elem H" % cof)
    c.color(DARK, "resn %s" % cof)
    c.color("atomic", "resn %s and not elem C" % cof)
    c.orient("polymer")
    c.turn("x", -20)
    c.zoom("polymer", 2)
    c.ray(WIDTH, HEIGHT)
    c.png(out, dpi=400)
    return os.path.isfile(out)


def kabsch(mobile, target):
    import numpy as np
    mc, tc = mobile.mean(0), target.mean(0)
    H = (mobile - mc).T @ (target - tc)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, mc, tc


def read_pdb(path):
    import numpy as np
    keys, xyz, lines = [], [], []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            keys.append((line[21], int(line[22:26]), line[12:16].strip(),
                         line[17:21].strip()))
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        lines.append(line)
    return keys, np.asarray(xyz), lines


def scored_copy(build):
    import json
    p = os.path.join(MD, "%s_report.json" % build)
    if not os.path.isfile(p):
        raise SystemExit("STOP: no %s, so which cofactor copy was measured is not "
                         "known and the picture would be of an arbitrary one." % p)
    ch, resid = json.load(open(p))["roles"]["scored_copy"]
    return ch, int(resid)


def superposed_last(build, work):
    import numpy as np
    ch, resid = scored_copy(build)
    k0, x0, _ = read_pdb(frame(build, "first"))
    k1, x1, lines = read_pdb(frame(build, "last"))
    if k0 != k1:
        raise SystemExit("STOP: the two frames of %s do not hold the same atoms in the "
                         "same order, so one cannot be put on the other." % build)
    sub = [i for i, (c, _r, _a, _n) in enumerate(k0) if c == ch]
    if not sub:
        raise SystemExit("STOP: %s has no chain %s, which is the subunit the measured "
                         "cofactor sits on." % (frame(build, "first"), ch))
    R, mc, tc = kabsch(x1[sub], x0[sub])
    moved = (x1 - mc) @ R.T + tc
    pick = [i for i, (c, r, a, _n) in enumerate(k0)
            if c == ch and r == resid and not a.startswith("H")]
    if not pick:
        raise SystemExit("STOP: no heavy atoms of residue %d in chain %s of %s."
                         % (resid, ch, frame(build, "first")))
    rmsd = float(np.sqrt(((moved[pick] - x0[pick]) ** 2).sum(1).mean()))
    out = os.path.join(work, "%s_last_fitted.pdb" % build)
    with open(out, "w") as fh:
        j = 0
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = moved[j]
                fh.write("%s%8.3f%8.3f%8.3f%s" % (line[:30], x, y, z, line[54:]))
                j += 1
            else:
                fh.write(line)
    return out, rmsd


def render_cofactor(p, out):
    import json
    import tempfile
    c = p.cmd
    need = [(b, w) for b in BUILDS for w in ("first", "last")]
    missing = [frame(b, w) for b, w in need if not os.path.isfile(frame(b, w))]
    if missing:
        say("skipped %s: run kynu_MD_3c_frames.ipynb first, no %s"
            % (os.path.basename(out),
               ", ".join(os.path.join(os.path.basename(os.path.dirname(m)),
                                      os.path.basename(m)) for m in missing)))
        return False
    drift = json.load(open(os.path.join(MD, "cofactor_drift.json")))
    work = tempfile.mkdtemp(prefix="kynu_fit_")
    tiles = []
    for build in BUILDS:
        ch, resid = scored_copy(build)
        fitted, rmsd = superposed_last(build, work)
        reported = float(drift[build]["last_tenth"])
        if not 0.5 * reported <= rmsd <= 2.0 * reported:
            raise SystemExit("STOP: the %s cofactor is %.2f A from its start in the "
                             "last frame and part 3b reports %.2f A over the last "
                             "tenth of the run. The two conventions differ, but not "
                             "by this much, so the picture and the series are not of "
                             "the same run." % (build, rmsd, reported))
        say("   %-10s cofactor %.2f A from its start at the end, against %.2f A "
            "reported over the last tenth" % (build, rmsd, reported))
        c.delete("all")
        c.load(frame(build, "first"), "pose_start")
        c.load(fitted, "pose_end")
        sel = "chain %s and resi %d and not elem H" % (ch, resid)
        c.hide("everything")
        c.show("cartoon", "pose_start and chain %s and polymer" % ch)
        c.color(GREY_A, "pose_start and chain %s" % ch)
        c.set("cartoon_transparency", 0.45, "pose_start")
        pocket = ("byres (pose_start and polymer within 5 of (pose_start and %s)) "
                  "and not (%s) and not elem H" % (sel, sel))
        c.show("sticks", pocket)
        c.color(POCKET, pocket)
        c.set("stick_radius", 0.045, "pose_start")
        c.show("sticks", "pose_start and " + sel)
        c.set_bond("stick_radius", 0.17, "pose_start and " + sel)
        c.show("sticks", "pose_end and " + sel)
        c.color(PALE, "pose_start and " + sel)
        c.color(RED, "pose_end and " + sel)
        c.set("stick_radius", 0.24, "pose_end")
        c.orient("(pose_start or pose_end) and " + sel)
        c.zoom("(pose_start or pose_end) and " + sel, 4.5)
        tile = os.path.join(work, "%s.png" % build)
        c.ray(int(WIDTH * 0.66), HEIGHT)
        c.png(tile, dpi=400)
        tiles.append(tile)
    try:
        from PIL import Image
    except ImportError:
        say("skipped %s: the two halves need PIL to be put side by side"
            % os.path.basename(out))
        return False
    ims = []
    for t in tiles:
        im = Image.open(t).convert("RGBA")
        box = im.getbbox()
        ims.append(im.crop(box) if box else im)
    h = max(i.height for i in ims)
    pad = int(0.03 * sum(i.width for i in ims))
    canvas = Image.new("RGBA", (sum(i.width for i in ims) + pad, h), (0, 0, 0, 0))
    x = 0
    for im in ims:
        canvas.paste(im, (x, (h - im.height) // 2), im)
        x += im.width + pad
    canvas.save(out)
    return os.path.isfile(out)


def main():
    os.makedirs(MD, exist_ok=True)
    p = start_pymol()
    if p is None:
        say("pymol is not available here, so neither picture was made. "
            "pip install pymol-open-source, where a wheel exists for the platform.")
        return 0
    made = []
    for fn, name in ((render_systems, "render_systems.png"),
                     (render_cofactor, "render_cofactor.png")):
        out = os.path.join(MD, name)
        if fn(p, out):
            made.append(out)
            gone = strip_credentials(out)
            say("wrote %s%s" % (out, ", dropped %s" % ", ".join(gone) if gone else ""))
    p.stop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
