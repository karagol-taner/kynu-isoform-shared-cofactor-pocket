import time, urllib.error, urllib.request

ACCESSION = "Q16719"
FAMILY = "IPR010111"
BLOCKS = (("retained", 1, 301), ("substituted", 302, 307), ("deleted", 308, 465))
INTERFACE_IN_LOST_SPAN = (303, 305, 306, 310, 313, 314, 316, 318, 332, 333)
CROSS_SUBUNIT_COFACTOR = (305, 332, 333)
RETRIES = 4
BACKOFF = 5

GAP_OPEN = -11
GAP_EXTEND = -1
MATRIX = "BLOSUM62"
MODE = "global"

CHAIN_LENGTH = BLOCKS[-1][2]


def get(url, accept=None):
    req = urllib.request.Request(url, headers={"Accept": accept} if accept else {})
    last = None
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            code = getattr(exc, "code", None)
            if code is not None and code not in (429, 500, 502, 503, 504):
                raise
            last = exc
            time.sleep(BACKOFF * (attempt + 1))
    raise SystemExit("%s unreachable after %d attempts: %s" % (url, RETRIES, last))


def block_of(position):
    for name, lo, hi in BLOCKS:
        if lo <= position <= hi:
            return name
    raise SystemExit("position %d is outside the sequence" % position)


def aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load(MATRIX)
    a.open_gap_score = GAP_OPEN
    a.extend_gap_score = GAP_EXTEND
    a.mode = MODE
    return a


def human_to_orthologue(human, seq, a=None):
    alignment = (a or aligner()).align(human, seq)[0]
    mapping = {}
    for (hs, he), (os_, oe) in zip(alignment.aligned[0], alignment.aligned[1]):
        for offset in range(he - hs):
            mapping[hs + offset + 1] = seq[os_ + offset]
    return mapping
