#!/usr/bin/env python3
"""Build a barcode -> consensus-genotype -> protein-consequence table from
Pacybara extraction output.

Inputs
  1. the Pacybara parameter file   (amplicon sequence + ORFSTART/ORFEND/BARCODE)
  2. bcExtract_combo.fastq.gz      (one barcode record per read)
  3. genoExtract.csv.gz            (readID, "var:qual;var:qual;..." per read)

Output
  barcode_genotype_table.csv - one row per real barcode:
    barcode, reads, geno, frac_exact_genotype, min_variant_support,
    aa_construct, aa_fl_precursor, aa_legacy, kind

Method
  Variant coordinates in genoExtract are ORF-relative (coordinate 1 == ORFSTART),
  and span the whole amplicon - negative coordinates upstream of the ORF and
  coordinates past the ORF end downstream. Only variants inside the ORF window
  can change the protein, so the genotype is restricted to that window before the
  per-barcode mode is taken. Two fractions are reported because they answer
  different questions:

    min_variant_support - lowest frequency of any consensus variant among the
      barcode's reads. This is the purity measure.
    frac_exact_genotype - share of reads whose ORF genotype matches the consensus
      exactly. Lower, because per-read sequencing noise ADDS spurious variants.

  Sanity-check both: a low support means a genuinely mixed or miscalled clone,
  whereas a low exact fraction alone just means noisy reads.

Numbering
  The same change is annotated three times, differing only in which residue
  numbering it is expressed in:

    aa_construct    the supplied ORF, residue 1 = the initiator Met
    aa_fl_precursor the full-length precursor, i.e. with the construct's
                    internal deletion restored (+--bdd-offset for residues at
                    or past --bdd-junction-aa). This is HGVS p. numbering.
    aa_legacy       the mature protein, --signal-peptide residues lower again.
                    For F8 this is the clinical/literature convention.

  Worked example, FVIII BDD defaults: Y1018F -> Y1699F -> Y1680F.
  Pass --no-legacy for a target with no internal deletion or signal peptide;
  all three columns then collapse to construct numbering.

Example
  python3 pacybara_barcode_genotypes.py \\
      --params 20260831_factor8_multistep_stds.txt \\
      --bc bcExtract_combo.fastq.gz \\
      --geno genoExtract.csv.gz \\
      --out barcode_genotype_table.csv --summary
"""
import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict

# --------------------------------------------------------------------------
# genetic code
# --------------------------------------------------------------------------
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[i] for i, (a, b, c) in enumerate(
    (a, b, c) for a in _B for b in _B for c in _B)}

POS_RE = re.compile(r"^(-?\d+)")
SNV_RE = re.compile(r"^(\d+)([ACGT])>([ACGT])$")
DEL_RE = re.compile(r"^(\d+)del$")
INS_RE = re.compile(r"^(\d+)ins([ACGT]+)$")


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


# --------------------------------------------------------------------------
# parameter file
# --------------------------------------------------------------------------
def read_params(path):
    """Return dict with amplicon, orf_start, orf_end, barcode_len, title."""
    args, seq, in_amp = {}, [], False
    with open(path) as fh:
        for line in fh:
            t = line.strip()
            if t.startswith("#BEGIN AMPLICON"):
                in_amp = True
                continue
            if t.startswith("#END AMPLICON"):
                in_amp = False
                continue
            if in_amp:
                if not t.startswith(">") and t:
                    seq.append(t)
                continue
            if t.startswith("#") or "=" not in t:
                continue
            k, _, v = t.partition("=")
            args[k.strip()] = v.strip()
    missing = [k for k in ("ORFSTART", "ORFEND") if k not in args]
    if missing:
        sys.exit(f"ERROR: {path} is missing {', '.join(missing)}")
    if not seq:
        sys.exit(f"ERROR: no amplicon sequence found in {path}")
    amp = "".join(seq)
    start, end = int(args["ORFSTART"]), int(args["ORFEND"])
    if not 1 <= start < end <= len(amp):
        sys.exit(f"ERROR: ORFSTART/ORFEND ({start}-{end}) outside the "
                 f"{len(amp)} bp amplicon")
    bc = args.get("BARCODE", "")
    return dict(amplicon=amp, orf_start=start, orf_end=end,
                orf=amp[start - 1:end].upper(),
                barcode_len=len(bc) if bc else None,
                title=args.get("TITLE", "(untitled)"))


# --------------------------------------------------------------------------
# variant application
# --------------------------------------------------------------------------
def apply_variants(orf, variants):
    """Apply a variant list to the ORF. Insertions are placed after the
    given position. Raises ValueError on an unparsable variant or a
    reference-base mismatch."""
    sub, dele, ins = {}, set(), defaultdict(str)
    for v in variants:
        m = SNV_RE.match(v)
        if m:
            p, ref, alt = int(m.group(1)), m.group(2), m.group(3)
            if orf[p - 1] != ref:
                raise ValueError(
                    f"{v}: ORF position {p} is {orf[p-1]}, not {ref}")
            sub[p] = alt
            continue
        m = DEL_RE.match(v)
        if m:
            dele.add(int(m.group(1)))
            continue
        m = INS_RE.match(v)
        if m:
            ins[int(m.group(1))] += m.group(2)
            continue
        raise ValueError(f"unparsable variant {v!r}")
    out = []
    for p in range(1, len(orf) + 1):
        if p not in dele:
            out.append(sub.get(p, orf[p - 1]))
        if p in ins:
            out.append(ins[p])
    return "".join(out)


def describe(wt, mut):
    """Compare two proteins and return (template, positions, kind).

    `template` is a format string whose {0}, {1}... placeholders stand in for
    residue positions, and `positions` are those positions in construct
    numbering. Keeping them separate is what lets the same change be rendered
    in a different numbering scheme: substituting into a pre-formatted string
    would only catch the first number and quietly corrupt any annotation
    carrying two (a range deletion, an insertion). Render with `render()`."""
    if wt == mut:
        return "(synonymous)", [], "synonymous"
    p = 0
    while p < min(len(wt), len(mut)) and wt[p] == mut[p]:
        p += 1
    s = 0
    while s < min(len(wt), len(mut)) - p and wt[-1 - s] == mut[-1 - s]:
        s += 1
    w, m = wt[p:len(wt) - s], mut[p:len(mut) - s]
    aa = p + 1
    if len(w) == len(m) == 1:
        return f"{w}{{0}}{m}", [aa], ("nonsense" if m == "*" else "missense")
    if not m and w:
        if len(w) == 1:
            return f"{w}{{0}}del", [aa], "in-frame deletion"
        return (f"{w[0]}{{0}}_{w[-1]}{{1}}del", [aa, aa + len(w) - 1],
                "in-frame deletion")
    if not w and m:
        return (f"{wt[p-1]}{{0}}_{wt[p]}{{1}}ins{m}", [p, p + 1],
                "in-frame insertion")
    if len(wt) != len(mut):
        return f"{(w or '?')[0]}{{0}}fs", [aa], "frameshift"
    return f"{w}{{0}}{m}", [aa], f"multi-residue ({len(w)}aa)"


def render(template, positions, renumber=None):
    """Fill a describe() template, mapping every position through `renumber`."""
    if not positions:
        return template
    f = renumber or (lambda x: x)
    return template.format(*[f(x) for x in positions])


# --------------------------------------------------------------------------
# input readers
# --------------------------------------------------------------------------
def read_barcodes(path, bc_len):
    """Return (Counter of barcode->reads, dict readID->barcode)."""
    counts, rid2bc = Counter(), {}
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            r = i & 3
            if r == 0:
                rid = line[1:].split(None, 1)[0]
            elif r == 1:
                s = line.strip()
                if len(s) == bc_len:
                    counts[s] += 1
                    rid2bc[rid] = s
    return counts, rid2bc


def pick_min_reads(counts):
    """Largest multiplicative gap in the reads-per-barcode spectrum.

    A clean library is bimodal - a spike of 1-read extraction artifacts, an
    empty band, then the real clones - so the cutoff is the first abundance
    above the widest gap. Returns (cutoff, gap_lo, gap_hi)."""
    vals = sorted(set(counts.values()))
    if len(vals) < 2:
        return 1, None, None
    best, lo, hi = 0.0, None, None
    for a, b in zip(vals, vals[1:]):
        if b / a > best:
            best, lo, hi = b / a, a, b
    return hi, lo, hi


def read_genotypes(path, rid2bc, orf_len):
    """Per barcode: Counter of ORF-restricted genotypes, read total, and
    per-variant counts."""
    genos = defaultdict(Counter)
    varfreq = defaultdict(Counter)
    nreads = Counter()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            rid, _, rest = line.rstrip("\n").partition(",")
            bc = rid2bc.get(rid.strip('"'))
            if bc is None:
                continue
            keep = []
            g = rest.strip().strip('"')
            if g and g != "=":
                for mut in g.split(";"):
                    name = mut.rsplit(":", 1)[0]
                    m = POS_RE.match(name)
                    if m and 1 <= int(m.group(1)) <= orf_len:
                        keep.append(name)
            keep.sort(key=lambda x: int(POS_RE.match(x).group(1)))
            nreads[bc] += 1
            genos[bc][tuple(keep)] += 1
            for name in keep:
                varfreq[bc][name] += 1
    return genos, varfreq, nreads


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default="parameters.txt",
                    help="Pacybara parameter file (amplicon + ORF coords)")
    ap.add_argument("--bc", default="bcExtract_combo.fastq.gz",
                    help="barcode FASTQ from Pacybara extraction")
    ap.add_argument("--geno", default="genoExtract.csv.gz",
                    help="genotype CSV from Pacybara extraction")
    ap.add_argument("--out", default="barcode_genotype_table.csv",
                    help="output CSV")
    ap.add_argument("--min-reads", default="31",
                    help="reads required to call a barcode real, or 'auto' to "
                         "take the widest gap in the abundance spectrum "
                         "(default: 31)")
    ap.add_argument("--bc-len", type=int, default=None,
                    help="barcode length (default: from BARCODE in --params)")
    ap.add_argument("--signal-peptide", type=int, default=19,
                    help="residues removed for mature numbering (default: 19)")
    ap.add_argument("--bdd-junction-aa", type=int, default=987,
                    help="construct residue where the internal deletion ends "
                         "(default: 987, FVIII BDD)")
    ap.add_argument("--bdd-offset", type=int, default=681,
                    help="residues deleted, added to positions at or past the "
                         "junction (default: 681, FVIII BDD)")
    ap.add_argument("--no-legacy", action="store_true",
                    help="skip reference-protein renumbering; all three aa_* "
                         "columns then repeat construct numbering")
    ap.add_argument("--summary", action="store_true",
                    help="also print one line per distinct consensus genotype")
    a = ap.parse_args()

    p = read_params(a.params)
    orf, orf_len = p["orf"], len(p["orf"])
    if orf_len % 3:
        print(f"WARNING: ORF is {orf_len} nt, not a multiple of 3",
              file=sys.stderr)
    wt_prot = translate(orf)
    log = lambda *m: print(*m, file=sys.stderr)
    log(f"{p['title']}: amplicon {len(p['amplicon'])} bp, "
        f"ORF {p['orf_start']}-{p['orf_end']} "
        f"({orf_len} nt / {len(wt_prot)} aa)")
    if not wt_prot.startswith("M"):
        log("WARNING: ORF does not begin with ATG")
    if "*" in wt_prot[:-1]:
        log(f"WARNING: {wt_prot[:-1].count('*')} internal stop codon(s) in ORF")

    bc_len = a.bc_len or p["barcode_len"]
    if not bc_len:
        sys.exit("ERROR: no BARCODE in the parameter file; pass --bc-len")
    counts, rid2bc = read_barcodes(a.bc, bc_len)
    log(f"{len(counts)} distinct {bc_len}bp barcodes over "
        f"{sum(counts.values())} reads")

    if a.min_reads == "auto":
        cutoff, lo, hi = pick_min_reads(counts)
        log(f"--min-reads auto: widest abundance gap {lo} -> {hi}; "
            f"using {cutoff}")
    else:
        cutoff = int(a.min_reads)
        _, lo, hi = pick_min_reads(counts)
        if lo is not None and not (lo < cutoff <= hi):
            log(f"NOTE: widest abundance gap is {lo} -> {hi}; --min-reads "
                f"{cutoff} does not fall in it. Check the distribution.")
    real = {b for b, n in counts.items() if n >= cutoff}
    dropped = sum(n for b, n in counts.items() if n < cutoff)
    log(f"{len(real)} barcodes with >={cutoff} reads "
        f"({len(counts) - len(real)} dropped, {dropped} reads)")
    if not real:
        sys.exit("ERROR: no barcodes passed the cutoff")
    rid2bc = {r: b for r, b in rid2bc.items() if b in real}

    genos, varfreq, nreads = read_genotypes(a.geno, rid2bc, orf_len)
    log(f"genotypes aggregated for {len(genos)} barcodes")

    # --no-legacy is just the identity mapping: no internal deletion restored,
    # no signal peptide removed, so all three columns collapse to construct
    # numbering without any special-casing downstream.
    offset = 0 if a.no_legacy else a.bdd_offset
    sigpep = 0 if a.no_legacy else a.signal_peptide

    def to_fl(aa):
        """Construct residue -> full-length precursor residue."""
        return aa + offset if aa >= a.bdd_junction_aa else aa

    def to_legacy(aa):
        """Construct residue -> mature reference-protein residue."""
        return to_fl(aa) - sigpep

    # resolve each distinct consensus genotype once
    cache = {}
    for bc in genos:
        gt = genos[bc].most_common(1)[0][0]
        if gt in cache:
            continue
        if not gt:
            cache[gt] = ("WT", "WT", "WT", "wild type")
            continue
        try:
            tpl, pos, kind = describe(wt_prot,
                                      translate(apply_variants(orf, gt)))
        except ValueError as e:
            log(f"WARNING: {bc} {';'.join(gt)}: {e}")
            cache[gt] = ("(error)", "(error)", "(error)", f"unparsed: {e}")
            continue
        cache[gt] = (render(tpl, pos), render(tpl, pos, to_fl),
                     render(tpl, pos, to_legacy), kind)

    rows = []
    for bc in sorted(genos, key=lambda b: (-nreads[b], b)):
        gt, hits = genos[bc].most_common(1)[0]
        ac, afl, al, kind = cache[gt]
        support = (min(varfreq[bc][m] / nreads[bc] for m in gt) if gt else None)
        rows.append((bc, nreads[bc], ";".join(gt) or "(none)",
                     hits / nreads[bc], support, ac, afl, al, kind))

    with open(a.out, "w") as fh:
        fh.write("barcode,reads,geno,frac_exact_genotype,min_variant_support,"
                 "aa_construct,aa_fl_precursor,aa_legacy,kind\n")
        for bc, n, g, frac, sup, ac, afl, al, kind in rows:
            fh.write("%s,%d,%s,%.4f,%s,%s,%s,%s,%s\n" % (
                bc, n, g, frac, "" if sup is None else "%.4f" % sup,
                ac, afl, al, kind))
    log(f"wrote {a.out} ({len(rows)} barcodes)")

    if a.summary:
        agg = defaultdict(lambda: [0, 0])
        for bc, n, g, frac, sup, ac, afl, al, kind in rows:
            e = agg[(g, ac, afl, al, kind)]
            e[0] += 1
            e[1] += n
        print("%-28s %-12s %-14s %-12s %6s %10s  %s" % (
            "nt genotype", "construct", "fl precursor", "legacy", "n_bc",
            "reads", "kind"))
        for (g, ac, afl, al, kind), (nbc, rd) in sorted(
                agg.items(), key=lambda x: -x[1][0]):
            print("%-28s %-12s %-14s %-12s %6d %10d  %s" % (
                g, ac, afl, al, nbc, rd, kind))


if __name__ == "__main__":
    main()
