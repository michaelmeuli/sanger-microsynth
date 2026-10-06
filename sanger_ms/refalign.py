"""Species identification of a Sanger read by plain alignment to one reference
genome per species (no MLSA, no locus extraction from the references).

Each read is seeded against every reference genome with exact 15-mers (both
strands), the best-supported genomic window is aligned end to end with
Biopython's PairwiseAligner, and the reference that differs at the fewest
positions of the read is the closest species. The read and all references are
then drawn as one read-anchored alignment plus a pairwise difference matrix.

Self-contained (numpy, Biopython, matplotlib). sanger-microsynth carries a copy
as sanger_ms/refalign.py; keep the two in sync.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

K = 15
PAD = 40  # bases of reference added on each side of the seeded window
MIN_SEEDS = 5
_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _LUT[_b + 32] = _i
_SENTINEL = np.uint32(0xFFFFFFFF)


def _kmer_codes(seq: str | bytes) -> npt.NDArray[np.uint32]:
    """2-bit code of every K-mer (uint32; sentinel where it holds a non-ACGT base)."""
    raw = seq.encode() if isinstance(seq, str) else seq
    b = _LUT[np.frombuffer(raw, dtype=np.uint8)]
    n = len(b) - K + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint32)
    bad = np.cumsum(np.concatenate([[0], b == 255]))
    invalid = (bad[K:] - bad[:-K]) > 0
    code = np.zeros(n, dtype=np.uint32)
    for j in range(K):
        code = (code << np.uint32(2)) | (b[j:n + j] & np.uint8(3)).astype(np.uint32)
    code[invalid] = _SENTINEL
    return code


_COMP = str.maketrans("ACGTN-", "TGCAN-")


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


_NAME_RE = re.compile(
    r"^\S+ Mycobacterium(?:\s+(?P<species>[a-z]+))(?:\s+(?:strain|subsp\.\s+\S+)?\s*(?P<strain>.+?))?"
    r"(?:\s+(?:isolate\b|plasmid\b|chromosome\b|NODE|contig|Scaffold|scaffold).*|,.*)?$")


@dataclass
class Reference:
    accession: str
    species: str
    strain: str
    path: Path
    contigs: dict[str, str] = field(default_factory=dict, repr=False)
    _ids: list[str] = field(default_factory=list, repr=False)
    # offset of each contig in _codes
    _starts: npt.NDArray[np.intp] = field(default_factory=lambda: np.empty(0, dtype=np.intp), repr=False)
    # all contigs' k-mer codes, sentinel-separated
    _codes: npt.NDArray[np.uint32] = field(default_factory=lambda: np.empty(0, dtype=np.uint32), repr=False)

    @property
    def label(self) -> str:
        return f"Mycobacterium {self.species} {self.strain}".strip()

    def load(self) -> None:
        if not self.contigs:
            self.contigs = {r.id: str(r.seq).upper() for r in SeqIO.parse(self.path, "fasta")}
            gap = np.full(K, _SENTINEL, dtype=np.uint32)
            parts: list[npt.NDArray[np.uint32]] = []
            starts: list[int] = []
            pos = 0
            self._ids = []
            for cid, seq in self.contigs.items():
                c = _kmer_codes(seq)
                parts += [c, gap]
                self._ids.append(cid)
                starts.append(pos)
                pos += len(c) + K
            self._starts = np.array(starts, dtype=np.intp)
            self._codes = np.concatenate(parts)


def load_references(ref_dir: Path, species_order: list[str] | None = None) -> list[Reference]:
    """One Reference per *.fasta/*.fna/*.fa in ref_dir, named from the first header."""
    refs = []
    for path in sorted(p for p in ref_dir.iterdir() if p.suffix in {".fasta", ".fna", ".fa"}):
        with open(path) as fh:
            header = fh.readline().strip().lstrip(">")
        accession = re.match(r"(GC[AF]_\d+\.\d)", path.name)
        m = _NAME_RE.match(header)
        refs.append(Reference(
            accession=accession.group(1) if accession else path.stem,
            species=m.group("species") if m else "unknown",
            strain=(m.group("strain") or "").strip() if m else header[:30],
            path=path))
    if species_order:
        refs.sort(key=lambda r: (species_order.index(r.species) if r.species in species_order else 99, r.species))
    return refs


_aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                           open_gap_score=-5, extend_gap_score=-2)
_aligner.end_deletion_score = 0  # the reference window may overhang the read freely


@dataclass
class RefHit:
    """The read laid over one reference: ref_cols[i] is the reference base at read position i
    ('-' = deletion in the reference or not covered, 'N' = ambiguous)."""
    ref: Reference
    ref_cols: str
    contig: str
    ref_start: int
    strand: str
    seeds: int


def _seed_clusters(ref: Reference, read: str) -> list[tuple[int, str, int, int]]:
    """(seed count, contig, diagonal, strand +1/-1) of the best-supported diagonals, best first."""
    out = []
    for strand, q in ((1, read), (-1, revcomp(read))):
        qc = _kmer_codes(q)
        valid = np.flatnonzero(qc != _SENTINEL)
        if not len(valid):
            continue
        uniq, idx, counts = np.unique(qc[valid], return_index=True, return_counts=True)
        keep = counts == 1
        qcodes, qoff = uniq[keep], valid[idx[keep]]
        gc = ref._codes
        pos = np.searchsorted(qcodes, gc)
        hit = np.flatnonzero(qcodes[np.minimum(pos, len(qcodes) - 1)] == gc)
        if len(hit) < MIN_SEEDS:
            continue
        ci = np.searchsorted(ref._starts, hit, side="right") - 1
        rel = hit - ref._starts[ci]  # position within the contig
        offs = qoff[pos[hit]]
        key = ci.astype(np.int64) * (1 << 32) + (rel - offs + (1 << 31))  # (contig, diagonal), sortable
        order = np.argsort(key)
        key = key[order]
        cuts = np.flatnonzero(np.diff(key) > 20) + 1
        for grp in np.split(key, cuts):
            if len(grp) >= MIN_SEEDS:
                med = int(np.median(grp))
                out.append((len(grp), ref._ids[med >> 32], (med & 0xFFFFFFFF) - (1 << 31), strand))
    return sorted(out, reverse=True)[:3]


def align_to_reference(ref: Reference, read: str) -> RefHit | None:
    """Best end-to-end placement of the read on the reference, None if no seeds."""
    ref.load()
    best = None
    for seeds, cid, diag, strand in _seed_clusters(ref, read):
        contig = ref.contigs[cid]
        lo, hi = max(0, diag - PAD), min(len(contig), diag + len(read) + PAD)
        q = read if strand == 1 else revcomp(read)
        aln = _aligner.align(contig[lo:hi], q)[0]
        if best is None or aln.score > best[0]:
            best = (aln.score, aln, cid, lo, strand, seeds)
    if best is None:
        return None
    _, aln, cid, lo, strand, seeds = best
    cols = ["-"] * len(read)  # in orientation of q
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for k in range(q1 - q0):
            cols[q0 + k] = ref.contigs[cid][lo + t0 + k]
    ref_cols = "".join(cols)
    ref_cols = "".join(c if c in "ACGT-" else "N" for c in ref_cols)
    if strand == -1:  # report in the orientation of the read as given
        ref_cols = ref_cols[::-1].translate(_COMP)
    return RefHit(ref, ref_cols, cid, lo, "+" if strand == 1 else "-", seeds)


def _diff(a: str, b: str) -> tuple[int, int]:
    """(differences, compared columns), over columns where both have an ACGT base."""
    m = d = 0
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            m += 1
            d += x != y
    return d, m


@dataclass
class ReadResult:
    name: str
    locus: str | None
    read: str
    hits: list[RefHit]
    diffs: list[int]     # per hit: differences vs the read
    compared: list[int]  # per hit: columns compared
    min_identity: float = 0.99
    min_margin: int = 2

    @property
    def order(self) -> list[int]:
        return sorted(range(len(self.hits)), key=lambda i: (self.diffs[i] / max(self.compared[i], 1), -self.compared[i]))

    @property
    def best(self) -> int | None:
        return self.order[0] if self.hits else None

    def identity(self, i: int) -> float:
        return 1 - self.diffs[i] / self.compared[i] if self.compared[i] else 0.0

    @property
    def status(self) -> str:
        """ok: >= min_identity and clearly closest; ambiguous: >= min_identity but second
        reference within min_margin differences; divergent: closest is < min_identity."""
        if not self.hits:
            return "no_hit"
        o = self.order
        if self.identity(o[0]) < self.min_identity:
            return "divergent"
        if len(o) > 1 and self.diffs[o[1]] - self.diffs[o[0]] < self.min_margin:
            return "ambiguous"
        return "ok"

    @property
    def closest_species(self) -> str:
        best = self.best
        return self.hits[best].ref.species if best is not None else "NA"

    @property
    def runner_up(self) -> int | None:
        return self.order[1] if len(self.hits) > 1 else None


def identify_read(refs: list[Reference], name: str, locus: str | None, seq: str,
                  min_identity: float = 0.99, min_margin: int = 2) -> ReadResult:
    hits = [h for h in (align_to_reference(r, seq) for r in refs) if h is not None]
    dc = [_diff(seq, h.ref_cols) for h in hits]
    return ReadResult(name, locus, seq, hits, [d for d, _ in dc], [c for _, c in dc], min_identity, min_margin)


def pairwise_matrix(res: ReadResult) -> tuple[list[str], npt.NDArray[np.int_], npt.NDArray[np.float64]]:
    """Labels (references then read), difference counts and identity (%) over pairwise-complete columns."""
    rows = [h.ref_cols for h in res.hits] + [res.read]
    labels = [h.ref.label for h in res.hits] + [res.name]
    n = len(rows)
    diff, ident = np.zeros((n, n), int), np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i != j:
                d, m = _diff(rows[i], rows[j])
                diff[i, j], ident[i, j] = d, 100 * (1 - d / m) if m else np.nan
    return labels, diff, ident


# ---------------------------------------------------------------- PDF figure

_CW = 0.602  # advance width of DejaVu Sans Mono in em
_BLOCK = 80
_FS = 7.0


def write_pdf(res: ReadResult, out: Path, title: str | None = None) -> Path:
    """Page 1: closest reference + difference matrix; alignment blocks below, continued on further pages."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Rectangle

    W, H = 842, 595  # A4 landscape in points
    labels, diff, ident = pairwise_matrix(res)
    n = len(labels)
    best = res.best
    names = [h.ref.label for h in res.hits] + [res.name]
    seqs = [h.ref_cols for h in res.hits] + [res.read]
    cw, lh = _CW * _FS, _FS * 1.35
    label_w = max(len(s) for s in names) * cw + 10
    x0 = 30 + label_w

    def new_page() -> tuple[Figure, Axes]:
        fig = plt.figure(figsize=(W / 72, H / 72))
        ax = fig.add_axes((0, 0, 1, 1))
        ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
        return fig, ax

    pages: list[Figure] = []
    fig, ax = new_page()
    pages.append(fig)
    ax.text(30, 30, title or f"{res.name}: alignment to the kansasii-complex reference strains",
            fontsize=12, weight="bold", va="center")
    if res.hits:
        o = res.order
        verdict = (f"Closest: M. {res.closest_species} ({res.hits[o[0]].ref.accession}), "
                   f"{res.identity(o[0]) * 100:.2f}% identity ({res.diffs[o[0]]} of {res.compared[o[0]]} positions differ)")
        if res.runner_up is not None:
            r = res.runner_up
            verdict += f"; next: M. {res.hits[r].ref.species} {res.identity(r) * 100:.2f}% ({res.diffs[r]})"
        verdict += f". Call: {res.status} (needs >= {res.min_identity * 100:g}% and >= {res.min_margin} fewer differences than the next)."
    else:
        verdict = "No reference shares enough k-mers with the read."
    ax.text(30, 50, verdict, fontsize=8, va="center")

    # difference matrix: upper triangle = number of differences, lower = % identity
    cell, my = 38, 75
    mx = 30 + label_w + 18
    vmax_d = max(diff.max(), 1)
    lo_id = np.nanmin(ident) if n > 1 and not np.all(np.isnan(ident)) else 100
    for j in range(n):
        ax.add_patch(Rectangle((mx + j * cell, my), cell - 2, 14, color="#d0d5dd"))
        ax.text(mx + j * cell + cell / 2, my + 7, str(j + 1), fontsize=6.5, ha="center", va="center")
    for i in range(n):
        y: float = my + 16 + i * 14
        ax.text(mx - 24, y + 6, names[i], fontsize=6.5, ha="right", va="center",
                weight="bold" if i == best or i == n - 1 else "normal")
        ax.add_patch(Rectangle((mx - 20, y), 18, 12, color="#d0d5dd"))
        ax.text(mx - 11, y + 6, str(i + 1), fontsize=6.5, ha="center", va="center")
        for j in range(n):
            if i == j:
                continue
            color: tuple[float, float, float]
            txt: str
            if j > i:
                a = diff[i, j] / vmax_d
                color, txt = (1, 1 - 0.55 * a, 1 - 0.55 * a), str(diff[i, j])
            else:
                a = (100 - ident[i, j]) / max(100 - lo_id, 1e-9) if not np.isnan(ident[i, j]) else 0
                color, txt = (1 - 0.5 * a, 1 - 0.5 * a, 1), "" if np.isnan(ident[i, j]) else f"{ident[i, j]:.2f}"
            ax.add_patch(Rectangle((mx + j * cell, y), cell - 2, 12, color=color))
            ax.text(mx + j * cell + cell - 5, y + 6, txt, fontsize=6.5, ha="right", va="center")
    ax.text(mx + n * cell + 8, my + 16 + 6, "upper: differences\nlower: % identity", fontsize=6, va="top", color="#555")

    y = my + 16 + n * 14 + 28
    L = len(res.read)
    consensus = [Counter(s[p] for s in seqs if s[p] in "ACGT").most_common(1) for p in range(L)]
    for b0 in range(0, L, _BLOCK):
        need = lh * (n + 2)
        if y + need > H - 30:
            fig, ax = new_page()
            pages.append(fig)
            y = 30
        b1 = min(b0 + _BLOCK, L)
        for p in range(b0, b1):  # ruler: tick + number every 20 bases
            if (p + 1) % 20 == 0:
                x = x0 + ((p - b0) + (p - b0) // 10) * cw + cw / 2
                ax.text(x, y, str(p + 1), fontsize=5.5, ha="center", va="bottom")
                ax.plot([x, x], [y + 1, y + 4], color="k", lw=0.5)
        y += 6
        for r, s in enumerate(seqs):
            yy = y + r * lh
            ax.text(x0 - 8, yy, names[r], fontsize=_FS, ha="right", va="top",
                    weight="bold" if r == best or r == n - 1 else "normal")
            seg = s[b0:b1]
            ax.text(x0, yy, " ".join(seg[i:i + 10] for i in range(0, len(seg), 10)),
                    fontsize=_FS, family="DejaVu Sans Mono", va="top")
            for k, ch in enumerate(seg):
                p = b0 + k
                if ch in "ACGT" and consensus[p] and ch != consensus[p][0][0]:
                    ax.add_patch(Rectangle((x0 + (k + k // 10) * cw, yy - 0.5), cw, lh - 1, color="#f4a7a7", zorder=0))
            ax.text(x0 + (len(seg) + len(seg) // 10) * cw + 6, yy, str(b1), fontsize=_FS, family="DejaVu Sans Mono", va="top")
        y += lh * (n + 1.4)

    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        for f in pages:
            pdf.savefig(f)
            plt.close(f)
    return out
