"""AB1 parsing and quality trimming (same trimming as mlsa-kansasii)."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

from Bio import SeqIO

_16S_RE = re.compile(r"16s|27f|1492r|mbak[-_ ]?14|(?<![a-z0-9])r?(?:259|264)r?(?![a-z0-9])", re.I)
_HSP65_RE = re.compile(r"hsp|hps65|65\s*kda|tb[-_ ]?1[12]w?", re.I)


def guess_locus(filename: str) -> str | None:
    if _HSP65_RE.search(filename):
        return "hsp65"
    if _16S_RE.search(filename):
        return "16S"
    return None


def mott_trim(qualities: list[int], error_threshold: float = 0.05) -> tuple[int, int]:
    """Modified-Mott trimming: the window maximising sum(threshold - 10**(-q/10)).
    Returns half-open (start, end); (0, 0) if no usable window."""
    best_sum = cur_sum = 0.0
    best_start = best_end = cur_start = 0
    for i, q in enumerate(qualities):
        val = error_threshold - 10 ** (-q / 10)
        if cur_sum <= 0:
            cur_start, cur_sum = i, val
        else:
            cur_sum += val
        if cur_sum > best_sum:
            best_sum, best_start, best_end = cur_sum, cur_start, i + 1
    return best_start, best_end


def mixed_peak_fraction(record: Any, start: int, end: int, min_ratio: float = 0.25) -> float | None:
    """Same as mlsa-kansasii mlsa.sanger_io.mixed_peak_fraction (keep in sync): a second template in the
    PCR (mixed culture, cross-contamination) shows as secondary peaks at the base-call positions.
    None if the trace lacks the channels or call positions."""
    raw = cast("dict[str, Any]", record.annotations.get("abif_raw", {}))
    order = raw.get("FWO_1")
    locs = raw.get("PLOC2")
    if order is None or locs is None or end <= start:
        return None
    if isinstance(order, bytes):
        order = order.decode()
    try:
        channels = {base: raw[f"DATA{i}"] for base, i in zip(order, (9, 10, 11, 12))}
    except KeyError:
        return None
    n_mixed = 0
    for i in range(start, end):
        pos = min(locs[i], len(channels["A"]) - 1)
        top, second = sorted((channels[b][pos] for b in "ACGT"), reverse=True)[:2]
        n_mixed += top > 0 and second / top >= min_ratio
    return n_mixed / (end - start)


@dataclass
class Read:
    name: str
    locus: str | None
    raw_length: int
    trim_start: int
    trim_end: int
    trimmed_seq: str
    mean_q: float
    frac_q20: float
    n_count: int
    qualities: list[int]
    mixed_fraction: float | None = None  # share of base calls with a secondary peak >= 25% of the primary

    @property
    def trimmed_length(self) -> int:
        return self.trim_end - self.trim_start


def load_read(path: Path, min_length: int = 100) -> Read | None:
    """Parse and trim an .ab1 trace; None if unparsable or without qualities."""
    try:
        rec = SeqIO.read(path, "abi")
    except Exception:
        return None
    quals = rec.letter_annotations.get("phred_quality")
    if not quals:
        return None
    start, end = mott_trim(quals)
    seq = str(rec.seq[start:end]).upper()
    tq = quals[start:end]
    return Read(
        name=path.name, locus=guess_locus(path.name), raw_length=len(quals),
        trim_start=start, trim_end=end, trimmed_seq=seq,
        mean_q=sum(tq) / len(tq) if tq else 0.0,
        frac_q20=sum(q >= 20 for q in tq) / len(tq) if tq else 0.0,
        n_count=seq.count("N"), qualities=list(quals),
        mixed_fraction=mixed_peak_fraction(rec, start, end),
    )
