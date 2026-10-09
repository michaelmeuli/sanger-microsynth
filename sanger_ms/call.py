"""Per-read action and per-sample species call (forward and reverse reads paired).

Rules validated against whole-genome species in mlsa-kansasii main4 (hsp65: 78/79 calls correct):
  report   closest >= 98% identical, >= 2 differences fewer than the next species, mixed peaks < 10%
  confirm  closest is the atypical-hsp65 kansasii lineage, or the margin is < 2: add gyrA
  repeat   no reference >= 98% identical, or the trimmed read is < 200 bp
  mixed    >= 10% of base calls carry a secondary peak: repeat from a pure colony
  16S      not species-informative for this complex
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from .refalign import ReadResult
from .sanger_io import Read

MIN_LENGTH = 200
MAX_MIXED = 0.10

REPORT, CONFIRM, REPEAT, MIXED, NOT_INFORMATIVE = (
    "report", "confirm with gyrA", "repeat sequencing", "possible mixed culture: repeat from pure colony",
    "16S: not species-informative for this complex")

_TNR_RE = re.compile(r"(?<!\d)(\d{10})(?!\d)")
_REV_RE = re.compile(r"tb[-_ ]?12|(?<![a-z0-9])12w?(?![a-z0-9])", re.I)
_FWD_RE = re.compile(r"tb[-_ ]?11|(?<![a-z0-9])11(?![a-z0-9])", re.I)


@dataclass
class ReadCall:
    species: str  # "" if none
    action: str
    note: str = ""


def read_call(read: Read, res: ReadResult | None) -> ReadCall:
    if read.locus == "16S":
        return ReadCall("", NOT_INFORMATIVE)
    if res is None or res.best is None or read.trimmed_length < MIN_LENGTH \
            or res.identity(res.best) < res.min_identity:
        return ReadCall("", REPEAT)
    species = res.closest_species
    if read.mixed_fraction is not None and read.mixed_fraction >= MAX_MIXED:
        return ReadCall(species, MIXED)
    if res.atypical:
        return ReadCall(species, CONFIRM, "persicum-like hsp65 type: matches the atypical kansasii genomes")
    if res.status != "ok":
        return ReadCall(species, CONFIRM, "closest species not clearly separated from the next")
    return ReadCall(species, REPORT)


def sample_key(name: str, locus: str | None) -> tuple[str, str]:
    m = _TNR_RE.search(name)
    return (m.group(1) if m else name, locus or "?")


def direction(name: str) -> str:
    """TB-11 = forward, TB-12/12w = reverse primer, from the file name."""
    tail = re.sub(r"^.*?(\d{10})", "", name)
    if _REV_RE.search(tail):
        return "rev"
    return "fwd" if _FWD_RE.search(tail) else "?"


@dataclass
class SampleCall:
    sample: str
    locus: str
    reads: list[str] = field(default_factory=list)
    species: str = ""
    action: str = ""
    note: str = ""


def sample_calls(rows: list[tuple[Read, ReadCall]]) -> list[SampleCall]:
    """One call per (sample, locus) from its forward/reverse/repeat reads: all usable reads must agree."""
    groups: dict[tuple[str, str], list[tuple[Read, ReadCall]]] = {}
    for read, call in rows:
        groups.setdefault(sample_key(read.name, read.locus), []).append((read, call))
    out: list[SampleCall] = []
    for (sample, locus), items in groups.items():
        sc = SampleCall(sample, locus, [f"{r.name} ({direction(r.name)})" for r, _ in items])
        if locus == "16S":
            sc.action = NOT_INFORMATIVE
        else:
            usable = [c for _, c in items if c.species and c.action in (REPORT, CONFIRM)]
            species = sorted({c.species for c in usable})
            if any(c.action == MIXED for _, c in items):
                sc.species = ", ".join(sorted({c.species for _, c in items if c.species}))
                sc.action = MIXED
            elif len(species) > 1:
                sc.species, sc.action = " / ".join(species), REPEAT
                sc.note = "reads disagree"
            elif usable:
                sc.species = species[0]
                sc.action = CONFIRM if any(c.action == CONFIRM for c in usable) else REPORT
                notes = [c.note for c in usable if c.note]
                sc.note = notes[0] if notes else (f"{len(usable)} reads agree" if len(usable) > 1 else "")
            else:
                sc.action = REPEAT
        out.append(sc)
    return out
