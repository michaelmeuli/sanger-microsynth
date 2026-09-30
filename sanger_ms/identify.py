"""Species identification of a trimmed read by BLAST against NCBI.

Only the quality-trimmed sequence is sent (never trace files or filenames).
"""
from __future__ import annotations

from dataclasses import dataclass

from Bio.Blast import NCBIWWW, NCBIXML

# Locus -> NCBI database. 16S has a curated RefSeq database; everything else
# (hsp65, unknown) goes against core_nt, which is slower.
BLAST_DB = {"16S": "16S_ribosomal_RNA", "hsp65": "core_nt", None: "core_nt"}
MIN_BLAST_LENGTH = 100


@dataclass
class Hit:
    title: str
    identity: float  # percent over the aligned region
    coverage: float  # percent of the query covered by the alignment
    evalue: float

    @property
    def confidence(self) -> str:
        if self.identity >= 99.0 and self.coverage >= 90:
            return "high"
        if self.identity >= 97.0 and self.coverage >= 80:
            return "medium"
        return "low"


def blast_read(seq: str, locus: str | None, n_hits: int = 3) -> list[Hit]:
    if len(seq) < MIN_BLAST_LENGTH:
        return []
    handle = NCBIWWW.qblast("blastn", BLAST_DB.get(locus, "core_nt"), seq, hitlist_size=n_hits)
    record = NCBIXML.read(handle)
    hits = []
    for aln in record.alignments[:n_hits]:
        hsp = aln.hsps[0]
        hits.append(Hit(
            title=aln.title.split(" ", 1)[-1][:110],
            identity=100 * hsp.identities / hsp.align_length,
            coverage=100 * (hsp.query_end - hsp.query_start + 1) / len(seq),
            evalue=hsp.expect,
        ))
    return hits
