"""Species identification of a trimmed read by BLAST against NCBI.

Only the quality-trimmed sequence is sent (never trace files or filenames).
"""
from __future__ import annotations

import threading
from dataclasses import dataclass

from Bio.Blast import NCBIWWW, NCBIXML

# Locus -> NCBI database. 16S: refseq_rna (curated NR_ 16S records; NCBI's 16S_ribosomal_RNA
# database hangs on remote BLAST). hsp65 and unknown loci: core_nt.
BLAST_DB = {"16S": "refseq_rna", "hsp65": "core_nt", None: "core_nt"}
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


BLAST_TIMEOUT = 300  # seconds; NCBI occasionally leaves a request queued for many minutes


def blast_read(seq: str, locus: str | None, n_hits: int = 3, timeout: int = BLAST_TIMEOUT) -> list[Hit]:
    """Top hits for a read; [] if too short, or if NCBI does not answer within
    timeout seconds (the stuck request is abandoned in a daemon thread)."""
    if len(seq) < MIN_BLAST_LENGTH:
        return []
    result: dict = {}

    def work():
        try:
            result["hits"] = _blast(seq, locus, n_hits)
        except Exception as exc:  # network/NCBI errors must not kill the batch
            result["error"] = exc

    t = threading.Thread(target=work, daemon=True)
    t.start()
    t.join(timeout)
    if "hits" not in result:
        print(f"    BLAST failed or timed out: {result.get('error', f'no answer in {timeout}s')}", flush=True)
        return []
    return result["hits"]


def _blast(seq: str, locus: str | None, n_hits: int) -> list[Hit]:
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
