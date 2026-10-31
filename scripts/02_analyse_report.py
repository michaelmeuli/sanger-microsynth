"""QC + species ID for Microsynth batches, PDF report, optional email.

Default: every batch folder in data/sanger/microsynth_mail/ without a report yet.
  --batch DIR   analyse one folder of .ab1 files (report is not emailed unless --send)
  --no-blast    skip the NCBI BLAST species identification
  --no-refalign skip the alignment to the kansasii-complex references (the 7 representative genomes)
  --no-send     build the PDF but don't email it
  --gtdb-representatives DIR  reference genomes for all reads (default GTDB_REPRESENTATIVES, the 7 representatives)
  --hsp65-amplicons [DIR]     opt in: align hsp65 reads to the per-genome hsp65 amplicons of the 72 genomes instead
                              (default DIR = HSP65_AMPLICONS). Flags the persicum-like (Korean) kansasii lineage for gyrA
                              confirmation, but needs gyrA in the lab; with the 7 genomes that lineage reads as persicum.
"""
from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the sanger_ms package

from sanger_ms import MAIL_DIR, OUTPUT  # noqa: E402
from sanger_ms.config import REQUIRED_SMTP, load_mail_config  # noqa: E402
from sanger_ms.identify import Hit, blast_read  # noqa: E402
from sanger_ms.kansasii_paths import GTDB_REPRESENTATIVES, HSP65_AMPLICONS, SPECIES  # noqa: E402
from sanger_ms.refalign import ReadResult, Reference, identify_read, load_references, write_pdf  # noqa: E402
from sanger_ms.report import build_pdf  # noqa: E402
from sanger_ms.sanger_io import Read, load_read  # noqa: E402
from sanger_ms.send import send_report  # noqa: E402


MIN_ALIGNED_IDENTITY = 0.90  # below this the read is not in the complex: no alignment PDF


def analyse_batch(batch: Path, blast: bool, refs: list[Reference] | None,
                  pdf_dir: Path, refs_hsp65: list[Reference] | None = None) -> list[tuple[Read, list[Hit], ReadResult | None]]:
    rows: list[tuple[Read, list[Hit], ReadResult | None]] = []
    for ab1 in sorted(batch.rglob("*.ab1")):
        read = load_read(ab1)
        if read is None:
            print(f"  skipped (unparsable or no qualities): {ab1.name}")
            continue
        t0 = time.time()
        hits = blast_read(read.trimmed_seq, read.locus) if blast else []
        top = f"{hits[0].title[:50]} ({hits[0].identity:.1f}%)" if hits else "-"
        print(f"  {ab1.name}: {read.locus}, {read.trimmed_length} bp, {top} [{time.time() - t0:.0f}s]", flush=True)
        ref = None
        if refs and read.trimmed_length >= 100:
            stem = Path(read.name).stem.replace(" ", "_")
            # default: the 7 representative genomes for every locus; with --hsp65-amplicons the hsp65 reads use one
            # amplicon per genome (diversity within species, atypical kansasii hsp65 explicit)
            use = refs_hsp65 if read.locus == "hsp65" and refs_hsp65 else refs
            ref = identify_read(use, Path(read.name).stem, read.locus, read.trimmed_seq)
            if ref.best is not None and ref.identity(ref.best) >= MIN_ALIGNED_IDENTITY:
                write_pdf(ref, pdf_dir / f"{stem}_alignment.pdf")
            print(f"    ref. alignment: {ref.status}, closest {ref.closest_species}"
                  + (f" {100 * ref.identity(ref.best):.2f}%" if ref.best is not None else ""), flush=True)
        rows.append((read, hits, ref))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=Path)
    ap.add_argument("--no-blast", action="store_true")
    ap.add_argument("--no-refalign", action="store_true")
    ap.add_argument("--no-send", action="store_true")
    ap.add_argument("--send", action="store_true", help="email even with --batch")
    ap.add_argument("--gtdb-representatives", type=Path, default=GTDB_REPRESENTATIVES)
    ap.add_argument("--hsp65-amplicons", type=Path, nargs="?", const=HSP65_AMPLICONS, default=None)
    args = ap.parse_args()

    OUTPUT.mkdir(parents=True, exist_ok=True)
    if args.batch:
        batches = [args.batch]
    else:
        batches = [d for d in sorted(MAIL_DIR.glob("*_*")) if d.is_dir() and not (OUTPUT / f"{d.name}.pdf").exists()]
    refs = None if args.no_refalign else load_references(args.gtdb_representatives, SPECIES)
    refs_hsp65 = None
    if refs and args.hsp65_amplicons is not None:
        if not args.hsp65_amplicons.is_dir():
            sys.exit(f"--hsp65-amplicons: not a directory: {args.hsp65_amplicons}")
        refs_hsp65 = load_references(args.hsp65_amplicons, SPECIES)
    results = {}
    attachments = []
    for b in batches:
        print(f"Analysing {b.name}")
        pdf_dir = OUTPUT / b.name
        rows = analyse_batch(b, not args.no_blast, refs, pdf_dir, refs_hsp65)
        attachments += sorted(pdf_dir.glob("*_alignment.pdf")) if pdf_dir.is_dir() else []
        if rows:
            results[b.name] = rows
    if not results:
        print("Nothing to report.")
        return
    pdf = build_pdf(results, OUTPUT / f"report_{datetime.now():%Y%m%d_%H%M}.pdf", hsp65_all_genomes=refs_hsp65 is not None)
    print(f"Wrote {pdf}")
    if not args.no_send and (args.send or not args.batch):
        send_report(load_mail_config(REQUIRED_SMTP), pdf, "Microsynth Sanger report",
                    f"Attached: QC and species identification for {len(results)} batch(es), "
                    "plus one alignment PDF per read against the kansasii-complex reference strains.",
                    extra=attachments)
        print("Report emailed.")
    if not args.batch:  # mark mail batches as done only after the mail went out (a failed send is retried)
        for name in results:
            marker = OUTPUT / f"{name}.pdf"
            if not marker.exists():
                marker.symlink_to(pdf.name)


if __name__ == "__main__":
    main()
