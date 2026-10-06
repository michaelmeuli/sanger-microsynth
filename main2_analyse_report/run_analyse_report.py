"""QC + species ID for Microsynth batches, PDF report, optional email.

Default: every batch folder in data/sanger/microsynth_mail/ without a report yet.
  --batch DIR   analyse one folder of .ab1 files (report is not emailed unless --send)
  --no-blast    skip the NCBI BLAST species identification
  --no-refalign skip the alignment to the 7 kansasii-complex reference strains
  --no-send     build the PDF but don't email it
"""
from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sanger_ms import MAIL_DIR, OUTPUT, REFERENCES, SPECIES
from sanger_ms.config import REQUIRED_SMTP, load_mail_config
from sanger_ms.identify import Hit, blast_read
from sanger_ms.refalign import ReadResult, Reference, identify_read, load_references, write_pdf
from sanger_ms.report import build_pdf
from sanger_ms.sanger_io import Read, load_read
from sanger_ms.send import send_report


MIN_ALIGNED_IDENTITY = 0.90  # below this the read is not in the complex: no alignment PDF


def analyse_batch(batch: Path, blast: bool, refs: list[Reference] | None,
                  pdf_dir: Path) -> list[tuple[Read, list[Hit], ReadResult | None]]:
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
            ref = identify_read(refs, Path(read.name).stem, read.locus, read.trimmed_seq)
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
    args = ap.parse_args()

    OUTPUT.mkdir(parents=True, exist_ok=True)
    if args.batch:
        batches = [args.batch]
    else:
        batches = [d for d in sorted(MAIL_DIR.glob("*_*")) if d.is_dir() and not (OUTPUT / f"{d.name}.pdf").exists()]
    refs = None if args.no_refalign else load_references(REFERENCES, SPECIES)
    results = {}
    attachments = []
    for b in batches:
        print(f"Analysing {b.name}")
        pdf_dir = OUTPUT / b.name
        rows = analyse_batch(b, not args.no_blast, refs, pdf_dir)
        attachments += sorted(pdf_dir.glob("*_alignment.pdf")) if pdf_dir.is_dir() else []
        if rows:
            results[b.name] = rows
    if not results:
        print("Nothing to report.")
        return
    pdf = build_pdf(results, OUTPUT / f"report_{datetime.now():%Y%m%d_%H%M}.pdf")
    if not args.batch:  # mark mail batches as done
        for name in results:
            marker = OUTPUT / f"{name}.pdf"
            if not marker.exists():
                marker.symlink_to(pdf.name)
    print(f"Wrote {pdf}")
    if not args.no_send and (args.send or not args.batch):
        send_report(load_mail_config(REQUIRED_SMTP), pdf, "Microsynth Sanger report",
                    f"Attached: QC and species identification for {len(results)} batch(es), "
                    "plus one alignment PDF per read against the kansasii-complex reference strains.",
                    extra=attachments)
        print("Report emailed.")


if __name__ == "__main__":
    main()
