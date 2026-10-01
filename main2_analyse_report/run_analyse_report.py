"""QC + species ID for Microsynth batches, PDF report, optional email.

Default: every batch folder in data/sanger/microsynth_mail/ without a report yet.
  --batch DIR   analyse one folder of .ab1 files (report is not emailed unless --send)
  --no-blast    skip the NCBI BLAST species identification
  --no-send     build the PDF but don't email it
"""
import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sanger_ms import MAIL_DIR, OUTPUT
from sanger_ms.config import REQUIRED_SMTP, load_mail_config
from sanger_ms.identify import blast_read
from sanger_ms.report import build_pdf
from sanger_ms.sanger_io import load_read
from sanger_ms.send import send_report


def analyse_batch(batch: Path, blast: bool):
    rows = []
    for ab1 in sorted(batch.rglob("*.ab1")):
        read = load_read(ab1)
        if read is None:
            print(f"  skipped (unparsable or no qualities): {ab1.name}")
            continue
        t0 = time.time()
        hits = blast_read(read.trimmed_seq, read.locus) if blast else []
        top = f"{hits[0].title[:50]} ({hits[0].identity:.1f}%)" if hits else "-"
        print(f"  {ab1.name}: {read.locus}, {read.trimmed_length} bp, {top} [{time.time() - t0:.0f}s]", flush=True)
        rows.append((read, hits))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=Path)
    ap.add_argument("--no-blast", action="store_true")
    ap.add_argument("--no-send", action="store_true")
    ap.add_argument("--send", action="store_true", help="email even with --batch")
    args = ap.parse_args()

    OUTPUT.mkdir(parents=True, exist_ok=True)
    if args.batch:
        batches = [args.batch]
    else:
        batches = [d for d in sorted(MAIL_DIR.glob("*_*")) if d.is_dir() and not (OUTPUT / f"{d.name}.pdf").exists()]
    results = {}
    for b in batches:
        print(f"Analysing {b.name}")
        rows = analyse_batch(b, blast=not args.no_blast)
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
                    f"Attached: QC and species identification for {len(results)} batch(es).")
        print("Report emailed.")


if __name__ == "__main__":
    main()
