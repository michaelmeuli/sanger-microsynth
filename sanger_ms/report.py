"""PDF report: summary table, per-read quality plots, species hits."""
from __future__ import annotations

import io
from datetime import date
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from .identify import Hit
from .call import ReadCall, SampleCall, read_call, sample_calls
from .refalign import ReadResult
from .sanger_io import Read

QUALITY_OK = 30.0  # mean post-trim Phred considered a good read
OK_COLOR, LOW_COLOR = "#1b7f3b", "#b3261e"


def _quality_plot(read: Read) -> Image:
    fig, ax = plt.subplots(figsize=(9, 1.6))
    ax.plot(read.qualities, lw=0.6, color="#444444")
    ax.axvspan(read.trim_start, read.trim_end, color="#9ecae1", alpha=0.4, label="kept after trimming")
    ax.axhline(20, color=LOW_COLOR, lw=0.6, ls="--")
    ax.set_xlim(0, len(read.qualities))
    ax.set_ylabel("Phred")
    ax.legend(loc="lower right", fontsize=7, frameon=False)
    ax.tick_params(labelsize=7)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return Image(buf, width=22 * cm, height=3.9 * cm)


def _ref_text(ref: ReadResult | None) -> str:
    b = ref.best if ref is not None else None
    if ref is None or b is None:
        return "no reference hit"
    return f"M. {ref.closest_species} {100 * ref.identity(b):.2f}% ({ref.status})"


def _sample_table(samples: list[SampleCall], small: Any) -> Table:
    """One row per sample and locus: forward/reverse reads combined into a single call."""
    data: list[list[Any]] = [["Sample", "Locus", "Reads", "Species", "Action", "Note"]]
    for sc in samples:
        data.append([sc.sample, sc.locus, Paragraph(escape("; ".join(sc.reads)), small),
                     f"M. {sc.species}" if sc.species else "", Paragraph(escape(sc.action), small),
                     Paragraph(escape(sc.note), small)])
    t = Table(data, repeatRows=1, colWidths=[2.4 * cm, 1.4 * cm, 9 * cm, 3.6 * cm, 5 * cm, 5 * cm])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eaed")),
                           ("FONTSIZE", (0, 0), (-1, -1), 7.5), ("VALIGN", (0, 0), (-1, -1), "TOP"),
                           ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#c4c7c5"))]))
    return t


def build_pdf(batches: dict[str, list[tuple[Read, list[Hit], ReadResult | None]]], out: Path) -> Path:
    styles = getSampleStyleSheet()
    small = styles["BodyText"].clone("small", fontSize=7.5, leading=9)
    doc = SimpleDocTemplate(str(out), pagesize=landscape(A4), leftMargin=1.5 * cm,
                            rightMargin=1.5 * cm, topMargin=1.5 * cm, bottomMargin=1.5 * cm,
                            title="Microsynth Sanger report")
    story = [Paragraph("Microsynth Sanger sequencing report", styles["Title"]),
             Paragraph(f"Generated {date.today().isoformat()}. Quality = Mott-trimmed Phred; "
                       "species = top NCBI BLAST hit of the trimmed read (confidence: high "
                       "&ge;99% identity and &ge;90% coverage, medium &ge;97%/80%). "
                       "Ref. alignment = closest of the 7 kansasii-complex reference strains by plain alignment "
                       "(call: ok = &ge;98% identity and &ge;2 fewer differences than the next species; hsp65 reads are aligned "
                       "to the hsp65 amplicon of 72 genomes, other reads to the 7 representative genomes). "
                       "Action: report = species can be reported; confirm = add gyrA; repeat = resequence; "
                       "mixed = &ge;10% secondary peaks. Forward and reverse reads of a sample are combined in "
                       "the sample table; "
                       "details in the attached &lt;read&gt;_alignment.pdf files.", styles["BodyText"]),
             Spacer(1, 0.4 * cm)]
    for batch, rows in batches.items():
        story.append(Paragraph(f"Batch {escape(batch)}", styles["Heading2"]))
        calls = [read_call(read, ref) for read, _, ref in rows]
        samples = sample_calls([(read, c) for (read, _, _), c in zip(rows, calls)])
        story += [_sample_table(samples, small), Spacer(1, 0.4 * cm)]
        header = ["Read", "Locus", "Raw bp", "Kept bp", "Mean Q", "% Q20", "N", "Top hit", "Ident %", "Cov %", "Conf.", "Ref. alignment", "Action"]
        data: list[list[Any]] = [header]
        style: list[tuple[Any, ...]] = [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eaed")),
                 ("FONTSIZE", (0, 0), (-1, -1), 7.5), ("VALIGN", (0, 0), (-1, -1), "TOP"),
                 ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#c4c7c5"))]
        for i, (read, hits, ref) in enumerate(rows, start=1):
            top = hits[0] if hits else None
            data.append([
                Paragraph(escape(read.name), small), read.locus or "?", read.raw_length,
                read.trimmed_length, f"{read.mean_q:.1f}", f"{100 * read.frac_q20:.0f}", read.n_count,
                Paragraph(escape(top.title), small) if top
                else ("failed read (<100 bp after trimming)" if read.trimmed_length < 100 else "not identified"),
                f"{top.identity:.1f}" if top else "", f"{top.coverage:.0f}" if top else "",
                top.confidence if top else "",
                Paragraph(escape(_ref_text(ref)), small),
                Paragraph(escape(calls[i - 1].action), small),
            ])
            ok = read.mean_q >= QUALITY_OK
            style.append(("TEXTCOLOR", (4, i), (4, i), colors.HexColor(OK_COLOR if ok else LOW_COLOR)))
        table = Table(data, repeatRows=1, colWidths=[4.3 * cm, 1.2 * cm, 1.2 * cm, 1.2 * cm, 1.2 * cm,
                                                      1.1 * cm, 0.7 * cm, 4.8 * cm, 1.2 * cm, 1.1 * cm, 1.3 * cm, 3.4 * cm, 3.4 * cm])
        table.setStyle(TableStyle(style))
        story += [table, PageBreak()]
        for read, _, _ in rows:
            story += [Paragraph(escape(read.name), styles["Heading4"]), _quality_plot(read), Spacer(1, 0.2 * cm)]
        story.append(PageBreak())
    doc.build(story[:-1])
    return out
