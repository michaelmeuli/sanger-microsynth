# sanger-microsynth

Fetch the Sanger results that [Microsynth](https://www.microsynth.com) mails
us, QC and identify the reads, and email a PDF report.

Microsynth has no public API (none is documented; results are emailed and kept
for at least 3 months in the webshop under "Status/Download"). Results come as
FASTA plus `.ab1` chromatograms. In the webshop under "Options & Preferences"
the sequences can be delivered trimmed or untrimmed; untrimmed is preferable
here because we trim ourselves.

## Mail access

UZH mail is Office 365 and only allows OAuth2 (no passwords, no app passwords),
so scripts can't log in directly. Instead an Outlook rule forwards Microsynth
mails (with attachments) to the Gmail account, and this pipeline reads and
sends through Gmail with an app password. (UZH may restrict auto-forwarding;
if the rule is blocked, save attachments into `data/sanger/microsynth_mail/<batch>/`
by hand and run `scripts/02_analyse_report.py` on them.)

## Layout

- `sanger_ms/` — `mail.py` (IMAP fetch, zip-slip-safe unzip, read-only via
  `BODY.PEEK`), `sanger_io.py` (AB1 + Mott trimming, same as mlsa-kansasii),
  `identify.py` (NCBI BLAST of the trimmed read), `refalign.py` (alignment of the read to the 7
  kansasii-complex reference strains, PDF per read; copy of mlsa-kansasii's `mlsa/refalign.py`), `report.py` (PDF, reportlab),
  `send.py` (SMTP), `config.py`.
- `scripts/01_fetch_mail.py` — downloads new result mails' attachments to
  `data/sanger/microsynth_mail/<date>_<msgid>/`.
- `scripts/02_analyse_report.py` — QC + species ID per batch, PDF to
  `output/sanger-microsynth/`, then emails it. References are flags:
  `--gtdb-representatives` (default `GTDB_REPRESENTATIVES`, the 7 representative genomes, for all reads) and the
  opt-in `--hsp65-amplicons [DIR]` (default dir `HSP65_AMPLICONS`, 72 amplicons, for hsp65 reads; see below).
- `sanger_ms/kansasii_paths.py` — shared project constants (vendored; same names as in every repo). Repo-specific
  paths (`MAIL_DIR`, `STATE_FILE`, `OUTPUT`, `MAIL_ENV`) are in `sanger_ms/__init__.py`.
- `submit_fetch_analyse.sbatch` — both steps as a job.

## Setup

```bash
conda env create -f environment.yml -p /home/mimeul/data/conda/envs/env_sanger_ms
mkdir -p ~/.config/sanger-microsynth && cp mail.env.example ~/.config/sanger-microsynth/mail.env
chmod 600 ~/.config/sanger-microsynth/mail.env   # then fill in the Gmail address and app password
```

## Run

```bash
conda activate env_sanger_ms
python scripts/01_fetch_mail.py
python scripts/02_analyse_report.py            # all new batches, emails the PDF
python scripts/02_analyse_report.py --batch DIR --no-blast   # local test, no mail
sbatch submit_fetch_analyse.sbatch                           # both steps
```

Batches that already have a `<batch>.pdf` marker in the output folder are skipped.

## Reference alignment and calls

Besides BLAST, every read is aligned to references of the *M. kansasii* complex:

- **Default, all reads: the 7 representative genomes** (`GTDB_REPRESENTATIVES` = `kansasii_complex_gtdb_representatives/`).
  This is the set validated against whole-genome species (hsp65-refalign: 110 of 136 isolates correct, 1 wrong, a mixed
  culture; 25 no call). One consequence: the persicum-like hsp65 of the Korean *M. kansasii* lineage is 1-3 differences from
  persicum and 10-11 from typical kansasii, so a read of such a strain would be called **persicum with `ok`**. It was not
  found among the lab isolates (0 of 121 kansasii) and in 1 of 880 public GTDB-kansasii assemblies. Confirm persicum calls
  with gyrA where this matters.
- **Option `--hsp65-amplicons`** (use it only if the lab sequences gyrA): hsp65 reads are aligned to the in-silico hsp65
  amplicon of each of the 72 GTDB genomes (`HSP65_AMPLICONS` = `data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_hsp65_amplicons/`,
  made by mlsa-kansasii `scripts/06_make_hsp65_references.py`). The three Korean genomes are marked "atypical hsp65", so
  a Korean-type read is `ambiguous` and gets the action *confirm with gyrA* instead of an `ok` persicum call. Cost: the
  persicum margin falls from a median of 7 to 3 differences, and 2 of 28 persicum reads in the lab data become
  `ambiguous` (also *confirm with gyrA*). Other loci (16S, unknown) always use the 7.

Call = reference with the fewest differences over the read; `ok` = &ge;98% identity and &ge;2 differences fewer than
the closest other species. `ambiguous` = tie, `divergent` = closest &lt;98%. One `<read>_alignment.pdf` per read
(the 7 reference strains: ATCC 12478, AFPC-000227, MK142, MK13, MK41, 241/15, DSM 43505; plus the Korean
atypical-hsp65 genome only when it is the closest hit, which can only happen with `--hsp65-amplicons`) goes to `output/sanger-microsynth/<batch>/` and is attached to the mail. No PDF below
90% identity. `--no-refalign` skips this.

The **Action** column and the per-sample table apply these rules (`sanger_ms/call.py`), validated in
mlsa-kansasii sanger_vs_wgs (hsp65: 78/79 calls correct against whole-genome species):

| Situation | Action |
|---|---|
| &ge;98% identity, &ge;2 differences ahead, mixed peaks &lt;10%, read &ge;200 bp | report |
| species not separated (margin < 2), or, with `--hsp65-amplicons`, closest = atypical-hsp65 kansasii genome | confirm with gyrA |
| no reference &ge;98% identical, or read &lt;200 bp | repeat sequencing |
| &ge;10% of base calls with a secondary peak | possible mixed culture: repeat from pure colony |
| 16S read | not species-informative for this complex |

**Forward and reverse reads** (TB-11 / TB-12(w)) are paired by the 10-digit sample number in the file name and
locus. Every read is still called on its own; the sample call needs all usable reads to agree (in the lab data,
51 samples have both reads and none disagree), otherwise the sample is "repeat sequencing, reads disagree".
No consensus sequence is built: each read covers about 85% of the 441 bp amplicon and the reads agree.

## Notes

- Only the trimmed sequence is sent to NCBI BLAST (16S: `refseq_rna`, since `16S_ribosomal_RNA` hangs remotely;
  otherwise `core_nt`). Use `--no-blast` to keep data local.
- Species confidence: high = ≥99% identity and ≥90% query coverage, medium =
  ≥97% / ≥80%. For mycobacteria, a BLAST top hit on 16S cannot separate
  close species (e.g. the *M. kansasii* complex); see mlsa-kansasii.
- Mail credentials are in `~/.config/sanger-microsynth/mail.env`, never in the repo.

## Windows checkout

- `.gitattributes` forces LF line endings, so a Windows checkout (even with
  `core.autocrlf=true`) keeps scripts runnable. Recommended: `git config core.autocrlf false`
  and `git config core.longpaths true`.
- Scripts default to the cluster root `/shares/sander.imm.uzh/MM/kansasii`; set the
  `KANSASII_ROOT` environment variable to point elsewhere (e.g. a mapped drive).
- Nextflow, Singularity and sbatch steps only run on the cluster (or WSL).

