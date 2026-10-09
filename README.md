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
by hand and run `main2` on them.)

## Layout

- `sanger_ms/` — `mail.py` (IMAP fetch, zip-slip-safe unzip, read-only via
  `BODY.PEEK`), `sanger_io.py` (AB1 + Mott trimming, same as mlsa-kansasii),
  `identify.py` (NCBI BLAST of the trimmed read), `refalign.py` (alignment of the read to the 7
  kansasii-complex reference strains, PDF per read; copy of mlsa-kansasii's `mlsa/refalign.py`), `report.py` (PDF, reportlab),
  `send.py` (SMTP), `config.py`.
- `main1_fetch_mail/` — downloads new result mails' attachments to
  `/shares/sander.imm.uzh/MM/kansasii/data/sanger/microsynth_mail/<date>_<msgid>/`.
- `main2_analyse_report/` — QC + species ID per batch, PDF to
  `/shares/sander.imm.uzh/MM/kansasii/output/sanger-microsynth/`, then emails it.

## Setup

```bash
conda env create -f environment.yml -p /home/mimeul/data/conda/envs/env_sanger_ms
mkdir -p ~/.config/sanger-microsynth && cp mail.env.example ~/.config/sanger-microsynth/mail.env
chmod 600 ~/.config/sanger-microsynth/mail.env   # then fill in the Gmail address and app password
```

## Run

```bash
conda activate env_sanger_ms
python main1_fetch_mail/run_fetch_mail.py
python main2_analyse_report/run_analyse_report.py            # all new batches, emails the PDF
python main2_analyse_report/run_analyse_report.py --batch DIR --no-blast   # local test, no mail
sbatch submit_fetch_analyse.sbatch                           # both steps
```

Batches that already have a `<batch>.pdf` marker in the output folder are skipped.

## Reference alignment and calls

Besides BLAST, every read is aligned to references of the *M. kansasii* complex:

- **hsp65 reads**: the in-silico hsp65 amplicon of each of the 72 GTDB genomes of the 7 species
  (`data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_hsp65_amplicons/`, made by mlsa-kansasii
  `scripts/make_hsp65_references.py`). This covers within-species diversity, and the three Korean *M. kansasii*
  genomes with a persicum-like hsp65 are marked as "atypical hsp65".
- **other reads (16S, unknown locus)**: the 7 representative genomes
  (`kansasii_complex_gtdb_representatives/`).

Call = reference with the fewest differences over the read; `ok` = &ge;98% identity and &ge;2 differences fewer than
the closest other species. `ambiguous` = tie, `divergent` = closest &lt;98%. One `<read>_alignment.pdf` per read
(best genome per species) goes to `output/sanger-microsynth/<batch>/` and is attached to the mail. No PDF below
90% identity. `--no-refalign` skips this.

The **Action** column and the per-sample table apply these rules (`sanger_ms/call.py`), validated in
mlsa-kansasii main4 (hsp65: 78/79 calls correct against whole-genome species):

| Situation | Action |
|---|---|
| &ge;98% identity, &ge;2 differences ahead, mixed peaks &lt;10%, read &ge;200 bp | report |
| closest = atypical-hsp65 kansasii genome, or species not separated | confirm with gyrA |
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

## Windows checkout

- `.gitattributes` forces LF line endings, so a Windows checkout (even with
  `core.autocrlf=true`) keeps scripts runnable. Recommended: `git config core.autocrlf false`
  and `git config core.longpaths true`.
- Scripts default to the cluster root `/shares/sander.imm.uzh/MM/kansasii`; set the
  `KANSASII_ROOT` environment variable to point elsewhere (e.g. a mapped drive).
- Nextflow, Singularity and sbatch steps only run on the cluster (or WSL).
