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

## Reference alignment

Besides BLAST, every read is aligned to the GTDB representative genome of each
*M. kansasii*-complex species
(`data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_gtdb_representatives/`, 7 genomes).
The closest reference (fewest differences over the read) is shown in the report table
("Ref. alignment"), and one `<read>_alignment.pdf` per read (pairwise difference matrix +
alignment of the read with all references, mismatches highlighted) goes to
`output/sanger-microsynth/<batch>/` and is attached to the mail. Call: `ok` = closest reference
&ge;99% identical and &ge;2 differences fewer than the next; `ambiguous` = tie; `divergent` = closest
&lt;99% (outside the references' diversity, or a bad read). No PDF below 90% identity. `--no-refalign` skips this.
Only 1 genome per species: intra-species diversity is not covered (see mlsa-kansasii README, main3).

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
