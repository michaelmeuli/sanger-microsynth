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
  `identify.py` (NCBI BLAST of the trimmed read), `report.py` (PDF, reportlab),
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

## Notes

- Only the trimmed sequence is sent to NCBI BLAST (16S: `refseq_rna`, since `16S_ribosomal_RNA` hangs remotely;
  otherwise `core_nt`). Use `--no-blast` to keep data local.
- Species confidence: high = ≥99% identity and ≥90% query coverage, medium =
  ≥97% / ≥80%. For mycobacteria, a BLAST top hit on 16S cannot separate
  close species (e.g. the *M. kansasii* complex); see mlsa-kansasii.
- Mail credentials are in `~/.config/sanger-microsynth/mail.env`, never in the repo.
