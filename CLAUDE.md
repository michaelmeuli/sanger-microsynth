# CLAUDE.md

See README.md. Easy to get wrong:

- **Environment:** `env_sanger_ms` (system python has no pandas/Biopython/reportlab).
- **Credentials** are in `~/.config/sanger-microsynth/mail.env`. Never read them
  into the conversation, log them, or commit them.
- **Sending mail is outward-facing.** `scripts/02_analyse_report.py` emails by default;
  use `--no-send` (or `--batch`) when testing.
- **Mailbox is never modified:** IMAP is opened read-only with `BODY.PEEK`;
  processed mails are tracked by Message-ID in
  `data/sanger/microsynth_mail/processed_messages.json`.
- **Paths:** shared constants (`DATA`, `OUTPUT_ROOT`, `GTDB_REPRESENTATIVES`, `HSP65_AMPLICONS`, `SPECIES`) come from
  `sanger_ms/kansasii_paths.py`, a vendored copy of `kansasii-docs/templates/kansasii_paths.py` (never edit it here;
  `bash repos/check_sync.sh --fix`); mail paths are in `sanger_ms/__init__.py`.
- The `.ab1` trimming in `sanger_io.py` is copied from mlsa-kansasii; keep them in sync.
- `sanger_ms/refalign.py` is copied from mlsa-kansasii (`mlsa/refalign.py`); keep them in sync.
