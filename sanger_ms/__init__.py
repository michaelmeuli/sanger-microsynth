"""Repo-specific constants; project-wide paths and species are in `kansasii_paths` (shared, vendored)."""
from pathlib import Path

from .kansasii_paths import DATA, OUTPUT_ROOT

# Attachments of the fetched Microsynth mails, one subfolder per mail.
MAIL_DIR = DATA / "sanger" / "microsynth_mail"
STATE_FILE = MAIL_DIR / "processed_messages.json"
OUTPUT = OUTPUT_ROOT / "sanger-microsynth"
# Credentials live outside the repo and the shares.
MAIL_ENV = Path.home() / ".config" / "sanger-microsynth" / "mail.env"
