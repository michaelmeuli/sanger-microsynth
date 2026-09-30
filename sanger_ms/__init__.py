from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SHARES = REPO_ROOT.parent.parent
# Attachments of the fetched Microsynth mails, one subfolder per mail.
MAIL_DIR = SHARES / "data" / "sanger" / "microsynth_mail"
OUTPUT = SHARES / "output" / "sanger-microsynth"
STATE_FILE = MAIL_DIR / "processed_messages.json"
# Credentials live outside the repo and the shares.
MAIL_ENV = Path.home() / ".config" / "sanger-microsynth" / "mail.env"
