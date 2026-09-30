"""Download attachments of new Microsynth result mails into data/sanger/microsynth_mail/."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sanger_ms.config import REQUIRED_IMAP, load_mail_config
from sanger_ms.mail import fetch_new

if __name__ == "__main__":
    new = fetch_new(load_mail_config(REQUIRED_IMAP))
    print(f"{len(new)} new batch(es) with .ab1 files")
    for d in new:
        print(d)
