"""Minimal KEY=VALUE parser for the mail credentials file."""
from __future__ import annotations

from . import MAIL_ENV

REQUIRED_IMAP = ("IMAP_HOST", "IMAP_USER", "IMAP_PASSWORD")
REQUIRED_SMTP = ("SMTP_HOST", "SMTP_USER", "SMTP_PASSWORD", "MAIL_FROM", "MAIL_TO")


def load_mail_config(required: tuple[str, ...]) -> dict[str, str]:
    if not MAIL_ENV.exists():
        raise SystemExit(f"Missing {MAIL_ENV}; copy mail.env.example there and fill it in.")
    cfg = {"IMAP_FOLDER": "INBOX", "SENDER_FILTER": "microsynth", "SMTP_PORT": "587"}
    for line in MAIL_ENV.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip()
    missing = [k for k in required if not cfg.get(k)]
    if missing:
        raise SystemExit(f"{MAIL_ENV} is missing: {', '.join(missing)}")
    return cfg
